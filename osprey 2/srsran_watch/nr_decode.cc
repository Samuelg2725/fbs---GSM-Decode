/*
 * nr_decode - 5G NR cell identity decoding on top of the srsRAN 4G libraries (used by nr_decode.py).
 *
 * The Python side does everything radio and timing related (search, PCI detection, timing, frequency offset,
 * per-symbol phase, OFDM demodulation with the correct NR cyclic prefixes, combining of the two antennas) and
 * hands this program resource grids. This program only decodes:
 *
 *   nr_decode mib GRIDS LMAX SCS PCI [PCI ...]      (SCS: SS/PBCH block subcarrier spacing, 15 or 30 kHz)
 *       GRIDS: complex float32 file, one SS/PBCH block per PCI given (4 symbols x 240 subcarriers, lowest
 *       subcarrier first). Every SSB index (0..LMAX-1) and half frame is tried; the PBCH CRC decides.
 *       Prints per grid:  MIB grid=.. pci=.. ssb_idx=.. hrf=.. sfn=.. scs_common=.. kssb=.. coreset0=.. ss0=..
 *                         typeA=.. barred=.. intra=.. bch=<BCCH-BCH message hex>     (or NOMIB grid=..)
 *
 *   nr_decode sib1 GRIDS PCI SCS NRB DURATION TYPEA SLOTS
 *       GRIDS: complex float32 file of slot grids (14 symbols x NRB*12 subcarriers) whose lowest resource block
 *       is the lowest resource block of CORESET 0 (so CORESET 0 and SIB1 lie inside, NRB = CORESET 0 size).
 *       SCS: 15 or 30 (kHz). DURATION: CORESET 0 symbols. TYPEA: 2 or 3 (dmrs-TypeA-Position from the MIB).
 *       SLOTS: comma list, one entry per grid: slot_index[:ssb_rb_start:ssb_rb_count:ssb_symbol_mask_hex]
 *       (the SS/PBCH block REs are not used by SIB1's PDSCH: given so they are rate-matched around).
 *       For each slot: PDCCH blind search for SI-RNTI in CORESET 0 / search space 0, then the PDSCH.
 *       Prints: DCI slot=.. L=.. ncce=.. mcs=.. rbs=..  /  SIB1 slot=.. hex=..  /  CELL nci=.. tac=.. plmn=..
 *
 * Built by print_sib1.py into srsRAN_4G/build/lib/examples/nr_decode.
 */
#include <complex.h>
extern "C" {
#include "srsran/phy/ch_estimation/dmrs_pbch.h"
#include "srsran/phy/phch/pbch_msg_nr.h"
#include "srsran/phy/phch/pbch_nr.h"
#include "srsran/phy/phch/ra_nr.h"
#include "srsran/phy/ue/ue_dl_nr.h"
#include "srsran/phy/utils/debug.h"
#include "srsran/phy/utils/re_pattern.h"
#include "srsran/phy/utils/vector.h"
}
#include "srsran/asn1/rrc_nr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

static std::vector<cf_t> read_file(const char* path)
{
  std::vector<cf_t> out;
  FILE*             f = fopen(path, "rb");
  if (f == NULL) {
    fprintf(stderr, "cannot open %s\n", path);
    exit(2);
  }
  fseek(f, 0, SEEK_END);
  long n = ftell(f) / (long)sizeof(cf_t);
  fseek(f, 0, SEEK_SET);
  out.resize((size_t)n);
  if (fread(out.data(), sizeof(cf_t), (size_t)n, f) != (size_t)n) {
    fprintf(stderr, "short read %s\n", path);
    exit(2);
  }
  fclose(f);
  return out;
}

static std::string hex(const uint8_t* b, uint32_t n)
{
  std::string s;
  char        t[3];
  for (uint32_t i = 0; i < n; i++) {
    snprintf(t, sizeof(t), "%02x", b[i]);
    s += t;
  }
  return s;
}

// ------------------------------------------------------------------------------------------------- MIB
static int cmd_mib(int argc, char** argv)
{
  if (argc < 6) {
    fprintf(stderr, "usage: nr_decode mib GRIDS LMAX SCS PCI [PCI ...]\n");
    return 1;
  }
  std::vector<cf_t> grids = read_file(argv[2]);
  uint32_t          Lmax  = (uint32_t)atoi(argv[3]);
  srsran_subcarrier_spacing_t ssb_scs =
      atoi(argv[4]) == 30 ? srsran_subcarrier_spacing_30kHz : srsran_subcarrier_spacing_15kHz;
  uint32_t n        = (uint32_t)(argc - 5);
  float    min_corr = getenv("NR_MIB_MIN_CORR") ? (float)atof(getenv("NR_MIB_MIN_CORR")) : 0.1f;
  if (grids.size() < (size_t)n * SRSRAN_SSB_NOF_RE) {
    fprintf(stderr, "grid file holds %zu values, need %u\n", grids.size(), n * SRSRAN_SSB_NOF_RE);
    return 1;
  }
  srsran_pbch_nr_t      pbch = {};
  srsran_pbch_nr_args_t pa   = {};
  pa.enable_decode           = true;
  if (srsran_pbch_nr_init(&pbch, &pa) < SRSRAN_SUCCESS) {
    fprintf(stderr, "pbch init failed\n");
    return 1;
  }
  for (uint32_t g = 0; g < n; g++) {
    uint32_t    pci  = (uint32_t)atoi(argv[5 + g]);
    const cf_t* grid = &grids[(size_t)g * SRSRAN_SSB_NOF_RE];
    bool        done = false;
    for (uint32_t ssb_idx = 0; ssb_idx < Lmax && !done; ssb_idx++) {
      for (uint32_t hf = 0; hf < 2 && !done; hf++) {
        srsran_dmrs_pbch_cfg_t dc = {};
        dc.N_id                   = pci;
        dc.n_hf                   = hf;
        dc.ssb_idx                = ssb_idx;
        dc.L_max                  = Lmax;
        dc.scs                    = ssb_scs;
        cf_t                    ce[SRSRAN_SSB_NOF_RE] = {};
        srsran_dmrs_pbch_meas_t meas                  = {};
        srsran_dmrs_pbch_measure(&dc, grid, &meas);
        if (srsran_dmrs_pbch_estimate(&dc, grid, ce) < SRSRAN_SUCCESS) {
          continue;
        }
        srsran_pbch_nr_cfg_t pc = {};
        pc.N_id                 = pci;
        pc.n_hf                 = hf;
        pc.ssb_idx              = ssb_idx;
        pc.Lmax                 = Lmax;
        srsran_pbch_msg_nr_t msg = {};
        if (srsran_pbch_nr_decode(&pbch, &pc, grid, ce, &msg) < SRSRAN_SUCCESS || !msg.crc) {
          continue;
        }
        // A weak or wrong block can decode to all zeros, which passes the CRC: a real PBCH also has its DMRS
        // (pilots) clearly present (normalised correlation >= 0.3 even at -2 dB SNR; false decodes <= 0.02)
        if (meas.corr < min_corr) {
          continue;
        }
        srsran_mib_nr_t mib = {};
        if (srsran_pbch_msg_nr_mib_unpack(&msg, &mib) < SRSRAN_SUCCESS) {
          printf("NOTMIB grid=%u pci=%u\n", g, pci);
          done = true;
          break;
        }
        uint8_t bch[3] = {};
        for (uint32_t i = 0; i < 24; i++) {
          bch[i / 8] |= (uint8_t)((msg.payload[i] & 1U) << (7 - i % 8));
        }
        printf("MIB grid=%u pci=%u ssb_idx=%u hrf=%u sfn=%u scs_common=%u kssb=%u coreset0=%u ss0=%u typeA=%u "
               "barred=%u intra=%u corr=%.3f cfo=%.1f bch=%s\n",
               g,
               pci,
               (uint32_t)msg.ssb_idx,
               (uint32_t)msg.hrf,
               mib.sfn,
               mib.scs_common == srsran_subcarrier_spacing_30kHz ? 30U : 15U,
               mib.ssb_offset,
               mib.coreset0_idx,
               mib.ss0_idx,
               mib.dmrs_typeA_pos == srsran_dmrs_sch_typeA_pos_3 ? 3U : 2U,
               mib.cell_barred ? 1U : 0U,
               mib.intra_freq_reselection ? 1U : 0U,
               meas.corr,
               meas.cfo_hz,
               hex(bch, 3).c_str());
        done = true;
      }
    }
    if (!done) {
      printf("NOMIB grid=%u pci=%u\n", g, pci);
    }
    fflush(stdout);
  }
  srsran_pbch_nr_free(&pbch);
  return 0;
}

// ------------------------------------------------------------------------------------------------- SIB1
template <class B>
static std::string bits_u(const B& b)
{
  // the bit string as text ("0101..."), then to a number (to_number() sign-extends past 32 bits)
  return std::to_string(std::stoull(b.to_string(), nullptr, 2));
}

static void print_sib1(const uint8_t* data, uint32_t len, uint32_t slot)
{
  asn1::rrc_nr::bcch_dl_sch_msg_s msg;
  asn1::cbit_ref                  bref(data, len);
  bool                            ok = msg.unpack(bref) == asn1::SRSASN_SUCCESS;
  uint32_t                        n  = ok ? SRSRAN_MIN((uint32_t)bref.distance_bytes(), len) : len;
  printf("SIB1 slot=%u hex=%s\n", slot, hex(data, n).c_str());
  if (!ok ||
      msg.msg.type().value != asn1::rrc_nr::bcch_dl_sch_msg_type_c::types_opts::c1) {
    printf("SIB1PARSE slot=%u failed\n", slot);
    return;
  }
  auto& c1 = msg.msg.c1();
  if (c1.type().value != asn1::rrc_nr::bcch_dl_sch_msg_type_c::c1_c_::types_opts::sib_type1) {
    printf("SIB1PARSE slot=%u other-SI\n", slot);
    return;
  }
  auto& sib1 = c1.sib_type1();
  for (uint32_t i = 0; i < sib1.cell_access_related_info.plmn_id_list.size(); i++) {
    auto&       info = sib1.cell_access_related_info.plmn_id_list[i];
    std::string plmns;
    for (uint32_t j = 0; j < info.plmn_id_list.size(); j++) {
      auto&       p = info.plmn_id_list[j];
      std::string mcc, mnc;
      if (p.mcc_present) {
        for (uint32_t k = 0; k < p.mcc.size(); k++) {
          mcc += (char)('0' + p.mcc[k]);
        }
      }
      for (uint32_t k = 0; k < p.mnc.size(); k++) {
        mnc += (char)('0' + p.mnc[k]);
      }
      plmns += (plmns.empty() ? "" : ",") + (mcc.empty() ? std::string("?") : mcc) + "-" + mnc;
    }
    printf("CELL slot=%u nci=%s tac=%s plmn=%s reserved=%u\n",
           slot,
           bits_u(info.cell_id).c_str(),
           info.tac_present ? bits_u(info.tac).c_str() : "-",
           plmns.c_str(),
           info.cell_reserved_for_oper.value == asn1::rrc_nr::plmn_id_info_s::cell_reserved_for_oper_opts::reserved
               ? 1U
               : 0U);
  }
  fflush(stdout);
}

static int cmd_sib1(int argc, char** argv)
{
  if (argc < 9) {
    fprintf(stderr, "usage: nr_decode sib1 GRIDS PCI SCS NRB DURATION TYPEA SLOTS\n");
    return 1;
  }
  std::vector<cf_t> grids = read_file(argv[2]);
  uint32_t          pci   = (uint32_t)atoi(argv[3]);
  uint32_t          scs   = (uint32_t)atoi(argv[4]);
  uint32_t          nrb   = (uint32_t)atoi(argv[5]);
  uint32_t          dur   = (uint32_t)atoi(argv[6]);
  uint32_t          typeA = (uint32_t)atoi(argv[7]);
  uint32_t          nre   = nrb * SRSRAN_NRE;
  size_t            slot_len = (size_t)nre * SRSRAN_NSYMB_PER_SLOT_NR;

  // slot list
  struct slot_t {
    uint32_t idx, rb0, nrb, mask;
  };
  std::vector<slot_t> slots;
  {
    std::string s(argv[8]);
    size_t      pos = 0;
    while (pos <= s.size()) {
      size_t      end  = s.find(',', pos);
      std::string item = s.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
      if (!item.empty()) {
        slot_t t = {};
        sscanf(item.c_str(), "%u:%u:%u:%x", &t.idx, &t.rb0, &t.nrb, &t.mask);
        slots.push_back(t);
      }
      if (end == std::string::npos) {
        break;
      }
      pos = end + 1;
    }
  }
  if (grids.size() < slots.size() * slot_len) {
    fprintf(stderr, "grid file holds %zu values, need %zu\n", grids.size(), slots.size() * slot_len);
    return 1;
  }

  srsran_carrier_nr_t carrier    = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci                    = pci;
  carrier.scs                    = scs == 30 ? srsran_subcarrier_spacing_30kHz : srsran_subcarrier_spacing_15kHz;
  carrier.nof_prb                = nrb;
  carrier.dl_center_frequency_hz = 0.0; // phase compensation is done by the caller
  carrier.ssb_center_freq_hz     = 0.0;

  cf_t*                  buffer[SRSRAN_MAX_PORTS] = {};
  buffer[0]                                       = srsran_vec_cf_malloc(SRSRAN_SF_LEN_PRB(nrb) * 2);
  srsran_ue_dl_nr_args_t args                     = {};
  args.nof_rx_antennas                            = 1;
  args.nof_max_prb                                = nrb;
  args.pdsch.sch.disable_simd                     = false;
  args.pdsch.sch.decoder_use_flooded              = false;
  args.pdsch.measure_evm                          = true;
  args.pdcch.measure_evm                          = true;
  srsran_ue_dl_nr_t ue_dl                         = {};
  if (srsran_ue_dl_nr_init(&ue_dl, buffer, &args) < SRSRAN_SUCCESS ||
      srsran_ue_dl_nr_set_carrier(&ue_dl, &carrier) < SRSRAN_SUCCESS) {
    fprintf(stderr, "ue_dl init failed\n");
    return 1;
  }

  // CORESET 0 = the whole grid (TS 38.211 7.3.2.2: interleaved, bundle 6, interleaver size 2, shift = PCI)
  srsran_pdcch_cfg_nr_t pdcch_cfg = {};
  srsran_coreset_t*     cs        = &pdcch_cfg.coreset[0];
  pdcch_cfg.coreset_present[0]    = true;
  cs->id                          = 0;
  cs->duration                    = dur;
  cs->offset_rb                   = 0;
  cs->mapping_type                = srsran_coreset_mapping_type_interleaved;
  cs->reg_bundle_size             = srsran_coreset_bundle_size_n6;
  cs->interleaver_size            = srsran_coreset_bundle_size_n2;
  cs->precoder_granularity        = srsran_coreset_precoder_granularity_reg_bundle;
  cs->shift_index                 = pci;
  for (uint32_t i = 0; i < SRSRAN_CORESET_FREQ_DOMAIN_RES_SIZE; i++) {
    cs->freq_resources[i] = i < nrb / 6;
  }
  srsran_search_space_t* ss       = &pdcch_cfg.search_space[0];
  pdcch_cfg.search_space_present[0] = true;
  ss->id                            = 0;
  ss->coreset_id                    = 0;
  ss->type                          = srsran_search_space_type_common_0;
  ss->formats[0]                    = srsran_dci_format_nr_1_0;
  ss->nof_formats                   = 1;
  for (uint32_t L = 0; L < SRSRAN_SEARCH_SPACE_NOF_AGGREGATION_LEVELS_NR; L++) {
    ss->nof_candidates[L] = srsran_pdcch_nr_max_candidates_coreset(cs, L);
  }
  srsran_dci_cfg_nr_t dci_cfg = {};
  dci_cfg.bwp_dl_initial_bw   = nrb;
  dci_cfg.bwp_ul_initial_bw   = nrb;
  dci_cfg.bwp_dl_active_bw    = nrb;
  dci_cfg.bwp_ul_active_bw    = nrb;
  dci_cfg.coreset0_bw         = nrb;
  dci_cfg.monitor_common_0_0  = true;
  dci_cfg.monitor_0_0_and_1_0 = true;
  if (srsran_ue_dl_nr_set_pdcch_config(&ue_dl, &pdcch_cfg, &dci_cfg) < SRSRAN_SUCCESS) {
    fprintf(stderr, "pdcch config failed\n");
    return 1;
  }
  srsran_sch_hl_cfg_nr_t hl = {};
  hl.typeA_pos              = typeA == 3 ? srsran_dmrs_sch_typeA_pos_3 : srsran_dmrs_sch_typeA_pos_2;

  srsran_softbuffer_rx_t softbuffer = {};
  srsran_softbuffer_rx_init_guru(&softbuffer, SRSRAN_SCH_NR_MAX_NOF_CB_LDPC, SRSRAN_LDPC_MAX_LEN_ENCODED_CB);
  uint8_t* data = srsran_vec_u8_malloc(SRSRAN_SLOT_MAX_NOF_BITS_NR);

  int found = 0;
  for (size_t i = 0; i < slots.size(); i++) {
    srsran_slot_cfg_t slot = {};
    slot.idx               = slots[i].idx;
    srsran_vec_cf_copy(ue_dl.sf_symbols[0], &grids[i * slot_len], (uint32_t)slot_len);
    srsran_dmrs_pdcch_estimate(&ue_dl.dmrs_pdcch[0], &slot, ue_dl.sf_symbols[0]);
    srsran_dci_dl_nr_t dci = {};
    int                n   = srsran_ue_dl_nr_find_dl_dci(&ue_dl, &slot, SRSRAN_SIRNTI, srsran_rnti_type_si, &dci, 1);
    if (n < 1) {
      continue;
    }
    srsran_sch_cfg_nr_t pdsch_cfg = {};
    // SS/PBCH block REs are not used by the PDSCH (TS 38.214 5.1.4); set before the allocation is worked out
    if (slots[i].nrb > 0 && slots[i].mask != 0) {
      srsran_re_pattern_t p = {};
      p.rb_begin            = slots[i].rb0;
      p.rb_end              = SRSRAN_MIN(slots[i].rb0 + slots[i].nrb, nrb);
      p.rb_stride           = 1;
      for (uint32_t k = 0; k < SRSRAN_NRE; k++) {
        p.sc[k] = true;
      }
      for (uint32_t l = 0; l < SRSRAN_NSYMB_PER_SLOT_NR; l++) {
        p.symbol[l] = (slots[i].mask >> l) & 1U;
      }
      srsran_re_pattern_merge(&pdsch_cfg.rvd_re, &p);
    }
    if (srsran_ra_dl_dci_to_grant_nr(&carrier, &slot, &hl, &dci, &pdsch_cfg, &pdsch_cfg.grant) < SRSRAN_SUCCESS) {
      printf("DCI slot=%u grant-error\n", slot.idx);
      continue;
    }
    uint32_t nof_rb = 0;
    for (uint32_t k = 0; k < SRSRAN_MAX_PRB_NR; k++) {
      nof_rb += pdsch_cfg.grant.prb_idx[k] ? 1 : 0;
    }
    printf("DCI slot=%u L=%u ncce=%u mcs=%u rv=%u rbs=%u symbols=%u+%u tbs=%u\n",
           slot.idx,
           dci.ctx.location.L,
           dci.ctx.location.ncce,
           dci.mcs,
           dci.rv,
           nof_rb,
           pdsch_cfg.grant.S,
           pdsch_cfg.grant.L,
           pdsch_cfg.grant.tb[0].tbs);
    srsran_softbuffer_rx_reset(&softbuffer);
    pdsch_cfg.grant.tb[0].softbuffer.rx = &softbuffer;
    srsran_pdsch_res_nr_t res           = {};
    res.tb[0].payload                   = data;
    if (srsran_ue_dl_nr_decode_pdsch(&ue_dl, &slot, &pdsch_cfg, &res) < SRSRAN_SUCCESS || !res.tb[0].crc) {
      printf("PDSCH slot=%u crc=KO\n", slot.idx);
      fflush(stdout);
      continue;
    }
    uint32_t nbytes = pdsch_cfg.grant.tb[0].tbs / 8;
    print_sib1(data, nbytes, slot.idx);
    found++;
    fflush(stdout);
  }
  srsran_softbuffer_rx_free(&softbuffer);
  free(data);
  srsran_ue_dl_nr_free(&ue_dl);
  free(buffer[0]);
  return found > 0 ? 0 : 3;
}

int main(int argc, char** argv)
{
  set_srsran_verbose_level(SRSRAN_VERBOSE_NONE);
  if (argc >= 2 && strcmp(argv[1], "mib") == 0) {
    return cmd_mib(argc, argv);
  }
  if (argc >= 2 && strcmp(argv[1], "sib1") == 0) {
    return cmd_sib1(argc, argv);
  }
  fprintf(stderr, "usage: nr_decode mib|sib1 ...\n");
  return 1;
}
