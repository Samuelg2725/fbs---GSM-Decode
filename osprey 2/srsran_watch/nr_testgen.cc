/*
 * nr_testgen - 5G NR test signal pieces made with srsRAN's own gNB encoders, to check nr_decode.py end to end.
 *
 *   nr_testgen OUT PCI SSB_SCS LMAX SFN KSSB CORESET0 SS0 TYPEA NRB DURATION SLOT NCI TAC MCC MNC [RB0 NRB MASK]
 *   (RB0 NRB MASK: SS/PBCH block resource blocks and symbol mask (hex) in the slot: the PDSCH goes around them,
 *    as a gNB does)
 *
 * Writes:
 *   OUT_ssb.cf32   LMAX SS/PBCH blocks (4 x 240 subcarriers each, lowest subcarrier first) for SSB index 0..LMAX-1,
 *                  first half frame, carrying a MIB with the given SFN, ssb-SubcarrierOffset (KSSB), CORESET 0 and
 *                  search space 0 indices and dmrs-TypeA-Position
 *   OUT_slot.cf32  one slot (14 symbols x NRB*12 subcarriers): CORESET 0 (NRB resource blocks, DURATION symbols)
 *                  with a DCI 1_0 for SI-RNTI and the PDSCH carrying a SIB1 with the given NCI, TAC and PLMN
 * nr_decode.py --selftest places these in a time signal the way a gNB transmits (cyclic prefixes, upconversion
 * phase, frequency positions), adds noise and decodes it.
 */
#include <complex.h>
extern "C" {
#include "srsran/phy/gnb/gnb_dl.h"
#include "srsran/phy/phch/pbch_msg_nr.h"
#include "srsran/phy/phch/ra_dl_nr.h"
#include "srsran/phy/phch/ra_nr.h"
#include "srsran/phy/sync/ssb.h"
#include "srsran/phy/utils/debug.h"
#include "srsran/phy/utils/re_pattern.h"
#include "srsran/phy/utils/vector.h"
}
#include "srsran/asn1/rrc_nr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <algorithm>

static void write_file(const std::string& path, const cf_t* x, size_t n)
{
  FILE* f = fopen(path.c_str(), "wb");
  if (f == NULL || fwrite(x, sizeof(cf_t), n, f) != n) {
    fprintf(stderr, "cannot write %s\n", path.c_str());
    exit(2);
  }
  fclose(f);
}

int main(int argc, char** argv)
{
  if (argc < 17) {
    fprintf(stderr, "usage: nr_testgen OUT PCI SSB_SCS LMAX SFN KSSB CORESET0 SS0 TYPEA NRB DURATION SLOT NCI TAC MCC MNC\n");
    return 1;
  }
  std::string out   = argv[1];
  uint32_t    pci   = (uint32_t)atoi(argv[2]);
  uint32_t    sscs  = (uint32_t)atoi(argv[3]);
  uint32_t    Lmax  = (uint32_t)atoi(argv[4]);
  uint32_t    sfn   = (uint32_t)atoi(argv[5]);
  uint32_t    kssb  = (uint32_t)atoi(argv[6]);
  uint32_t    c0    = (uint32_t)atoi(argv[7]);
  uint32_t    ss0   = (uint32_t)atoi(argv[8]);
  uint32_t    typeA = (uint32_t)atoi(argv[9]);
  uint32_t    nrb   = (uint32_t)atoi(argv[10]);
  uint32_t    dur   = (uint32_t)atoi(argv[11]);
  uint32_t    slotn = (uint32_t)atoi(argv[12]);
  uint64_t    nci   = strtoull(argv[13], NULL, 10);
  uint32_t    tac   = (uint32_t)atoi(argv[14]);
  std::string mcc   = argv[15];
  std::string mnc   = argv[16];
  srsran_subcarrier_spacing_t scs = sscs == 30 ? srsran_subcarrier_spacing_30kHz : srsran_subcarrier_spacing_15kHz;

  // ---------------------------------------------------------------- SS/PBCH blocks
  srsran_ssb_t      ssb  = {};
  srsran_ssb_args_t sa   = {};
  sa.enable_encode       = true;
  sa.max_srate_hz        = 61.44e6;
  sa.min_scs             = srsran_subcarrier_spacing_15kHz;
  if (srsran_ssb_init(&ssb, &sa) < SRSRAN_SUCCESS) {
    fprintf(stderr, "ssb init\n");
    return 1;
  }
  srsran_ssb_cfg_t sc = {};
  sc.srate_hz         = 61.44e6;
  sc.center_freq_hz   = Lmax == 8 ? 3.5e9 : 1.8e9;
  sc.ssb_freq_hz      = sc.center_freq_hz;
  sc.scs              = scs;
  sc.pattern          = sscs == 30 ? SRSRAN_SSB_PATTERN_C : SRSRAN_SSB_PATTERN_A;
  sc.duplex_mode      = sscs == 30 ? SRSRAN_DUPLEX_MODE_TDD : SRSRAN_DUPLEX_MODE_FDD;
  sc.periodicity_ms   = 20;
  if (srsran_ssb_set_cfg(&ssb, &sc) < SRSRAN_SUCCESS) {
    fprintf(stderr, "ssb cfg\n");
    return 1;
  }
  uint32_t          nsym = 2 * 56;
  std::vector<cf_t> big((size_t)nsym * SRSRAN_SSB_BW_SUBC);
  std::vector<cf_t> blocks((size_t)Lmax * SRSRAN_SSB_NOF_RE);
  for (uint32_t i = 0; i < Lmax; i++) {
    srsran_mib_nr_t mib = {};
    mib.sfn             = sfn;
    mib.ssb_idx         = (uint8_t)i;
    mib.hrf             = false;
    mib.scs_common      = scs;
    mib.ssb_offset      = kssb;
    mib.dmrs_typeA_pos  = typeA == 3 ? srsran_dmrs_sch_typeA_pos_3 : srsran_dmrs_sch_typeA_pos_2;
    mib.coreset0_idx    = c0;
    mib.ss0_idx         = ss0;
    mib.cell_barred     = false;
    mib.intra_freq_reselection = true;
    srsran_pbch_msg_nr_t msg   = {};
    if (srsran_pbch_msg_nr_mib_pack(&mib, &msg) < SRSRAN_SUCCESS) {
      fprintf(stderr, "mib pack\n");
      return 1;
    }
    std::fill(big.begin(), big.end(), 0);
    if (srsran_ssb_put_grid(&ssb, pci, &msg, big.data(), SRSRAN_SSB_BW_SUBC) < SRSRAN_SUCCESS) {
      fprintf(stderr, "ssb put\n");
      return 1;
    }
    uint32_t l_first = ssb.l_first[i];
    for (uint32_t l = 0; l < 4; l++) {
      srsran_vec_cf_copy(&blocks[(size_t)i * SRSRAN_SSB_NOF_RE + l * SRSRAN_SSB_BW_SUBC],
                         &big[(size_t)(l_first + l) * SRSRAN_SSB_BW_SUBC],
                         SRSRAN_SSB_BW_SUBC);
    }
    printf("SSB idx=%u l_first=%u\n", i, l_first);
  }
  write_file(out + "_ssb.cf32", blocks.data(), blocks.size());

  // ---------------------------------------------------------------- SIB1 slot
  srsran_carrier_nr_t carrier    = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci                    = pci;
  carrier.scs                    = scs;
  carrier.nof_prb                = nrb;
  carrier.dl_center_frequency_hz = 0.0;
  carrier.ssb_center_freq_hz     = 0.0;
  cf_t*                buf[SRSRAN_MAX_PORTS] = {};
  buf[0]                                     = srsran_vec_cf_malloc(SRSRAN_SF_LEN_PRB(nrb) * 2);
  srsran_gnb_dl_args_t ga                    = {};
  ga.nof_tx_antennas                         = 1;
  ga.nof_max_prb                             = nrb;
  ga.scs                                     = scs;
  ga.srate_hz = (double)srsran_min_symbol_sz_rb(nrb) * (sscs == 30 ? 30e3 : 15e3);
  ga.pdsch.sch.disable_simd                  = false;
  srsran_gnb_dl_t gnb                        = {};
  if (srsran_gnb_dl_init(&gnb, buf, &ga) < SRSRAN_SUCCESS || srsran_gnb_dl_set_carrier(&gnb, &carrier) < SRSRAN_SUCCESS) {
    fprintf(stderr, "gnb init\n");
    return 1;
  }
  srsran_pdcch_cfg_nr_t pdcch_cfg = {};
  srsran_coreset_t*     cs        = &pdcch_cfg.coreset[0];
  pdcch_cfg.coreset_present[0]    = true;
  cs->id                          = 0;
  cs->duration                    = dur;
  cs->mapping_type                = srsran_coreset_mapping_type_interleaved;
  cs->reg_bundle_size             = srsran_coreset_bundle_size_n6;
  cs->interleaver_size            = srsran_coreset_bundle_size_n2;
  cs->precoder_granularity        = srsran_coreset_precoder_granularity_reg_bundle;
  cs->shift_index                 = pci;
  for (uint32_t i = 0; i < SRSRAN_CORESET_FREQ_DOMAIN_RES_SIZE; i++) {
    cs->freq_resources[i] = i < nrb / 6;
  }
  srsran_search_space_t* ss         = &pdcch_cfg.search_space[0];
  pdcch_cfg.search_space_present[0] = true;
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
  if (srsran_gnb_dl_set_pdcch_config(&gnb, &pdcch_cfg, &dci_cfg) < SRSRAN_SUCCESS) {
    fprintf(stderr, "pdcch cfg\n");
    return 1;
  }

  // SIB1 contents
  uint8_t sib1[512] = {};
  uint32_t sib1_len = 0;
  {
    asn1::rrc_nr::bcch_dl_sch_msg_s msg;
    auto&                           s1 = msg.msg.set_c1().set_sib_type1();
    s1.cell_access_related_info.plmn_id_list.resize(1);
    auto& info = s1.cell_access_related_info.plmn_id_list[0];
    info.plmn_id_list.resize(1);
    info.plmn_id_list[0].mcc_present = true;
    for (int i = 0; i < 3; i++) {
      info.plmn_id_list[0].mcc[i] = (uint8_t)(mcc[i] - '0');
    }
    info.plmn_id_list[0].mnc.resize((uint32_t)mnc.size());
    for (size_t i = 0; i < mnc.size(); i++) {
      info.plmn_id_list[0].mnc[(uint32_t)i] = (uint8_t)(mnc[i] - '0');
    }
    info.tac_present = true;
    info.tac.from_number(tac);
    info.cell_id.from_number(nci & 0xFFFFFFFFFULL);
    info.cell_reserved_for_oper.value = asn1::rrc_nr::plmn_id_info_s::cell_reserved_for_oper_opts::not_reserved;
    asn1::bit_ref bref(sib1, sizeof(sib1));
    if (msg.pack(bref) != asn1::SRSASN_SUCCESS) {
      fprintf(stderr, "sib1 pack\n");
      return 1;
    }
    sib1_len = (uint32_t)bref.distance_bytes();
  }

  srsran_slot_cfg_t slot = {};
  slot.idx               = slotn;
  srsran_gnb_dl_base_zero(&gnb);
  srsran_dci_dl_nr_t dci     = {};
  dci.ctx.location.L         = getenv("NR_TESTGEN_AL") ? (uint32_t)atoi(getenv("NR_TESTGEN_AL")) : 2; // 2 = level 4
  dci.ctx.location.ncce      = 0;
  dci.ctx.ss_type            = srsran_search_space_type_common_0;
  dci.ctx.coreset_id         = 0;
  dci.ctx.coreset_start_rb   = 0;
  dci.ctx.rnti_type          = srsran_rnti_type_si;
  dci.ctx.rnti               = SRSRAN_SIRNTI;
  dci.ctx.format             = srsran_dci_format_nr_1_0;
  dci.coreset0_bw            = nrb;
  dci.time_domain_assigment  = 0;
  dci.mcs                    = 2;
  dci.rv                     = 0;
  dci.sii                    = 0;
  srsran_sch_hl_cfg_nr_t hl  = {};
  hl.typeA_pos               = typeA == 3 ? srsran_dmrs_sch_typeA_pos_3 : srsran_dmrs_sch_typeA_pos_2;
  srsran_re_pattern_t ssb_rvd = {};
  bool                has_rvd = argc >= 20;
  if (has_rvd) {
    srsran_re_pattern_t& p = ssb_rvd;
    p.rb_begin            = (uint32_t)atoi(argv[17]);
    p.rb_end              = SRSRAN_MIN(p.rb_begin + (uint32_t)atoi(argv[18]), nrb);
    p.rb_stride           = 1;
    uint32_t mask         = (uint32_t)strtoul(argv[19], NULL, 16);
    for (uint32_t k = 0; k < SRSRAN_NRE; k++) {
      p.sc[k] = true;
    }
    for (uint32_t l = 0; l < SRSRAN_NSYMB_PER_SLOT_NR; l++) {
      p.symbol[l] = (mask >> l) & 1U;
    }
  }
  srsran_sch_cfg_nr_t pdsch_cfg = {};
  uint32_t            len_rb    = has_rvd ? 4 : nrb / 2;
  bool                ok        = false;
  for (; len_rb <= nrb && !ok; len_rb++) {    // smallest allocation whose transport block holds the SIB1
    // with sync blocks in the slot the SIB1 goes beside them (below their resource blocks), as gNBs do
    if (has_rvd && len_rb > ssb_rvd.rb_begin) {
      break;
    }
    dci.freq_domain_assigment = srsran_ra_nr_type1_riv(nrb, 0, len_rb);
    pdsch_cfg                 = {};
    if (srsran_ra_dl_dci_to_grant_nr(&carrier, &slot, &hl, &dci, &pdsch_cfg, &pdsch_cfg.grant) < SRSRAN_SUCCESS) {
      fprintf(stderr, "grant\n");
      return 1;
    }
    ok = (uint32_t)(pdsch_cfg.grant.tb[0].tbs / 8) >= sib1_len;
  }
  if (!ok) {
    fprintf(stderr, "SIB1 (%u bytes) does not fit\n", sib1_len);
    return 1;
  }
  uint8_t* data = srsran_vec_u8_malloc(SRSRAN_SLOT_MAX_NOF_BITS_NR);
  memset(data, 0, SRSRAN_SLOT_MAX_NOF_BITS_NR);
  memcpy(data, sib1, sib1_len);
  srsran_softbuffer_tx_t sb = {};
  srsran_softbuffer_tx_init_guru(&sb, SRSRAN_SCH_NR_MAX_NOF_CB_LDPC, SRSRAN_LDPC_MAX_LEN_ENCODED_CB);
  pdsch_cfg.grant.tb[0].softbuffer.tx = &sb;
  uint8_t* datas[SRSRAN_MAX_TB]       = {data, NULL};
  if (srsran_gnb_dl_pdcch_put_dl(&gnb, &slot, &dci) < SRSRAN_SUCCESS ||
      srsran_gnb_dl_pdsch_put(&gnb, &slot, &pdsch_cfg, datas) < SRSRAN_SUCCESS) {
    fprintf(stderr, "put\n");
    return 1;
  }
  write_file(out + "_slot.cf32", gnb.sf_symbols[0], (size_t)nrb * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR);
  printf("SIB1 bytes=%u tbs=%u rbs=%u symbols=%u+%u\n",
         sib1_len,
         pdsch_cfg.grant.tb[0].tbs / 8,
         len_rb - 1,
         pdsch_cfg.grant.S,
         pdsch_cfg.grant.L);
  return 0;
}
