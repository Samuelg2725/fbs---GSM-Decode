#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "srsran/asn1/rrc.h"
#include "srsran/srsran.h"

static void save_line(const char* line)
{
  char   ts[32];
  time_t now = time(NULL);
  struct tm when = *localtime(&now);
  strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", &when);

  FILE* out = fopen("/home/revector/sib1_cell.txt", "a");
  if (out) {
    fprintf(out, "%s  %s\n", ts, line);
    fclose(out);
  }
  fprintf(stderr, "\n%s  %s\n", ts, line);
  fflush(stderr);
}

/* The PDSCH payload is packed bytes. nbits is the transport-block size.
 * The result is appended to /home/revector/sib1_cell.txt because the pdsch_ue
 * status display moves the cursor and paints over a line printed on screen. */
extern "C" void srsran_sib1_print(const uint8_t* payload, uint32_t nbits)
{
  static int reported_failure = 0;
  char       line[256];

  if (payload == NULL || nbits < 16) {
    return;
  }

  /* raw bytes of every system-information block that passed its CRC, for the pcap file */
  {
    char     hex[2 * 2048 + 1];
    uint32_t nbytes = (nbits + 7) / 8;
    if (nbytes > 2048) {
      nbytes = 2048;
    }
    for (uint32_t i = 0; i < nbytes; i++) {
      snprintf(hex + 2 * i, 3, "%02x", payload[i]);
    }
    fprintf(stderr, "\nBCCH-DL-SCH hex %s\n", hex);
    fflush(stderr);
  }

  asn1::rrc::bcch_dl_sch_msg_s msg;
  asn1::cbit_ref               bref(payload, (nbits + 7) / 8);
  if (msg.unpack(bref) != asn1::SRSASN_SUCCESS) {
    if (!reported_failure) {
      snprintf(line, sizeof(line), "CRC passed, parse failed, transport block %u bits", nbits);
      save_line(line);
      reported_failure = 1;
    }
    return;
  }
  if (msg.msg.type() != asn1::rrc::bcch_dl_sch_msg_type_c::types::c1) {
    return;
  }
  if (msg.msg.c1().type() != asn1::rrc::bcch_dl_sch_msg_type_c::c1_c_::types::sib_type1) {
    return;
  }

  const asn1::rrc::sib_type1_s& sib1 = msg.msg.c1().sib_type1();
  const uint32_t                cell = (uint32_t)sib1.cell_access_related_info.cell_id.to_number();
  const uint32_t                tac  = (uint32_t)sib1.cell_access_related_info.tac.to_number();
  int                           used = snprintf(line,
                      sizeof(line),
                      "SIB1 cell identity %u   eNB %u   sector %u   TAC %u",
                      cell,
                      cell >> 8,
                      cell & 0xffu,
                      tac);
  const auto&                   list = sib1.cell_access_related_info.plmn_id_list;
  for (uint32_t i = 0; i < list.size() && used > 0 && (uint32_t)used < sizeof(line); i++) {
    const asn1::rrc::plmn_id_s& plmn = list[i].plmn_id;
    unsigned                    mcc  = plmn.mcc[0] * 100u + plmn.mcc[1] * 10u + plmn.mcc[2];
    unsigned                    mnc  = 0;
    for (uint32_t d = 0; d < plmn.mnc.size(); d++) {
      mnc = mnc * 10u + plmn.mnc[d];
    }
    int wrote;
    if (plmn.mnc.size() == 2) {
      wrote = snprintf(line + used, sizeof(line) - (uint32_t)used, "   PLMN %03u-%02u", mcc, mnc);
    } else {
      wrote = snprintf(line + used, sizeof(line) - (uint32_t)used, "   PLMN %03u-%03u", mcc, mnc);
    }
    if (wrote < 0) {
      break;
    }
    used += wrote;
  }
  save_line(line);
}


/* The decoded MIB: 24 bits, one bit per byte as srsRAN keeps it. Printed packed (3 bytes) for the pcap file. */
extern "C" void srsran_mib_print(const uint8_t* bits)
{
  uint8_t packed[3] = {0, 0, 0};
  if (bits == NULL) {
    return;
  }
  for (int i = 0; i < 24; i++) {
    if (bits[i]) {
      packed[i / 8] |= (uint8_t)(0x80u >> (i % 8));
    }
  }
  fprintf(stderr, "\nBCCH-BCH hex %02x%02x%02x\n", packed[0], packed[1], packed[2]);
  fflush(stderr);
}


/* ------------------------------------------------------------------------------------------------------------
 * SIB1 decoding with combining of repeated transmissions.
 *
 * A cell sends the same SIB1 every 20 ms (subframe 5 of even frames) in four redundancy versions per 80 ms.
 * Stock pdsch_ue decodes every transmission on its own: it clears the soft buffer each time, and gives up on
 * a subframe whose control message (PDCCH, format 1A for SI-RNTI) does not decode. For a weak cell this needs
 * several dB more signal than necessary. Here:
 *   1. PDCCH: the soft bits of the control region are added up over the SIB1 occasions with the same
 *      redundancy version (every 80 ms, so the control message is bit-identical), one sum per control-region
 *      size; when a single subframe gives no control message, the search is repeated on the sums.
 *   2. PDSCH: the soft buffer is kept between transmissions of the same allocation, so the redundancy versions
 *      combine the way LTE intends (as HARQ does). At most PDSCH_MAX_COMBINE transmissions are added up
 *      Soft bits are divided by 4 first: srsRAN's 16-bit turbo decoder loses accuracy with large inputs, and
 *      8 transmissions measured best (16 was worse).
 *   3. When a subframe's control message is missed, the last decoded SIB1 allocation is used (it does not
 *      change between transmissions) with the redundancy version for this frame.
 * The CRC still decides: a block is only reported when it passes. Set SIB1_NO_COMBINE=1 to switch this off.
 * ---------------------------------------------------------------------------------------------------------- */
#define PDCCH_MAX_COMBINE 16
#define PDSCH_MAX_COMBINE 8

static float*               pdcch_acc[4][3] = {};      /* [redundancy version slot k][CFI - 1] */
static uint32_t             pdcch_n[4][3]   = {};
static srsran_pdsch_grant_t sib1_grant;
static uint32_t             sib1_grant_cfi = 0;
static bool                 have_grant     = false;
static uint32_t             pdsch_n        = PDSCH_MAX_COMBINE;
static uint32_t             reported       = 0;
static srsran_dci_dl_t      sib1_dci;                /* last SIB1 control message found */
static bool                 have_dci       = false;
static bool                 dci_printed    = false;

extern "C" int srsran_sib1_combining_enabled(void)
{
  const char* off = getenv("SIB1_NO_COMBINE");
  return !(off && *off && *off != '0');
}

static bool same_allocation(const srsran_pdsch_grant_t* a, const srsran_pdsch_grant_t* b)
{
  return a->tb[0].tbs == b->tb[0].tbs && a->nof_prb == b->nof_prb && a->nof_re == b->nof_re &&
         memcmp(a->prb_idx, b->prb_idx, sizeof(a->prb_idx)) == 0;
}

extern "C" int srsran_sib1_decode(srsran_ue_dl_t*     q,
                                  srsran_dl_sf_cfg_t* sf,
                                  srsran_ue_dl_cfg_t* cfg,
                                  srsran_pdsch_cfg_t* pdsch_cfg,
                                  uint8_t*            data[SRSRAN_MAX_CODEWORDS],
                                  bool                acks[SRSRAN_MAX_CODEWORDS])
{
  srsran_dci_dl_t dci[SRSRAN_MAX_DCI_MSG];
  memset(dci, 0, sizeof(dci));

  srsran_ue_dl_set_mi_auto(q);
  if (srsran_ue_dl_decode_fft_estimate(q, sf, cfg) < 0) {
    return SRSRAN_ERROR;
  }
  int n = srsran_ue_dl_find_dl_dci(q, sf, cfg, SRSRAN_SIRNTI, dci);

  /* 1. control message: add up the control-region soft bits over occasions with the same redundancy version.
   *    One sum per control-region size (CFI): at low SNR the CFI itself is sometimes misread, and a subframe read
   *    with the wrong CFI must not spoil the sum for the right one. The sum is scaled by 1/sqrt(n): the noise
   *    stays at the level of one subframe while the signal grows (srsRAN skips candidates whose soft bits are
   *    small on average, so a plain average would never be tried). */
  uint32_t k        = (sf->tti / 10 / 2) % 4;
  bool     from_avg = false;
  uint32_t cfi      = sf->cfi;
  if (cfi >= 1 && cfi <= 3) {
    uint32_t nbits = 72 * q->pdcch.nof_cce[cfi - 1];
    float**  acc   = &pdcch_acc[k][cfi - 1];
    if (*acc == NULL) {
      *acc = (float*)calloc(q->pdcch.max_bits, sizeof(float));
    }
    if (*acc != NULL && nbits <= q->pdcch.max_bits) {
      if (pdcch_n[k][cfi - 1] >= PDCCH_MAX_COMBINE) {
        memset(*acc, 0, q->pdcch.max_bits * sizeof(float));
        pdcch_n[k][cfi - 1] = 0;
      }
      for (uint32_t i = 0; i < nbits; i++) {
        (*acc)[i] += q->pdcch.llr[i];
      }
      pdcch_n[k][cfi - 1]++;
    }
  }
  if (n <= 0) {
    /* search the sum for the CFI just read first, then the other two */
    uint32_t order[3] = {cfi, cfi == 1 ? 2u : 1u, cfi == 3 ? 2u : 3u};
    for (int o = 0; o < 3 && n <= 0; o++) {
      uint32_t c = order[o];
      if (c < 1 || c > 3 || pdcch_acc[k][c - 1] == NULL || pdcch_n[k][c - 1] < 2) {
        continue;
      }
      uint32_t nbits = 72 * q->pdcch.nof_cce[c - 1];
      float    scale = 1.0f / sqrtf((float)pdcch_n[k][c - 1]);
      for (uint32_t i = 0; i < nbits && i < q->pdcch.max_bits; i++) {
        q->pdcch.llr[i] = pdcch_acc[k][c - 1][i] * scale;
      }
      sf->cfi = c;
      n       = srsran_ue_dl_find_dl_dci(q, sf, cfg, SRSRAN_SIRNTI, dci);
      if (n > 0) {
        from_avg = true;
        cfi      = c;
      } else {
        sf->cfi = cfi;
      }
    }
  }

  /* grant: from this subframe's control message, else the last decoded SIB1 allocation */
  bool from_saved = false;
  if (n > 0) {
    if (srsran_ue_dl_dci_to_pdsch_grant(q, sf, cfg, &dci[0], &pdsch_cfg->grant)) {
      return SRSRAN_ERROR;
    }
    sib1_dci = dci[0];
    have_dci = true;
  } else if (have_grant) {
    pdsch_cfg->grant = sib1_grant;
    sf->cfi          = sib1_grant_cfi; /* SIB1's control region size does not change */
    cfi              = sib1_grant_cfi;
    from_saved       = true;
  } else {
    if (getenv("SIB1_DEBUG")) {
      fprintf(stderr, "\nSIB1dbg tti=%u k=%u cfi=%u dci=none\n", sf->tti, k, cfi);
    }
    return 0;
  }
  for (int i = 0; i < SRSRAN_MAX_CODEWORDS; i++) {
    if (pdsch_cfg->grant.tb[i].enabled && (from_saved || pdsch_cfg->grant.tb[i].rv < 0)) {
      pdsch_cfg->grant.tb[i].rv = ((uint32_t)ceilf(1.5f * (float)k)) % 4;   /* 36.321 5.3.1 for SI */
    }
  }

  /* 2. data: keep the soft buffer while the allocation is the same */
  uint32_t max_combine = PDSCH_MAX_COMBINE;
  if (getenv("SIB1_MAX_COMBINE")) {
    max_combine = (uint32_t)atoi(getenv("SIB1_MAX_COMBINE"));
  }
  if (!have_grant || !same_allocation(&pdsch_cfg->grant, &sib1_grant) || pdsch_n >= max_combine) {
    for (int i = 0; i < SRSRAN_MAX_CODEWORDS; i++) {
      if (pdsch_cfg->grant.tb[i].enabled) {
        srsran_softbuffer_rx_reset_tbs(pdsch_cfg->softbuffers.rx[i], (uint32_t)pdsch_cfg->grant.tb[i].tbs);
      }
    }
    pdsch_n = 0;
  }
  pdsch_n++;
  if (!from_saved) {
    sib1_grant     = pdsch_cfg->grant;
    sib1_grant_cfi = cfi;
    have_grant     = true;
  }

  srsran_pdsch_res_t res[SRSRAN_MAX_CODEWORDS];
  memset(res, 0, sizeof(res));
  for (int i = 0; i < SRSRAN_MAX_CODEWORDS; i++) {
    if (pdsch_cfg->grant.tb[i].enabled) {
      res[i].payload = data[i];
      res[i].crc     = false;
    }
  }
  /* Soft bits are divided by 4 before they are added up: srsRAN's 16-bit turbo decoder loses accuracy when its
   * input grows, which otherwise caps the gain from combining (measured: 33 SIB1 decodes instead of 17 at the
   * same noise level). SIB1_LLR_SCALE overrides the factor (1 = srsRAN's own scale). */
  float llr_scale = getenv("SIB1_LLR_SCALE") ? (float)atof(getenv("SIB1_LLR_SCALE")) : 0.25f;
  if (llr_scale > 0.0f && llr_scale != 1.0f) {
    /* equalised symbols (and so soft bits) shrink by this factor when the channel estimate grows by it */
    for (uint32_t p = 0; p < q->cell.nof_ports; p++) {
      for (uint32_t a = 0; a < q->nof_rx_antennas; a++) {
        srsran_vec_sc_prod_cfc(q->chest_res.ce[p][a], 1.0f / llr_scale, q->chest_res.ce[p][a], q->chest_res.nof_re);
      }
    }
    q->chest_res.noise_estimate /= llr_scale * llr_scale;
  }
  if (srsran_ue_dl_decode_pdsch(q, sf, pdsch_cfg, res)) {
    return SRSRAN_ERROR;
  }
  for (int i = 0; i < SRSRAN_MAX_CODEWORDS; i++) {
    if (pdsch_cfg->grant.tb[i].enabled) {
      acks[i] = res[i].crc;
    }
  }
  if (getenv("SIB1_DEBUG")) {
    fprintf(stderr,
            "\nSIB1dbg tti=%u k=%u cfi=%u dci=%s combined=%u crc=%d rv=%d tbs=%u prb=%u\n",
            sf->tti,
            k,
            cfi,
            from_saved ? "saved" : (from_avg ? "averaged" : "single"),
            pdsch_n,
            (int)acks[0],
            pdsch_cfg->grant.tb[0].rv,
            pdsch_cfg->grant.tb[0].tbs,
            pdsch_cfg->grant.nof_prb);
  }
  if (acks[0] && have_dci && !dci_printed) {
    /* what is needed to rebuild this cell's SIB1 subframe (for interference cancellation) */
    fprintf(stderr,
            "\nSIB1dci cfi=%u L=%u ncce=%u riv=%u mode=%d nprb1a=%d ngap=%d mcs=%u\n",
            cfi,
            sib1_dci.location.L,
            sib1_dci.location.ncce,
            sib1_dci.type2_alloc.riv,
            (int)sib1_dci.type2_alloc.mode,
            (int)sib1_dci.type2_alloc.n_prb1a,
            (int)sib1_dci.type2_alloc.n_gap,
            sib1_dci.tb[0].mcs_idx);
    fflush(stderr);
    dci_printed = true;
  }
  if (acks[0]) {
    if ((pdsch_n > 1 || from_avg || from_saved) && reported < 3) {
      fprintf(stderr,
              "\nSIB1 combining: decoded from %u transmission(s)%s%s\n",
              pdsch_n,
              from_avg ? ", control message from averaged repeats" : "",
              from_saved ? ", allocation from an earlier control message" : "");
      fflush(stderr);
      reported++;
    }
    pdsch_n = PDSCH_MAX_COMBINE; /* start afresh after a success */
  }
  return 1;
}
