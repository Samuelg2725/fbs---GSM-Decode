/*
 * sib1_regen: rebuild a decoded cell's SIB1 subframe (subframe 5: PCFICH, the SIB1 control message on PDCCH and
 * SIB1 on PDSCH) with srsRAN's own transmitter code, for interference cancellation (lte_cancel.py).
 *
 * The control message and allocation come from the decoder ("SIB1dci" line printed by pdsch_ue with the patch),
 * the payload is the decoded SIB1 transport block. SIB1 is sent in four redundancy versions per 80 ms
 * (36.321 5.3.1: RV = ceil(3k/2) mod 4 for k = (SFN/2) mod 4), so four subframes are written, k = 0..3.
 *
 * Output (-o): for k = 0..3, for each antenna port: the subframe's resource grid as complex float32, symbol-major
 * (14 OFDM symbols x 12*nof_prb subcarriers, lowest subcarrier first). Only PCFICH, PDCCH and PDSCH are in it.
 *
 *   sib1_regen -c PCI -p PRB -P ports -g phich_res(0..3: 1/6,1/2,1,2) -d phich_extended(0/1) -f CFI
 *              -L agg_level_log2 -N ncce -r riv -M mode -a nprb1a -G ngap -m mcs -x payload_hex -o out
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "srsran/srsran.h"

int main(int argc, char** argv)
{
  srsran_cell_t   cell = {};
  srsran_dci_dl_t dci  = {};
  uint32_t        cfi = 2, L = 2, ncce = 0, phich_res = 2, phich_ext = 0;
  const char*     hex = NULL;
  const char*     out = NULL;
  cell.nof_prb        = 50;
  cell.nof_ports      = 2;
  cell.cp             = SRSRAN_CP_NORM;
  int opt;
  while ((opt = getopt(argc, argv, "c:p:P:g:d:f:L:N:r:M:a:G:m:x:o:")) != -1) {
    switch (opt) {
      case 'c': cell.id = atoi(optarg); break;
      case 'p': cell.nof_prb = atoi(optarg); break;
      case 'P': cell.nof_ports = atoi(optarg); break;
      case 'g': phich_res = atoi(optarg); break;
      case 'd': phich_ext = atoi(optarg); break;
      case 'f': cfi = atoi(optarg); break;
      case 'L': L = atoi(optarg); break;
      case 'N': ncce = atoi(optarg); break;
      case 'r': dci.type2_alloc.riv = atoi(optarg); break;
      case 'M': dci.type2_alloc.mode = atoi(optarg); break;
      case 'a': dci.type2_alloc.n_prb1a = atoi(optarg); break;
      case 'G': dci.type2_alloc.n_gap = atoi(optarg); break;
      case 'm': dci.tb[0].mcs_idx = atoi(optarg); break;
      case 'x': hex = optarg; break;
      case 'o': out = optarg; break;
      default: fprintf(stderr, "bad option\n"); return 1;
    }
  }
  if (!hex || !out) {
    fprintf(stderr, "need -x payload_hex and -o output\n");
    return 1;
  }
  const srsran_phich_r_t res_map[4] = {SRSRAN_PHICH_R_1_6, SRSRAN_PHICH_R_1_2, SRSRAN_PHICH_R_1, SRSRAN_PHICH_R_2};
  cell.phich_resources = res_map[phich_res & 3];
  cell.phich_length    = phich_ext ? SRSRAN_PHICH_EXT : SRSRAN_PHICH_NORM;

  uint8_t  payload[2048] = {};
  uint32_t nbytes        = strlen(hex) / 2;
  for (uint32_t i = 0; i < nbytes && i < sizeof(payload); i++) {
    unsigned v;
    sscanf(hex + 2 * i, "%2x", &v);
    payload[i] = (uint8_t)v;
  }

  srsran_regs_t   regs;
  srsran_pcfich_t pcfich;
  srsran_pdcch_t  pdcch;
  srsran_pdsch_t  pdsch;
  if (srsran_regs_init(&regs, cell) || srsran_pcfich_init(&pcfich, 1) || srsran_pcfich_set_cell(&pcfich, &regs, cell) ||
      srsran_pdcch_init_enb(&pdcch, cell.nof_prb) || srsran_pdcch_set_cell(&pdcch, &regs, cell) ||
      srsran_pdsch_init_enb(&pdsch, cell.nof_prb) || srsran_pdsch_set_cell(&pdsch, cell)) {
    fprintf(stderr, "init failed\n");
    return 1;
  }
  srsran_softbuffer_tx_t softbuffer;
  srsran_softbuffer_tx_init(&softbuffer, cell.nof_prb);

  uint32_t nof_re = SRSRAN_SF_LEN_RE(cell.nof_prb, cell.cp);
  cf_t*    sf_symbols[SRSRAN_MAX_PORTS];
  for (uint32_t p = 0; p < SRSRAN_MAX_PORTS; p++) {
    sf_symbols[p] = srsran_vec_cf_malloc(nof_re);
  }
  uint8_t* data[SRSRAN_MAX_CODEWORDS] = {payload, NULL};

  FILE* f = fopen(out, "wb");
  if (!f) {
    perror(out);
    return 1;
  }
  srsran_dl_sf_cfg_t dl_sf = {};
  dl_sf.tti                = 5; /* subframe 5: scrambling of PCFICH, PDCCH and PDSCH depends on the subframe */
  dl_sf.cfi                = cfi;
  dl_sf.sf_type            = SRSRAN_SF_NORM;

  dci.rnti        = SRSRAN_SIRNTI;
  dci.format      = SRSRAN_DCI_FORMAT1A;
  dci.alloc_type  = SRSRAN_RA_ALLOC_TYPE2;
  dci.tb[0].ndi   = 0;
  dci.tb[1].mcs_idx = 0;
  srsran_tm_t tm  = cell.nof_ports == 1 ? SRSRAN_TM1 : SRSRAN_TM2;

  for (uint32_t k = 0; k < 4; k++) {
    for (uint32_t p = 0; p < SRSRAN_MAX_PORTS; p++) {
      srsran_vec_cf_zero(sf_symbols[p], nof_re);
    }
    dci.tb[0].rv = ((uint32_t)ceilf(1.5f * (float)k)) % 4;
    srsran_pcfich_encode(&pcfich, &dl_sf, sf_symbols);

    srsran_dci_msg_t dci_msg = {};
    srsran_dci_msg_pack_pdsch(&cell, &dl_sf, NULL, &dci, &dci_msg);
    dci_msg.location.L    = L;
    dci_msg.location.ncce = ncce;
    if (srsran_pdcch_encode(&pdcch, &dl_sf, &dci_msg, sf_symbols)) {
      fprintf(stderr, "PDCCH encode failed\n");
      return 1;
    }
    srsran_pdsch_cfg_t pdsch_cfg = {};
    pdsch_cfg.rnti               = SRSRAN_SIRNTI;
    pdsch_cfg.softbuffers.tx[0]  = &softbuffer;
    if (srsran_ra_dl_dci_to_grant(&cell, &dl_sf, tm, false, &dci, &pdsch_cfg.grant)) {
      fprintf(stderr, "grant failed\n");
      return 1;
    }
    if (k == 0) {
      srsran_softbuffer_tx_reset(&softbuffer);
    }
    if (srsran_pdsch_encode(&pdsch, &dl_sf, &pdsch_cfg, data, sf_symbols)) {
      fprintf(stderr, "PDSCH encode failed\n");
      return 1;
    }
    for (uint32_t p = 0; p < cell.nof_ports; p++) {
      fwrite(sf_symbols[p], sizeof(cf_t), nof_re, f);
    }
  }
  fclose(f);
  printf("ok\n");
  return 0;
}
