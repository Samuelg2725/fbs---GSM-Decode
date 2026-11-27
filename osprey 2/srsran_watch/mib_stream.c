/*
 * mib_stream - decode the LTE MIB of every radio frame of one cell (used by lte_mib_stream.py).
 *
 *   mib_stream FILE PCI NPRB PORTS NFRAMES
 *
 * FILE: complex float32 samples at the cell's standard rate (1.92/3.84/7.68/15.36/23.04/30.72 Msps for 6..100
 * resource blocks), sample 0 = the start of a subframe 0. PORTS: the cell's transmit antennas (1, 2 or 4).
 * Every frame (every 10 subframes) is decoded on its own (the PBCH is not combined across frames), so each
 * line is one MIB actually received in that frame:
 *   FRAME k found=1 sfn=N ports=P snr=S hex=XXXXXX       (k = frame index in the file, N = full SFN)
 *   FRAME k found=0 snr=S
 */
#include "srsran/phy/common/phy_common.h"
#include "srsran/phy/ue/ue_mib.h"
#include "srsran/phy/utils/debug.h"
#include "srsran/phy/utils/vector.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
  if (argc < 6) {
    fprintf(stderr, "usage: mib_stream FILE PCI NPRB PORTS NFRAMES\n");
    return 1;
  }
  uint32_t pci     = (uint32_t)atoi(argv[2]);
  uint32_t nprb    = (uint32_t)atoi(argv[3]);
  uint32_t ports   = (uint32_t)atoi(argv[4]);
  uint32_t nframes = (uint32_t)atoi(argv[5]);
  srsran_use_standard_symbol_size(true); // the recordings are at the standard LTE rates
  FILE*    f       = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  uint32_t sf_len = SRSRAN_SF_LEN_PRB(nprb);
  cf_t*    buf    = srsran_vec_cf_malloc(sf_len);
  srsran_ue_mib_t mib;
  if (srsran_ue_mib_init(&mib, buf, nprb) < SRSRAN_SUCCESS) {
    fprintf(stderr, "ue_mib init\n");
    return 1;
  }
  srsran_cell_t cell = {};
  cell.id            = pci;
  cell.nof_prb       = nprb;
  cell.nof_ports     = ports;
  cell.cp            = SRSRAN_CP_NORM;
  if (srsran_ue_mib_set_cell(&mib, cell) < SRSRAN_SUCCESS) {
    fprintf(stderr, "ue_mib cell\n");
    return 1;
  }
  for (uint32_t k = 0; k < nframes; k++) {
    if (fseek(f, (long)((size_t)k * 10 * sf_len * sizeof(cf_t)), SEEK_SET) != 0 ||
        fread(buf, sizeof(cf_t), sf_len, f) != sf_len) {
      break;
    }
    srsran_ue_mib_reset(&mib); // this frame on its own
    uint8_t  bch[SRSRAN_BCH_PAYLOAD_LEN] = {};
    uint32_t nports                       = 0;
    int      sfn_offset                   = 0;
    int      ret                          = srsran_ue_mib_decode(&mib, bch, &nports, &sfn_offset);
    if (ret == SRSRAN_UE_MIB_FOUND) {
      srsran_cell_t c   = cell;
      uint32_t      sfn = 0;
      srsran_pbch_mib_unpack(bch, &c, &sfn);
      uint8_t hex[3] = {};
      for (uint32_t i = 0; i < 24; i++) {
        hex[i / 8] |= (uint8_t)((bch[i] & 1U) << (7 - i % 8));
      }
      printf("FRAME %u found=1 sfn=%u ports=%u snr=%.1f hex=%02x%02x%02x\n",
             k,
             (sfn + (uint32_t)sfn_offset) % 1024,
             nports,
             mib.chest_res.snr_db,
             hex[0],
             hex[1],
             hex[2]);
    } else {
      printf("FRAME %u found=0 snr=%.1f\n", k, mib.chest_res.snr_db);
    }
    fflush(stdout);
  }
  srsran_ue_mib_free(&mib);
  free(buf);
  fclose(f);
  return 0;
}
