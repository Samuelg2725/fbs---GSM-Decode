#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "srsran/asn1/rrc.h"

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
