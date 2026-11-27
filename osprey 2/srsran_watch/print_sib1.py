#!/usr/bin/env python3
"""Make stock pdsch_ue print the 28-bit cell identity from a passing SIB1 block, and the raw bytes of the
MIB and of every system-information block (used for the pcap files).

This file is not a decoder. It copies sib1_print.cc into the srsRAN examples
directory and adds one call in lib/examples/pdsch_ue.c. srsRAN still does the
receive, the CRC check, and the ASN.1 parse. A block is printed only after
that CRC has passed.
"""

from pathlib import Path
import shutil
import sys

EXAMPLES = Path(sys.argv[1] if len(sys.argv) > 1 else "~/srsRAN_4G/lib/examples").expanduser()
SOURCE = Path(__file__).resolve().parent / "sib1_print.cc"

CALL = """nof_detected++;
        if (prog_args.rnti == SRSRAN_SIRNTI) {
          for (uint32_t sib_tb = 0; sib_tb < SRSRAN_MAX_CODEWORDS; sib_tb++) {
            if (acks[sib_tb]) {
              srsran_sib1_print(data[sib_tb], pdsch_cfg.grant.tb[sib_tb].tbs);
            }
          }
        }"""

PROTO = "void srsran_sib1_print(const uint8_t* payload, uint32_t nbits);\n"
MIB_PROTO = "void srsran_mib_print(const uint8_t* bits);\n"
MIB_ANCHOR = "srsran_pbch_mib_unpack(bch_payload, &cell, &sfn);"
MIB_CALL = MIB_ANCHOR + "\n              srsran_mib_print(bch_payload);"
# In file mode pdsch_ue sets up its decoder once with PHICH resource "1" and never applies the PHICH settings
# the MIB reports, so a cell configured with another PHICH resource (e.g. 1/2) decodes its MIB but never SIB1.
PHICH_ANCHOR = "srsran_mib_print(bch_payload);"
# SIB1 decoding with combining of repeated transmissions (srsran_sib1_decode in sib1_print.cc)
COMBINE_PROTO = ("int srsran_sib1_combining_enabled(void);\n"
                 "int srsran_sib1_decode(srsran_ue_dl_t* q, srsran_dl_sf_cfg_t* sf, srsran_ue_dl_cfg_t* cfg,\n"
                 "                       srsran_pdsch_cfg_t* pdsch_cfg, uint8_t* data[SRSRAN_MAX_CODEWORDS],\n"
                 "                       bool acks[SRSRAN_MAX_CODEWORDS]);\n")
COMBINE_ANCHOR = "n = srsran_ue_dl_find_and_decode(&ue_dl, &dl_sf, &ue_dl_cfg, &pdsch_cfg, data, acks);"
COMBINE_CALL = """if (prog_args.rnti == SRSRAN_SIRNTI && srsran_sib1_combining_enabled()) {
                  /* SIB1 uses one port or transmit diversity: decode once per subframe, combining repeats */
                  n = (tm == (cell.nof_ports == 1 ? 0u : 1u))
                          ? srsran_sib1_decode(&ue_dl, &dl_sf, &ue_dl_cfg, &pdsch_cfg, data, acks)
                          : 0;
                } else {
                  n = srsran_ue_dl_find_and_decode(&ue_dl, &dl_sf, &ue_dl_cfg, &pdsch_cfg, data, acks);
                }"""
PHICH_CALL = PHICH_ANCHOR + """
              {
                /* file mode: use the PHICH settings the MIB reports (srsRAN starts from a fixed guess) */
                srsran_cell_t mib_cell   = ue_dl.cell;
                mib_cell.phich_resources = cell.phich_resources;
                mib_cell.phich_length    = cell.phich_length;
                srsran_ue_dl_set_cell(&ue_dl, mib_cell);
              }"""


def main():
    if not SOURCE.is_file():
        sys.exit("sib1_print.cc was not found next to this script")
    if not EXAMPLES.is_dir():
        sys.exit("srsRAN examples were not found at %s" % EXAMPLES)
    shutil.copy(SOURCE, EXAMPLES / "sib1_print.cc")
    shutil.copy(SOURCE.parent / "sib1_regen.c", EXAMPLES / "sib1_regen.c")
    shutil.copy(SOURCE.parent / "nr_decode.cc", EXAMPLES / "nr_decode.cc")
    shutil.copy(SOURCE.parent / "mib_stream.c", EXAMPLES / "mib_stream.c")
    if (SOURCE.parent / "nr_testgen.cc").is_file():
        shutil.copy(SOURCE.parent / "nr_testgen.cc", EXAMPLES / "nr_testgen.cc")

    pdsch = EXAMPLES / "pdsch_ue.c"
    text = pdsch.read_text()
    if "srsran_sib1_print" not in text:
        anchor = "nof_detected++;"
        if anchor not in text:
            sys.exit("pdsch_ue.c has no decode-success point to hook")
        text = text.replace(anchor, CALL, 1)
        main_at = text.find("int main(")
        if main_at < 0:
            sys.exit("pdsch_ue.c has no main")
        text = text[:main_at] + PROTO + "\n" + text[main_at:]
        pdsch.write_text(text)
        print("patched %s" % pdsch)
    else:
        print("pdsch_ue.c already prints SIB1")
    # second hook: the raw MIB, for the pcap file
    if "srsran_mib_print" not in text:
        if MIB_ANCHOR not in text:
            sys.exit("pdsch_ue.c has no MIB decode point to hook")
        text = text.replace(MIB_ANCHOR, MIB_CALL, 1)
        main_at = text.find("int main(")
        text = text[:main_at] + MIB_PROTO + "\n" + text[main_at:]
        pdsch.write_text(text)
        print("patched %s (MIB bytes)" % pdsch)
    else:
        print("pdsch_ue.c already prints the MIB bytes")
    # fourth hook: SIB1 decoding with combining of repeats
    if "srsran_sib1_decode" not in text:
        if COMBINE_ANCHOR not in text:
            sys.exit("pdsch_ue.c has no PDSCH decode call to hook")
        text = text.replace(COMBINE_ANCHOR, COMBINE_CALL, 1)
        main_at = text.find("int main(")
        text = text[:main_at] + COMBINE_PROTO + "\n" + text[main_at:]
        pdsch.write_text(text)
        print("patched %s (SIB1 combining)" % pdsch)
    else:
        print("pdsch_ue.c already combines SIB1 repeats")
    # third hook: apply the MIB's PHICH settings
    if "use the PHICH settings the MIB reports" not in text:
        text = text.replace(PHICH_ANCHOR, PHICH_CALL, 1)
        pdsch.write_text(text)
        print("patched %s (PHICH settings from the MIB)" % pdsch)
    else:
        print("pdsch_ue.c already uses the MIB's PHICH settings")

    cmake = EXAMPLES / "CMakeLists.txt"
    cmake_text = cmake.read_text()
    changed = False
    old_exe = "add_executable(pdsch_ue pdsch_ue.c)"
    new_exe = "add_executable(pdsch_ue pdsch_ue.c sib1_print.cc)"
    if old_exe in cmake_text:
        cmake_text = cmake_text.replace(old_exe, new_exe)
        changed = True
    # A previous edit left "rrc_asn1" after the closing parenthesis, which
    # CMake rejects. Put it back inside the link command.
    if ") rrc_asn1" in cmake_text:
        cmake_text = cmake_text.replace(") rrc_asn1", " rrc_asn1)")
        changed = True
    for line in (
        "target_link_libraries(pdsch_ue srsran_phy srsran_common srsran_rf pthread)",
        "target_link_libraries(pdsch_ue srsran_common srsran_phy pthread)",
    ):
        if line in cmake_text:
            cmake_text = cmake_text.replace(line, line[:-1] + " rrc_asn1)")
            changed = True
    if "sib1_regen" not in cmake_text:
        cmake_text += ("\n# SIB1 subframe rebuild for interference cancellation (lte_cancel.py)\n"
                       "add_executable(sib1_regen sib1_regen.c)\n"
                       "target_link_libraries(sib1_regen srsran_phy srsran_common pthread)\n")
        changed = True
    if "nr_decode" not in cmake_text:
        cmake_text += ("\n# 5G NR MIB / SIB1 decoding for nr_decode.py\n"
                       "add_executable(nr_decode nr_decode.cc)\n"
                       "target_link_libraries(nr_decode srsran_phy srsran_common rrc_nr_asn1 asn1_utils pthread)\n")
        changed = True
    if "mib_stream" not in cmake_text:
        cmake_text += ("\n# LTE MIB of every frame of one cell for lte_mib_stream.py\n"
                       "add_executable(mib_stream mib_stream.c)\n"
                       "target_link_libraries(mib_stream srsran_phy srsran_common pthread)\n")
        changed = True
    if "nr_testgen" not in cmake_text and (EXAMPLES / "nr_testgen.cc").is_file():
        cmake_text += ("\n# 5G NR test signals (SS/PBCH block + SIB1) for checking nr_decode.py\n"
                       "add_executable(nr_testgen nr_testgen.cc)\n"
                       "target_link_libraries(nr_testgen srsran_phy srsran_common rrc_nr_asn1 asn1_utils pthread)\n")
        changed = True
    if changed:
        cmake.write_text(cmake_text)
        print("patched %s" % cmake)
    else:
        print("CMakeLists.txt already builds the SIB1 print")


if __name__ == "__main__":
    main()
