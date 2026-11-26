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
    if changed:
        cmake.write_text(cmake_text)
        print("patched %s" % cmake)
    else:
        print("CMakeLists.txt already builds the SIB1 print")


if __name__ == "__main__":
    main()
