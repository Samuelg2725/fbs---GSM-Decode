#!/usr/bin/env python3
"""Find a USB cellular modem (e.g. Teltonika CALYX) on its serial ports and ask it what it is.

  python3 modem_probe.py              # try every /dev/ttyUSB* and /dev/ttyACM* port
  python3 modem_probe.py /dev/ttyUSB2 # one port
  python3 modem_probe.py --cmd 'AT+COPS?' /dev/ttyUSB2

Sends standard AT commands (ATI, AT+CGMI, AT+CGMM, AT+CGMR) and prints the answers, so we know which cellular chip
is inside and which cell-scan commands it supports. No extra packages needed. Run with sudo if the ports are not
readable (or add your user to the 'dialout' group).
"""
import argparse
import glob
import os
import select
import sys
import termios
import time

IDENTIFY = ["ATE0", "ATI", "AT+CGMI", "AT+CGMM", "AT+CGMR", "AT+CPIN?", "AT+COPS?", "AT+CSQ"]


def open_port(path):
    fd = os.open(path, os.O_RDWR | os.O_NOCTTY | os.O_NONBLOCK)
    attrs = termios.tcgetattr(fd)
    attrs[0] = 0                                     # raw input
    attrs[1] = 0                                     # raw output
    attrs[2] = termios.CS8 | termios.CREAD | termios.CLOCAL
    attrs[3] = 0                                     # no echo / canonical mode
    attrs[4] = attrs[5] = termios.B115200
    termios.tcsetattr(fd, termios.TCSANOW, attrs)
    termios.tcflush(fd, termios.TCIOFLUSH)
    return fd


def command(fd, cmd, timeout=3.0):
    """Send one AT command, return the reply text (up to OK / ERROR or the timeout)."""
    os.write(fd, (cmd + "\r").encode())
    out, end = b"", time.time() + timeout
    while time.time() < end:
        r, _, _ = select.select([fd], [], [], 0.2)
        if r:
            try:
                chunk = os.read(fd, 4096)
            except BlockingIOError:
                continue
            out += chunk
            text = out.decode(errors="replace")
            if "\r\nOK\r\n" in text or "ERROR" in text:
                break
    return out.decode(errors="replace").strip()


def main():
    ap = argparse.ArgumentParser(description="Identify a USB cellular modem through AT commands.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("ports", nargs="*", help="serial ports (default: all /dev/ttyUSB* and /dev/ttyACM*)")
    ap.add_argument("--cmd", action="append", help="send this AT command instead of the identification set")
    ap.add_argument("--timeout", type=float, default=3.0, help="seconds to wait per command (default 3)")
    args = ap.parse_args()
    ports = args.ports or sorted(glob.glob("/dev/ttyUSB*") + glob.glob("/dev/ttyACM*"))
    if not ports:
        sys.exit("No /dev/ttyUSB* or /dev/ttyACM* ports. Is the modem plugged in and powered? Check: lsusb")
    for port in ports:
        print("== %s" % port)
        try:
            fd = open_port(port)
        except OSError as e:
            print("   cannot open: %s (try sudo)" % e)
            continue
        try:
            first = command(fd, "AT", 1.5)
            if "OK" not in first:
                print("   no answer to AT (not the AT port)")
                continue
            for cmd in (args.cmd or IDENTIFY):
                reply = command(fd, cmd, args.timeout)
                print("   %-10s -> %s" % (cmd, " | ".join(l for l in reply.splitlines() if l.strip()) or "(nothing)"))
        finally:
            os.close(fd)
    return 0


if __name__ == "__main__":
    sys.exit(main())
