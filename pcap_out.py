#!/usr/bin/env python3
"""Write decoded GSM and LTE messages to a .pcap file that Wireshark opens and decodes.

Every message is stored the standard way for radio captures: a GSMTAP header (the format gr-gsm, srsRAN and
osmocom use) inside UDP to port 4729, inside IPv4/Ethernet. Wireshark recognises port 4729 as GSMTAP and
decodes the message underneath (GSM RR System Information, LTE RRC MIB / SIB1 / SIBs).

  GSM: the messages gr-gsm decodes already carry a GSMTAP header - they are written as they are.
  LTE: the MIB (BCCH-BCH) and system information blocks (BCCH-DL-SCH) srsRAN decodes are given a GSMTAP
       header of type LTE RRC. Reading an LTE packet in Wireshark:
         time          when the message was on air (the recording's time + the cell's frame position in it;
                       each new MIB 40 ms apart, SIB1 in subframe 5) - not when it was decoded
         ARFCN         the EARFCN
         frame number  the LTE system frame number (SFN, 0-1023) from the cell's MIB
         source port   10000 + PCI (PCI 21 -> 10021, PCI 430 -> 10430), so the Info column shows the cell

Open a file with:  wireshark ~/field/pcap/<file>.pcap   (or tshark -r <file> -V for text)
"""
import os
import struct
import threading
import time

GSMTAP_PORT = 4729
GSMTAP_TYPE_LTE_RRC = 0x0d
LTE_RRC_BCCH_BCH = 4           # MIB
LTE_RRC_BCCH_DL_SCH = 5        # SIB1 and the other system information blocks


def default_path(tool):
    folder = os.path.expanduser("~/field/pcap")
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, "%s_%s.pcapng" % (tool, time.strftime("%Y%m%d_%H%M%S")))


def utc_text(ts):
    """'2026-11-27 13:15:51.123456 UTC' for a packet comment."""
    return time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(ts)) + (".%06d UTC" % int(round((ts % 1) * 1e6) % 1000000))


class PcapWriter(object):
    """Writes packets to a pcapng file (Ethernet link type), in time order, each with a comment that spells out
    when it was on air and what it is (shown by Wireshark under "Packet comments", and as a column with
    right-click > Apply as Column). Safe to call from several threads.

    Packets are kept until flush() (or close()), then written sorted by their on-air time: decoding finishes
    cell by cell, so without this a later-decoded cell's earlier messages would appear after later ones."""

    def __init__(self, path, buffered=True):
        self.path = os.path.expanduser(path)
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        self.fh = open(self.path, "wb")
        # Section Header Block (byte order magic, version 1.0, section length unknown)
        shb = struct.pack("<IHHq", 0x1A2B3C4D, 1, 0, -1)
        self.fh.write(struct.pack("<II", 0x0A0D0D0A, 12 + len(shb)) + shb + struct.pack("<I", 12 + len(shb)))
        # Interface Description Block: Ethernet, timestamps in microseconds (the default)
        idb = struct.pack("<HHI", 1, 0, 65535)
        self.fh.write(struct.pack("<II", 1, 12 + len(idb)) + idb + struct.pack("<I", 12 + len(idb)))
        self.fh.flush()
        self.count = 0
        self.lock = threading.Lock()
        self.buffered = buffered
        self.pending = []

    def _packet(self, gsmtap, ts=None, sport=GSMTAP_PORT, comment=None):
        ts = time.time() if ts is None else ts
        udp = struct.pack(">HHHH", sport, GSMTAP_PORT, 8 + len(gsmtap), 0) + gsmtap
        ip = struct.pack(">BBHHHBBH4s4s", 0x45, 0, 20 + len(udp), 0, 0x4000, 64, 17, 0,
                         bytes([127, 0, 0, 1]), bytes([127, 0, 0, 1]))
        ip = ip[:10] + struct.pack(">H", _ip_checksum(ip)) + ip[12:]
        frame = b"\x00" * 12 + b"\x08\x00" + ip + udp
        text = utc_text(ts) + (" | " + comment if comment else "")
        with self.lock:
            self.pending.append((ts, len(self.pending) + self.count, frame, text))
            if not self.buffered or len(self.pending) >= 2000:
                self._write_pending()

    def _write_pending(self):
        for ts, _n, frame, text in sorted(self.pending, key=lambda p: (p[0], p[1])):
            us = int(round(ts * 1e6))
            data = frame + b"\0" * (-len(frame) % 4)
            c = text.encode("utf-8")
            opts = struct.pack("<HH", 1, len(c)) + c + b"\0" * (-len(c) % 4) + struct.pack("<HH", 0, 0)
            body = struct.pack("<IIIII", 0, us >> 32, us & 0xffffffff, len(frame), len(frame)) + data + opts
            self.fh.write(struct.pack("<II", 6, 12 + len(body)) + body + struct.pack("<I", 12 + len(body)))
            self.count += 1
        self.pending = []
        self.fh.flush()

    def flush(self):
        """Write the packets kept so far, in time order (call when a recording has been fully decoded)."""
        with self.lock:
            self._write_pending()

    def gsm(self, blob, antenna=0, ts=None):
        """A gr-gsm message (GSMTAP header + GSM frame). antenna: 0 = RX1, 1 = RX2."""
        if not blob or len(blob) < 16:
            return
        blob = bytearray(blob)
        if antenna:
            blob[13] = antenna & 0xff                   # GSMTAP antenna number
        arfcn = ((blob[4] << 8) | blob[5]) & 0x3fff
        fn = int.from_bytes(bytes(blob[8:12]), "big")
        self._packet(bytes(blob), ts, comment="GSM ARFCN %d | frame %d%s" % (arfcn, fn,
                                                                           " | RX2" if antenna else ""))

    def lte_rrc(self, subtype, payload, earfcn=0, pci=0, ts=None, sfn=None, note=None):
        """An LTE RRC message: subtype LTE_RRC_BCCH_BCH (MIB) or LTE_RRC_BCCH_DL_SCH (SIBs).
        ts: when it was on air; sfn: its system frame number (for a MIB it is read from the MIB itself)."""
        if not payload:
            return
        if sfn is None and subtype == LTE_RRC_BCCH_BCH:
            sfn = mib_sfn(payload)
        # GSMTAP v2 header (16 bytes): version, header length in 32-bit words, type, timeslot, ARFCN,
        # signal dBm, SNR dB, frame number, sub type, antenna, sub slot, reserved
        hdr = struct.pack(">BBBBHbbIBBBB", 2, 4, GSMTAP_TYPE_LTE_RRC, 0, int(earfcn or 0) & 0x3fff, 0, 0,
                          int(sfn or 0) & 0x3ff, subtype, 0, 0, 0)
        what = "MIB" if subtype == LTE_RRC_BCCH_BCH else "SIB"
        comment = "LTE EARFCN %d | PCI %d | %s | SFN %s%s" % (int(earfcn or 0), int(pci or 0), what,
                                                            "?" if sfn is None else sfn, " | " + note if note else "")
        self._packet(hdr + bytes(payload), ts, sport=10000 + int(pci or 0), comment=comment)

    def close(self):
        try:
            with self.lock:
                self._write_pending()
            self.fh.close()
        except OSError:
            pass


class NRPcapWriter(object):
    """5G NR RRC messages (MIB, SIB1) for Wireshark. GSMTAP has no NR type, so this uses Wireshark's own
    "exported PDU" link type (252): each packet names the NR RRC dissector to use (nr-rrc.bcch.bch,
    nr-rrc.bcch.dl.sch) and Wireshark decodes it with that, as if it had captured it itself."""

    def __init__(self, path):
        self.path = os.path.expanduser(path)
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        self.fh = open(self.path, "wb")
        self.fh.write(struct.pack("<IHHiIII", 0xa1b2c3d4, 2, 4, 0, 0, 65535, 252))
        self.fh.flush()
        self.count = 0
        self.lock = threading.Lock()

    def nr_rrc(self, dissector, payload, ts=None):
        if not payload:
            return
        name = dissector.encode()
        # tag 12 = protocol (dissector) name (Wireshark reads exactly its length: no padding), tag 0 = end of tags
        frame = struct.pack(">HH", 12, len(name)) + name + struct.pack(">HH", 0, 0) + bytes(payload)
        ts = time.time() if ts is None else ts
        sec = int(ts)
        with self.lock:
            self.fh.write(struct.pack("<IIII", sec, int((ts - sec) * 1e6), len(frame), len(frame)) + frame)
            self.fh.flush()
            self.count += 1

    def close(self):
        try:
            self.fh.close()
        except OSError:
            pass


def mib_sfn(payload):
    """System frame number in a MIB (3 bytes): bandwidth 3 bits, PHICH 3 bits, then the SFN's top 8 bits."""
    if len(payload) < 3:
        return None
    v = int.from_bytes(bytes(payload[:3]), "big")
    return ((v >> 10) & 0xff) * 4


def add_args(ap):
    ap.add_argument("--pcap", default=None, metavar="FILE",
                    help="pcap file for the decoded messages (default ~/field/pcap/<tool>_<time>.pcap)")
    ap.add_argument("--no-pcap", action="store_true", help="do not write a pcap file")


def open_from_args(args, tool, nr=False):
    """PcapWriter (NRPcapWriter with nr=True) for the tool's --pcap / --no-pcap options, or None."""
    if getattr(args, "no_pcap", False):
        return None
    return (NRPcapWriter if nr else PcapWriter)(args.pcap or default_path(tool))


def report(writer):
    if writer is not None:
        import common as C
        writer.close()
        C.info("%d decoded message(s) saved to %s (open with Wireshark)" % (writer.count, writer.path))
        if writer.count:
            C.info("Wireshark: each packet's comment gives its UTC time, channel and cell. For clock times in the "
                   "Time column: View > Time Display Format > UTC Date and Time of Day")


def _ip_checksum(header):
    total = sum(struct.unpack(">10H", header))
    while total >> 16:
        total = (total & 0xffff) + (total >> 16)
    return (~total) & 0xffff


if __name__ == "__main__":
    # self-check: a file with one LTE MIB and one GSM-looking packet
    w = PcapWriter("/tmp/pcap_out_selftest.pcapng")
    w.lte_rrc(LTE_RRC_BCCH_BCH, bytes.fromhex("a8f080"), earfcn=6300, pci=430, ts=time.time() + 1)
    w.lte_rrc(LTE_RRC_BCCH_BCH, bytes.fromhex("a8f480"), earfcn=6300, pci=21, ts=time.time())
    w.close()
    print("wrote %d packet(s) to %s" % (w.count, w.path))
