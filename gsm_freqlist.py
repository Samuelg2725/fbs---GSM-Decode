#!/usr/bin/env python3
"""Decoder for GSM frequency lists (3GPP TS 44.018 sections 10.5.2.1b, 10.5.2.13, 10.5.2.22, 10.5.2.22a).

Used by gsm_neigh.py for the neighbour lists a cell broadcasts:
  SI2    (type 0x1a)  Neighbour Cell Description          -> BCCH frequencies of neighbour cells
  SI2bis (type 0x02)  Neighbour Cell Description (ext.)   -> more of them (when SI2 has EXT-IND set)
  SI2ter (type 0x03)  Neighbour Cell Description 2        -> other-band neighbours (e.g. 1800 from 900)
Formats: bit map 0, range 1024/512/256/128, variable bit map. Pure Python, no radio needed.
"""


class Bits(object):
    def __init__(self, data, pos=0):
        self.data = bytes(data)
        self.pos = pos

    def get(self, n):
        v = 0
        for _ in range(n):
            byte = self.pos >> 3
            if byte >= len(self.data):
                raise EOFError
            bit = 7 - (self.pos & 7)
            v = (v << 1) | ((self.data[byte] >> bit) & 1)
            self.pos += 1
        return v


def _gp2(n):
    """Greatest power of 2 <= n."""
    p = 1
    while p * 2 <= n:
        p *= 2
    return p


def _range_values(w, rng):
    """44.018 Annex J decoding: W(1..k) -> F(1..k) for a range of size rng (1024, 512, 256, 128)."""
    out = []
    for k in range(1, len(w)):
        if w[k] == 0:
            break
        index = k
        j = _gp2(index)
        n = w[index]
        while index > 1:
            if 2 * index < 3 * j:                       # left child
                index -= j // 2
                n = (n + w[index] - rng // j - 1) % (2 * rng // j - 1) + 1
            else:                                        # right child
                index -= j
                n = (n + w[index] - 1) % (2 * rng // j - 1) + 1
            j //= 2
        out.append(n)
    return out


# W(k) field sizes per range (index 1..)
_WSIZES = {
    1024: [10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 6],
    512: [9, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5],
    256: [8, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4],
    128: [7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4] + [3] * 13,
}


def _read_w(bits, rng):
    w = [0]
    for size in _WSIZES[rng]:
        try:
            w.append(bits.get(size))
        except EOFError:
            break
    return w


def decode(octets, mask=0xce):
    """16-octet frequency list -> (format name, sorted list of ARFCNs).
    mask: which octet-1 bits are format bits (0xce for SI2/SI2bis, 0x8e for SI2ter)."""
    cd = bytes(octets)
    if len(cd) < 2:
        return "short", []
    b0 = cd[0]
    found = set()
    if (b0 & 0xc0 & mask) == 0x00:                       # bit map 0: ARFCN 1..124
        if len(cd) < 16:
            return "bitmap0 (short)", []
        for i in range(1, 125):
            if cd[15 - ((i - 1) >> 3)] & (1 << ((i - 1) & 7)):
                found.add(i)
        return "bitmap0", sorted(found)
    if (b0 & 0xc8 & mask) == 0x80:                       # range 1024
        if b0 & 0x04:
            found.add(0)                                 # F0 = ARFCN 0
        bits = Bits(cd, 6)
        w = _read_w(bits, 1024)
        found.update(_range_values(w, 1024))
        return "range1024", sorted(found)
    fmt = b0 & 0xce & mask
    rng = {0x88: 512, 0x8a: 256, 0x8c: 128}.get(fmt)
    bits = Bits(cd, 7)
    try:
        orig = bits.get(10)
    except EOFError:
        return "short", []
    if rng:
        found.add(orig)
        w = _read_w(bits, rng)
        for f in _range_values(w, rng):
            found.add((orig + f) % 1024)
        return "range%d" % rng, sorted(found)
    if fmt == 0x8e:                                      # variable bit map
        found.add(orig)
        for k in range(1, 112):
            try:
                if bits.get(1):
                    found.add((orig + k) % 1024)
            except EOFError:
                break
        return "varbitmap", sorted(found)
    return "unknown format 0x%02x" % b0, []


SI_NAMES = {0x1a: "SI2", 0x02: "SI2bis", 0x03: "SI2ter", 0x07: "SI2quater", 0x1b: "SI3", 0x1c: "SI4"}


def rr_body(raw, types):
    """Find an RR message of one of `types` in a gr-gsm decoded frame (GSMTAP header + L2 pseudo length
    + PD 0x06 + type). Returns (type, body) with body[2] == type, or (None, None)."""
    if not raw:
        return None, None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)
    starts.extend((16, 0))
    for off in starts:
        if off + 19 <= len(raw) and (raw[off + 1] & 0x0f) == 0x06 and raw[off + 2] in types:
            return raw[off + 2], raw[off:]
    return None, None


def neighbours(raw):
    """(message name, format, ARFCN list, ext_ind) for SI2/SI2bis/SI2ter frames, else None."""
    t, body = rr_body(raw, (0x1a, 0x02, 0x03))
    if t is None:
        return None
    octets = body[3:19]
    mask = 0x8e if t == 0x03 else 0xce
    fmt, arfcns = decode(octets, mask)
    ext = bool(octets[0] & 0x20) if t in (0x1a, 0x02) else False
    return SI_NAMES[t], fmt, arfcns, ext


if __name__ == "__main__":
    # self-check with a bit map 0 list containing ARFCNs 1, 6, 13, 22 and 124
    data = bytearray(16)
    for a in (1, 6, 13, 22, 124):
        data[15 - ((a - 1) >> 3)] |= 1 << ((a - 1) & 7)
    print(decode(data))
