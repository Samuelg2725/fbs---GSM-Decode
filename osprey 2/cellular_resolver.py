"""Turn a reported cell line into an exact downlink frequency.

The channel number is an EARFCN, UARFCN, or GSM ARFCN. This module only
converts that number. It does not transmit, and it does not read a subscriber
identity. Frequencies follow 3GPP TS 36.101, TS 25.101, and TS 45.005.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field

# Tenths of a megahertz keep the 100 kHz LTE raster exact.
# F_MHz = tenths / 10.  F_Hz = tenths * 100_000.


@dataclass(frozen=True)
class HardwareProfile:
    name: str
    min_mhz: float
    max_mhz: float


KRAKEN_POC = HardwareProfile("KrakenSDR PoC", 24.0, 1766.0)
WIDEBAND = HardwareProfile("Wideband DF", 24.0, 6000.0)
PROFILES = {
    "kraken_poc": KRAKEN_POC,
    "wideband": WIDEBAND,
}

# band, duplex, dl_low_tenths, dl_noffs, dl_min, dl_max,
# ul_low_tenths, ul_noffs, ul_min, ul_max
# SDL rows leave the uplink fields as None. TDD uses the downlink numbers on both links.
_LTE_ROWS = (
    (1, "FDD", 21100, 0, 0, 599, 19200, 18000, 18000, 18599),
    (2, "FDD", 19300, 600, 600, 1199, 18500, 18600, 18600, 19199),
    (3, "FDD", 18050, 1200, 1200, 1949, 17100, 19200, 19200, 19949),
    (4, "FDD", 21100, 1950, 1950, 2399, 17100, 19950, 19950, 20399),
    (5, "FDD", 8690, 2400, 2400, 2649, 8240, 20400, 20400, 20649),
    (6, "FDD", 8750, 2650, 2650, 2749, 8300, 20650, 20650, 20749),
    (7, "FDD", 26200, 2750, 2750, 3449, 25000, 20750, 20750, 21449),
    (8, "FDD", 9250, 3450, 3450, 3799, 8800, 21450, 21450, 21799),
    (9, "FDD", 18449, 3800, 3800, 4149, 17499, 21800, 21800, 22149),
    (10, "FDD", 21100, 4150, 4150, 4749, 17100, 22150, 22150, 22749),
    (11, "FDD", 14759, 4750, 4750, 4949, 14279, 22750, 22750, 22949),
    (12, "FDD", 7290, 5010, 5010, 5179, 6990, 23010, 23010, 23179),
    (13, "FDD", 7460, 5180, 5180, 5279, 7770, 23180, 23180, 23279),
    (14, "FDD", 7580, 5280, 5280, 5379, 7880, 23280, 23280, 23379),
    (17, "FDD", 7340, 5730, 5730, 5849, 7040, 23730, 23730, 23849),
    (18, "FDD", 8600, 5850, 5850, 5999, 8150, 23850, 23850, 23999),
    (19, "FDD", 8750, 6000, 6000, 6149, 8300, 24000, 24000, 24149),
    (20, "FDD", 7910, 6150, 6150, 6449, 8320, 24150, 24150, 24449),
    (21, "FDD", 14959, 6450, 6450, 6599, 14479, 24450, 24450, 24599),
    (22, "FDD", 35100, 6600, 6600, 7399, 34100, 24600, 24600, 25399),
    (23, "FDD", 21800, 7500, 7500, 7699, 20000, 25500, 25500, 25699),
    (24, "FDD", 15250, 7700, 7700, 8039, 16265, 25700, 25700, 26039),
    (25, "FDD", 19300, 8040, 8040, 8689, 18500, 26040, 26040, 26689),
    (26, "FDD", 8590, 8690, 8690, 9039, 8140, 26690, 26690, 27039),
    (27, "FDD", 8520, 9040, 9040, 9209, 8070, 27040, 27040, 27209),
    (28, "FDD", 7580, 9210, 9210, 9659, 7030, 27210, 27210, 27659),
    (29, "SDL", 7170, 9660, 9660, 9769, None, None, None, None),
    (30, "FDD", 23500, 9770, 9770, 9869, 23050, 27660, 27660, 27759),
    (31, "FDD", 4625, 9870, 9870, 9919, 4525, 27760, 27760, 27809),
    (32, "SDL", 14520, 9920, 9920, 10359, None, None, None, None),
    (33, "TDD", 19000, 36000, 36000, 36199, None, None, None, None),
    (34, "TDD", 20100, 36200, 36200, 36349, None, None, None, None),
    (35, "TDD", 18500, 36350, 36350, 36949, None, None, None, None),
    (36, "TDD", 19300, 36950, 36950, 37549, None, None, None, None),
    (37, "TDD", 19100, 37550, 37550, 37749, None, None, None, None),
    (38, "TDD", 25700, 37750, 37750, 38249, None, None, None, None),
    (39, "TDD", 18800, 38250, 38250, 38649, None, None, None, None),
    (40, "TDD", 23000, 38650, 38650, 39649, None, None, None, None),
    (41, "TDD", 24960, 39650, 39650, 41589, None, None, None, None),
    (42, "TDD", 34000, 41590, 41590, 43589, None, None, None, None),
    (43, "TDD", 36000, 43590, 43590, 45589, None, None, None, None),
    (44, "TDD", 7030, 45590, 45590, 46589, None, None, None, None),
    (45, "TDD", 14470, 46590, 46590, 46789, None, None, None, None),
    (46, "TDD", 51500, 46790, 46790, 54539, None, None, None, None),
    (47, "TDD", 58550, 54540, 54540, 55239, None, None, None, None),
    (48, "TDD", 35500, 55240, 55240, 56739, None, None, None, None),
    (49, "TDD", 35500, 56740, 56740, 58239, None, None, None, None),
    (50, "TDD", 14320, 58240, 58240, 59089, None, None, None, None),
    (51, "TDD", 14270, 59090, 59090, 59139, None, None, None, None),
    (52, "TDD", 33000, 59140, 59140, 60139, None, None, None, None),
    (53, "TDD", 24835, 60140, 60140, 60254, None, None, None, None),
    (54, "TDD", 16700, 60255, 60255, 60304, None, None, None, None),
    (65, "FDD", 21100, 65536, 65536, 66435, 19200, 131072, 131072, 131971),
    (66, "FDD", 21100, 66436, 66436, 67335, 17100, 131972, 131972, 132671),
    (67, "SDL", 7380, 67336, 67336, 67535, None, None, None, None),
    (68, "FDD", 7530, 67536, 67536, 67835, 6980, 132672, 132672, 132971),
    (69, "SDL", 25700, 67836, 67836, 68335, None, None, None, None),
    (70, "FDD", 19950, 68336, 68336, 68585, 16950, 132972, 132972, 133121),
    (71, "FDD", 6170, 68586, 68586, 68935, 6630, 133122, 133122, 133471),
    (72, "FDD", 4610, 68936, 68936, 68985, 4510, 133472, 133472, 133521),
    (73, "FDD", 4600, 68986, 68986, 69035, 4500, 133522, 133522, 133571),
    (74, "FDD", 14750, 69036, 69036, 69465, 14270, 133572, 133572, 134001),
    (75, "SDL", 14320, 69466, 69466, 70315, None, None, None, None),
    (76, "SDL", 14270, 70316, 70316, 70365, None, None, None, None),
    (85, "FDD", 7280, 70366, 70366, 70545, 6980, 134002, 134002, 134181),
    (87, "FDD", 4200, 70546, 70546, 70595, 4100, 134182, 134182, 134231),
    (88, "FDD", 4220, 70596, 70596, 70645, 4120, 134232, 134232, 134281),
    (103, "FDD", 7570, 70646, 70646, 70655, 7870, 134282, 134282, 134291),
    (106, "FDD", 9350, 70656, 70656, 70705, 8960, 134292, 134292, 134341),
    (111, "FDD", 18200, 73386, 73386, 73485, 18000, 134342, 134342, 134441),
)

def _utra_rows():
    """General FDD UARFCN ranges. Offset is in tenths of a megahertz."""
    rows = [
        (1, 10562, 10838, 0, 9612, 9888, 0),
        (2, 9662, 9938, 0, 9262, 9538, 0),
        (3, 1162, 1513, 15750, 937, 1288, 15250),
        (4, 1537, 1738, 18050, 1312, 1513, 14500),
        (5, 4357, 4458, 0, 4132, 4233, 0),
        (6, 4387, 4413, 0, 4162, 4188, 0),
        (7, 2237, 2563, 21750, 2012, 2338, 21000),
        (8, 2937, 3088, 3400, 2712, 2863, 3400),
        (9, 9237, 9387, 0, 8762, 8912, 0),
        (10, 3112, 3388, 14900, 2887, 3163, 11350),
        (11, 3712, 3787, 7360, 3487, 3562, 7330),
        (12, 3842, 3903, -370, 3617, 3678, -220),
        (13, 4017, 4043, -550, 3792, 3818, 210),
        (14, 4117, 4143, -630, 3892, 3918, 120),
        (19, 712, 763, 7350, 312, 363, 7700),
        (20, 4512, 4638, -1090, 4287, 4413, -230),
        (21, 862, 912, 13260, 462, 512, 13580),
        (22, 4662, 5038, 25800, 4437, 4813, 25250),
        # XXV DL 1930 MHz. 5112/5 = 1022.4, so offset = 1932.4 - 1022.4 = 910.0 MHz.
        # UL 1850 MHz. 4887/5 = 977.4, first carrier 1852.4, offset = 875.0 MHz.
        (25, 5112, 5413, 9100, 4887, 5188, 8750),
        # XXVI DL 859 MHz. 5762/5 = 1152.4, first 861.4, offset = -291.0 MHz.
        # UL 814 MHz. 5537/5 = 1107.4, first 816.4, offset = -291.0 MHz.
        (26, 5762, 5913, -2910, 5537, 5688, -2910),
    ]
    return rows


# Additional UARFCNs sit on a 100 kHz raster: F = offset_mhz + 0.2 * N.
# offset stored in tenths, so F_tenths = offset_tenths + 2 * N.
_UTRA_EXTRA_DL = (
    (2, (412, 437, 462, 487, 512, 537, 562, 587, 612, 637, 662, 687), 18501),
    (5, (1007, 1012, 1032, 1037, 1062, 1087), 6701),
    (6, (1037, 1062), 6701),
)
_UTRA_EXTRA_UL = (
    (2, (12, 37, 62, 87, 112, 137, 162, 187, 212, 237, 262, 287), 18501),
    (5, (782, 787, 807, 812, 837, 862), 6701),
    (6, (812, 837), 6701),
)


@dataclass
class Carrier:
    rat: str
    band: int
    band_name: str
    duplex: str
    channel: int
    downlink_tenths: int
    uplink_tenths: int | None

    @property
    def downlink_mhz(self) -> float:
        return self.downlink_tenths / 10.0

    @property
    def downlink_hz(self) -> int:
        return self.downlink_tenths * 100_000

    @property
    def uplink_mhz(self) -> float | None:
        if self.uplink_tenths is None:
            return None
        return self.uplink_tenths / 10.0

    @property
    def uplink_hz(self) -> int | None:
        if self.uplink_tenths is None:
            return None
        return self.uplink_tenths * 100_000


@dataclass
class ParsedTarget:
    channel: int
    channel_key: str
    pci: int | None = None
    cell_id: int | None = None
    lat: float | None = None
    lon: float | None = None


@dataclass
class Resolution:
    parsed: ParsedTarget
    primary: Carrier
    alternates: list[Carrier] = field(default_factory=list)
    kraken_poc_compatible: bool = False
    poc_priority_band: bool = False
    status_message: str = ""

    def to_dict(self) -> dict:
        data = asdict(self.parsed)
        data.update(_carrier_dict(self.primary))
        data["kraken_poc_compatible"] = self.kraken_poc_compatible
        data["poc_priority_band"] = self.poc_priority_band
        data["status_message"] = self.status_message
        data["alternates"] = [_carrier_dict(item) for item in self.alternates]
        return data


def _carrier_dict(carrier: Carrier) -> dict:
    return {
        "rat": carrier.rat,
        "band": carrier.band,
        "band_name": carrier.band_name,
        "duplex": carrier.duplex,
        "downlink_mhz": carrier.downlink_mhz,
        "downlink_hz": carrier.downlink_hz,
        "uplink_mhz": carrier.uplink_mhz,
        "uplink_hz": carrier.uplink_hz,
    }


def parse_target(text: str) -> ParsedTarget:
    """Read one bracketed target line. Colons, missing colons, and spaces are all accepted."""
    raw = text.strip()
    if raw.startswith("[") and raw.endswith("]"):
        raw = raw[1:-1]
    channel = None
    channel_key = None
    pci = None
    cell_id = None
    lat = None
    lon = None
    for piece in raw.split(","):
        piece = piece.strip()
        if not piece:
            continue
        match = re.fullmatch(r"([A-Za-z]+)\s*:?\s*([+-]?\d+(?::\d+)?(?:\.\d+)?)", piece.replace(" ", ""))
        if not match:
            raise ValueError("unrecognised field %r" % piece)
        key = match.group(1).upper()
        value = match.group(2)
        if key in ("EARFCN", "ARFCN", "UARFCN"):
            channel = int(value)
            channel_key = key
        elif key == "PCI":
            pci = int(value)
        elif key in ("CELLID", "CELL", "ECI", "CID"):
            cell_id = int(value)
        elif key in ("LAT", "LATITUDE"):
            lat = _coordinate(value)
        elif key in ("LONG", "LON", "LONGITUDE"):
            lon = _coordinate(value)
        else:
            raise ValueError("unknown field %s" % key)
    if channel is None or channel_key is None:
        raise ValueError("target has no EARFCN, ARFCN, or UARFCN")
    if pci is not None and not 0 <= pci <= 503:
        raise ValueError("PCI %s is outside 0-503" % pci)
    return ParsedTarget(channel, channel_key, pci, cell_id, lat, lon)


def _coordinate(value: str) -> float:
    """3:21654 is 3.21654. 53.4084 and 123 stay ordinary numbers."""
    dotted = re.fullmatch(r"([+-]?)(\d+):(\d+)", value)
    if dotted:
        sign, whole, frac = dotted.groups()
        return float("%s%s.%s" % (sign, whole, frac))
    return float(value)


def lte_carrier(earfcn: int) -> Carrier | None:
    for band, duplex, dl_low, dl_off, dl_min, dl_max, ul_low, ul_off, ul_min, ul_max in _LTE_ROWS:
        if not dl_min <= earfcn <= dl_max:
            continue
        dl = dl_low + (earfcn - dl_off)
        ul = None
        if duplex == "TDD":
            ul = dl
        elif duplex == "FDD" and ul_low is not None:
            paired = ul_off + (earfcn - dl_off)
            if ul_min <= paired <= ul_max:
                ul = ul_low + (paired - ul_off)
        return Carrier("LTE", band, "Band %d" % band, duplex, earfcn, dl, ul)
    return None


def utra_carrier(uarfcn: int) -> Carrier | None:
    for band, dl_min, dl_max, dl_off, ul_min, ul_max, ul_off in _utra_rows():
        if dl_min <= uarfcn <= dl_max:
            dl = dl_off + 2 * uarfcn
            span = uarfcn - dl_min
            paired = ul_min + span
            ul = ul_off + 2 * paired if ul_min <= paired <= ul_max else None
            return Carrier("UMTS", band, "Band %d" % band, "FDD", uarfcn, dl, ul)
    for band, channels, offset in _UTRA_EXTRA_DL:
        if uarfcn in channels:
            dl = offset + 2 * uarfcn
            ul = _paired_extra_uplink(band, channels.index(uarfcn))
            return Carrier("UMTS", band, "Band %d" % band, "FDD", uarfcn, dl, ul)
    return None


def _paired_extra_uplink(band: int, index: int) -> int | None:
    for item_band, channels, offset in _UTRA_EXTRA_UL:
        if item_band == band and index < len(channels):
            return offset + 2 * channels[index]
    return None


def gsm_carriers(arfcn: int) -> list[Carrier]:
    """Every GSM plan that owns this ARFCN. 512-810 belongs to both DCS-1800 and PCS-1900."""
    found = []
    if 1 <= arfcn <= 124:
        found.append(_gsm("P-GSM 900", "GSM-900", arfcn, 9350 + 2 * arfcn, 8900 + 2 * arfcn))
    if arfcn == 0 or 975 <= arfcn <= 1023:
        shift = arfcn - 1024
        found.append(_gsm("E-GSM 900", "GSM-900", arfcn, 9350 + 2 * shift, 8900 + 2 * shift))
    if 128 <= arfcn <= 251:
        uplink = 8242 + 2 * (arfcn - 128)
        found.append(_gsm("GSM-850", "GSM-850", arfcn, uplink + 450, uplink))
    if 512 <= arfcn <= 885:
        uplink = 17102 + 2 * (arfcn - 512)
        found.append(_gsm("DCS-1800", "DCS-1800", arfcn, uplink + 950, uplink))
    if 512 <= arfcn <= 810:
        uplink = 18502 + 2 * (arfcn - 512)
        found.append(_gsm("PCS-1900", "PCS-1900", arfcn, uplink + 800, uplink))
    return found


def _gsm(name: str, family: str, arfcn: int, downlink: int, uplink: int) -> Carrier:
    return Carrier("GSM", 0, name, "FDD", arfcn, downlink, uplink)


def resolve(text: str) -> Resolution:
    parsed = parse_target(text)
    primary, alternates = _select(parsed)
    if primary is None:
        raise ValueError("channel %s %s is not a known LTE, UMTS, or GSM channel" % (parsed.channel_key, parsed.channel))
    compatible = KRAKEN_POC.min_mhz <= primary.downlink_mhz <= KRAKEN_POC.max_mhz
    priority = _priority(primary)
    return Resolution(
        parsed=parsed,
        primary=primary,
        alternates=alternates,
        kraken_poc_compatible=compatible,
        poc_priority_band=priority and compatible,
        status_message=_status(primary, compatible, priority),
    )


def _select(parsed: ParsedTarget):
    channel = parsed.channel
    key = parsed.channel_key
    lte = lte_carrier(channel)
    umts = utra_carrier(channel)
    gsm = gsm_carriers(channel) if channel <= 1023 else []
    pci_marks_lte = parsed.pci is not None and 0 <= parsed.pci <= 503

    if key == "UARFCN":
        return umts, _others(umts, lte, gsm)
    if key == "EARFCN" or pci_marks_lte or channel > 1023:
        extra_gsm = gsm if channel <= 1023 else []
        return lte, _others(lte, None, extra_gsm)
    # Ambiguous ARFCN at or below 1023. Show GSM first and LTE beside it.
    primary = gsm[0] if gsm else lte
    return primary, _others(primary, lte, gsm[1:] if gsm else [])


def _others(primary, lte, gsm):
    pool = []
    if lte is not None:
        pool.append(lte)
    pool.extend(gsm or [])
    if primary is None:
        return pool
    return [item for item in pool if (item.rat, item.band_name, item.downlink_tenths) != (primary.rat, primary.band_name, primary.downlink_tenths)]


def _priority(carrier: Carrier) -> bool:
    if carrier.rat == "LTE" and carrier.band in (8, 20, 28):
        return True
    if carrier.rat == "GSM" and carrier.band_name in ("P-GSM 900", "E-GSM 900", "GSM-850"):
        return True
    return False


def _status(carrier: Carrier, compatible: bool, priority: bool) -> str:
    mhz = "%.1f" % carrier.downlink_mhz
    if not compatible:
        if carrier.downlink_mhz > KRAKEN_POC.max_mhz:
            return "OUT OF RANGE FOR KRAKEN POC (%s MHz > 1766 MHz limit; requires Wideband DF hardware)" % mhz
        return "OUT OF RANGE FOR KRAKEN POC (%s MHz < 24 MHz limit)" % mhz
    if carrier.rat == "LTE" and carrier.band == 20:
        label = "Band 20 800MHz"
    elif carrier.rat == "LTE" and carrier.band == 28:
        label = "Band 28 700MHz"
    elif carrier.rat == "LTE" and carrier.band == 8:
        label = "Band 8 900MHz"
    elif carrier.rat == "GSM" and carrier.band_name in ("P-GSM 900", "E-GSM 900"):
        label = "GSM-900"
    elif carrier.rat == "GSM" and carrier.band_name == "GSM-850":
        label = "GSM-850"
    elif priority:
        label = carrier.band_name
    else:
        label = "%s %s" % (carrier.band_name, mhz + "MHz")
    return "READY (KrakenSDR PoC - %s)" % label
