"""LTE bit operations used by the cell-identity decoder.

The polynomials and the interleavers are the ones in 3GPP TS 36.212.
This module does not touch the radio.
"""

from __future__ import annotations

# Tail-biting convolutional code, constraint length 7.
# The spec writes the generators as octal 133, 171, 165, with the newest bit
# on the left of that number. This encoder shifts the newest bit in at the
# bottom, so the taps are the bit-reversed values 0x6D, 0x4F, 0x57.
_CONV_GEN = (0x6D, 0x4F, 0x57)


def _ones(value: int) -> int:
    """Popcount. int.bit_count() is missing on the Pi's Python 3.9."""
    return bin(value).count("1")
_CONV_COLUMNS = (
    0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,
    1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31,
)

# Quadratic-permutation interleaver constants, TS 36.212 Table 5.1.3-3.
# Only the block sizes a SIB1-sized transport block can land on.
_QPP = {
    40: (3, 10), 48: (7, 12), 56: (19, 42), 64: (7, 16), 72: (7, 18),
    80: (11, 20), 88: (5, 22), 96: (11, 24), 104: (7, 26), 112: (41, 84),
    120: (103, 90), 128: (15, 32), 136: (9, 34), 144: (17, 108), 152: (9, 38),
    160: (21, 120), 168: (101, 84), 176: (21, 44), 184: (57, 46), 192: (23, 48),
    200: (13, 50), 208: (27, 52), 216: (11, 36), 224: (27, 56), 232: (85, 58),
    240: (29, 60), 248: (33, 62), 256: (15, 32), 264: (17, 198), 272: (33, 68),
    280: (103, 210), 288: (19, 36), 296: (19, 74), 304: (37, 76), 312: (19, 78),
    320: (21, 120), 328: (21, 82), 336: (115, 84), 344: (193, 86), 352: (21, 44),
    360: (133, 90), 368: (81, 46), 376: (45, 94), 384: (23, 48), 392: (243, 98),
    400: (151, 40), 408: (155, 102), 416: (25, 52), 424: (51, 106), 432: (47, 72),
    440: (91, 110), 448: (29, 168), 456: (29, 114), 464: (247, 58), 472: (29, 118),
    480: (89, 180), 488: (91, 122), 496: (157, 62), 504: (55, 84), 512: (31, 64),
}


def bits_of(value: int, width: int) -> list[int]:
    return [(value >> (width - 1 - i)) & 1 for i in range(width)]


def bits_to_int(bits) -> int:
    value = 0
    for bit in bits:
        value = (value << 1) | (int(bit) & 1)
    return value


def crc_bits(data, polynomial: int, width: int) -> list[int]:
    """Remainder of the polynomial division, highest-degree bit first.

    The message bits are the high-order coefficients. polynomial is the
    generator without its degree-width term, so CRC-16 is 0x1021.
    """
    register = 0
    mask = (1 << width) - 1
    top = 1 << (width - 1)
    for bit in data:
        feedback = ((register & top) >> (width - 1)) ^ (int(bit) & 1)
        register = ((register << 1) & mask)
        if feedback:
            register ^= polynomial
    return bits_of(register, width)


def crc16(data) -> list[int]:
    # D^16 + D^12 + D^5 + 1.
    return crc_bits(data, 0x1021, 16)


def crc24a(data) -> list[int]:
    # D^24 + D^23 + D^18 + D^17 + D^14 + D^11 + D^10 + D^7 + D^6 + D^5 + D^4 + D^3 + D + 1.
    return crc_bits(data, 0x864CFB, 24)


def gold(c_init: int, length: int) -> list[int]:
    """The LTE pseudo-random sequence c(n), TS 36.211 7.2."""
    total = length + 1600
    x1 = [0] * (total + 31)
    x2 = [0] * (total + 31)
    x1[0] = 1
    for index in range(31):
        x2[index] = (int(c_init) >> index) & 1
    for index in range(total):
        x1[index + 31] = x1[index + 3] ^ x1[index]
        x2[index + 31] = x2[index + 3] ^ x2[index + 2] ^ x2[index + 1] ^ x2[index]
    return [x1[index + 1600] ^ x2[index + 1600] for index in range(length)]


def qpsk_mod(bits) -> list[complex]:
    symbols = []
    scale = 1.0 / (2 ** 0.5)
    for index in range(0, len(bits) - 1, 2):
        symbols.append(scale * ((1 - 2 * bits[index]) + 1j * (1 - 2 * bits[index + 1])))
    return symbols


def qpsk_soft(symbol: complex) -> tuple[float, float]:
    """Positive means the bit is more likely 0."""
    return (float(symbol.real), float(symbol.imag))


def _conv_state_after(bits) -> int:
    state = 0
    for bit in bits:
        state = ((state << 1) | (int(bit) & 1)) & 0x3F
    return state


def conv_encode(bits) -> list[list[int]]:
    """Three tail-biting streams, each as long as the message."""
    state = _conv_state_after(bits[-6:])
    streams = [[], [], []]
    for bit in bits:
        state = ((state << 1) | (int(bit) & 1)) & 0x7F
        for gen_index, polynomial in enumerate(_CONV_GEN):
            parity = _ones(state & polynomial) & 1
            streams[gen_index].append(parity)
        state &= 0x3F
    return streams


def _transition_table(generators):
    """dest, parity for every register state and input bit. Built once."""
    table = []
    for state in range(64):
        row = []
        for bit in (0, 1):
            entered = ((state << 1) | bit) & 0x7F
            parity = tuple(_ones(entered & polynomial) & 1 for polynomial in generators)
            row.append((entered & 0x3F, parity))
        table.append(row)
    return table


# Live cells use the bit-reversed taps. The octal table is the other convention.
_TRANS = _transition_table(_CONV_GEN)
_TRANS_OCTAL = _transition_table((0o133, 0o171, 0o165))
VITERBI_TABLES = (_TRANS, _TRANS_OCTAL)


def viterbi_tailbiting(soft_streams, trans=None) -> list[int]:
    """Tail-biting decode. Two laps of the block, which is enough at the MIB's code rate.

    The old search tried all 64 starting states. That is exact and far too slow
    on the Pi. Repeating the block once lets the register settle, then the
    second lap is the message.
    """
    if trans is None:
        trans = _TRANS
    length = len(soft_streams[0])
    score = [-1e300] * 64
    score[0] = 0.0
    history = []
    steps = length * 2
    for step in range(steps):
        metrics = (
            soft_streams[0][step % length],
            soft_streams[1][step % length],
            soft_streams[2][step % length],
        )
        nxt = [-1e300] * 64
        prev = [0] * 64
        choice = [0] * 64
        for state in range(64):
            current = score[state]
            for bit in (0, 1):
                dest, parity = trans[state][bit]
                branch = (
                    (metrics[0] if parity[0] == 0 else -metrics[0])
                    + (metrics[1] if parity[1] == 0 else -metrics[1])
                    + (metrics[2] if parity[2] == 0 else -metrics[2])
                )
                total = current + branch
                if total > nxt[dest]:
                    nxt[dest] = total
                    prev[dest] = state
                    choice[dest] = bit
        history.append((prev, choice))
        score = nxt
    state = max(range(64), key=score.__getitem__)
    bits = [0] * length
    for step in range(steps - 1, steps - 1 - length, -1):
        prev, choice = history[step]
        bits[step % length] = choice[state]
        state = prev[state]
    return bits


def viterbi_exact(soft_streams, trans=None) -> list[int]:
    """Maximum-likelihood tail-biting decode. Used when the fast pass is one bit off."""
    if trans is None:
        trans = _TRANS
    length = len(soft_streams[0])
    best_bits = [0] * length
    best_score = -1e300
    for start in range(64):
        score = [-1e300] * 64
        score[start] = 0.0
        history = []
        for step in range(length):
            metrics = (soft_streams[0][step], soft_streams[1][step], soft_streams[2][step])
            nxt = [-1e300] * 64
            prev = [0] * 64
            choice = [0] * 64
            for state in range(64):
                current = score[state]
                for bit in (0, 1):
                    dest, parity = trans[state][bit]
                    branch = (
                        (metrics[0] if parity[0] == 0 else -metrics[0])
                        + (metrics[1] if parity[1] == 0 else -metrics[1])
                        + (metrics[2] if parity[2] == 0 else -metrics[2])
                    )
                    total = current + branch
                    if total > nxt[dest]:
                        nxt[dest] = total
                        prev[dest] = state
                        choice[dest] = bit
            history.append((prev, choice))
            score = nxt
        if score[start] <= best_score:
            continue
        best_score = score[start]
        state = start
        bits = [0] * length
        for step in range(length - 1, -1, -1):
            prev, choice = history[step]
            bits[step] = choice[state]
            state = prev[state]
        best_bits = bits
    return best_bits


def _interleave_stream(bits):
    """Column interleaver used by convolutional rate matching. None marks a hole."""
    length = len(bits)
    rows = (length + 31) // 32
    padded = rows * 32
    holes = padded - length
    matrix = [None] * padded
    for index, bit in enumerate(bits):
        matrix[holes + index] = bit
    output = []
    for column in _CONV_COLUMNS:
        for row in range(rows):
            output.append(matrix[row * 32 + column])
    return output


def conv_rate_map(length: int):
    """Where each transmitted convolutional bit sits in the three streams.

    The PBCH reads this map repeatedly until it has filled 1920 bits.
    """
    holes = []
    for stream in range(3):
        interleaved = _interleave_stream([0] * length)
        holes.append(interleaved)
    mapping = []
    # Rebuild with indices instead of the dummy zeros.
    indexed = []
    for stream in range(3):
        interleaved = _interleave_stream(list(range(length)))
        indexed.append(interleaved)
    collected = indexed[0] + indexed[1] + indexed[2]
    for item in collected:
        if item is not None:
            # item is the bit index; stream is implied by which third we are in.
            pass
    # The dummy path above loses the stream number. Build it directly.
    mapping = []
    per_stream = []
    for stream in range(3):
        per_stream.append(_interleave_stream([(stream, index) for index in range(length)]))
    for item in per_stream[0] + per_stream[1] + per_stream[2]:
        if item is not None:
            mapping.append(item)
    return mapping


def conv_rate_match(streams, count: int) -> list[int]:
    mapping = conv_rate_map(len(streams[0]))
    output = []
    cursor = 0
    while len(output) < count:
        stream, index = mapping[cursor % len(mapping)]
        output.append(streams[stream][index])
        cursor += 1
    return output


def conv_rate_recover(soft, length: int):
    """Add repeated soft bits back onto the three convolutional streams."""
    mapping = conv_rate_map(length)
    streams = [[0.0] * length for _ in range(3)]
    for index, metric in enumerate(soft):
        stream, bit_index = mapping[index % len(mapping)]
        streams[stream][bit_index] += metric
    return streams


def code_block_size(bit_count: int) -> int:
    """Smallest turbo block that can hold bit_count bits."""
    for size in sorted(_QPP):
        if size >= bit_count:
            return size
    raise ValueError("transport block is larger than this decoder handles")


def _rsc_tables():
    parity = [[0, 0] for _ in range(8)]
    nxt = [[0, 0] for _ in range(8)]
    for state in range(8):
        s0 = state & 1
        s1 = (state >> 1) & 1
        s2 = (state >> 2) & 1
        for bit in (0, 1):
            feedback = bit ^ s1 ^ s2
            parity[state][bit] = feedback ^ s0 ^ s2
            nxt[state][bit] = feedback | (s0 << 1) | (s1 << 2)
    return nxt, parity


_RSC_NEXT, _RSC_PARITY = _rsc_tables()


def _encode_rsc(bits, state=0):
    parity = []
    for bit in bits:
        parity.append(_RSC_PARITY[state][bit])
        state = _RSC_NEXT[state][bit]
    return parity, state


def _termination_bits(state):
    bits = []
    parity = []
    for _ in range(3):
        s1 = (state >> 1) & 1
        s2 = (state >> 2) & 1
        bit = s1 ^ s2
        parity.append(_RSC_PARITY[state][bit])
        state = _RSC_NEXT[state][bit]
        bits.append(bit)
    return bits, parity


def turbo_encode(bits) -> list[list[int]]:
    """Three streams of length K+4, including the LTE tail-bit arrangement."""
    size = len(bits)
    f1, f2 = _QPP[size]
    order = [((f1 * index + f2 * index * index) % size) for index in range(size)]
    interleaved = [bits[index] for index in order]
    par_a, state_a = _encode_rsc(bits, 0)
    par_b, state_b = _encode_rsc(interleaved, 0)
    tail_a, tail_par_a = _termination_bits(state_a)
    tail_b, tail_par_b = _termination_bits(state_b)
    systematic = list(bits)
    first = list(par_a)
    second = list(par_b)
    systematic += [tail_a[0], tail_par_a[1], tail_b[0], tail_par_b[1]]
    first += [tail_par_a[0], tail_a[2], tail_par_b[0], tail_b[2]]
    second += [tail_a[1], tail_par_a[2], tail_b[1], tail_par_b[2]]
    return [systematic, first, second]


def _branch(systematic, parity, prior, state, bit):
    sign = 1.0 if bit == 0 else -1.0
    metric = 0.5 * sign * (systematic + prior)
    if _RSC_PARITY[state][bit] == 0:
        metric += 0.5 * parity
    else:
        metric -= 0.5 * parity
    return metric


def _max_log_map(systematic, parity, prior):
    """Extrinsic LLRs. Positive favours bit 0. The trellis starts and ends at state 0."""
    length = len(systematic)
    forward = [[None] * 8 for _ in range(length + 1)]
    forward[0][0] = 0.0
    for index in range(length):
        for state in range(8):
            if forward[index][state] is None:
                continue
            for bit in (0, 1):
                score = forward[index][state] + _branch(systematic[index], parity[index], prior[index], state, bit)
                dest = _RSC_NEXT[state][bit]
                if forward[index + 1][dest] is None or score > forward[index + 1][dest]:
                    forward[index + 1][dest] = score
        peak = max(value for value in forward[index + 1] if value is not None)
        forward[index + 1] = [None if value is None else value - peak for value in forward[index + 1]]
    # The information bits do not end in state 0. The tail bits do that
    # afterwards, so the backward pass starts with every state equally likely.
    backward = [0.0] * 8
    extrinsic = [0.0] * length
    for index in range(length - 1, -1, -1):
        best = [None, None]
        nxt = [None] * 8
        for state in range(8):
            if forward[index][state] is None:
                continue
            for bit in (0, 1):
                dest = _RSC_NEXT[state][bit]
                if backward[dest] is None:
                    continue
                gamma = _branch(systematic[index], parity[index], prior[index], state, bit)
                total = forward[index][state] + gamma + backward[dest]
                if best[bit] is None or total > best[bit]:
                    best[bit] = total
                if nxt[state] is None or gamma + backward[dest] > nxt[state]:
                    nxt[state] = gamma + backward[dest]
        if best[0] is not None and best[1] is not None:
            extrinsic[index] = (best[0] - best[1]) - systematic[index] - prior[index]
        peak = max((value for value in nxt if value is not None), default=0.0)
        backward = [None if value is None else value - peak for value in nxt]
    return extrinsic


def turbo_decode(soft_streams, iterations: int = 4) -> list[int]:
    """Decode one code block. soft_streams are length K+4, positive favours bit 0."""
    systematic, first, second = soft_streams
    size = len(systematic) - 4
    f1, f2 = _QPP[size]
    order = [((f1 * index + f2 * index * index) % size) for index in range(size)]
    inverse = [0] * size
    for index, source in enumerate(order):
        inverse[source] = index
    prior_a = [0.0] * size
    decoded = [0] * size
    for _ in range(iterations):
        extrinsic_a = _max_log_map(systematic[:size], first[:size], prior_a)
        prior_b = [extrinsic_a[order[index]] for index in range(size)]
        extrinsic_b = _max_log_map(
            [systematic[order[index]] for index in range(size)],
            second[:size],
            prior_b,
        )
        prior_a = [extrinsic_b[inverse[index]] for index in range(size)]
        posterior = [systematic[index] + prior_a[index] + extrinsic_a[index] for index in range(size)]
        decoded = [0 if value >= 0 else 1 for value in posterior]
    return decoded


def _turbo_interleave(values):
    length = len(values)
    rows = (length + 31) // 32
    holes = rows * 32 - length
    matrix = [None] * (rows * 32)
    for index, value in enumerate(values):
        matrix[holes + index] = value
    output = []
    for column in _CONV_COLUMNS:
        for row in range(rows):
            output.append(matrix[row * 32 + column])
    return output


def _turbo_interleave_second(values):
    """The parity-2 interleaver, which is rotated by one position."""
    plain = _turbo_interleave(values)
    count = len(plain)
    return [plain[(index + 1) % count] for index in range(count)]


def turbo_rate_match(streams, count: int, redundancy: int = 0) -> list:
    """Circular-buffer selection. None entries are filler and are not transmitted."""
    first = _turbo_interleave(streams[0])
    second = _turbo_interleave(streams[1])
    third = _turbo_interleave_second(streams[2])
    buffer = first + second + third
    rows = (len(streams[0]) + 31) // 32
    start = rows * (2 * ((len(buffer) + 8 * rows - 1) // (8 * rows)) * redundancy + 2)
    output = []
    cursor = start
    while len(output) < count:
        item = buffer[cursor % len(buffer)]
        if item is not None:
            output.append(item)
        cursor += 1
        if cursor > start + len(buffer) * 4 and len(output) < count:
            break
    return output


def turbo_rate_recover(soft, block_size: int, redundancy: int = 0, filler: int = 0):
    """Scatter received soft bits back onto the three streams of length K+4."""
    length = block_size + 4
    systematic_ids = [None if index < filler else index for index in range(block_size)]
    systematic_ids += list(range(block_size, length))
    first = _turbo_interleave(systematic_ids)
    second = _turbo_interleave(list(range(length)))
    third = _turbo_interleave_second(list(range(length)))
    buffer = [("s", item) for item in first] + [("p", item) for item in second] + [("q", item) for item in third]
    rows = (length + 31) // 32
    start = rows * (2 * ((len(buffer) + 8 * rows - 1) // (8 * rows)) * redundancy + 2)
    streams = [[0.0] * length for _ in range(3)]
    names = {"s": 0, "p": 1, "q": 2}
    produced = 0
    cursor = start
    guard = start + len(buffer) * 4
    while produced < len(soft) and cursor < guard:
        kind, index = buffer[cursor % len(buffer)]
        if index is not None:
            streams[names[kind]][index] += soft[produced]
            produced += 1
        cursor += 1
    return streams
