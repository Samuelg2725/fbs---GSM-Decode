"""Shared helpers for the GSM and LTE tools: colours, result tables, operator names, timestamps.

Colour is used when writing to a terminal or a `tee` pipe; set the environment variable NO_COLOR=1
(or pass --no-color to a tool) to turn it off.
"""
import os
import re
import sys
import time
from datetime import datetime, timezone

# ------------------------------------------------------------------------------------------- colours
_ANSI = {"red": "31", "green": "32", "yellow": "33", "blue": "34", "magenta": "35", "cyan": "36",
         "grey": "90", "bold": "1"}
_USE_COLOUR = [os.environ.get("NO_COLOR") is None]


def set_colour(enabled):
    _USE_COLOUR[0] = bool(enabled)


def c(text, colour):
    """Wrap text in an ANSI colour (or return it unchanged when colour is off)."""
    if not _USE_COLOUR[0] or not colour:
        return str(text)
    codes = ";".join(_ANSI[p] for p in colour.split("+") if p in _ANSI)
    return "\x1b[%sm%s\x1b[0m" % (codes, text)


def strip(text):
    return re.sub("\x1b\\[[0-9;]*m", "", str(text))


_STATUS = [False]


def status(text):
    """Show a one-line progress message that the next printed line replaces (terminal only)."""
    if sys.stdout.isatty():
        sys.stdout.write("\r\x1b[K" + c(text, "grey"))
        sys.stdout.flush()
        _STATUS[0] = True


def clear_status():
    """Remove the progress line shown by status()."""
    if _STATUS[0]:
        sys.stdout.write("\r\x1b[K")
        sys.stdout.flush()
        _STATUS[0] = False


def heading(title):
    """Section heading line."""
    clear_status()
    print("", flush=True)
    print(c("== %s " % title + "=" * max(0, 70 - len(title)), "cyan+bold"), flush=True)


def info(text):
    clear_status()
    print(c(text, "grey"), flush=True)


def warn(text):
    clear_status()
    print(c(text, "yellow"), flush=True)


# --------------------------------------------------------------------------------------------- tables
class Table(object):
    """Fixed-width text table with optional per-cell colour.

        t = Table(["ARFCN", "Cell ID", "Operator"])
        t.add([13, 942, ("Vodafone", "red")])     # a cell can be a value or (value, colour)
        t.show()
    """

    def __init__(self, columns, align=None):
        self.columns = list(columns)
        self.align = align or {}
        self.rows = []

    def add(self, cells, colour=None):
        """Add a row. colour (optional) is applied to every cell without its own colour."""
        row = []
        for cell in cells:
            if isinstance(cell, tuple):
                row.append((str(cell[0]), cell[1]))
            else:
                row.append(("" if cell is None else str(cell), colour))
        self.rows.append(row)

    def show(self, title=None):
        if title:
            heading(title)
        if not self.rows:
            print(c("  (nothing to show)", "grey"), flush=True)
            return
        widths = [len(h) for h in self.columns]
        for row in self.rows:
            for i, (val, _col) in enumerate(row):
                widths[i] = max(widths[i], len(val))

        def fmt(i, val):
            return val.rjust(widths[i]) if self.align.get(self.columns[i]) == "right" else val.ljust(widths[i])
        line = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
        print(c(line, "grey"), flush=True)
        print(c("|", "grey") + c("|", "grey").join(" " + c(fmt(i, h), "bold") + " "
                                                   for i, h in enumerate(self.columns)) + c("|", "grey"), flush=True)
        print(c(line, "grey"), flush=True)
        for row in self.rows:
            print(c("|", "grey") + c("|", "grey").join(" " + c(fmt(i, val), col) + " "
                                                       for i, (val, col) in enumerate(row)) + c("|", "grey"),
                  flush=True)
        print(c(line, "grey"), flush=True)


# -------------------------------------------------------------------------------------------- networks
OPERATORS = {
    (234, 10): "O2", (234, 15): "Vodafone", (234, 20): "Three", (234, 30): "EE", (234, 31): "EE",
    (234, 32): "EE", (234, 33): "EE", (234, 50): "JT", (234, 55): "Sure", (234, 58): "Manx Telecom",
}
OPERATOR_COLOUR = {"Vodafone": "red", "O2": "blue", "EE": "cyan", "Three": "magenta"}


def operator(mcc, mnc):
    try:
        return OPERATORS.get((int(mcc), int(mnc)), "")
    except (TypeError, ValueError):
        return ""


def op_cell(name):
    """Table cell for an operator name, coloured per network."""
    return (name or "unknown", OPERATOR_COLOUR.get(name, ""))


# -------------------------------------------------------------------------------------------- utilities
def now_stamp():
    t = datetime.now(timezone.utc)
    return t.strftime("%H:%M:%S.") + "%03d" % (t.microsecond // 1000)


def utc_iso(t=None):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t if t is not None else time.time()))


class quiet(object):
    """Context manager that hides console output from C/C++ libraries (radio drivers print a lot).
    Ctrl-C is held back while output is redirected, so the terminal can never be left silenced."""

    def __enter__(self):
        import signal
        try:
            signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT})
            self._masked = True
        except (AttributeError, ValueError):
            self._masked = False
        sys.stdout.flush()
        sys.stderr.flush()
        self._null = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
        self._save = (os.dup(1), os.dup(2))
        os.dup2(self._null[0], 1)
        os.dup2(self._null[1], 2)
        return self

    def __exit__(self, *exc):
        import signal
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(self._save[0], 1)
        os.dup2(self._save[1], 2)
        for fd in list(self._null) + list(self._save):
            os.close(fd)
        if self._masked:
            signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGINT})
        return False


def stop_on_first_ctrl_c():
    """First Ctrl-C stops the tool cleanly; further Ctrl-C (e.g. forwarded by `timeout`) are ignored so the
    summary is always printed."""
    import signal

    def handler(signum, frame):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        raise KeyboardInterrupt
    signal.signal(signal.SIGINT, handler)


def ignore_ctrl_c():
    import signal
    try:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    except (KeyboardInterrupt, ValueError):
        pass


def stop_leftovers(script_names=(), program_names=()):
    """Stop processes left from earlier runs: python scripts with one of script_names in their command line, and
    programs named in program_names (e.g. srsRAN's pdsch_ue left running after its tool was killed - they keep
    the CPU at 100% and slow every new run down). This process and its parents are never touched.
    Returns the number stopped."""
    import signal
    keep, pid = set(), os.getpid()
    while pid > 1:                                   # this process and its ancestors
        keep.add(pid)
        try:
            pid = int(open("/proc/%d/stat" % pid).read().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            break
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit() or int(name) in keep:
            continue
        try:
            cmd = open("/proc/%s/cmdline" % name, "rb").read().replace(b"\x00", b" ").decode(errors="replace")
        except OSError:
            continue
        if not cmd:
            continue
        exe = os.path.basename(cmd.split()[0])
        if exe in program_names or ("python" in exe and any(t in cmd for t in script_names)):
            found.append(int(name))
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for p in found:
            try:
                os.kill(p, sig)
            except OSError:
                pass
        if found and sig == signal.SIGTERM:
            time.sleep(1.0)
    return len(found)


# Heavy processing leaves one CPU core free and runs at lower priority, so the Pi (and an SSH session to it)
# stays responsive while a tool is working.
WORKERS = max(1, (os.cpu_count() or 4) - 1)


def be_nice():
    """Lower this tool's priority (srsRAN and other child programs inherit it)."""
    try:
        os.nice(10)
    except (AttributeError, OSError):
        pass


def rate_colour(pct):
    return "green" if pct >= 95 else ("yellow" if pct >= 70 else "red")
