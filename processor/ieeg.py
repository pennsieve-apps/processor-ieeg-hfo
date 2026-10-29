"""
Finding and reading the recordings of a BIDS-iEEG derivative (the `ieeg` port of
ieeg-prep), a plain BIDS-iEEG dataset, or loose EDF/BDF/BrainVision/FIF/EEGLAB files.

The contract is docs/conventions.md of workflow-epilepsy-ieeg:

- `recordings.csv` (anywhere under INPUT_DIR) lists the recordings, with `path`
  relative to the folder the CSV is in, plus sfreq, onset_s and the baseline interval;
- `<stem>[_desc-preproc]_channels.tsv` next to the signal file says which channels are
  good, their type and labels;
- `<stem>[_desc-preproc]_events.tsv` holds `seizure`, `baseline` and `artifact` rows.

Signal files that no recordings.csv lists are picked up as well (a plain BIDS dataset or
loose files), with whatever sidecars they have. Nothing here loads signal data; the
callers read one time window at a time.

This module is shared, unchanged, by processor-ieeg-spikes and processor-ieeg-hfo.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

SIGNAL_EXTENSIONS = (".edf", ".bdf", ".vhdr", ".fif", ".set")
# Channel types analysed (BIDS `type`, upper case). Everything else (ECG, EMG, EOG, MISC,
# TRIG, REF, ...) is ignored.
NEURAL_TYPES = {"SEEG", "ECOG", "EEG", "DBS", "IEEG"}
# Used only when there is no _channels.tsv: names that are clearly not brain signals.
NON_NEURAL_NAME = re.compile(r"^(ECG|EKG|EMG|EOG|RESP|PHOTIC|EVENT|STATUS|TRIG|DC\d*|SPO2|PULSE|PR|OSAT|ANNOT)",
                             re.IGNORECASE)
SEIZURE_TYPES = {"seizure", "ictal", "sz"}


def walk_files(root: Path):
    """Every file under root, following symlinked folders (INPUT_DIR entries are symlinks), once."""
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for name in sorted(filenames):
            if not name.startswith("."):
                yield Path(dirpath) / name


def _float(value) -> float | None:
    try:
        v = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return [{(k or "").strip(): (v or "").strip() for k, v in row.items()}
                for row in csv.DictReader(fh, delimiter="\t")]


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return [{(k or "").strip(): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(fh)]


@dataclass
class Recording:
    recording_id: str
    path: Path                         # the signal file (.edf, .vhdr, ...)
    rel: str                           # path relative to INPUT_DIR, for logs and provenance
    source: str                        # "recordings.csv" or "scan"
    sfreq: float | None = None         # from recordings.csv (the file is authoritative)
    kind: str = ""                     # ictal | interictal | unknown
    onset_s: float | None = None       # seizure onset from recordings.csv
    baseline: tuple[float, float] | None = None
    channels_tsv: Path | None = None
    events_tsv: Path | None = None
    sidecar_json: Path | None = None
    notes: list[str] = field(default_factory=list)


def signal_prefix(path: Path) -> str:
    """`sub-01_task-x_desc-preproc_ieeg.edf` -> `sub-01_task-x_desc-preproc`; `foo.edf` -> `foo`."""
    stem = path.name[: -len(path.suffix)]
    return stem[: -len("_ieeg")] if stem.endswith("_ieeg") else stem


def recording_id_from_prefix(prefix: str) -> str:
    """The BIDS entity chain of the source recording: drop a `_desc-<label>` entity."""
    return re.sub(r"_desc-[A-Za-z0-9]+", "", prefix)


def _sidecars(rec: Recording) -> None:
    prefix = signal_prefix(rec.path)
    folder = rec.path.parent
    for attr, suffix in (("channels_tsv", "_channels.tsv"), ("events_tsv", "_events.tsv"),
                         ("sidecar_json", "_ieeg.json")):
        if getattr(rec, attr) is None:
            p = folder / f"{prefix}{suffix}"
            if p.is_file():
                setattr(rec, attr, p)


def find_recordings(input_dir: Path) -> tuple[list[Recording], list[str], list[dict]]:
    """
    (recordings, every input file seen (relative), problems). A problem is a listed
    recording whose file is missing, or a duplicate recording id; each is a dict with
    recording_id, input and reason.
    """
    files = list(walk_files(input_dir))
    rel = lambda p: os.path.relpath(p, input_dir)  # noqa: E731
    seen = [rel(p) for p in files]
    problems: list[dict] = []
    recs: list[Recording] = []
    listed: set[str] = set()

    for csv_path in (p for p in files if p.name == "recordings.csv"):
        try:
            rows = read_csv(csv_path)
        except (OSError, UnicodeDecodeError, csv.Error) as exc:
            problems.append({"recording_id": "", "input": rel(csv_path), "reason": f"unreadable recordings.csv: {exc}"})
            continue
        for row in rows:
            rid, path = row.get("recording_id", ""), row.get("path", "")
            if not rid or not path:
                continue
            sig = csv_path.parent / path
            if not sig.is_file():
                problems.append({"recording_id": rid, "input": rel(csv_path),
                                 "reason": f"recordings.csv lists {path}, which is not in the input"})
                continue
            listed.add(os.path.realpath(sig))
            bs, be = _float(row.get("baseline_start_s")), _float(row.get("baseline_end_s"))
            recs.append(Recording(
                recording_id=rid, path=sig, rel=rel(sig), source="recordings.csv",
                sfreq=_float(row.get("sfreq")), kind=row.get("kind", ""), onset_s=_float(row.get("onset_s")),
                baseline=(bs, be) if bs is not None and be is not None and be > bs else None))

    for p in files:
        if p.suffix.lower() not in SIGNAL_EXTENSIONS or os.path.realpath(p) in listed:
            continue
        prefix = signal_prefix(p)
        recs.append(Recording(recording_id=recording_id_from_prefix(prefix), path=p, rel=rel(p), source="scan"))

    for r in recs:
        _sidecars(r)

    # one recording per id: recordings.csv entries win over scanned files (a derivative
    # next to its source dataset), otherwise the first one found
    out: list[Recording] = []
    by_id: dict[str, Recording] = {}
    for r in sorted(recs, key=lambda r: (r.source != "recordings.csv", r.rel)):
        if r.recording_id in by_id:
            keep = by_id[r.recording_id]
            problems.append({"recording_id": r.recording_id, "input": r.rel, "status": "skipped",
                             "reason": f"duplicate recording id; {keep.rel} is used"})
            continue
        by_id[r.recording_id] = r
        out.append(r)
    out.sort(key=lambda r: r.recording_id)
    return out, seen, problems


# --- channels -------------------------------------------------------------------

@dataclass
class Channel:
    name: str
    type: str          # upper-case BIDS type, or MNE's when there is no channels.tsv
    status: str        # good | bad
    reason: str = ""   # why it is not analysed


def channels(rec: Recording, raw_names: list[str], raw_types: list[str]) -> list[Channel]:
    """
    Every channel of the recording with good/bad status. With a channels.tsv, its type and
    status columns decide (a channel absent from the file is bad); without one, MNE's
    channel type and a name filter decide.
    """
    in_file = set(raw_names)
    out: list[Channel] = []
    if rec.channels_tsv is not None:
        for row in read_tsv(rec.channels_tsv):
            name = row.get("name", "")
            if not name:
                continue
            ctype = (row.get("type") or "").upper()
            status = (row.get("status") or "good").lower()
            reason = ""
            if name not in in_file:
                status, reason = "bad", "listed in channels.tsv but not in the signal file"
            elif status == "bad":
                reason = row.get("status_description") or "status bad in channels.tsv"
                if reason.lower() == "n/a":
                    reason = "status bad in channels.tsv"
            elif ctype and ctype not in NEURAL_TYPES:
                status, reason = "bad", f"type {ctype} is not analysed"
            out.append(Channel(name, ctype or "EEG", "bad" if reason else "good", reason))
        return out
    for name, mtype in zip(raw_names, raw_types):
        ctype = mtype.upper()
        if ctype not in NEURAL_TYPES:
            out.append(Channel(name, ctype, "bad", f"type {ctype} is not analysed"))
        elif NON_NEURAL_NAME.match(name):
            out.append(Channel(name, ctype, "bad", "name looks like a non-brain channel"))
        else:
            out.append(Channel(name, ctype, "good"))
    return out


# --- time intervals -------------------------------------------------------------

def events(rec: Recording) -> list[dict]:
    if rec.events_tsv is None:
        return []
    try:
        return read_tsv(rec.events_tsv)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        rec.notes.append(f"events.tsv unreadable: {exc}")
        return []


def seizure_intervals(rec: Recording, evs: list[dict], duration_s: float) -> list[tuple[float, float, bool]]:
    """(start, end, end_known) of every known seizure; end = recording end when its duration is unknown."""
    out = []
    for e in evs:
        if (e.get("trial_type") or "").strip().lower() in SEIZURE_TYPES:
            on, dur = _float(e.get("onset")), _float(e.get("duration"))
            if on is None:
                continue
            known = dur is not None and dur > 0
            out.append((on, on + dur if known else duration_s, known))
    if rec.onset_s is not None and not any(abs(s - rec.onset_s) < 1.0 for s, _, _ in out):
        out.append((rec.onset_s, duration_s, False))
    return sorted(out)


def baseline_interval(rec: Recording, evs: list[dict]) -> tuple[float, float] | None:
    if rec.baseline is not None:
        return rec.baseline
    for e in evs:
        if (e.get("trial_type") or "").strip().lower() == "baseline":
            on, dur = _float(e.get("onset")), _float(e.get("duration"))
            if on is not None and dur is not None and dur > 0:
                return on, on + dur
    return None


def artifact_intervals(evs: list[dict]) -> list[tuple[float, float]]:
    out = []
    for e in evs:
        if (e.get("trial_type") or "").strip().lower() == "artifact":
            on, dur = _float(e.get("onset")), _float(e.get("duration"))
            if on is not None and dur is not None and dur > 0:
                out.append((on, on + dur))
    return out


def subtract(intervals: list[tuple[float, float]], cut: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """intervals minus the union of cut."""
    out = list(intervals)
    for c0, c1 in cut:
        nxt = []
        for a, b in out:
            if c1 <= a or c0 >= b:
                nxt.append((a, b))
                continue
            if c0 > a:
                nxt.append((a, c0))
            if c1 < b:
                nxt.append((c1, b))
        out = nxt
    return [(a, b) for a, b in out if b > a]


@dataclass
class Plan:
    segments: list[tuple[float, float]]   # seconds from the recording start, analysed
    description: str                      # for the log and provenance
    excluded: list[dict]                  # what was cut and why


def analysis_plan(rec: Recording, duration_s: float, segment: str, exclude_seizure: bool,
                  seizure_margin_s: float) -> Plan:
    """
    Which parts of the recording to analyse.

    segment=all: the whole recording; segment=baseline: only the baseline interval
    (recordings without one are analysed whole, with a note). Then, when exclude_seizure,
    every known seizure is cut out with seizure_margin_s on both sides (to the end of the
    recording when the seizure's end is unknown), and annotated artifacts always are.
    """
    evs = events(rec)
    excluded: list[dict] = []
    whole = [(0.0, duration_s)]
    desc = "whole recording"
    if segment == "baseline":
        bl = baseline_interval(rec, evs)
        if bl is None:
            rec.notes.append("SEGMENT=baseline but the recording has no baseline interval; analysed it whole")
        else:
            whole = [(max(0.0, bl[0]), min(duration_s, bl[1]))]
            desc = f"baseline {whole[0][0]:.1f}-{whole[0][1]:.1f} s"
    cuts: list[tuple[float, float]] = []
    if exclude_seizure:
        for s, e, known in seizure_intervals(rec, evs, duration_s):
            a, b = max(0.0, s - seizure_margin_s), min(duration_s, e + seizure_margin_s if known else duration_s)
            cuts.append((a, b))
            excluded.append({"start_s": round(a, 3), "end_s": round(b, 3),
                             "reason": "seizure" + ("" if known else " (end unknown: cut to the end of the recording)")})
    for a, b in artifact_intervals(evs):
        cuts.append((a, b))
        excluded.append({"start_s": round(a, 3), "end_s": round(b, 3), "reason": "artifact"})
    segs = subtract(whole, cuts)
    if excluded:
        desc += ", minus " + ", ".join(f"{x['reason'].split(' ')[0]} {x['start_s']:.1f}-{x['end_s']:.1f} s"
                                       for x in excluded)
    return Plan(segs, desc, excluded)


def chunks(segments: list[tuple[float, float]], window_s: float, pad_s: float, min_s: float,
           edge_s: float = 0.0):
    """
    Split segments into analysis windows. Yields (core_start, core_end, read_start, read_end)
    in seconds: detections are kept when they fall in the core; the read interval adds
    pad_s on each side where the segment allows it (never across excluded time). edge_s is
    trimmed from both ends of each segment (a detector that cannot see its first and last
    second). Segments shorter than min_s are skipped (yielded with read_start = None so
    the caller can log them).
    """
    for s0, s1 in segments:
        c0, c1 = s0 + edge_s, s1 - edge_s
        if s1 - s0 < min_s or c1 <= c0:
            yield s0, s1, None, None
            continue
        n = max(1, math.ceil((c1 - c0) / window_s - 1e-9))
        step = (c1 - c0) / n
        for i in range(n):
            a, b = c0 + i * step, c0 + (i + 1) * step
            ra, rb = max(s0, a - pad_s), min(s1, b + pad_s)
            yield a, b, ra, rb


def sidecar_line_freq(rec: Recording) -> float | None:
    if rec.sidecar_json is None:
        return None
    try:
        return _float(json.loads(rec.sidecar_json.read_text()).get("PowerLineFrequency"))
    except (OSError, ValueError, AttributeError):
        return None


# --- signal ---------------------------------------------------------------------

def open_raw(path: Path):
    """The recording, lazily (nothing loaded). EDF/BDF/BrainVision/FIF/EEGLAB through MNE."""
    import mne

    return mne.io.read_raw(str(path), preload=False, verbose="error")


def read_window(raw, picks: list[str], start: int, stop: int, max_bytes: int = 256 << 20):
    """
    Samples [start, stop) of the picked channels as float32 microvolts, shape (channels,
    samples). Read in groups of channels so the float64 buffer MNE returns stays under
    max_bytes.
    """
    import numpy as np

    n = stop - start
    out = np.empty((len(picks), n), dtype=np.float32)
    group = max(1, int(max_bytes // max(1, n * 8)))
    for i in range(0, len(picks), group):
        block = raw.get_data(picks=picks[i:i + group], start=start, stop=stop)
        out[i:i + group] = block * 1e6  # MNE returns volts
        del block
    return out
