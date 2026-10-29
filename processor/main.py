"""
processor-ieeg-hfo: high-frequency oscillation (ripple, fast ripple) detection on iEEG.

Reads the BIDS-iEEG derivative of ieeg-prep (or a plain BIDS-iEEG dataset, or loose
EDF/BrainVision/FIF files) from INPUT_DIR and writes to OUTPUT_DIR (or OUTPUT_DIR/<LABEL>/):

  hfo.csv          recording_id, channel, start_s, end_s, band, detector, amplitude_uv
  hfo_rates.csv    recording_id, channel, band, n_hfo, minutes, rate_per_min, detector
  provenance.json

A recording sampled below MIN_SFREQ (default 1000 Hz) is not analysed: it is logged and
listed in provenance.json with status "skipped: sampling rate X Hz < 1000", and the run
still succeeds (empty tables when every recording is skipped). One recording at a time, one
time window (WINDOW_S) of the good channels at a time, channels in parallel (THREADS).
"""
from __future__ import annotations

import csv
import json
import logging
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

from processor import detect as D
from processor import ieeg as I

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-ieeg-hfo")
log.setLevel(logging.INFO)
for noisy in ("mne", "numba", "joblib"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-ieeg-hfo"
REPO_ROOT = Path(__file__).resolve().parent.parent
PAD_S = 2.0          # context read on each side of a window (never across excluded time)
MIN_SECONDS = 10.0   # shorter analysable segments are not analysed (too little for the threshold statistics)

HFO_HEADER = ["recording_id", "channel", "start_s", "end_s", "band", "detector", "amplitude_uv"]
RATES_HEADER = ["recording_id", "channel", "band", "n_hfo", "minutes", "rate_per_min", "detector"]


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for pkg in ("epycom", "mne", "numpy", "scipy", "joblib", "numba"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def cpu_count() -> int:
    try:
        from joblib import cpu_count as jc

        return max(1, jc(only_physical_cores=False))
    except Exception:  # noqa: BLE001
        return max(1, os.cpu_count() or 1)


def _bool(name: str, value: str) -> bool:
    v = value.strip().lower()
    if v in ("true", "1", "yes"):
        return True
    if v in ("false", "0", "no"):
        return False
    raise ValueError(f"{name} must be true or false (got {value!r})")


def env_params() -> dict:
    """Parameters from the environment; raises ValueError on a bad value."""
    e = lambda k, d: (os.environ.get(k) or d).strip()  # noqa: E731
    raw = {
        "DETECTOR": e("DETECTOR", "linelength").lower(),
        "BANDS": e("BANDS", "ripple,fast_ripple").lower(),
        "MIN_SFREQ": e("MIN_SFREQ", "1000"),
        "FAST_RIPPLE_MIN_SFREQ": e("FAST_RIPPLE_MIN_SFREQ", "2000"),
        "THRESHOLD": e("THRESHOLD", "auto").lower(),
        "WINDOW_MS": e("WINDOW_MS", "auto").lower(),
        "MIN_DURATION_MS": e("MIN_DURATION_MS", "6"),
        "SEGMENT": e("SEGMENT", "all").lower(),
        "EXCLUDE_SEIZURE": e("EXCLUDE_SEIZURE", "true"),
        "SEIZURE_MARGIN_S": e("SEIZURE_MARGIN_S", "10"),
        "WINDOW_S": e("WINDOW_S", "300"),
        "THREADS": e("THREADS", "0"),
        "LABEL": e("LABEL", ""),
    }
    p = dict(raw)
    det = p["DETECTOR"]
    if det not in D.DEFAULTS:
        raise ValueError(f"DETECTOR must be one of {', '.join(D.DEFAULTS)} (got {raw['DETECTOR']!r})")
    bands = [b.strip().replace("-", "_") for b in raw["BANDS"].split(",") if b.strip()]
    bad = [b for b in bands if b not in D.BANDS]
    if not bands or bad:
        raise ValueError(f"BANDS must be a comma-separated list of {', '.join(D.BANDS)} (got {raw['BANDS']!r})")
    p["bands"] = [b for b in D.BANDS if b in bands]
    p["min_sfreq"] = float(raw["MIN_SFREQ"])
    p["fr_min_sfreq"] = float(raw["FAST_RIPPLE_MIN_SFREQ"])
    if p["min_sfreq"] < 500 or p["fr_min_sfreq"] < 1000:
        raise ValueError("MIN_SFREQ must be >= 500 and FAST_RIPPLE_MIN_SFREQ >= 1000 (Hz)")
    p["threshold"] = D.DEFAULTS[det]["threshold"] if raw["THRESHOLD"] == "auto" else float(raw["THRESHOLD"])
    if p["threshold"] <= 0:
        raise ValueError("THRESHOLD must be > 0 or auto")
    if det == "cs":
        p["window_ms"] = None
    else:
        p["window_ms"] = D.DEFAULTS[det]["window_ms"] if raw["WINDOW_MS"] == "auto" else float(raw["WINDOW_MS"])
        if not 1 <= p["window_ms"] <= 100:
            raise ValueError("WINDOW_MS must be auto or 1-100")
    p["min_duration"] = float(raw["MIN_DURATION_MS"]) / 1000.0
    if p["min_duration"] < 0:
        raise ValueError("MIN_DURATION_MS must be >= 0")
    if p["SEGMENT"] not in ("all", "baseline"):
        raise ValueError(f"SEGMENT must be all or baseline (got {raw['SEGMENT']!r})")
    p["exclude_seizure"] = _bool("EXCLUDE_SEIZURE", raw["EXCLUDE_SEIZURE"])
    p["margin"] = float(raw["SEIZURE_MARGIN_S"])
    if p["margin"] < 0:
        raise ValueError("SEIZURE_MARGIN_S must be >= 0")
    p["window"] = float(raw["WINDOW_S"])
    if p["window"] < 30:
        raise ValueError("WINDOW_S must be >= 30")
    threads = int(raw["THREADS"])
    if threads < 0:
        raise ValueError("THREADS must be >= 0 (0 = every vCPU)")
    p["threads"] = threads or cpu_count()
    if p["LABEL"] and not all(c.isalnum() or c in "-_" for c in p["LABEL"]):
        raise ValueError("LABEL may hold letters, digits, - and _ only")
    return p


class CsvOut:
    def __init__(self, path: Path, header: list[str]):
        self.fh = open(path, "w", newline="", encoding="utf-8")
        self.w = csv.DictWriter(self.fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        self.w.writeheader()

    def rows(self, rows):
        self.w.writerows(rows)
        self.fh.flush()

    def close(self):
        self.fh.close()


class Skipped(Exception):
    """A designed reason not to analyse a recording (not a failure)."""


def process_recording(rec: I.Recording, params: dict, parallel) -> dict:
    """Detect HFOs in one recording. Returns hfo rows, rate rows and a status summary."""
    from joblib import delayed

    t0 = time.monotonic()
    raw = I.open_raw(rec.path)
    fs = float(raw.info["sfreq"])
    if fs < params["min_sfreq"]:
        raise Skipped(f"skipped: sampling rate {fs:g} Hz < {params['min_sfreq']:g}")
    n_times = int(raw.n_times)
    duration = n_times / fs
    bands, band_notes = D.bands_for(fs, params["bands"], params["fr_min_sfreq"])
    rec.notes.extend(band_notes)
    if not bands:
        raise Skipped("skipped: no requested band can be analysed (" + "; ".join(band_notes) + ")")
    chans = I.channels(rec, list(raw.ch_names), list(raw.get_channel_types()))
    good = [c for c in chans if c.status == "good"]
    if not good:
        raise ValueError(f"no good iEEG channel among {len(chans)}")
    index = {name: i for i, name in enumerate(raw.ch_names)}
    picks = [index[c.name] for c in good]
    plan = I.analysis_plan(rec, duration, params["SEGMENT"], params["exclude_seizure"], params["margin"])

    minutes = {c.name: 0.0 for c in good}
    found: dict[str, list] = {c.name: [] for c in good}
    errors: dict[str, str] = {}
    windows = 0
    t_read = t_detect = 0.0
    det = params["DETECTOR"]
    for c0, c1, r0, r1 in I.chunks(plan.segments, params["window"], PAD_S, MIN_SECONDS):
        if r0 is None:
            rec.notes.append(f"segment {c0:.1f}-{c1:.1f} s is shorter than {MIN_SECONDS:g} s; not analysed")
            continue
        s, e = max(0, int(round(r0 * fs))), min(n_times, int(round(r1 * fs)))
        ta = time.monotonic()
        data = I.read_window(raw, picks, s, e)
        tb = time.monotonic()
        results = parallel(delayed(D.detect_channel)(data[i], fs, det, bands, params["window_ms"],
                                                     params["threshold"], params["min_duration"], s / fs, (c0, c1))
                           for i in range(len(good)))
        del data
        t_read += tb - ta
        t_detect += time.monotonic() - tb
        windows += 1
        for c, (events, err) in zip(good, results):
            if err:
                errors.setdefault(c.name, f"{c0:.0f}-{c1:.0f} s: {err}")
                continue
            minutes[c.name] += (c1 - c0) / 60.0
            found[c.name].extend(events)

    analysed = [c for c in good if minutes[c.name] > 0]
    hfo = []
    for c in good:
        for t0_, t1_, band, amp in sorted(found[c.name]):
            hfo.append({"recording_id": rec.recording_id, "channel": c.name, "start_s": f"{t0_:.4f}",
                        "end_s": f"{t1_:.4f}", "band": band, "detector": det, "amplitude_uv": f"{amp:.2f}",
                        "_t": t0_})
    hfo.sort(key=lambda r: (r["_t"], r["channel"], r["band"]))
    rates = []
    counts = {b: 0 for b in bands}
    for c in analysed:
        m = minutes[c.name]
        for b in bands:
            n = sum(1 for ev in found[c.name] if ev[2] == b)
            counts[b] += n
            rates.append({"recording_id": rec.recording_id, "channel": c.name, "band": b, "n_hfo": n,
                          "minutes": f"{m:.3f}", "rate_per_min": f"{n / m:.4f}", "detector": det})
    summary = {
        "sfreq": fs, "duration_s": round(duration, 3), "bands": bands, "n_channels": len(chans), "n_good": len(good),
        "n_analysed": len(analysed), "analysed": plan.description,
        "segments": [[round(a, 3), round(b, 3)] for a, b in plan.segments], "excluded": plan.excluded,
        "minutes": round(max(minutes.values()) if minutes else 0.0, 3), "windows": windows,
        "n_hfo": counts, "channels_not_analysed": {c.name: c.reason for c in chans if c.status != "good"},
        "channel_errors": errors, "notes": rec.notes,
        "seconds": {"read": round(t_read, 2), "detect": round(t_detect, 2), "total": round(time.monotonic() - t0, 2)},
    }
    return {"hfo": hfo, "rates": rates, "summary": summary}


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    try:
        params = env_params()
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    out_dir = output_dir / params["LABEL"] if params["LABEL"] else output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    recs, seen, problems = I.find_recordings(input_dir)
    log.info("input: %d recording(s), %d file(s); detector %s (threshold %g%s), bands %s, segment %s, seizure %s, "
             "%d thread(s)", len(recs), len(seen), params["DETECTOR"], params["threshold"],
             f", window {params['window_ms']:g} ms" if params["window_ms"] else "", ",".join(params["bands"]),
             params["SEGMENT"],
             f"excluded (margin {params['margin']:g} s)" if params["exclude_seizure"] else "kept", params["threads"])
    for p in problems:
        log.warning("%s: %s", p["input"], p["reason"])

    hfo_out = CsvOut(out_dir / "hfo.csv", HFO_HEADER)
    rates_out = CsvOut(out_dir / "hfo_rates.csv", RATES_HEADER)
    prov_recs: list[dict] = [{"recording_id": p["recording_id"], "input": p["input"],
                              "status": p.get("status", "failed"), "error": p["reason"]} for p in problems]
    n_ok = n_skipped = 0
    n_failed = sum(1 for p in problems if p.get("status", "failed") == "failed")

    from joblib import Parallel

    with Parallel(n_jobs=params["threads"], backend="loky", max_nbytes=None) as parallel:
        for rec in recs:
            try:
                res = process_recording(rec, params, parallel)
            except Skipped as why:
                log.warning("%s: %s", rec.recording_id, why)
                prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": str(why),
                                  "notes": rec.notes})
                n_skipped += 1
                continue
            except Exception as exc:  # noqa: BLE001 - one recording's failure must not sink the run
                log.error("%s: FAILED: %s\n%s", rec.recording_id, exc, traceback.format_exc())
                prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": "failed",
                                  "error": f"{type(exc).__name__}: {exc}"})
                n_failed += 1
                continue
            s = res["summary"]
            if s["n_analysed"] == 0:
                status = "skipped: no analysable time" + (
                    f" ({'; '.join(s['notes'])})" if s["notes"] else f" ({s['analysed']})")
                log.warning("%s: %s", rec.recording_id, status)
                prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": status, **s})
                n_skipped += 1
                continue
            hfo_out.rows(res["hfo"])
            rates_out.rows(res["rates"])
            n_ok += 1
            for ch, err in s["channel_errors"].items():
                log.warning("%s: channel %s: %s", rec.recording_id, ch, err)
            for note in s["notes"]:
                log.info("%s: %s", rec.recording_id, note)
            per = ", ".join(f"{n} {b} ({n / max(s['minutes'], 1e-9) / max(s['n_analysed'], 1):.2f}/min/channel)"
                            for b, n in s["n_hfo"].items())
            log.info("%s: %d channels (%d good, %d analysed), %.1f Hz, %.2f min analysed (%s), %s, %.1f s "
                     "(read %.1f s, detect %.1f s)", rec.recording_id, s["n_channels"], s["n_good"], s["n_analysed"],
                     s["sfreq"], s["minutes"], s["analysed"], per, s["seconds"]["total"], s["seconds"]["read"],
                     s["seconds"]["detect"])
            prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": "ok", **s})
    hfo_out.close()
    rates_out.close()

    status = "ok"
    if not recs and not problems:
        log.error("no recording found under INPUT_DIR (looked for recordings.csv and %s files); it holds %d file(s): %s",
                  "/".join(I.SIGNAL_EXTENSIONS), len(seen), ", ".join(seen[:40]) + (" …" if len(seen) > 40 else ""))
        status = "failed: no input"
    elif n_ok == 0 and n_skipped == 0 and n_failed > 0:
        status = "failed: every recording failed"
    elif n_ok == 0:
        status = "ok: nothing analysed (" + (f"{n_skipped} skipped" + (f", {n_failed} failed" if n_failed else "")) + ")"
    prov = {
        "app": APP_NAME, "version": version_of(REPO_ROOT), "tools": tool_versions(),
        "parameters": {k: params[k] for k in ("DETECTOR", "BANDS", "MIN_SFREQ", "FAST_RIPPLE_MIN_SFREQ", "THRESHOLD",
                                              "WINDOW_MS", "MIN_DURATION_MS", "SEGMENT", "EXCLUDE_SEIZURE",
                                              "SEIZURE_MARGIN_S", "WINDOW_S", "THREADS", "LABEL")},
        "resolved": {"threshold": params["threshold"], "window_ms": params["window_ms"],
                     "threads": params["threads"]},
        "inputs": seen, "status": status,
        "started": started.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2),
        "recordings": prov_recs,
    }
    with open(out_dir / "provenance.json", "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2)
    log.info("done: %d recording(s) ok, %d skipped, %d failed, %.1f s -> %s", n_ok, n_skipped, n_failed,
             prov["seconds"], out_dir)
    return 1 if status.startswith("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
