"""
A tiny synthetic BIDS-iEEG derivative that follows workflow-epilepsy-ieeg/docs/conventions.md
(the `ieeg` output of ieeg-prep), for the tests and `make run`:

  <out>/ieeg/dataset_description.json
  <out>/ieeg/recordings.csv
  <out>/ieeg/sub-01/ieeg/sub-01_task-rest_run-01_desc-preproc_{ieeg.edf,ieeg.json,channels.tsv,events.tsv}   2048 Hz
  <out>/ieeg/sub-02/ieeg/sub-02_task-ictal_run-01_desc-preproc_{...}                                        1024 Hz
  <out>/ieeg/sub-03/ieeg/sub-03_task-rest_run-01_desc-preproc_{...}                                         512 Hz

sub-01 (2048 Hz, 90 s): LA01-LA02 carries ripples (120 Hz, 80 ms, 80 µV) every ~3 s (~20/min),
LA02-LA03 fast ripples (350 Hz, 40 ms, 50 µV) every ~4 s (~15/min), LA03-LA04 background only,
LB01-LB02 sharp spikes without oscillations (a classic source of false HFOs after
band-pass filtering), LB02-LB03 is marked bad.
sub-02 (1024 Hz, 120 s, a seizure at 90 s, baseline 0-80 s): ripples on RA01-RA02 before
the seizure, and bursts of both kinds on every channel during it (EXCLUDE_SEIZURE must keep
them out). Fast ripples are not analysed at 1024 Hz by default.
sub-03 (512 Hz): below MIN_SFREQ, skipped.

Usage: python -m tests.make_fixture <out-dir>
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

SEED = 42


def pink(n: int, rng: np.random.Generator, std: float) -> np.ndarray:
    """1/f background, `std` µV."""
    spec = np.fft.rfft(rng.standard_normal(n))
    f = np.fft.rfftfreq(n)
    f[0] = f[1]
    x = np.fft.irfft(spec / np.sqrt(f), n)
    return (x / x.std() * std).astype(np.float64)


def burst(fs: float, freq: float, dur: float, amp: float) -> np.ndarray:
    """A Hann-windowed oscillation of `dur` seconds at `freq` Hz, peak `amp` µV."""
    t = np.arange(int(dur * fs)) / fs
    return amp * np.hanning(len(t)) * np.sin(2 * np.pi * freq * t)


def sharp(fs: float, amp: float) -> np.ndarray:
    t = np.arange(int(0.1 * fs)) / fs
    return amp * np.exp(-np.abs(t - 0.03) / 0.004)


def times(duration: float, rate_per_min: float, rng, start: float = 1.0, stop: float | None = None) -> np.ndarray:
    stop = duration - 1.0 if stop is None else stop
    step = 60.0 / rate_per_min
    t = np.arange(start, stop, step)
    return t + rng.uniform(-0.2 * step, 0.2 * step, size=len(t))


def add_at(x: np.ndarray, fs: float, ts, shape: np.ndarray) -> None:
    for t in ts:
        i = int(round(t * fs))
        j = min(len(x), i + len(shape))
        if 0 <= i < len(x):
            x[i:j] += shape[: j - i]


def write_edf(path: Path, data: np.ndarray, names: list[str], fs: float) -> None:
    from edfio import Edf, EdfSignal

    sigs = []
    for x, name in zip(data, names):
        r = max(1.0, float(np.max(np.abs(x))) * 1.1)
        sigs.append(EdfSignal(x, sampling_frequency=fs, label=name, physical_dimension="uV",
                              physical_range=(-r, r)))
    Edf(sigs).write(path)


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


CHANNEL_COLS = ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status",
                "status_description", "contacts", "soz", "resected", "onset", "spread"]
REC_COLS = ["recording_id", "subject", "session", "task", "run", "path", "sfreq", "n_channels", "n_good",
            "duration_s", "montage", "kind", "onset_s", "baseline_start_s", "baseline_end_s", "source_format"]


def write_recording(root: Path, sub: str, task: str, fs: float, data: np.ndarray, names: list[str],
                    bad=(), events=(), onset="", baseline=("", "")) -> dict:
    stem = f"sub-{sub}_task-{task}_run-01"
    d = root / f"sub-{sub}" / "ieeg"
    d.mkdir(parents=True, exist_ok=True)
    write_edf(d / f"{stem}_desc-preproc_ieeg.edf", data, names, fs)
    write_tsv(d / f"{stem}_desc-preproc_channels.tsv", CHANNEL_COLS, [
        {"name": n, "type": "SEEG", "units": "µV", "low_cutoff": "0.5", "high_cutoff": f"{fs / 2:g}",
         "sampling_frequency": f"{fs:g}", "status": "bad" if n in bad else "good",
         "status_description": "noisy" if n in bad else "n/a", "contacts": n.replace("-", "|"), "soz": "",
         "resected": "", "onset": "", "spread": ""} for n in names])
    write_tsv(d / f"{stem}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], list(events))
    (d / f"{stem}_desc-preproc_ieeg.json").write_text(json.dumps({"SamplingFrequency": fs, "PowerLineFrequency": 60}))
    return {"recording_id": stem, "subject": sub, "session": "", "task": task, "run": "01",
            "path": f"sub-{sub}/ieeg/{stem}_desc-preproc_ieeg.edf", "sfreq": f"{fs:g}", "n_channels": len(names),
            "n_good": len(names) - len(bad), "duration_s": f"{data.shape[1] / fs:g}", "montage": "bipolar",
            "kind": "ictal" if onset else "interictal", "onset_s": onset, "baseline_start_s": baseline[0],
            "baseline_end_s": baseline[1], "source_format": "bids-edf"}


def make(out: Path) -> Path:
    rng = np.random.default_rng(SEED)
    root = out / "ieeg"
    root.mkdir(parents=True, exist_ok=True)
    (root / "dataset_description.json").write_text(json.dumps({
        "Name": "synthetic fixture", "BIDSVersion": "1.9.0", "DatasetType": "derivative",
        "GeneratedBy": [{"Name": "ieeg-prep", "Description": "tests/make_fixture.py of processor-ieeg-hfo"}]},
        indent=2))
    truth: dict = {}
    recs = []

    fs, dur = 2048.0, 90.0
    names = ["LA01-LA02", "LA02-LA03", "LA03-LA04", "LB01-LB02", "LB02-LB03"]
    data = np.vstack([pink(int(fs * dur), rng, 30.0) for _ in names])
    rip, fr, sp = times(dur, 20, rng), times(dur, 15, rng), times(dur, 20, rng)
    add_at(data[0], fs, rip, burst(fs, 120.0, 0.08, 80.0))
    add_at(data[1], fs, fr, burst(fs, 350.0, 0.04, 50.0))
    add_at(data[3], fs, sp, sharp(fs, 400.0))
    data[4] += 200 * rng.standard_normal(data.shape[1])
    recs.append(write_recording(root, "01", "rest", fs, data, names, bad=("LB02-LB03",)))
    truth[recs[-1]["recording_id"]] = {"ripple": rip.tolist(), "fast_ripple": fr.tolist()}

    fs, dur = 1024.0, 120.0
    names = ["RA01-RA02", "RA02-RA03", "RA03-RA04"]
    data = np.vstack([pink(int(fs * dur), rng, 30.0) for _ in names])
    rip2 = times(dur, 20, rng, stop=79.0)
    add_at(data[0], fs, rip2, burst(fs, 150.0, 0.08, 80.0))
    sz = np.arange(90.5, 105.0, 0.5)
    for i in range(3):
        add_at(data[i], fs, sz, burst(fs, 110.0, 0.08, 50.0))
    recs.append(write_recording(root, "02", "ictal", fs, data, names, events=[
        {"onset": "0", "duration": "80", "trial_type": "baseline", "source": "fixture"},
        {"onset": "90", "duration": "15", "trial_type": "seizure", "source": "fixture"}], onset="90",
        baseline=("0", "80")))
    truth[recs[-1]["recording_id"]] = {"ripple": rip2.tolist(), "seizure": [90, 105]}

    fs, dur = 512.0, 30.0
    names = ["RB01-RB02", "RB02-RB03"]
    data = np.vstack([pink(int(fs * dur), rng, 30.0) for _ in names])
    recs.append(write_recording(root, "03", "rest", fs, data, names))

    with open(root / "recordings.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REC_COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(recs)
    (out / "truth.json").write_text(json.dumps(truth))
    return out


if __name__ == "__main__":
    print(make(Path(sys.argv[1] if len(sys.argv) > 1 else "data/fixture")))
