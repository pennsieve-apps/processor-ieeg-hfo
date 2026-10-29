"""
Unit tests on a synthetic fixture (tests/make_fixture.py); no Docker, no real data.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from processor import detect as D
from processor import ieeg as I
from processor import main as M
from tests import make_fixture

REC1 = "sub-01_task-rest_run-01"    # 2048 Hz
REC2 = "sub-02_task-ictal_run-01"   # 1024 Hz, seizure at 90 s
REC3 = "sub-03_task-rest_run-01"    # 512 Hz
PARAMS = ("DETECTOR", "BANDS", "MIN_SFREQ", "FAST_RIPPLE_MIN_SFREQ", "THRESHOLD", "WINDOW_MS", "MIN_DURATION_MS",
          "SEGMENT", "EXCLUDE_SEIZURE", "SEIZURE_MARGIN_S", "WINDOW_S", "THREADS", "LABEL")


@pytest.fixture(scope="module")
def fixture_dir(tmp_path_factory) -> Path:
    return make_fixture.make(tmp_path_factory.mktemp("fixture"))


def run(monkeypatch, input_dir: Path, output_dir: Path, **env) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("THREADS", "1")
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    return M.main()


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def rates(out: Path) -> dict[tuple[str, str, str], dict]:
    return {(r["recording_id"], r["channel"], r["band"]): r for r in read_csv(out / "hfo_rates.csv")}


def statuses(out: Path) -> dict[str, str]:
    prov = json.loads((out / "provenance.json").read_text())
    return {r["recording_id"]: r["status"] for r in prov["recordings"]}


@pytest.fixture(scope="module")
def default_run(fixture_dir, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    mp = pytest.MonkeyPatch()
    try:
        assert run(mp, fixture_dir, out) == 0
    finally:
        mp.undo()
    return out


def test_outputs_and_provenance(default_run):
    for name, header in (("hfo.csv", M.HFO_HEADER), ("hfo_rates.csv", M.RATES_HEADER)):
        with open(default_run / name) as fh:
            assert next(csv.reader(fh)) == header
    prov = json.loads((default_run / "provenance.json").read_text())
    assert prov["app"] == "processor-ieeg-hfo" and prov["version"] == "0.1.0" and prov["status"] == "ok"
    assert prov["resolved"]["threshold"] == 3.0 and prov["resolved"]["window_ms"] == 5.0
    st = statuses(default_run)
    assert st[REC1] == "ok" and st[REC2] == "ok"
    assert st[REC3] == "skipped: sampling rate 512 Hz < 1000"


@pytest.mark.parametrize("detector", ["linelength", "rms", "cs"])
def test_detectors_find_the_planted_events(fixture_dir, tmp_path, monkeypatch, detector):
    assert run(monkeypatch, fixture_dir, tmp_path, DETECTOR=detector) == 0
    truth = json.loads((fixture_dir / "truth.json").read_text())
    r = rates(tmp_path)
    n_rip = int(r[(REC1, "LA01-LA02", "ripple")]["n_hfo"])
    n_fr = int(r[(REC1, "LA02-LA03", "fast_ripple")]["n_hfo"])
    assert 0.9 * len(truth[REC1]["ripple"]) <= n_rip <= 1.1 * len(truth[REC1]["ripple"])
    assert 0.7 * len(truth[REC1]["fast_ripple"]) <= n_fr <= 1.1 * len(truth[REC1]["fast_ripple"])
    # background and a spike-only channel: (almost) nothing
    for ch in ("LA03-LA04", "LB01-LB02"):
        for band in ("ripple", "fast_ripple"):
            assert int(r[(REC1, ch, band)]["n_hfo"]) <= 2
    # ripples are found where they were planted
    ev = [e for e in read_csv(tmp_path / "hfo.csv") if e["recording_id"] == REC1 and e["channel"] == "LA01-LA02"]
    planted = np.array(truth[REC1]["ripple"])
    hits = [np.any((planted >= float(e["start_s"]) - 0.08) & (planted <= float(e["end_s"]))) for e in ev]
    assert np.mean(hits) > 0.9
    assert all(e["detector"] == detector for e in ev)


def test_ripples_only_at_1024_hz_and_seizure_excluded(default_run):
    r = rates(default_run)
    assert (REC2, "RA01-RA02", "ripple") in r and (REC2, "RA01-RA02", "fast_ripple") not in r
    assert float(r[(REC2, "RA01-RA02", "ripple")]["minutes"]) == pytest.approx(80 / 60, abs=0.01)
    ev = [float(e["start_s"]) for e in read_csv(default_run / "hfo.csv") if e["recording_id"] == REC2]
    assert ev and not [t for t in ev if t >= 80]
    prov = json.loads((default_run / "provenance.json").read_text())
    rec2 = [x for x in prov["recordings"] if x["recording_id"] == REC2][0]
    assert rec2["bands"] == ["ripple"] and any("fast ripples not analysed" in n for n in rec2["notes"])


def test_fast_ripples_forced_at_1024_hz(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, FAST_RIPPLE_MIN_SFREQ="1000") == 0
    assert (REC2, "RA01-RA02", "fast_ripple") in rates(tmp_path)


def test_keep_seizure(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, EXCLUDE_SEIZURE="false", BANDS="ripple") == 0
    ev = [float(e["start_s"]) for e in read_csv(tmp_path / "hfo.csv")
          if e["recording_id"] == REC2 and e["channel"] == "RA03-RA04"]
    assert len([t for t in ev if 90 <= t <= 105]) > 10
    assert {e["band"] for e in read_csv(tmp_path / "hfo.csv")} == {"ripple"}


def test_every_recording_below_min_sfreq_is_a_success(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path / "out", MIN_SFREQ="5000") == 0
    assert read_csv(tmp_path / "out" / "hfo.csv") == [] and read_csv(tmp_path / "out" / "hfo_rates.csv") == []
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    assert prov["status"].startswith("ok: nothing analysed")
    assert set(statuses(tmp_path / "out").values()) == {
        "skipped: sampling rate 2048 Hz < 5000", "skipped: sampling rate 1024 Hz < 5000",
        "skipped: sampling rate 512 Hz < 5000"}


def test_segment_baseline_and_label(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, SEGMENT="baseline", LABEL="pre") == 0
    r = rates(tmp_path / "pre")
    assert float(r[(REC2, "RA01-RA02", "ripple")]["minutes"]) == pytest.approx(80 / 60, abs=0.01)


def test_plain_files_without_recordings_csv(fixture_dir, tmp_path, monkeypatch):
    edf = next((fixture_dir / "ieeg" / "sub-01").rglob("*.edf"))
    (tmp_path / "in").mkdir()
    shutil.copy(edf, tmp_path / "in" / "recording.edf")  # no sidecars at all
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out", BANDS="ripple") == 0
    r = rates(tmp_path / "out")
    assert ("recording", "LA01-LA02", "ripple") in r
    assert ("recording", "LB02-LB03", "ripple") in r  # without channels.tsv nothing says it is bad


def test_no_input_and_bad_parameters(tmp_path, fixture_dir, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 1
    assert run(monkeypatch, fixture_dir, tmp_path / "o2", DETECTOR="hilbert") == 2
    assert run(monkeypatch, fixture_dir, tmp_path / "o2", BANDS="gamma") == 2


def test_failed_recording_is_skipped(fixture_dir, tmp_path, monkeypatch):
    shutil.copytree(fixture_dir / "ieeg", tmp_path / "in")
    next((tmp_path / "in" / "sub-02").rglob("*.edf")).write_bytes(b"broken")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out", BANDS="ripple") == 0
    st = statuses(tmp_path / "out")
    assert st[REC1] == "ok" and st[REC2] == "failed"


# --- units ------------------------------------------------------------------------

def test_bands_for():
    assert D.bands_for(2048, ["ripple", "fast_ripple"], 2000)[0] == ["ripple", "fast_ripple"]
    assert D.bands_for(1024, ["ripple", "fast_ripple"], 2000)[0] == ["ripple"]
    assert D.bands_for(1000, ["ripple", "fast_ripple"], 1000)[0] == ["ripple"]  # 500 Hz is not below Nyquist
    assert D.bands_for(1024, ["ripple", "fast_ripple"], 1000)[0] == ["ripple", "fast_ripple"]


def test_energy_event_bounds_are_corrected():
    """epycom measures the first half of each window; the reported event must hug the burst."""
    fs = 2048.0
    rng = np.random.default_rng(0)
    x = rng.standard_normal(int(fs * 20)) * 5
    s = int(10 * fs)
    b = make_fixture.burst(fs, 150.0, 0.05, 80.0)
    x[s:s + len(b)] += b
    y = D.bandpass(x, fs, "ripple")
    for det, w, thr in (("linelength", 5.0, 3.0), ("rms", 3.0, 5.0)):
        ev = D.staba_criteria(y, D.energy_events(y, fs, det, w, thr), fs, 0.006)
        assert len(ev) == 1
        start, stop = ev[0][0] / fs, ev[0][1] / fs
        assert 10.0 <= start <= 10.025 and 10.025 <= stop <= 10.055, (det, start, stop)


def test_staba_merge_and_peaks():
    fs = 1000.0
    y = np.zeros(1000)
    y[100:160] = np.sin(2 * np.pi * 100 * np.arange(60) / fs) * 10   # 6 cycles: 12 rectified peaks
    y[500:515] = np.sin(2 * np.pi * 100 * np.arange(15) / fs) * 10   # 1.5 cycles: 3 peaks
    kept = D.staba_criteria(y, [[100, 130], [135, 160], [500, 515]], fs, 0.006)
    assert kept == [[100, 160]]
