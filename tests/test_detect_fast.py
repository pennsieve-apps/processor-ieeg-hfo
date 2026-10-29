"""
v0.1.1 computes epycom's windowed line length / RMS for every window at once instead of in
epycom's Python loop. These tests hold it to epycom 0.3 itself (still installed: its detectors
are the oracle) and to the hfo.csv v0.1.0 wrote for the synthetic fixture.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from epycom.event_detection.hfo.ll_detector import detect_hfo_ll
from epycom.event_detection.hfo.rms_detector import detect_hfo_rms

from processor import detect as D
from tests import make_fixture
from tests.test_processor import fixture_dir, run  # noqa: F401 - the module-scoped fixture

V010 = Path(__file__).parent / "data" / "v0.1.0"
EPYCOM = {"linelength": detect_hfo_ll, "rms": detect_hfo_rms}


def _windows(w: int) -> tuple[int, float, int]:
    """(window given to epycom, overlap, increment) as energy_events derives them."""
    ws = 2 * w
    overlap = max(1, w // 2) / ws
    return ws, overlap, int(math.ceil(ws * overlap))


@pytest.mark.parametrize("detector", ["linelength", "rms"])
def test_vectorised_windows_equal_epycom(detector):
    """Random lengths, windows and bursty signals: the same events as epycom's own detector."""
    rng = np.random.default_rng(7)
    for _ in range(150):
        n, w = int(rng.integers(40, 30000)), int(rng.integers(2, 40))
        ws, overlap, inc = _windows(w)
        if n <= ws:
            continue
        y = rng.standard_normal(n) * np.where(rng.random(n) < 0.01, 8.0, 1.0)
        y = np.convolve(y, np.ones(3), "same")
        want = [(int(s), int(e)) for s, e in EPYCOM[detector](y, fs=1024, threshold=3.0, window_size=ws,
                                                              window_overlap=overlap)]
        assert D._epycom_events(y, detector, 3.0, ws, inc) == want, (n, ws, inc)


@pytest.mark.parametrize("detector,w,thr", [("linelength", 5, 3.0), ("rms", 3, 5.0)])
def test_energy_events_equal_epycom_on_band_passed_bursts(detector, w, thr):
    """A band-passed signal with planted ripples, through energy_events and through epycom."""
    fs = 2048.0
    rng = np.random.default_rng(11)
    x = rng.standard_normal(int(fs * 60)) * 5
    for t in rng.uniform(1, 59, 40):
        b = make_fixture.burst(fs, float(rng.uniform(90, 230)), 0.05, 40.0)
        s = int(t * fs)
        x[s:s + len(b)] += b[: len(x) - s]
    y = D.bandpass(x, fs, "ripple")
    ws, overlap, inc = _windows(max(2, int(round(w / 1000 * fs))))
    raw = EPYCOM[detector](y, fs=int(fs), threshold=thr, window_size=ws, window_overlap=overlap)
    assert len(raw) > 30
    assert D._epycom_events(y, detector, thr, ws, inc) == [(int(s), int(e)) for s, e in raw]


@pytest.mark.parametrize("detector", ["linelength", "rms"])
def test_fixture_output_equals_v010(fixture_dir, tmp_path, monkeypatch, detector):  # noqa: F811
    """hfo.csv of the synthetic fixture (both bands) is byte for byte what v0.1.0 wrote."""
    assert run(monkeypatch, fixture_dir, tmp_path, DETECTOR=detector, FAST_RIPPLE_MIN_SFREQ="1000") == 0
    assert (tmp_path / "hfo.csv").read_text() == (V010 / f"hfo_{detector}.csv").read_text()
