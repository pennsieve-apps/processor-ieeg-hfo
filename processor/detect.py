"""
High-frequency oscillation (HFO) detectors, from epycom 0.3 (Apache-2.0, Cimbálník et al.,
FNUSA-ICRC Brno):

- linelength: epycom.detect_hfo_ll, the line-length energy detector of Gardner et al. 2007,
  on the band-passed signal;
- rms: epycom.detect_hfo_rms, the RMS detector of Staba et al. 2002, on the band-passed signal;
- cs: epycom.detect_hfo_cs_beta, the CS algorithm of Cimbálník et al. 2018, on the broadband
  signal (it filters into its own 17 bands from 44 Hz up and models each band's statistics).

For linelength and rms, the candidate events then go through the Staba et al. 2002 criteria:
events closer than 10 ms are merged, and an event is kept when the rectified band-passed
signal has at least 6 peaks above 3 SD of the band-passed signal (3 oscillations) and it
lasts at least MIN_DURATION_MS.

epycom's windowed line length / RMS value of a window is computed over the window's first
half (a 'same'-mode convolution read at index 0). We pass twice the effective window we want
and correct the event end for it, so that WINDOW_MS is the window that is really measured.
"""
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.signal import butter, find_peaks, sosfiltfilt

BANDS = {"ripple": (80.0, 250.0), "fast_ripple": (250.0, 500.0)}
MERGE_S = 0.010        # Staba 2002: events closer than 10 ms are one event
PEAK_SD = 3.0          # Staba 2002: peaks of the rectified signal > 3 SD ...
MIN_PEAKS = 6          # ... at least 6 of them
DEFAULTS = {           # per detector: effective window (ms), threshold (SD above the mean)
    "linelength": {"window_ms": 5.0, "threshold": 3.0},
    "rms": {"window_ms": 3.0, "threshold": 5.0},
    "cs": {"window_ms": None, "threshold": 0.1},
}


def bands_for(fs: float, wanted: list[str], fast_ripple_min_sfreq: float) -> tuple[list[str], list[str]]:
    """(bands analysed, notes on bands left out) for a sampling rate."""
    out, notes = [], []
    for b in wanted:
        lo, hi = BANDS[b]
        if b == "fast_ripple" and fs < fast_ripple_min_sfreq:
            notes.append(f"fast ripples not analysed: sampling rate {fs:g} Hz < FAST_RIPPLE_MIN_SFREQ "
                         f"{fast_ripple_min_sfreq:g} Hz")
            continue
        if hi >= fs / 2:
            notes.append(f"{b} not analysed: its upper edge {hi:g} Hz is not below Nyquist ({fs / 2:g} Hz)")
            continue
        out.append(b)
    return out, notes


def bandpass(x: np.ndarray, fs: float, band: str) -> np.ndarray:
    lo, hi = BANDS[band]
    sos = butter(4, [lo, hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x)


def _merge(events: list[list[int]], gap: int) -> list[list[int]]:
    out: list[list[int]] = []
    for s, e in sorted(events):
        if out and s - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def staba_criteria(y: np.ndarray, events: list[list[int]], fs: float, min_duration_s: float) -> list[list[int]]:
    """Merge events < 10 ms apart; keep those >= min_duration_s with >= 6 rectified peaks > 3 SD."""
    events = _merge(events, int(round(MERGE_S * fs)))
    thr = PEAK_SD * float(np.std(y))
    keep = []
    for s, e in events:
        if (e - s) / fs < min_duration_s:
            continue
        peaks, _ = find_peaks(np.abs(y[s:e]), height=thr)
        if len(peaks) >= MIN_PEAKS:
            keep.append([s, e])
    return keep


def energy_events(y: np.ndarray, fs: float, detector: str, window_ms: float, threshold: float) -> list[list[int]]:
    """epycom line-length / RMS detection on a band-passed signal; returns [start, stop) samples."""
    from epycom.event_detection.hfo.ll_detector import detect_hfo_ll
    from epycom.event_detection.hfo.rms_detector import detect_hfo_rms

    w = max(2, int(round(window_ms / 1000.0 * fs)))      # effective window, samples
    ws = 2 * w                                           # what epycom is given (it measures the first half)
    step = max(1, w // 2)
    overlap = step / ws                                  # epycom: increment = ceil(ws * overlap)
    inc = int(math.ceil(ws * overlap))
    fn = detect_hfo_ll if detector == "linelength" else detect_hfo_rms
    raw = fn(y, fs=int(round(fs)), threshold=threshold, window_size=ws, window_overlap=overlap)
    out = []
    for start, stop in raw:
        stop_eff = min(len(y), max(start + w, stop - w - inc))  # end of the last window above threshold
        out.append([int(start), int(stop_eff)])
    return out


def cs_events(x: np.ndarray, fs: float, bands: list[str], threshold: float) -> dict[str, list[list[int]]]:
    """epycom CS detector on the broadband signal; band detections assigned to ripple/fast ripple by their centre."""
    from epycom.event_detection.hfo.cs_detector import detect_hfo_cs_beta

    out: dict[str, list[list[int]]] = {b: [] for b in bands}
    dets = detect_hfo_cs_beta(np.asarray(x, dtype=np.float64), fs=int(round(fs)), threshold=threshold,
                              cycs_per_detect=4.0)
    for d in dets:
        start, stop, low, high, is_conglom = d[0], d[1], d[2], d[3], d[8]
        if is_conglom:
            continue
        centre = math.sqrt(float(low) * float(high))
        for b in bands:
            lo, hi = BANDS[b]
            if lo <= centre < hi:
                out[b].append([int(start), int(stop)])
    return {b: _merge(v, 0) for b, v in out.items()}


def detect_channel(x: np.ndarray, fs: float, detector: str, bands: list[str], window_ms: float | None,
                   threshold: float, min_duration_s: float, offset_s: float,
                   core: tuple[float, float]) -> tuple[list[tuple[float, float, str, float]], str | None]:
    """
    HFOs in one channel window. x starts at offset_s (seconds from the recording start); only
    events starting inside core are returned, as (start_s, end_s, band, amplitude_uv), where
    amplitude is the largest |band-passed signal| within the event. Returns (events, error).
    """
    try:
        x = np.asarray(x, dtype=np.float64)
        if not np.all(np.isfinite(x)):
            return [], "non-finite samples"
        if np.ptp(x) == 0:
            return [], None
        found: list[tuple[float, float, str, float]] = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            filtered = {b: bandpass(x, fs, b) for b in bands}
            if detector == "cs":
                per_band = cs_events(x, fs, bands, threshold)
            else:
                per_band = {b: staba_criteria(filtered[b], energy_events(filtered[b], fs, detector, window_ms,
                                                                         threshold), fs, min_duration_s)
                            for b in bands}
        for b in bands:
            y = filtered[b]
            for s, e in per_band[b]:
                t0 = offset_s + s / fs
                if core[0] <= t0 < core[1]:
                    found.append((round(t0, 4), round(offset_s + e / fs, 4), b, round(float(np.max(np.abs(y[s:e]))), 2)))
    except Exception as exc:  # noqa: BLE001 - one channel's failure must not sink the recording
        return [], f"{type(exc).__name__}: {exc}"
    return found, None
