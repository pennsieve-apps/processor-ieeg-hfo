# processor-ieeg-hfo

A step of the **Interictal Biomarker Mapping** workflow
(`ieeg-prep → {ieeg-spikes, ieeg-bandpower → ieeg-normative, ieeg-hfo} → ieeg-report`). It
finds high-frequency oscillations (HFOs) — ripples (80–250 Hz) and fast ripples (250–500 Hz)
— on every good channel of an intracranial EEG recording (SEEG or ECoG) and counts them per
channel per minute.

**For the epileptologist.** Interictal HFOs, fast ripples in particular, are more frequent in
the seizure onset zone, and removing tissue with high HFO rates has been associated with good
surgical outcome in several studies. They are short (tens of ms), small (µV to tens of µV) and
tedious to mark by eye on an expanded time scale; this step marks them automatically on every
channel over the same interictal time, seizures left out, so channels can be ranked by rate
next to spike rates and band-power abnormality. HFOs need a high sampling rate: below
1000 Hz the step skips the recording on purpose (and says so), and fast ripples are only
looked for at 2000 Hz or more by default.

**For the researcher.** Three detectors from epycom 0.3 (FNUSA-ICRC Brno): line length of the
band-passed signal (Gardner et al. 2007; default), its RMS (Staba et al. 2002), both followed by
the Staba et al. 2002 event criteria, and the CS algorithm (Cimbálník et al. 2018), which models
the statistics of each of its 17 frequency bands on the broadband signal. Zelmann et al. 2012
compared such detectors and found them sensitive to their thresholds and to the
background; run more than one (`LABEL`) when it matters. Outputs are long-format CSV keyed by
`recording_id`, `channel` and `band` (workflow-epilepsy-ieeg `docs/conventions.md`).

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `ieeg` | input | root | The BIDS-iEEG derivative of ieeg-prep: `recordings.csv`, `*_desc-preproc_ieeg.edf` with `_channels.tsv`, `_events.tsv`. Also any BIDS-iEEG dataset, or plain EDF/BDF/BrainVision/FIF/EEGLAB files, found anywhere under `INPUT_DIR` (folders and symlinks followed). |
| `hfo` | output | `OUTPUT_DIR/` (or `OUTPUT_DIR/<LABEL>/`) | `hfo.csv`, `hfo_rates.csv`, `provenance.json` |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `DETECTOR` | `linelength` | `linelength`, `rms` or `cs` (see *How it works*). |
| `BANDS` | `ripple,fast_ripple` | Bands to look in: `ripple` 80–250 Hz, `fast_ripple` 250–500 Hz. A band is analysed only when its upper edge is below Nyquist. |
| `MIN_SFREQ` | `1000` | Recordings sampled below this are **not analysed**: status `skipped: sampling rate X Hz < 1000` in `provenance.json` and the log, no rows, and the run still succeeds (exit 0), also when every recording is skipped. |
| `FAST_RIPPLE_MIN_SFREQ` | `2000` | Fast ripples are looked for only at this rate or above (see *Sampling rate*). Set `1000` to force them at 1024 Hz. |
| `THRESHOLD` | `auto` | linelength / rms: SD above the mean of the line length / RMS within the window (`auto`: 3 for line length, epycom's default; 5 for RMS, Staba 2002). cs: its detection threshold between 0 and 1 (`auto`: 0.1, epycom's default). |
| `WINDOW_MS` | `auto` | Window of the line length / RMS (`auto`: 5 ms line length, 3 ms RMS as in Staba 2002), slid in half-window steps. Not used by cs. |
| `MIN_DURATION_MS` | `6` | Shortest event kept by linelength and rms (Staba 2002: 6 ms). |
| `SEGMENT` | `all` | `all`: the whole recording. `baseline`: only the baseline interval of `recordings.csv` or a `baseline` event; recordings without one are analysed whole (noted in `provenance.json`). |
| `EXCLUDE_SEIZURE` | `true` | Leave known seizures out (`onset_s` of `recordings.csv`, `seizure` events; `ictal` too for plain BIDS), with `SEIZURE_MARGIN_S` on both sides; to the end of the recording when the seizure's end is unknown. `artifact` events are always left out. |
| `SEIZURE_MARGIN_S` | `10` | Seconds left out before the onset and after the end of each seizure. |
| `WINDOW_S` | `300` | Windows the recording is read and analysed in: all good channels at once, plus 2 s of context on each side. Detection thresholds are computed per channel and window. |
| `THREADS` | `0` | Channels analysed in parallel (processes); `0` = one per vCPU. |
| `LABEL` | (empty) | When set, write to `OUTPUT_DIR/<LABEL>/` (e.g. one run per detector). |

## Outputs

**`hfo.csv`**: one row per HFO.

| Column | |
|---|---|
| `recording_id` | the BIDS entity chain of the source recording |
| `channel` | the `_channels.tsv` name (`LA01-LA02`) |
| `start_s`, `end_s` | seconds from the recording start |
| `band` | `ripple` or `fast_ripple` |
| `detector` | `linelength`, `rms` or `cs` |
| `amplitude_uv` | largest absolute value of the band-passed signal within the event, µV |

**`hfo_rates.csv`**: one row per analysed channel and band, zeros included: `recording_id`,
`channel`, `band`, `n_hfo`, `minutes` (time of that channel actually analysed), `rate_per_min`,
`detector`. Bad channels, non-iEEG channels and bands not analysed at the recording's sampling
rate get no row; `provenance.json` says why.

**`provenance.json`**: app and version, tool versions, parameters (and the resolved threshold
and window), every input file seen, timings, and per recording: status (`ok`, `failed` with the
error, `skipped: sampling rate X Hz < 1000`, `skipped: no analysable time …`), sampling rate,
bands analysed, analysed segments and what was excluded and why, minutes, HFO counts per band,
channels not analysed with the reason, notes (e.g. why fast ripples were not analysed).

## How it works

1. **Finds the recordings** (every `recordings.csv`, plus signal files no CSV lists, with
   their sidecars) and handles them one at a time, opened lazily with MNE.
2. **Checks the sampling rate**: below `MIN_SFREQ`, the recording is skipped (a designed
   branch, not a failure). Bands are chosen from `BANDS` and the rate.
3. **Per window** of `WINDOW_S` (all good channels read at once as float32 µV), each channel in
   a worker process:
   - **linelength / rms**: zero-phase 4th-order Butterworth band-pass per band; epycom's
     `detect_hfo_ll` / `detect_hfo_rms` mark windows whose line length / RMS exceeds mean +
     `THRESHOLD` × SD; then the **Staba et al. 2002 criteria**: events closer than 10 ms are
     merged, and an event is kept when it lasts ≥ `MIN_DURATION_MS` and the rectified
     band-passed signal has ≥ 6 peaks above 3 SD of the band-passed signal (three oscillations).
   - **cs**: epycom's `detect_hfo_cs_beta` on the broadband signal; its per-band detections
     are assigned to ripple or fast ripple by their band's centre frequency and overlapping
     ones merged.
4. Keeps the events that start inside the window's core and writes the rows.

epycom computes each window's line length / RMS over the first half of the window it is
given (a `'same'`-mode convolution read at index 0). We pass twice the window we want and
correct the event end, so `WINDOW_MS` is the window actually measured
(`tests/test_processor.py::test_energy_event_bounds_are_corrected`).

## Sampling rate

Ripples need a rate well above 500 Hz and fast ripples well above 1000 Hz: a band must sit
below Nyquist *and* below the recording system's anti-aliasing filter, which typically starts
attenuating at a fraction of Nyquist. At 1024 Hz (the demo data), the 500 Hz edge of the
fast-ripple band is 12 Hz below Nyquist and inside that roll-off, so fast ripples would be
both attenuated and aliased: by default they are only looked for at ≥ 2000 Hz, and a 1024 Hz
recording gets **ripples only**, with a note in `provenance.json`. Between 1000 and 2000 Hz,
`FAST_RIPPLE_MIN_SFREQ=1000` forces them; treat those rates as lower bounds.

## Runtime

`computeTypes: [standard]`, **2 vCPU / 4 GiB**, 6 h timeout. Channels run in parallel
(`THREADS`, default one per vCPU); memory is bounded by the window, not the recording length.
The cost is CPU: epycom slides the line-length / RMS window in a Python loop.

Measured in Docker (linux/amd64 **emulated** on an M-series Mac, `--cpus 2 --memory 4g`) on 4
seizure clips of the demo dataset (Multi-Expert Seizure Annotation, 1024 Hz, so ripples only;
as-recorded referential montage: HUP224 230 channels × 338 s, HUP230 106 × 384 s, HUP250
170 × 296 s, HUP257 118 × 463 s), default `linelength`:

| Run | Analysed | CPU time | Wall time | Peak memory (cgroup `memory.peak`) |
|---|---|---|---|---|
| whole recordings (`EXCLUDE_SEIZURE=false`) | 24.7 min, 3,725 channel-minutes | 1,706 s | 895 s | 1.46 GB |
| 120 s pre-ictal baselines (`SEGMENT=baseline`) | 8.0 min, 1,248 channel-minutes | 481 s | 249 s | 1.50 GB |
| HUP atlas, 60 recordings at 200 Hz (all skipped) | — | — | 7 s | — |

`memory.peak` includes the page cache of the EDF files read. About 0.46 CPU-seconds per
channel-minute per band emulated, ~0.29 natively (M-series, measured per channel; `rms` and
`cs` cost about the same, 1.4–1.7 s per 5 channel-minutes natively). Expected on real amd64
hardware: a 1 h, 128-channel recording at 1024 Hz (ripples) ≈ 40 CPU-minutes ≈ 20 min on 2
vCPU; at 2048 Hz with both bands about four times that, well inside the 6 h timeout. More
vCPUs scale nearly linearly (`THREADS`). The image is ~200 MB compressed (~840 MB unpacked).

## Results on real data

**Seizure clips, 120 s pre-ictal baseline** (`SEGMENT=baseline SEIZURE_MARGIN_S=0`, the 4
clips above, 1024 Hz, **ripples only**: fast ripples are not analysed at 1024 Hz by default,
see *Sampling rate*). Ripple rate in the channels of the experts' consensus unequivocal
electrographic onset (`ueo_consensus`) vs all other channels (AUC: probability that a UEO
channel has a higher rate than a non-UEO channel):

| Recording | UEO channels | linelength: UEO vs others (median /min), AUC | rms AUC | cs AUC |
|---|---|---|---|---|
| HUP224 ictal339135 (230 ch.) | 6 | 1.25 vs 1.00, 0.78 | 0.50 (1 ripple in all) | 0.66 |
| HUP230 ictal296207 (106 ch.) | 3 | 9.00 vs 1.50, **0.89** | 0.89 | 0.73 |
| HUP250 ictal209153 (170 ch.) | 2 | 3.00 vs 3.00, 0.38 | 0.12 | 0.08 |
| HUP257 ictal387681 (118 ch.) | 5 | 3.50 vs 4.00, 0.44 | 0.80 | 0.39 |

Mean ripple rates per channel: linelength 0.9–5.0/min, rms 0.0–2.6/min (Staba's 5 SD
threshold is strict), cs 1.5–8.4/min. Two minutes of referential 1024 Hz data per seizure is
a hard test: the onset channels stand out in two of four clips. Expect better separation on
the bipolar montage of ieeg-prep and on longer interictal recordings, and do not read a
single clip's HFO ranking as a localisation.

**ieeg-prep output** (`data/work/prep/ieeg` of workflow-epilepsy-ieeg: all 24 demo clips,
bipolar montage, 1024 Hz so ripples only, `SEGMENT=baseline`: ieeg-prep's pre-ictal baseline
of up to 90 s). Ripple rate in channels labelled `onset=true` in `_channels.tsv` vs the others:

| Detector | AUC, median (IQR) | AUC > 0.5 | median rate, onset vs other | mean rate per channel |
|---|---|---|---|---|
| `linelength` | **0.78** (0.61–0.92) | 19 of 24 | 2.2 vs 0.7 /min | 2.3 /min |
| `rms` | 0.81 (0.60–0.88) | 19 of 24 | 1.3 vs 0.0 /min | 0.9 /min |
| `cs` | 0.65 (0.50–0.76) | 18 of 24 | 1.3 vs 0.0 /min | 6.7 /min |

On the bipolar montage the onset channels stand out clearly better than on the referential
clips above.

**HUP atlas** (200 Hz): all 60 recordings are skipped with `skipped: sampling rate 200 Hz <
1000`, empty tables, exit 0, as designed.

**Synthetic fixture** (`tests/`): all three detectors find the planted ripples (30/30 at
2048 Hz, 26/26 at 1024 Hz) and fast ripples (linelength and rms 22/22, cs 17/22), with at most
one false event on the background and spike-only channels.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python        # unit tests on a synthetic fixture
make run PYTHON=.venv/bin/python         # Docker over the fixture -> data/output/
make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg   # an ieeg-prep output
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-ieeg-hfo:dev .
docker run --rm --platform linux/amd64 --cpus 2 --memory 4g --env-file dev.env \
  -v /path/to/ieeg:/data/input:ro -v "$PWD/data/output":/data/output processor-ieeg-hfo:dev
```

## Limitations

- **Automated HFO detection has many false positives.** Sharp transients (spikes, artifacts)
  produce ringing in a band-pass filter that looks like an oscillation; the 6-peak criterion
  removes many of these but not all, and nothing here separates "HFOs on spikes" from
  isolated HFOs, or physiological from pathological ripples (both occur, e.g. in
  hippocampus and occipital cortex). Review high-rate channels on the raw signal.
- **Thresholds are relative** to each channel's own window (mean + k·SD, or CS's band
  statistics): a channel with continuous high-frequency activity raises its own threshold and
  reports fewer events; a very quiet channel reports small ones. Rates compare channels
  within a recording better than across patients.
- Gardner et al. 2007 tuned the line-length detector for 30–85 Hz events; here it is applied to
  the ripple and fast-ripple bands, as epycom and other toolboxes do.
- **Montage matters**: on a common or distant reference, an HFO at the reference shows on
  every channel and muscle activity leaks in; ieeg-prep's bipolar montage is the intended
  input. The demo results below are on the referential (as-recorded) clips.
- Short recordings give noisy rates: 2 minutes of baseline yields single-digit counts per
  channel.
- Line noise harmonics inside the bands (e.g. 120/180/240 Hz) must be removed upstream
  (ieeg-prep notches them); they would otherwise be counted as ripples.
- Stretches shorter than 10 s are not analysed; `minutes` says how much time each rate is
  based on.

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles epycom (Apache-2.0), MNE-Python
(BSD-3-Clause) and their dependencies; see `THIRD_PARTY.md`. Nothing GPL is linked in-process.

## Citations

HFO detectors and their evaluation:

> Gardner AB, Worrell GA, Marsh E, Dlugos D, Litt B. **Human and automated detection of
> high-frequency oscillations in clinical intracranial EEG recordings.** Clin Neurophysiol
> 118(5):1134–1143 (2007). doi:10.1016/j.clinph.2006.12.019, PMID 17382583.

> Staba RJ, Wilson CL, Bragin A, Fried I, Engel J Jr. **Quantitative Analysis of
> High-Frequency Oscillations (80–500 Hz) Recorded in Human Epileptic Hippocampus and
> Entorhinal Cortex.** J Neurophysiol 88(4):1743–1752 (2002). doi:10.1152/jn.2002.88.4.1743,
> PMID 12364503.

> Cimbálník J, Hewitt A, Worrell G, Stead M. **The CS algorithm: A novel method for high
> frequency oscillation detection in EEG.** J Neurosci Methods 293:6–16 (2018).
> doi:10.1016/j.jneumeth.2017.08.023, PMID 28860077.

> Zelmann R, Mari F, Jacobs J, Zijlmans M, Dubeau F, Gotman J. **A comparison between
> detectors of high frequency oscillations.** Clin Neurophysiol 123(1):106–116 (2012).
> doi:10.1016/j.clinph.2011.06.006, PMID 21763191.

Implementation: epycom 0.3, ElectroPhYsiology COmputational Module, FNUSA-ICRC Brno,
https://github.com/ICRC-BME/epycom (Apache-2.0).

Reading the recordings:

> Gramfort A, Luessi M, Larson E, Engemann DA, Strohmeier D, Brodbeck C, Goj R, Jas M,
> Brooks T, Parkkonen L, Hämäläinen M. **MEG and EEG data analysis with MNE-Python.** Front
> Neurosci 7:267 (2013). doi:10.3389/fnins.2013.00267, PMID 24431986.

> Holdgraf C, Appelhoff S, Bickel S, et al. (39 authors; all in `CITATION.cff`).
> **iEEG-BIDS, extending the Brain Imaging Data Structure specification to human
> intracranial electrophysiology.** Sci Data 6:102 (2019). doi:10.1038/s41597-019-0105-7.

Test data:

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. **Multi-Expert Seizure
> Annotation.** Pennsieve Discover (2026), doi:10.26275/dely-fzwx (CC BY 4.0).

Every DOI was checked against Crossref (articles) or DataCite (datasets): title, author list
and year; PMIDs against PubMed. Crossref lists only the first author of the MNE-Python paper;
its full author list here and in `CITATION.cff` is PubMed's.
