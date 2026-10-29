# sample-input

Real input for this app is the `ieeg` output of **ieeg-prep**: a BIDS-iEEG derivative with
`recordings.csv`, and per recording `*_desc-preproc_ieeg.edf`, `_channels.tsv`,
`_events.tsv` (see `docs/conventions.md` in workflow-epilepsy-ieeg). Any BIDS-iEEG dataset
or loose EDF/BrainVision/FIF files work too. Only recordings sampled at `MIN_SFREQ` (1000 Hz)
or faster are analysed; slower ones are skipped by design.

The workflows' demo data are public on Pennsieve Discover:

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. Multi-Expert Seizure
> Annotation. Pennsieve Discover 530 (2026), doi:10.26275/dely-fzwx (CC BY 4.0).
> BIDS-iEEG SEEG seizure clips at 1024 Hz with 120 s before the clinically marked onset.

At 1024 Hz only ripples (80–250 Hz) are analysed by default: fast ripples (250–500 Hz) need
`FAST_RIPPLE_MIN_SFREQ` (2000 Hz). The HUP interictal atlas (Discover 179) is 200 Hz and is
skipped.

```sh
cd ../workflow-epilepsy-ieeg && scripts/fetch-sample-data.sh        # -> data/seizures
# run ieeg-prep on data/seizures (-> data/work/prep/ieeg), then, for the pre-ictal baselines:
cd ../processor-ieeg-hfo && make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg
#   (with SEGMENT=baseline in dev.env)
```

For tests and `make run`, `tests/make_fixture.py` writes a small synthetic derivative
instead (`make fixture` -> `data/fixture/ieeg/`): a 2048 Hz clip with planted ripples
(120 Hz) and fast ripples (350 Hz) on two channels, a background channel, a channel with
sharp spikes only and a bad channel; a 1024 Hz clip with ripples, a baseline and a seizure;
and a 512 Hz clip that is skipped (`data/fixture/truth.json` holds the planted times).
