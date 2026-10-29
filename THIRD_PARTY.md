# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, installed from PyPI as binary wheels at the versions pinned in
`requirements.txt`. All are permissively licensed (BSD, MIT, Apache-2.0, PSF-based, MPL-2.0
for two pure-Python helpers); nothing GPL is linked in-process.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| epycom | 0.3 | HFO detectors: `detect_hfo_ll` (line length, Gardner et al. 2007), `detect_hfo_rms` (RMS, Staba et al. 2002), `detect_hfo_cs_beta` (CS, Cimbálník et al. 2018) | Apache-2.0 (its `LICENSE.txt`; the PyPI metadata says "BSD 3.0", the file and the GitHub repository say Apache-2.0, which is what we follow) | https://github.com/ICRC-BME/epycom |
| MNE-Python | 1.12.1 | Reads EDF/BDF/BrainVision/FIF/EEGLAB lazily, one window at a time | BSD-3-Clause | https://mne.tools |
| NumPy | 1.26.4 | Arrays (epycom 0.3 needs NumPy < 1.27). The wheel bundles OpenBLAS (BSD-3-Clause) and the GCC Fortran runtime (GPL-3.0 with the GCC Runtime Library Exception, which places no obligations on code linked against it) | BSD-3-Clause | https://numpy.org |
| SciPy | 1.17.1 | Filters, resampling, peak finding. Bundles OpenBLAS and the GCC Fortran runtime as above | BSD-3-Clause | https://scipy.org |
| numba, llvmlite | 0.68.1, 0.50.0 | JIT compiler of the CS detector | BSD-2-Clause; llvmlite BSD-2-Clause AND Apache-2.0 WITH LLVM-exception (bundled LLVM) | https://numba.pydata.org |
| pandas | 3.0.6 | Imported by epycom | BSD-3-Clause | https://pandas.pydata.org |
| joblib, cloudpickle | 1.6.0, 3.1.2 | Channels in parallel (process pool) | BSD-3-Clause | https://joblib.readthedocs.io |
| Matplotlib | 3.11.2 | MNE dependency (no plots are made) | Matplotlib License (PSF-based, BSD-compatible); bundled DejaVu/STIX fonts under their own free licences | https://matplotlib.org |
| contourpy, cycler, kiwisolver | 1.3.3, 0.12.1, 1.5.1 | Matplotlib dependencies | BSD-3-Clause | PyPI |
| fonttools, pyparsing | 4.68.0, 3.3.3 | Matplotlib dependencies | MIT | PyPI |
| pillow | 13.0.0 | Matplotlib dependency | MIT-CMU | https://python-pillow.org |
| decorator, lazy-loader, packaging | 5.3.1, 0.6, 26.3 | MNE dependencies | BSD-2-Clause; BSD-3-Clause; Apache-2.0 OR BSD-2-Clause | PyPI |
| Jinja2, MarkupSafe | 3.1.6, 3.0.4 | MNE dependencies (HTML reports, unused) | BSD-3-Clause | https://palletsprojects.com |
| pooch, platformdirs, requests, certifi, charset-normalizer, idna, urllib3 | 1.9.0, 4.14.0, 2.34.2, 2026.10.23, 3.5.2, 3.20, 2.8.0 | MNE dependencies (dataset downloads, never used: no network at run time) | BSD-3-Clause; MIT; Apache-2.0; MPL-2.0; MIT; BSD-3-Clause; MIT | PyPI |
| tqdm | 4.70.1 | MNE dependency (progress bars) | MPL-2.0 AND MIT | https://tqdm.github.io |
| python-dateutil, six, tzdata | 3.0.0, 1.17.0, 2026.6 | pandas/Matplotlib dependencies | Apache-2.0/BSD dual; MIT; Apache-2.0 | PyPI |
| python:3.12.16-slim-bookworm | — | Base image (Debian 12 + CPython 3.12, pip) | PSF-2.0 (Python); Debian packages under their own free licences; pip MIT | https://hub.docker.com/_/python |

scikit-learn, which epycom 0.3 declares as a dependency, is not installed: none of the
detectors imports it.

Data used for testing (not in the image or the repository): Multi-Expert Seizure Annotation (Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E; Pennsieve Discover 530, doi:10.26275/dely-fzwx, CC BY 4.0).
