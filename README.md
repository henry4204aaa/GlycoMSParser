# GlycoMSP

**GlycoMSP** (Glyco Mass Spectra Parser) is a Python (tkinter) desktop workflow for
glycan MS/MS data: RAW/mzML → MS2-indexed tables → MAS / CGA annotation →
selected-fragment-ion feature matrices → local Random Forest training and
scan-level prediction reports.

Current version: **v1.10**

---

## Overview

GlycoMSP connects several steps of a glycan LC–MS/MS analysis that are usually
handled by separate tools, while preserving traceability between source spectra,
annotations, and the datasets generated from them:

1. RAW / mzML conversion to MS2-indexed spectral tables
2. Manual annotation support (MAS) for trainable-dataset preparation
3. Constraint-based glycan annotation (CGA) of candidate compositions
4. Selected-fragment-ion feature-matrix construction
5. Local Random Forest training and scan-level prediction reports

GlycoMSP was originally developed for positive-ion-mode permethylated N-glycan
LC–MS/MS datasets.

---

## Key features

- **Two annotation routes** — Manual Annotation Support (MAS) and Constraint-based
  Glycan Annotation (CGA), both producing MS2-indexed, traceable tables.
- **Traceable MS2-indexed data** — links are preserved between source spectra,
  annotations, metadata JSON, and the generated trainable datasets.
- **Local machine learning** — selected-fragment-ion feature matrices feed a local
  Random Forest classifier; no data leaves the machine.
- **Scan-level prediction reports** — per-scan predictions with summary reporting.
- **Score B (optional, experimental)** — a motif-aware re-ranking module driven by a
  user-curated Excel workbook. Score B was *not* used in the main analyses and is
  provided as an optional advanced feature. See
  [`docs/score_b_design.md`](docs/score_b_design.md).

Composition notation follows the project convention:

```text
H = Hexose      N = HexNAc      S = Neu5Ac
G = Neu5Gc      KDN = KDN         F = Fucose
```

Example: `F1H5N4S1`

---

## Installation

GlycoMSP is currently used from source. Clone the repository and create a Python
environment:

```bash
git clone https://github.com/henry4204aaa/GlycoMSParser
cd GlycoMSParser/src
pip install -r requirements.txt          # Python 3.10–3.12
# Python 3.13+:
# pip install -r requirementspy313.txt
```

**Environment**

- Python 3.10–3.12 recommended. Python 3.13 has also been tested (H.-C. Chang);
use requirementspy313.txt, as a separate pin is needed because some dependency wheels are not yet available for 3.13.

- Core packages: `pandas`, `numpy`, `scikit-learn`, `openpyxl`, `pymzml`,
  `joblib` / `skops` (and `pymsfilereader` on Windows). See the repo
  `requirements*.txt` for the full list.
- **Windows is required for Thermo RAW conversion** (`pymsfilereader` + the Thermo
  libraries). **mzML input is cross-platform.**

---

## Quick start

Launch the GUI (the active, code-reviewed entry point):

```bash
python ./src/mspfileloaderv14.py
```

Typical workflow:

1. Convert RAW or mzML files into MS2-indexed spectral CSV / TSV tables.
2. Create or link metadata JSON files.
3. Prepare a MAS or CGA-derived annotation set.
4. Build a trainable dataset (selected-fragment-ion feature matrix).
5. Train or apply a Random Forest model.
6. Export prediction results and the scan-level report.

For the full GUI walkthrough, see the **user manual** (maintained separately) and the
documentation wiki under [`docs/`](docs/).

### Default ML configuration

The default training configuration matches the manuscript:

- Random Forest, **400 trees**, `class_weight="balanced"`, `random_state=42`,
  `min_samples_split=2`, `min_samples_leaf=1`
- Class-eligibility filter: **≥ 5 spectra** per class
- **80/20 stratified** train/test split (no validation fold by default)

---

## Data & reproducibility

- **Primary data (zebrafish RAW + original annotations):** GlycoPOST —
  <http://doi.org/10.50821/GLYCOPOST-GPST000224>
- **Reproducibility archive (code + reference workbook + U-937 reanalysis files
  where redistribution is permitted):** Zenodo —
  <https://doi.org/10.5281/zenodo.20823042>
- The reference Score B workbook for positive-ion-mode permethylated N-glycans is
  provided with GlycoMSP and in the Zenodo package.

---

## License

GlycoMSP is released under the **MIT License**.
Copyright (c) 2026 Huan-Chuan Tseng. See [`LICENSE.md`](LICENSE.md).

---

## Citation

> Tseng, H.-C. *GlycoMSP enables traceable glycan MS/MS annotation and local
> machine-learning-ready dataset construction.* Manuscript submitted to
> *Bioinformatics Advances* (2026).

A formal citation will be added once the reference is available. Until then, please
cite this repository and the Zenodo DOI above.

---

## Author

Developed by **Huan-Chuan Tseng** as a PhD research project at the University of Tokyo 
(glycomics / mass spectrometry / bioinformatics).
