# GlycoMSP

**GlycoMSP** (Glyco Mass Spectra Parser) is a Python-based workflow for glycomic LC–MS/MS data preprocessing, glycan annotation support, dataset preparation, and machine-learning-assisted glycan spectral analysis.

Current version: **v1.09**  
Last README update: **2026-05-09**

> This repository is under active development for PhD thesis and manuscript preparation. The codebase is being reorganized toward a cleaner publication-ready structure.

---

## Overview

GlycoMSP is designed to support glycan LC–MS/MS analysis by connecting several steps that are often handled separately:

1. Raw spectral data conversion and preprocessing
2. Manual annotation dataset preparation
3. Constraint-based glycan annotation (CGA)
4. Fragment-ion-aware scoring and label refinement
5. Machine-learning dataset construction
6. Random Forest-based glycan spectral classification and prediction

The project was originally developed for permethylated glycan LC–MS/MS datasets and is currently being extended and reorganized for broader usability.

---

## Main Features

### 1. Raw / mzML Data Preprocessing

GlycoMSP supports conversion and extraction of MS/MS spectral information from compatible mass spectrometry data files.

Current support includes:

- Thermo RAW file extraction through Thermo MSFileReader / pymsfilereader
- mzML reader support under active development
- Export of converted spectral tables for downstream annotation and ML processing
- Metadata generation for traceability

Thermo RAW extraction requires a Windows environment with the Thermo library installed.

---

### 2. Manual Annotation Support

The manual annotation workflow supports preparation of trainable datasets from:

- Converted MS/MS CSV or TSV files
- Annotation Excel workbooks
- Fragment ion lists
- Metadata JSON files

The workflow is intended to preserve links between source spectra, annotations, metadata, and generated datasets.

---

### 3. Constraint-Based Glycan Annotation

The CGA workflow generates candidate glycan compositions based on user-defined biological and structural constraints.

Supported concepts include:

- N-glycan and O-glycan composition generation
- Permethylated glycan mass calculation
- User-configurable glycan composition boundaries
- Optional glycan motif / glycotope constraints
- In-silico composition export
- Pseudo-label generation for ML training support

Composition notation currently follows the project convention:

```text
H = Hexose
N = HexNAc
S = Neu5Ac
G = Neu5Gc
K = KDN
F = Fucose
```

Example:

```text
F1H5N4S1
```

---

### 4. Score B / Motif-Aware Re-ranking

GlycoMSP includes an experimental motif-aware scoring workflow for improving CGA-derived candidate ranking.

The scoring system can use:

- Fragment ion evidence
- Motif-support rules
- Composition-consistency checks
- Unexpected-evidence penalties
- User-selected motif policies

This component is currently being refined for manuscript and supplementary documentation.

---

### 5. Machine Learning Workflow

GlycoMSP includes tools for building ML-compatible glycan spectral datasets and training Random Forest models.

Current ML-related functions include:

- Trainable dataset construction
- Optional Non-glycan class handling
- Class filtering by minimum sample count
- Majority-class balancing
- Train/test/validation splitting
- Random Forest model training
- Label encoding and prediction export
- Prediction summary report generation

The current ML implementation primarily uses `scikit-learn`.

---

## Repository Structure

The repository is currently being reorganized into the following structure:

```text
GlycoMSParser/
├── README.md
├── .gitignore
├── src/
│   └── active GlycoMSP source code
├── docs/
│   └── documentation, workflow notes, and method descriptions
├── templates/
│   └── reusable annotation, metadata, and workbook templates
├── tmpdata/
│   └── local-only temporary data and archived development files
├── logs/
│   └── local-only log files
└── .venv/
    └── local Python virtual environment
```

`tmpdata/`, `logs/`, `.venv/`, and editor-specific folders are intended to remain local and are excluded from normal repository tracking.

---

## Installation

A stable package installation workflow is still under preparation.

For current development use, clone the repository and create a Python environment manually:

```bash
git clone https://github.com/henry4204aaa/GlycoMSParser
cd GlycoMSParser
```

Install commonly required packages:

```bash
cd src
pip install requirements.txt
```

If you are using Python3.13+

```bash
cd src
pip install requirementspy313.txt
```

Additional packages may be required depending on the workflow being used.

---

## Requirements

### General Python Requirements

- Python 3.10+ recommended
- pandas
- numpy
- openpyxl
- scikit-learn
- joblib
- tkinter (built-in)

### Thermo RAW File Extraction

Thermo RAW support requires:

- Windows OS
- Thermo MSFileReader installed
- pymsfilereader installed and correctly linked

If these are not available, use mzML-based workflows where possible.

---

## Basic Usage

The current GUI entry point is under `src/`.
The current active code-reviewed GlycoMSP GUI is `mspfileloaderv14.py`
Please run the file 
```bash
python .\src\mspfileloaderv14.py
```

Typical workflow:

1. Launch the GlycoMSP GUI
2. Convert RAW or mzML files into spectral CSV / TSV format
3. Create or link metadata JSON files
4. Prepare manual annotation or CGA-derived datasets
5. Generate trainable datasets
6. Train or apply ML models
7. Export prediction results and reports

Example command structure may change as the repository is reorganized.

---

## Current Development Status

As of **2026-05-09**, GlycoMSP is in an active pre-publication development stage.

Current priorities include:

- Repository cleanup and restructuring
- README and documentation updates
- Template file organization
- GUI workflow stabilization
- JSON schema clarification
- Score B documentation
- Manuscript- and thesis-aligned software release preparation

---

## Notes on Data and Templates

This repository may include small template files or demonstration files only.

Large raw files, private experimental datasets, generated intermediate files, logs, and temporary JSON files should not be committed to the repository.

Recommended local-only folders:

```text
tmpdata/
logs/
archive/
```

---

## Citation

A formal citation will be added after manuscript submission or publication.

For now, please cite the repository or contact the author if using GlycoMSP in collaborative work.

---

## Author

Developed by **Huan-Chuan Tseng**  
PhD research project, glycomics / mass spectrometry / bioinformatics

---

## License

License information will be added before public release.

---

## Changelog

### v1.09 — 2026-05-09

- Repository cleanup and restructuring initiated
- Active code moved toward `src/`
- Temporary data and archived development files separated from active source code
- README updated for publication-preparation stage
- Documentation and template folders introduced

### Early Development

- 2023-01-02: Initial GitHub repository upload
- 2022-12: Fundamental core scripts and early pipeline prototypes created
