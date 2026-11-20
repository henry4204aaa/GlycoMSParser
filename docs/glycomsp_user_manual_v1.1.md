# GlycoMSP User Manual

**Version:** v1.1
**Tool version covered:** GlycoMSP 1.14.1 (`mspfileloaderv14.py`)
**Last updated:** 2026-11-20
**Audience:** Users (chemists and glycobiologists), thesis defense committee, manuscript reviewers

---

## Versions covered by this manual

Manual v1.0 described GlycoMSP v1.10. This manual, v1.1, describes GlycoMSP 1.14.1 and marks in the text which version added each feature that v1.10 did not have. The table lists what changed for you in each version. `CHANGELOG.md` in the repository holds the full release notes.

| Version | Date | What changed for you | Where to read |
|---|---|---|---|
| **v1.10** | 2026-06-26 | The published, frozen release. Everything in this manual that carries no later version mark. | Parts 1–8 |
| **1.11** (preview) | 2026-09-23 | An optional GlyTouCan / WURCS reference database. A checkbox in Prepare Dataset fills empty accession and WURCS cells from a local database, and a command-line tool builds that database. Off by default. | §3.7, §3.9 |
| **1.12** (preview) | 2026-10-01 | The Manage References window to select, build and inspect the reference database, a settings file that remembers your choice, and a showcase database snapshot. Score B reads lowercase `s` and `p` as sulfate and phosphate. | §3.7, §5.6 |
| **1.13** (preview) | 2026-10-31 | Conversion reports what it read: an outcome in the log, a status file and a diagnostics report. Spectra without a precursor charge are kept. A conversion options dialog offers two ways to supply a missing charge. Native (underivatized) glycan masses in the CGA route. Two command-line tools. Reader corrections. | §2.2, §3.8, §3.9, §4.9, §5.1 |
| **1.14** | 2026-11-01 | The publication build for the manuscript revision and the first version that can annotate low-resolution (ion-trap) MS2 data. CGA mass tolerances in ppm or Da, a visible precursor gate in CGA → Trainable, recorded tolerances that prefill prediction, and negative-mode native N-glycan ion lists. | §2.4, §2.6, §2.8, §3.8, §4.10 |
| **1.14.1** | 2026-11-18 | A test release. The application is identical to 1.14. It adds the automated test suite and the validation tools, and documents five defects that the new tests found. | §7.1, §7.2, §8.5 |

**Which version the manuscript cites.** v1.10 is the release cited by the manuscript and archived on Zenodo. 1.14 is the publication build for the manuscript revision, and 1.14.1 is the test release that accompanies it.

**Where to get each version.** All versions are in the GitHub repository `https://github.com/henry4204aaa/GlycoMSParser`. The repository has no release tags. A version is identified by the version line in `CHANGELOG.md` and by the title of the main window. The reproducibility archive is on Zenodo (see Appendix C).

**What stays the same.** GlycoMSP is developed additively. The default behaviour of v1.10 is preserved and new behaviour ships behind explicit options, so results produced with v1.10 remain reproducible. The exceptions are the reader corrections of 1.13, which are listed in §4.9.

**What the title bar shows.** The main window of GlycoMSP 1.14.1 reports `v1.14`, because 1.14.1 did not change the application.

---

## How to read this manual

This manual is built so different readers can enter at different points and not waste time.

If you're a **chemist or glycobiologist** running your own samples through the tool, start at **Part 2 — Quickstart** and follow it linearly. Part 1 gives you the conceptual picture in two minutes; Part 4 explains *why* the scoring works the way it does if you're curious. If your data are native glycans or come from an ion-trap instrument, read **§3.8** before you convert your files.

If you're on the **thesis committee** and want a structural picture of the system, read **Part 1** for the pipeline orientation, then jump to **Part 4 — Concepts** for the rationale behind the two annotation routes and the scoring design.

If you're a **manuscript reviewer** evaluating reproducibility, read **Part 1** for orientation, then go to **Part 8 — Reproducing the manuscript results**, which maps each result to its Zenodo package with exact commands and expected numbers. §8.5 covers the automated tests and validation tools, and §8.6 the external-dataset validation runs. **Part 5 — Reference** lists default parameters, JSON schemas, and output file specifications, and the **Concepts** section (Part 4, including §4.8 on the MAS-ML evaluation design) is where the manuscript-aligned formulas and evaluation rationale live.

A glossary of mass-spectrometry and machine-learning terms is in **Appendix A**. The manual assumes you already work with mass spectrometry; machine-learning terminology is defined in-prose on first use and again in the glossary. **Appendix D** lists what changed between manual v1.0 and v1.1, section by section.

---

# Part 1 — Introduction

## 1.1 What GlycoMSP does

GlycoMSP is a desktop application that converts glycomics tandem mass spectrometry (MS2) data into machine-learning-ready training matrices, trains classifiers on those matrices, and applies the trained models to predict glycan compositions for unlabeled spectra.

It supports two complementary paths to a labeled training set. The first, **Constraint-based Glycan Annotation (CGA)**, generates labels computationally from in-silico libraries and an ion-evidence scoring scheme; this path is appropriate when no expert-validated annotations exist. The second, **Manual annotation support (MAS)**, ingests labels from an expert-prepared annotation Excel sheet and validates them against the spectra; this path is appropriate when ground-truth annotations are available. The two paths converge on the same feature representation, so models can be trained on data from either path or a mixture of both.

Around the two annotation paths, the application provides supporting workflows: per-RAW-file metadata capture and conversion, optional negative-class (Non-glycan) sampling, optional motif-aware re-ranking via Score B, ion-mining for diagnostic-fragment discovery, and a prediction step that produces ranked predictions with confidence scores plus reviewer-friendly output reports.

GlycoMSP 1.11 to 1.14 add features that leave the v1.10 defaults unchanged. They are a local GlyTouCan / WURCS reference database that annotates compositions with accession and WURCS strings (§3.7), conversion reports and options for spectra without a precursor charge (§2.2 and §4.9), native (underivatized) glycan masses and annotation of ion-trap data (§3.8), and command-line tools (§3.9). The new options are off by default, and with all of them at their defaults, results equal those of v1.10 except for the reader corrections listed in §4.9.

## 1.2 Where GlycoMSP fits

```
┌──────────────────┐    ┌─────────────┐    ┌──────────────────┐    ┌────────────────┐
│  Mass            │ →  │  RAW or     │ →  │     GlycoMSP     │ →  │  Trained model │
│  spectrometer    │    │  mzML files │    │   (this manual)  │    │  + predictions │
└──────────────────┘    └─────────────┘    └──────────────────┘    └────────────────┘
```

Upstream, GlycoMSP is fed by any LC-MS/MS instrument that produces RAW (Thermo) or mzML files. On Windows with the Thermo MSFileReader library installed, GlycoMSP reads `.raw` directly; on macOS, Linux, or Windows machines without MSFileReader, you supply pre-converted `.mzML` files. Every conversion also writes a status file and a diagnostics report that state what was read (§2.2, §4.9).

Downstream, GlycoMSP produces three artifact families: (1) a trained classifier saved as `.skops` (with `.joblib` and `.pkl` fallbacks accepted on load), (2) prediction CSVs accompanied by a class-summary file and a JSON/Markdown run summary, and (3) review queues highlighting spectra that fell below confidence thresholds for human inspection. With the reference-database fill ticked, the prediction table also carries the columns `Predicted_GlyToucan_ID` and `Predicted_WURCS` (§3.7).

## 1.3 The two paths at a glance

| | **CGA pipeline (computed labels)** | **MAS pipeline (expert labels)** |
|---|---|---|
| **When to use** | You have raw spectra but no expert annotations | You have expert annotations validated against the spectra |
| **You provide** | RAW/mzML, sample metadata, configuration flags for the in-silico library generator | RAW/mzML, sample metadata, an annotation Excel (MSlist + ionlist sheets) |
| **GlycoMSP generates** | An in-silico glycan library, candidate compositions per spectrum, Score A (and optionally Score B), per-spectrum top-ranked composition labels | A merged annotation table linking expert labels to spectra, with optional automatic Non-glycan negatives |
| **Output** | ML-ready wide CSV with constraint-based labels | ML-ready wide CSV with expert-verified labels |
| **Mixable for training?** | Yes — both produce the same column schema | Yes |

The two paths converge after labeling: feature extraction (log-transformed intensity at PPM-matched fragment masses) is identical, and the ML training and prediction steps make no distinction between CGA-derived and MAS-derived rows. Manuscript Methods uses the names **CGA-ML** and **MAS-ML** for the end-to-end pipelines built around each labeling path.

## 1.4 What you'll need

Before starting, gather the following:

- **One or more `.raw` or `.mzML` files** with MS2 fragmentation data. CGA and MAS both require MS2 spectra. To avoid confusion, each conversion batch must be **either all RAW or all mzML** — a single batch does not mix the two formats. Before each conversion in the GUI, the conversion options dialog offers two ways to supply a precursor charge for spectra whose file reports none, and both are off by default (§2.2).
- **Sample metadata** for each RAW or mzML file: *Experiment Title*, *Experiment Description*, *Author running this analysis*, *Raw data acquired date*, *Glycan Type* (`N`, `O` or `N+O`), *Mass Analyzer charge mode* (`+` or `-`) and *Derivatization Type*. Sulfate and phosphate options for negative mode are set later in the CGA Setup window. The *Derivatization Type* list of the Metadata Editor includes `Native` and `Native(Reduced)` for underivatized glycans (§3.8). GlycoMSP will prompt you for these the first time it sees a file; subsequent files in the same batch can reuse the previous metadata as a starting point.
- **An ion list (Excel or CSV)** with reference fragment m/z values and glycotope/structural identifiers. This is required for Score A computation in the CGA pipeline and for fragment-feature extraction. The same ion list is reusable across samples sharing the same derivatization and ionization mode.
- **An expert annotation Excel** — only for the MAS pipeline. Must contain an `MSlist` sheet (scan IDs and glycan labels) and an `ionlist` sheet (reference fragment ions). See Part 5.1 for the exact column requirements.
- **A Score B workbook** (optional). An Excel file (`.xlsx`, `.xlsm`, `.xltx` or `.xltm`) with six required sheets that encode glycotopic motif rules and penalty parameters. It is only needed if you wish to enable motif-aware re-ranking. By default Score B is disabled, and Score A alone drives candidate selection. See Part 3.3 for an overview and `score_b_config_instructions.md` (in the manuscript supplementary materials) for full configuration.

You start the application from the repository root with `python src/mspfileloaderv14.py` (§2.1). Two further notes on the operating environment. First, file paths inside saved JSON method files are stored as POSIX-style for the experiment file and native (OS-dependent) for the per-sample method files; the manual conversion is handled internally and you do not need to think about it. Second, the GUI title bar shows the build version (for example, *GlycoMSP File Manager GUI v1.14 (20261101 publication build)*). GlycoMSP 1.14.1 has the same application code as 1.14, so its title bar also reads v1.14. Cite the title when reporting issues.

---

# Part 2 — Quickstart: your first end-to-end run via the CGA pipeline

This Quickstart walks one realistic end-to-end run through the **CGA pipeline** (constraint-based glycan annotation): from raw instrument files all the way to a trained model and predictions on unlabeled data. We chose CGA because it does not require pre-existing expert annotations; everything you need to label your data is generated by the tool. The MAS pipeline (Part 3.1) and Score B enrichment (Part 3.3) build on the same backbone and only diverge at specific steps.

> **Conventions used here.** Italicised words are screen labels (button text, field labels) reproduced verbatim. `Code-style` text refers to file names, column names, or values you should literally type. Curly braces like `{this}` mark placeholders you replace with your own values. Glycan composition labels are written **fucose-first** in the fixed order `F` (Fuc) → `H` (Hex) → `N` (HexNAc) → `S` (Neu5Ac) → `G` (Neu5Gc) → `KDN`, with a count after each present monosaccharide (e.g., `F1H5N4` is the bi-antennary composition with 1 fucose, 5 hexoses, and 4 N-acetylhexosamines; `F1H5N5KDN2` carries 2 KDN). KDN is written as the literal token `KDN`, not `K`.

> **What is new since v1.10.** Everything added in GlycoMSP 1.11 to 1.14.1 is off by default or additive. With all new options at their defaults, results equal those of v1.10, except for the corrections listed under *Fixed* in the changelog entry of GlycoMSP 1.13 (for example negative-mode precursor masses, ion-trap RAW peak lists and MS3 rows). Sections of this Quickstart that describe a later addition say which version added it.

## 2.1 Launch GlycoMSP and meet the main window

Launch GlycoMSP by running `python src/mspfileloaderv14.py` from the project directory. The title bar reads *GlycoMSP File Manager GUI v1.14 (20261101 publication build)*. GlycoMSP 1.14.1 has the same application code as 1.14, so the title bar still says v1.14. The first thing you see is the **main window**:

> **[SCREENSHOT: `GlycoMSP_screenshots/v1.14main.jpg` — Main launcher window.]**
> Title bar: *GlycoMSP File Manager GUI v1.14 (20261101 publication build)*. Top button row (left to right): *Add Raw File*, *Add mzML File*, *Clear*, *About*, *Save Log*. Centre is a large white text panel — this is the **log** where the application prints status messages as it works. Below the log is a status indicator (showing *Idle* in blue with a small green progress bar). Bottom of the window is a labelled section titled *Analysis Tools* containing three primary action buttons: *Convert Raw to CSV*, *Prepare Dataset*, *Run ML Analysis*.

The three *Analysis Tools* buttons map directly to the three stages of a GlycoMSP run. You'll work them left to right: **convert** raw instrument files to scan CSVs, **prepare** a labeled training dataset from those CSVs, then **run ML analysis** to train a model and use it to make predictions.

The top-row file-add buttons (*Add Raw File* / *Add mzML File*) populate a working list of files for batch conversion. *Clear* empties this list. *About* opens a window with the version and last-update date of each module. *Save Log* asks for a file name and writes the application log, including its debug lines, to a `.log` file.

## 2.2 Add raw files and convert to CSV

Click **Add Raw File** (Windows) or **Add mzML File** (macOS / Linux / pre-converted) and select one or more files. Selected paths are appended to the application's working list and shown in the log panel.

Now click **Convert Raw to CSV** in the *Analysis Tools* section. GlycoMSP first shows the **conversion options dialog** (added in GlycoMSP 1.13). The **Metadata Editor** window opens only after you leave this dialog with *Convert*. *Cancel* returns you to the main window and converts nothing.

> **[SCREENSHOT: `GlycoMSP_screenshots/v1.14_convertoptions.jpg` — Conversion options dialog.]**
> Dialog titled *Conversion options*. First line of text: *Precursor charge for MS2 spectra whose file does not report one.* followed by *Leave both off for the usual conversion.* An unticked checkbox *Calculate the charge from the MS1 isotope pattern (charge 1–3)* with a grey explanation below it. An entry field, empty, labelled *Dataset charge for the rest (blank = off):* with a grey explanation below it. A closing line *Applies to every file in this batch; saved in each project's metadata file.* Buttons *Convert* and *Cancel* at the bottom right.

Both options are off every time the dialog opens, and with both off the conversion is the same as in earlier versions. The dialog's own explanations read as follows:

- *Calculate the charge from the MS1 isotope pattern (charge 1–3)*: "Reads the spacing of the precursor's isotope peaks in its parent MS1 scan. Used only where the file gives no charge; a charge the file reports is never changed. Needs profile MS1 (ion-trap RAW; for mzML, MS1 converted without peak picking)."
- *Dataset charge for the rest (blank = off):* "One charge for every MS2 spectrum still without one after the step above, e.g. 1. The polarity comes from each scan, so enter 1, not −1. Leave blank if your files report charges."

An entry that is not a whole number other than 0 stops with a dialog titled *Conversion options* that reads "Dataset charge must be a whole number other than 0 (for example 1 or 2), or blank." The choice is saved in the project metadata file as `Conversion charge comb` and `Conversion dataset charge`. §4.9 explains what the two options do and what the readers do and do not do. §3.8 covers ion-trap data.

After you click *Convert*, the Metadata Editor window opens for the first file in your list. **The conversion of that file starts at once.** The form stays locked while the conversion runs in a background thread.

> **[SCREENSHOT: `GlycoMSP_screenshots/convert_metadatawindow.jpg` — Metadata Editor window while a RAW file is converting (zebrafish brain example).]**
> Window titled *Converting: {file name}*. Status text *Converting raw file... please wait.* in blue. The form has seven fields, all greyed out while the conversion runs: *Experiment Title*, *Experiment Description*, *Author running this analysis* and *Raw data acquired date* (text entries), then *Glycan Type* (list: *N*, *O*, *N+O*), *Mass Analyzer charge mode* (list: *+*, *-*) and *Derivatization Type* (list: *PerMe(Freeend)*, *PerMe(Reduced)*, *Native*, *Native(Reduced)*, *Others*). Below the form is a button row: *Import Metadata*, *Use Last Metadata*, *Clear*, *Generate Dataset*, *Cancel*. The last row holds the entry *Raw File (optional):* with a *Browse* button. When the conversion has finished, the window title changes to *Metadata for: {file name}*, the status text changes (see the outcomes below) and the form unlocks.

Fill in the metadata for this first file. GlycoMSP does not check the fields before it saves them. *Glycan Type*, *Mass Analyzer charge mode* and *Derivatization Type* matter later because they prefill the CGA Setup window (§2.4). For underivatized glycans choose *Native* or *Native(Reduced)* as *Derivatization Type* (§3.8). If you have a previous metadata JSON saved from an earlier run, click **Import Metadata** to populate the fields from it. For the second and later files in a batch the fields are filled with the previous file's metadata when the conversion ends, and **Use Last Metadata** does the same on demand.

When the form is complete, click **Generate Dataset**. The application:

1. Asks you for an output folder (*Select folder to save metadata and output files*). This question is asked once per batch, for the first file only.
2. Asks for a project name (dialog *Project Name*, prompt *Enter a name for this batch/project:*). The save name for the file is `{project name}_{YYYYmmdd_HHMMSS}_{raw file name without extension}`.
3. Writes a JSON metadata file `{save name}.json` in the output folder.
4. Promotes the temporary scan CSV and MS3 CSV to their final names `ms2_{save name}.csv` and `ms3_{save name}.csv`, moves the diagnostics files described below next to them, and saves the application log of the session as `{save name}.log`.
5. Advances to the next file of the batch in the same window (a new conversion starts at once), or shows the message *All metadata have been completed.* and closes the editor when the batch is complete.

The resulting **scan CSV** is the input to everything else. It is a tab-separated text file with the extension `.csv`. One row per MS2 spectrum, with nine columns in this order: `entry_no`, `MS1scan_no`, `MS1_isolationmass`, `MS1_monoisolationmass`, `chargeState`, `protonatedmass`, `MS2scan_no`, `peaklist` and `peakintensity`. The key columns for the rest of the workflow are `MS2scan_no`, `protonatedmass`, `peaklist` and `peakintensity`. Spectra without peaks are not exported. A spectrum whose precursor charge is unknown is kept with `chargeState` 0 and `protonatedmass` 0 (§4.9).

The second file, `ms3_{save name}.csv`, holds MS3 spectra in the same tab-separated format with seven columns: `entry no` (written with a space), `MS2scan_no`, `MS2_isolationmass`, `MS2_monoisomass`, `MS3scan_no`, `peaklist` and `peakintensity`. The file is written for every conversion, and the MS2-only workflow of this Quickstart does not use it.

> **Note.** GlycoMSP focuses on MS2 scan-level traceability tied to the peak list and intensities; retention time and other LC/EIC-related parameters are not retained in the current design.

**Conversion outcome (added in GlycoMSP 1.13).** Every conversion ends with one of four outcomes. GlycoMSP writes the outcome to the log panel in a line that starts with `[reader]`:

```
[reader] outcome: complete (exit code 0)
[reader] outcome: complete with warnings (exit code 3)
[reader] outcome: partial (exit code 4)
[reader] outcome: failed (exit code 5)
```

The outcome line is followed by one line per warning, `[reader] warning: {text}`, and for a failed conversion by `[reader] reason: {text}`. The last line gives the diagnostics report, `[reader] diagnostics: {path of the .reader_diag.txt file}`. If no usable report exists, the line reads `[reader] diagnostics file not used ({state})` instead. The exit code is the same number the command-line converter returns (§3.9). The status text of the Metadata Editor shows the outcome too:

| Outcome | What it means | Status text in the Metadata Editor | What GlycoMSP does |
|---|---|---|---|
| complete | All spectra were read and nothing needs checking. | *Conversion complete. Ready for metadata.* (green) | Saves the files as usual. |
| complete with warnings | The conversion worked and no spectrum was lost to a decoding error, but the reader lists something you should read. The usual case is MS2 spectra without a precursor charge. They are kept with their peaks and precursor m/z, with `chargeState` 0 and `protonatedmass` 0. The reader also warns when no spectrum has a usable precursor mass and when MS2 spectra are profile data. | *Conversion complete with warnings (see log). Ready for metadata.* (orange) | Saves the files as usual. Read the `[reader] warning:` lines. |
| partial | Some spectra could not be decoded and were not exported. The CSV files hold the rest. | *Conversion partial (see log). Ready for metadata.* (orange) | Shows a dialog titled *Partial conversion* that lists the warnings. The files are saved with this warning when you click *Generate Dataset*. |
| failed | The input could not be opened, the conversion stopped, or no MS2 spectrum with peaks was exported. | See the error state below. | Never saves the result as if the conversion had succeeded. Nothing is moved to the output folder. |

A typical warning for unknown charges reads: "{n} of {m} exported MS2 spectra have no charge (written as chargeState 0, protonatedmass 0); charge-dependent steps cannot use them." Such rows pass through the later steps with mass 0 and are not excluded there (§7.2).

**Files written next to the CSV files (added in GlycoMSP 1.13).** Every conversion also produces a **status file** and a **diagnostics report**. After you click *Generate Dataset* they sit in the output folder under the same save name as the CSV files:

- `{save name}.scan_status.tsv`: the status file with one row per scan.
- `{save name}.reader_diag.json` and `{save name}.reader_diag.txt`: the diagnostics report, with the outcome, the counts of exported and skipped spectra, the counts of known and unknown precursor charges, and the reason for any failure.

The reader first writes these files next to the input file under the input file's name. The first conversion of a batch has no output folder yet, so the reader uses the folder of the input file. *Generate Dataset* then moves the files to the output folder, renames them to the save name and updates the file names written inside the JSON report. §5.1 gives the formats.

If a conversion fails (for example a file the reader cannot open), the editor enters an error state. The status text turns red and reads *An error occurred. Please check the log. You can close this editor.* The window title is prefixed with "Error:", the rest of the batch is halted, and the main window shows *Idle (error)*. A dialog titled *Conversion failed* gives the reason and the exit code, and names the diagnostics report. When the reader stopped with an unexpected exception, a dialog titled *Thermo Library Error* appears instead, with the exception text. The editor then closes by itself after a short time. No scan CSV is saved for the failed file, and the files the reader wrote stay where it wrote them. Use **Save Log** on the main window to capture the diagnostic output.

> **Reporting conversion errors.** GlycoMSP's MS2-indexed conversion has so far been tested specifically on Thermo RAW and their corresponding mzML spectra. If you hit an error, we would appreciate a report so we can improve the reliability of the initial MS2-indexed dataset preparation. Please include: the **error log**; **screenshots** (optional); the **file(s) that failed** (optional); and a short description of the spectral data — **vendor, instrument, and acquisition mode (DDA / DIA / SRM or other specific settings)**.

## 2.3 Open Prepare Dataset and add your scan CSV

Back at the main window, click **Prepare Dataset**. This opens the **Prepare Dataset** window, which is where the CGA labeling work happens.

> **[SCREENSHOT: `GlycoMSP_screenshots/Prepare_dataset_withrefDB.jpg` — Prepare Dataset window.]**
> Top-left: a *Check integrity* button (greyed out, a placeholder for a future feature). Below it: the **Dataset Explorer** tree, a hierarchical view that will hold your experiments and samples. Centre: the action button grid, organised into rows. Row 1 — *Add MS2 CSV(s)*, *Add Excels (ion list/man annotation)*, *Add Sample using metadata*. Row 2 — *Add Sample*, *Clean up empty unassigned sample tags*, *Link Sample*, *Merge Sample* (the last two start greyed out). Row 3 — *Negative Sampling Setting*, *Ion Mining*, *Browse Ion Mining result*, *Load Method*. Row 4 — *CGA manager*, *CGA → Trainable*, and an unticked checkbox *Fill GlyToucan ID / WURCS from reference DB (see Manage References)*. Row 5 — *Manage References…* and the label *Active DB:* followed by the database name, version and entry count. Below the action grid: a *Close* button. Status footer: *Experiment method file: None* and *Sample method will be saved at: None*. Bottom row: *Set Experiment Method File*, *Set Sample Method Folder*, *Load .exp.json*, *Save .exp.json*, *Export Method v1 (MAS)*, *Export Method v1 (CGA)*.

Click **Add MS2 CSV(s)** and select the scan CSV(s) you produced in Step 2.2. They appear in the Dataset Explorer under a top-level node called *Unassigned*, with one *Sample* sub-node per file. The sample names default to the CSV filename stem.

If you'd rather group samples under a meaningful experiment name from the start, you can instead click **Add Sample using metadata** and select your metadata JSON files. GlycoMSP reads the *Experiment Title* and *Raw filename* fields and creates the experiment node and sample node accordingly. (Internally these two paths just call different file-handler routines; the resulting tree structure is identical.)

**Reference database controls (added in GlycoMSP 1.11 and 1.12).** Three controls in the last rows belong to the optional reference database. *Manage References…* opens the Manage References window, where you build, select and inspect the database. The label *Active DB:* shows the database that is in use, with its name, version and entry count. The checkbox *Fill GlyToucan ID / WURCS from reference DB (see Manage References)* is off by default. When you tick it, empty GlyToucan ID and WURCS cells are filled from the database in the later steps. With the box unticked, results equal those of v1.10. §3.7 explains the database.

The Dataset Explorer is your working surface for the rest of the workflow. Right-clicking a sample exposes per-sample actions, and drag-and-drop is available, but the tree view is not fully reliable — not every click is guaranteed to register, so treat it as a convenience rather than a dependable control surface.

> **Note.** The validation status indicator (✅ linked-and-validated, ⚠️ / ❌ issues) is only updated in the **MAS** workflow, where it gates the *Link/Validate* and *Merge* buttons. The CGA workflow does not refresh this indicator in the tree; this is a minor display-only issue and is safe to ignore, as the indicator does not guard the CGA path.

## 2.4 Run CGA Manager — set up the in-silico library and Score A inputs

Select your sample (single-click in the tree), then click **CGA manager**. This opens the **Constraint-based Glycan Annotation Setup** window (called CGA Setup in this manual). The sample needs a linked scan CSV and a metadata JSON. Without them GlycoMSP shows a message and does not open the window.

> **[SCREENSHOT: `GlycoMSP_screenshots/CGA_Setup_v1.14.jpg` — CGA Setup window (N-glycan sample, *Derivatization Type* list open).]**
> Window titled *Constraint-based Glycan Annotation Setup* (about 900 × 800 pixels, resizable). From the top: the frame *Metadata (can override here)* with the lists *Glycan Type:* (*N* or *O*), *Mass Analyzer charge mode:* (*+* or *-*) and *Derivatization Type:* (open in the screenshot, with the entries *PerMe*, *PerMe(Freeend)*, *PerMe(Reduced)*, *Native*, *Native(Reduced)* and *None*), and a *Reset from file* button. The frame *In-Silico Generation Flags* with the checkbox *Show advanced options (debug /future options)* and the sub-frames *glyco-feature restraints* (checkboxes such as `allow5ac`, `allow5gc`, `allowfuc` and a spinbox `allowpsa`), *substituent(neg)* (`SO3`, `SO3_range`, `PO3H`, `PO3H_range`), *internal extension/elongation* (`arm_count (N-glycans):`, `internal_minrep:`, `internal_maxrep:`), *NG Core Flags* (`corefuc`, `bicorefuc`, `highman`, `hybrid`) and *composition restraints* (`HexA_range`, `compcheck`, `Hex_range`, `HexNAc_range`, `Neu5Ac_range`, `Neu5Gc_range`, `KDN_range`, `Fucose_range`). The frame *Files selected for CGA analysis* with three read-only rows: *In-silico glycan list (compositions):*, *Fragmentation ion list (features):* and *Score B reference:*. The frame *Mass tolerances (CGA matching)* with the rows *Precursor match:* and *Score A fragment match:*, each with a number entry, a unit list (*ppm* in the screenshot) and the grey text *default 20 ppm* or *default 10 ppm*, and a grey note about ion-trap data. At the bottom, a button row *Load CGA settings*, *Save CGA settings*, a separator, *Generate In-Silico CSV*, *Link in silico glycan list*, *Add Ion List*, *Add Score B Workbook*, *Validate Score B Workbook*, *Clear Score B Workbook*, and below it *Start CGA analysis* and *Close*.

For an *O* sample the window shows the frame *O-glycan core types (select 1–4)* (*Tn Antigen*, *Core 1* to *Core 4*) and hides the N-glycan fields `allowldnf`, `arm_count` and the *NG Core Flags*. The three values in the metadata frame come from the sample's metadata file and you can override them for this analysis. *Reset from file* restores the values read from the file.

This window has three tasks: (a) generate or link an in-silico glycan library, (b) attach an ion list, and (c) optionally attach and validate a Score B workbook. Since GlycoMSP 1.14 it also sets (d) the CGA mass tolerances.

**(a) In-silico library.** If you don't have a library file already, configure the flag panel for your glycan class (the fields shown depend on whether *Glycan Type* is *N* or *O*) and click **Generate In-Silico CSV**. GlycoMSP asks for an output folder (*Select output folder for in-silico CSV*). The application calls the composition generator (`compnewv4.NGlaunch` for N-glycans or `compnewv4.OGlaunch` for O-glycans), produces a CSV `{sample name}_insilico_{YYYYmmdd_HHMMSS}.csv` listing every theoretical composition consistent with your flags and modifier counts, and shows its path in *In-silico glycan list (compositions):*. The library is **reusable across samples** that share derivatization and ionization mode; if you already generated one, click **Link in silico glycan list** instead and pick the file directly. The entries *Native* and *Native(Reduced)* of *Derivatization Type* give underivatized precursor masses (§3.8, added in GlycoMSP 1.13).

**(b) Ion list.** Click **Add Ion List** and pick your ion-list file (`.csv`, `.xlsx` or `.xls`). Its path appears in *Fragmentation ion list (features):*. Score A and the feature extraction read one column from it, the fragment masses. The header of that column must be `mass`, `mz`, `ion_mz` or `m/z`, in any capitalisation. A CSV file must be comma-separated. In an Excel file GlycoMSP uses the sheet named `ionlist`, `ions` or `ion_list`, or otherwise the first sheet that has such a column. All other columns are ignored at this step.

**(c) Score B workbook (optional).** If you want motif-aware re-ranking, click **Add Score B Workbook** and pick your six-sheet Excel file (`.xlsx`, `.xlsm`, `.xltx`, or `.xltm`). Then click **Validate Score B Workbook** to run a quick load test. A dialog titled *Score B Workbook* reports either *Workbook validated successfully.* or the reason why the validation failed. *Clear Score B Workbook* removes the attachment. **Skip this for your first run** — Score B is disabled by default, and the CGA pipeline produces useful results from Score A alone. Part 3.3 covers Score B in depth.

**(d) Mass tolerances (added in GlycoMSP 1.14, defaults reproduce v1.10).** The *Mass tolerances (CGA matching)* frame sets the two windows that CGA analysis uses. Each is a number plus a unit, *ppm* or *Da*:

- **Precursor match** (default **20 ppm**): how far a library mass may be from the scan's precursor mass. A ppm window is relative to the observed precursor mass; a Da window is ± that many Da around it.
- **Score A fragment match** (default **10 ppm**): how far a peak may be from an ion-list mass to count as a hit (§4.1).

The fields show the defaults each time the window opens, and with the defaults CGA analysis writes the same TSV as v1.10, where these two values were fixed. An empty, zero, negative or non-numeric value stops the run with an error dialog titled *CGA tolerance*. A non-default tolerance is marked *NON-DEFAULT* in the log panel and recorded in the sample's method JSON (§5.2).

Keep the defaults for Orbitrap-class data with a monoisotopic precursor m/z. Widen them only when your instrument's mass errors are larger, for example on ion-trap data whose precursor mass is the isolation m/z (§3.8). When you widen them, do two more things:

- **Widen the CGA → Trainable fields too.** They are separate ppm fields and do not follow these settings (§2.6).
- **Check the wider window against chance (decoy check).** Make copies of the in-silico library and of the ion list with every mass increased by 3.37 Da. That offset avoids residue, isotope (≈ 1.003 Da) and charge spacings, so no decoy mass is a glycan or fragment of your sample. Run CGA analysis on the same sample with the same tolerances against the decoy copies, and compare the number of matched scans, and of scans passing Score A, with the real run. Every decoy match is a chance match: real − decoy estimates the matches above chance, and decoy ÷ real the fraction expected to be chance. A wider window is worth it while real − decoy keeps growing. On a dense library with accurate-mass data, a wide window can give more decoy matches than real ones, so decide per dataset.

## 2.5 Start CGA analysis and produce the pseudolabel TSV

With the in-silico library and ion list attached (Score B is optional), click **Start CGA analysis**. The application:

1. Loads the in-silico library and normalises composition columns.
2. For each MS2 scan, performs binary-search precursor mass matching against the library at the precursor match tolerance (default 20 ppm; set in CGA Setup from GlycoMSP 1.14, fixed at 20 ppm before), producing candidate compositions per scan. At most three candidates are kept per scan, the smallest `ppm_error` first.
3. For each scan with at least one candidate, attaches an ion-fitting score (Score A) by counting how many reference fragments match peaks in the spectrum within the Score A fragment match tolerance (default 10 ppm). §4.1 explains how hits are counted. Score A is the number of hits divided by the number of masses in the ion list, so it lies between 0 and 1. It is computed against the whole ion list, so all candidates of one scan carry the same Score A.
4. If a Score B workbook was attached and validated, runs the Score B enrichment pass: each candidate is re-evaluated against motif-evidence rules and assigned a Score B value `b ∈ [0, 1]`. Score B failure is non-fatal — the original Score A TSV is preserved.
5. Writes the result to a tab-separated file `{sample name}_CGA_{YYYYmmdd_HHMMSS}.tsv` in the folder of the scan CSV (the **pseudolabel TSV**). It holds all columns of the scan CSV and adds `composition`, `theoretical_mass`, `ppm_error`, `ion score`, `ion hit count`, `ion hits m/z`, and the label columns `Glycanannotation2`, `GlyToucan ID` and `WURCS`. The three label columns are empty unless you fill them (§3.7). When Score B is active the file also has `score_b`, `score_b_rank`, `score_b_selected` and `score_b_motif_summary`.

When the run ends, a dialog titled *CGA analysis complete* shows the path of the TSV. If no scan has a precursor match, the run stops with a message titled *CGA Aborted* and writes no TSV.

The pseudolabel TSV is your raw labeling output: one row per (scan, candidate composition) pair, plus one row for every scan without a candidate. Its size therefore follows the raw file, about 1.7 GB for a run of 100,000 MS2 scans. A file with matched rows only is planned. The next step collapses the TSV to one row per scan.

## 2.6 Build the trainable CSV (CGA → Trainable)

Back in the Prepare Dataset window with your sample selected, click **CGA → Trainable**. This opens the *CGA → Trainable (one-pass)* modal.

> **[SCREENSHOT: `GlycoMSP_screenshots/CGA_Trainable_v1.14.jpg` — CGA → Trainable modal.]**
> A modal dialog titled *CGA → Trainable (one-pass)* with seven numbered parts. Three path rows, each with a *Choose…* button: *1) CGA analysis report TSV/CSV (long):*, *2) Ion sheet (CSV/XLSX) or Manual Excel:* (with an optional *Sheet name (optional):* field below it) and *3) Salvage composition file (optional):*. The frame *4) Thresholds / selection* with the fields *min ion score* (0.07), *max |ppm_error|* (20.0), *Top N per scan* (1) and *ion list ppm* (10.0), and a second row with the checkbox *precursor gate* (ticked), a ppm entry (10) and the grey text *ppm, measured vs theoretical precursor mass (untick to skip this gate)*. The frame *5) Negative sampling* with the ticked checkbox *Add Non-glycan entries (easy negatives)*, *Max ratio (neg:pos)* (3.0), *Keep scans with < min hits to ion list* (3), the ticked checkbox *Shuffle negatives (random)* and *Seed:* (42). The frame *6) Feature building* with two radio buttons, *Extract features from long-form peaks and build feature matrix (log10(I) + 1)* (selected) and *Reuse existing trainable CSV features*, and a path entry with a *Choose…* button for the second option. The frame *7) Include mass as feature* with the unticked checkbox *Include precursor mass as feature*. Below: the text *Output: (auto-named)* and the button *Build Trainable CSV*.

Boxes 1) and 2) fill in from the sample's CGA result and ion list. Confirm that *1) CGA analysis report TSV/CSV (long):* points at the file from Step 2.5 and that *2) Ion sheet (CSV/XLSX) or Manual Excel:* points at the file from Step 2.4. Leave *5) Negative sampling* enabled and box *6) Feature building* at *Extract features from long-form peaks and build feature matrix (log10(I) + 1)*. For CGA-derived data tick *Include precursor mass as feature* in box *7)*. The mass-feature option of the training step (§2.7) can only use the mass if the trainable CSV carries the `protonatedmass` column. Otherwise accept the defaults for your first run.

Click **Build Trainable CSV**. The application:

1. Reads the pseudolabel TSV and keeps the rows that pass the *precursor gate*, whose Score A is at least *min ion score* and whose `ppm_error` is within *max |ppm_error|*. From these it keeps *Top N per scan* candidates per scan (default 1). The candidates of a scan share one Score A, so the first one in the TSV is kept, which is the one with the smallest precursor error (§4.1). If the TSV carries Score B columns, the row that Score B selected (lowest `score_b_rank`) is used for its scan, and Score A ranking is used for the scans Score B did not select.
2. Optionally adds low-evidence scans as Non-glycan negatives. A scan becomes a candidate when it was not selected as a glycan and has fewer ion-list hits than the *Keep scans with < min hits to ion list* value (default 3). The candidates are shuffled with the seed and capped at *Max ratio (neg:pos)* times the number of glycan rows.
3. Builds the feature matrix. Each spectrum becomes one row; each reference fragment in the ion list becomes one column. Cell values are `i_feature,n = log10(I_n) + 1` where `I_n` is the maximum intensity within the ion-matching tolerance of fragment `n`'s mass (default 10 ppm for training-set feature encoding), or `1.0` if no peak matched. The transformation compresses the dynamic range while keeping zero-intensity hits distinguishable from missing matches.
4. Writes the **trainable CSV** `{sample name}_trainable_fromCGA_{YYYYmmdd_HHMMSS}_withMass.csv` (or `_noMass.csv` when box *7)* is unticked) in the folder of the pseudolabel TSV. Its columns are `MS2scan_no`, `Structure` (the composition label, or the literal string `Non-glycan` for negatives), `Glycanannotation2`, `GlyToucan ID`, `WURCS`, then `protonatedmass` when box *7)* is ticked, then one column per ion-list mass, named by the mass. The three label columns are empty unless you fill them. With the fill option on (§2.3), GlyToucan ID and WURCS cells of the trainable CSV are filled from the reference database (§3.7).

**The thresholds box: two precursor filters and the feature window.** Box *4) Thresholds / selection* holds every tolerance this step uses. They are all in ppm, and none of them follows the CGA Setup tolerances:

| Field | Default | What it does |
|---|---|---|
| *min ion score* | 0.07 | Score A pass mark (§4.1). |
| *max \|ppm_error\|* | 20 | **Precursor filter 1.** Drops candidates whose `ppm_error` in the pseudolabel TSV is larger. `ppm_error` is the measured-vs-theoretical precursor error, relative to the theoretical mass. |
| *precursor gate* checkbox and ppm entry | on, 10 | **Precursor filter 2.** Recomputes the measured-vs-theoretical precursor error from the scan CSV, relative to the theoretical mass, and drops rows beyond the value. Untick it to skip this filter. **Field added in GlycoMSP 1.14:** v1.10 applied this gate at 10 ppm, with no field to change it. |
| *Top N per scan* | 1 | Candidates kept per scan. |
| *ion list ppm* | 10 | **MS2 feature window**, not a precursor filter. A peak within this many ppm of an ion-list mass gives that ion's feature value (step 3). This modal's Non-glycan sampling uses the same window. |

A row must pass both precursor filters, so the stricter one decides. With the defaults that is the 10 ppm gate: CGA analysis matches precursors at 20 ppm, but CGA → Trainable keeps only rows within 10 ppm. Unticking the gate leaves *max |ppm_error|* active.

The *5) Negative sampling* box also has its own settings (on, max ratio 3.0, scans with fewer than 3 ion-list hits). The shared Negative Sampling Setting dialog (§3.4) does not apply to this step.

**If you widened the CGA tolerances** (§2.4 (d)), widen these fields to match. Otherwise CGA → Trainable drops what the wider CGA run found:

- Set *max |ppm_error|* and the *precursor gate* at least as wide as the precursor match tolerance. For a Da window, convert at the lowest precursor mass you expect: 0.5 Da at m/z 1,000 is 500 ppm.
- The *ion list ppm* feature window is in ppm only. A 0.1 Da fragment error is 200 ppm at m/z 500 but 100 ppm at m/z 1,000, so choose a value that covers your fragment errors over the ion list's mass range.
- Keep *ion list ppm* at 500 ppm or less. That is the largest value the unlabeled-dataset dialog accepts (§2.8), so prediction features can then be built with the same window.

The values used are recorded with the trainable CSV in the method JSON (§5.2). At training they are copied into `training_run.json` when GlycoMSP can find that record (§4.5). At prediction, values that differ from the defaults prefill the precursor gate and the feature window (§2.8).

When the run completes, a dialog titled *Done* shows the path of the trainable CSV with its number of rows and columns and the count of each class. The path is also shown next to *Output:*, recorded in the sample, and written into the per-sample method JSON.

You now have one trainable CSV per processed sample, ready for ML training. If you have multiple samples, you can train per-sample or combine them; combining is covered in Step 2.7.

## 2.7 Train a Random Forest model

Close the Prepare Dataset window (your work is auto-saved into the per-sample method JSON if a method folder is set, or held in the in-memory experiment state otherwise; see Part 5.3 for details). Back at the main window, click **Run ML Analysis**. The window opens on the **Train Model** tab.

> **[SCREENSHOT: `GlycoMSP_screenshots/Run_ML_analysis.jpg` — ML Analysis window, Train Model tab.]**
> Title bar: *ML Analysis*. Two tabs at the top: *Train Model* (active) and *Predict*. Step-numbered form below: *Step 1: Load Trainable Dataset (.csv)* with a *Select CSV File* button, plus two checkboxes — *CGA dataset?* and *Use protonated mass as a model feature (Non-glycan = 0)*. *Step 2: Select Label Column* combobox (default `Structure`). *Step 3: Choose Classifier* combobox (default *Random Forest*). *Step 4: Train/Test Parameters* — *Set Parameters / Train the Model* button + *Params:* status line (default `(using built-ins)`). Below: a large blue *Train Model* button. Below that: a *Trainable File Info (Origin Tracking)* read-only panel and a *Combine Datasets for Training* row with *Select Datasets* button. *Close* button at the bottom.

**Step 1 — Load.** Click **Select CSV File** and pick the trainable CSV from Step 2.6. The *Trainable File Info* panel shows the file name, the path and, when the folder holds an `.exp.json` file, the name of that experiment file. Leave **CGA dataset?** as it is. The training code of this version does not read that checkbox. Tick **Use protonated mass as a model feature (Non-glycan = 0)** to add the precursor mass as an input feature (recommended for CGA-derived datasets — mass agreement is part of how the labels were assigned, so the model benefits from seeing it explicitly; Non-glycan rows have their precursor mass set to zero). The trainable CSV must contain the `protonatedmass` column for this option to take effect, which is the case when box *7) Include mass as feature* was ticked in Step 2.6.

**Step 2 — Label column.** Leave at `Structure` unless your CSV uses a different column name. The list also offers `IUPACname(optional)`, `Glycanannotation2` and `GlyToucan ID`.

**Step 3 — Classifier.** Leave at *Random Forest*. The only other entry, *Coming soon*, is a placeholder. GlycoMSP trains a single Random Forest pipeline.

**Step 4 — Parameters.** Click **Set Parameters / Train the Model**. This opens the **Train/Test Parameters** dialog, where you set the train/validation/test split, class-balancing options, and Random Forest hyperparameters.

> **[SCREENSHOT: `GlycoMSP_screenshots/train_parameters.jpg` — Train/Test Parameters dialog. (Optional companion: `trainlabel.jpg` for label/classifier selection.)]**
> Dialog titled *Train/Test Parameters*. Sliders *Test split* (default 0.20) and *Validation split* (default 0.00) with the note *Note: Test + Validation must be < 1.0* and a live readout (*Train 80% | Val 0% | Test 20%*). A spinbox *Minimum samples per class* (default 5). A checkbox *Enable class balancing (cap majority)* (ticked) with the fields *Majority label:* (`Non-glycan`) and *Majority factor (×max minor):* (3). A checkbox *Stratify by label* (ticked). The block *Random Forest Options* with *n_estimators (trees):* (400), a checkbox *Use class_weight = "balanced"* (ticked) and a checkbox *Real-world test (cap training only)* (ticked). The block *Prediction Thresholding* with the checkbox *Enable confidence threshold → fallback to Majority label* (**unticked by default**), *Threshold τ (0.00–0.99):* (0.65) and *Majority support margin δ (0.00–0.20):* (0.05). Buttons *Edit ML Parameters…*, *Force current params* and *Close*.

For your first run, accept the defaults. They are unchanged since v1.10: an 80/20 train/test split (`test_split = 0.20`, `val_split = 0.00`; no validation fold by default), stratified split, balanced class weights, 400 trees in the Random Forest (`max_depth=None`, `min_samples_split=2`, `min_samples_leaf=1`, `random_state=42`), a class-eligibility filter that drops compositions with fewer than 5 spectra, and the Non-glycan class capped at 3× the largest glycan class on the training fold only (the held-out test is left uncapped). The optional confidence-threshold step (`τ = 0.65`, `margin = 0.05`) is **disabled by default**. When enabled, it relabels weak predictions as *Non-glycan* at prediction time (§2.8). Close the dialog with **Close**. The settings apply as soon as you change them, and the *Params:* line in the Train Model tab changes from `(using built-ins)` to a short summary such as `Params: RF Trees: 400, Split test=0.2, val=0.0, min_class_size=5, real_world=True, balance=True`.

> **Note — fixed (non-configurable) RF parameters.** The dialog exposes only `n_estimators` and `class_weight`. `max_depth=None`, `min_samples_split=2`, `min_samples_leaf=1` and `random_state=42` are the defaults of the pipeline and have no field in the dialog. The button *Edit ML Parameters…* opens a separate parameter editor that the Quickstart does not use.

**Train.** Click the large blue **Train Model** button. The pipeline:

1. Reads the trainable CSV and normalises the composition labels of the label column.
2. Drops classes with fewer than *Minimum samples per class* rows, caps the Non-glycan class according to your balancing settings, and splits into train / validation / test with stratification.
3. If *Use protonated mass as a model feature* was ticked, sets `protonatedmass = 0` for Non-glycan rows (string-label match) before label encoding. Otherwise it removes the column.
4. Drops administrative columns (`uid`, `origin_*`, `predicted_label` if present) to produce the numeric feature matrix.
5. Fits a `sklearn.RandomForestClassifier` on the training fold.
6. Predicts the test fold, applies the confidence-threshold rule when it is enabled, and logs a summary of glycan versus Non-glycan classification (precision, recall, F1 and accuracy). A per-class report is written to a file (see item 7).
7. Saves the trained model as `{dataset stem}_rf_model.skops` next to the trainable CSV, or as `{dataset stem}_withMass_rf_model.skops` when the mass feature was ticked. Sidecar files with the same stem hold the label encoder (`_labelencoder.joblib`), the feature column order (`_features.json`) and the per-class report (`_rf_performance.txt`). The reproducibility manifest `training_run.json` in the same folder captures the effective parameters, column order, class list, input file hashes, and library versions.

When training finishes, a dialog titled *Training Complete* summarises the run: the split mode, the number of trees and class weight, whether thresholding is on, the number of classes kept, the majority cap that was applied and the number of classes dropped for being too small. If you provided multiple trainable CSVs via *Combine Datasets for Training*, the application concatenates them with row-level UIDs preserved before training (useful for tracing which spectrum came from which sample after the fact).

## 2.8 Run prediction on unlabeled data

Switch to the **Predict** tab.

> **[SCREENSHOT: `GlycoMSP_screenshots/Run_ML_analysis_predict.jpg` — ML Analysis window, Predict tab.]**
> Same window, *Predict* tab active. *Step 1: Load Trained Model (.skops/.joblib/.pkl)* with *Select Model File* button. *Step 2: Prepare Unlabeled Input Dataset* with two buttons stacked vertically — *Create unlabeled dataset of certain experiment* (uses an experiment JSON to assemble inputs) and *Load unlabeled dataset* (loads a pre-built unlabeled CSV directly). Below: a labelled panel *Precursor (MS1) Gate for Predictions* containing — *Apply precursor gate* checkbox (enabled by default), *Tolerance (ppm):* entry (default `10`), *In-silico CSV (composition ↔ theoretical_mass):* entry + *Select in-silico CSV file* button. Below the gate panel: a green *Run Prediction* button. *Close* button at the bottom.

**Step 1 — Model.** Click **Select Model File** and pick the `.skops` file from Step 2.7. The application loads the model, the label encoder, and the feature column list. (For older models saved as `.joblib` or `.pkl`, those formats also load; the path of the sidecar files is derived from the model path with the appropriate extension.)

**Prefill from the model (added in GlycoMSP 1.14).** GlycoMSP looks for the `cga_tolerances` record in the `training_run.json` next to the model (§4.5). If the CGA → Trainable fields recorded there differ from the defaults, loading the model prefills the matching prediction settings, and the dialog titled *Model Loaded* lists what was prefilled:

- A non-default CGA → Trainable precursor gate becomes the *Tolerance (ppm)* below. If the gate was unticked at CGA → Trainable, *Apply precursor gate* is unticked too.
- A non-default *ion list ppm* becomes the suggested value in the unlabeled-dataset dialog (Step 2).

A model trained at the defaults, a MAS model, or a model without the record leaves everything as in v1.10. A precursor-gate field you have already changed yourself is kept (the dialog says so), also when you load another model.

**Step 2 — Unlabeled input.** You have two options. **Load unlabeled dataset** lets you pick a pre-built unlabeled CSV directly — useful if you already prepared the input. **Create unlabeled dataset of certain experiment** is the integrated workflow: pick an experiment JSON, and the application traverses each sample's method JSON, calls `extract_fragment_masses` to build the column order, calls `extract_ion_intensities` to compute features for each spectrum, and writes a `<csv_stem>_unlabeled_ppm{ppm}.csv` next to each sample's converted CSV. A `UID` column of the form `{experiment hash}:{sample hash}:{six-digit scan number}` is added so predictions remain traceable to source spectra (§3.6).

The integrated workflow first asks for the feature window in ppm (dialog *PPM Tolerance*, prompt *Enter PPM tolerance (e.g., 20):*), the prediction-side counterpart of the trainable CSV's *ion list ppm* (§2.6). The field has no fixed value. The 20 in the prompt is only an example, so enter the value you want. GlycoMSP 1.14 accepts 1–500 ppm with decimals; v1.10 accepted whole numbers from 1 to 200. When the model was trained with a non-default feature window, the dialog is prefilled with it. If you enter a different value, a notice titled *Tolerance differs from the model* appears and is written to the log panel.

**Precursor gate.** The *Apply precursor gate* checkbox (enabled by default) tells the predictor to compare each spectrum's measured precursor mass against the theoretical mass of the predicted composition. The comparison covers glycan predictions only. If the deviation exceeds the *Tolerance (ppm)* setting (default 10 ppm), the prediction does not pass the gate. It is listed in `*_gate_debug.csv` and marked `pred_ok` = False in the full table (§2.9). Provide an in-silico CSV via the *Select in-silico CSV file* button; this file maps composition labels to theoretical masses (the same in-silico CSV from Step 2.4 works). If the gate is on and no in-silico CSV is selected, the gate is skipped and the log panel says so.

When the model has a non-default precursor-gate record and the gate you run with differs from it, a *Tolerance differs from the model* notice appears and is written to the log panel. `run_summary.json` records what was requested, what was applied and whether it matches the model (§5.4). The output folder name shows the gate as a whole number of ppm, rounded down (`_PG12ppm` for 12.5 ppm, as before).

**Run.** Click the green **Run Prediction** button. The pipeline:

1. Loads the unlabeled feature matrix and aligns its columns to the model's training column order. Missing training-feature columns are filled with `0` and extra columns are dropped. The predictor does not encode features again. The features of a CSV made with *Create unlabeled dataset of certain experiment* were encoded at the ppm value you entered in the *PPM Tolerance* dialog, and a CSV you load with *Load unlabeled dataset* keeps the window it was built with. For a model trained at the default 10 ppm feature window, entering 10 gives features encoded like the training set.
2. Calls `model.predict_proba` to obtain class probabilities and takes the model's predicted class as the prediction. The prediction reporter later computes the confidence as the top class probability.
3. Applies the confidence threshold, if it is enabled (it is off by default, see the *Prediction Thresholding* block of the Train/Test Parameters dialog in §2.7). A prediction of a glycan class whose top probability is below `τ`, and whose *Non-glycan* probability is within `margin` of the top probability, is relabelled as the majority label, which is *Non-glycan* by default.
4. Applies the precursor gate (if enabled): glycan predictions whose theoretical precursor mass deviates beyond the tolerance are marked as not passed.
5. Calls the prediction reporter, which writes `predictions_with_conf.csv`, `class_summary.csv`, `by_sample_summary.csv` (when the input has sample columns), `run_summary.json`, `run_summary.md` and the review queues `to_review_*.csv`. A window titled *Prediction summary* then shows the total spectra predicted, the accepted and rejected counts and a table of the top classes. A dialog titled *Prediction Complete* shows the output folder and the path of the main table.

## 2.9 Read your outputs

Open the output folder named for the prediction run. GlycoMSP creates it next to your unlabeled input CSV. Its name is `{method base}` plus `_MS1feat` or `_noMS1` (whether the model uses the precursor mass) plus `_PG{n}ppm` or `_noPG` (the precursor gate). You will find:

**Recommended for reviewing results:**

- **`class_summary.csv`** — one row per predicted class, with spectra counts per composition (`pred_label`, `pred_count`). This is your headline result for profiling.
- **`predictions_with_conf.csv`** — one row per predicted spectrum, with the columns of the unlabeled input (`MS2scan_no`, `protonatedmass`, `UID` when the input has it, and the features), the predicted composition `pred_label`, and the probability columns (`max_proba`, `second_label`, `second_proba`, `margin`, `entropy`, `accepted`, `reject_reason`, `topk_labels`, `topk_probas`). It also has the columns `ppm_precursor` and `pred_ok`, which hold values only when the gate is applied. Use it to trace any individual composition back to its scans.
- **`{input stem}_predicted{suffix}.csv`** — the complete table of all predicted spectra, with `Predicted_Label` and the columns `ppm_precursor` and `pred_ok`.
- **Precursor-gate outputs** (`*_gate_passed.csv`, `*_gate_debug.csv`, `*_PG_counts.csv`, `*_PG{n}ppm_only.csv`) — show which predictions agree with the measured precursor mass. They are the practical quality filter for composition-level results. `*_gate_debug.csv` lists the glycan predictions that did not pass.
- **`run_summary.json` / `run_summary.md`** — reproducibility metadata (model path, input path, settings, version, output inventory) plus per-run statistics.

When the precursor gate is applied, the reporter files (`predictions_with_conf.csv`, `class_summary.csv`, `by_sample_summary.csv`, `run_summary.*`, `to_review_*.csv`) are computed from the glycan predictions that passed the gate. *Non-glycan* predictions have no theoretical mass to compare and are not in these files. They are in the complete table.

With the fill option on (§2.3), the prediction tables also carry the columns `Predicted_GlyToucan_ID` and `Predicted_WURCS` next to `Predicted_Label`, filled from the reference database (added in GlycoMSP 1.11, §3.7).

**Not recommended at laboratory scale** (`to_review_low_conf.csv`, `to_review_low_margin.csv`, `to_review_ambiguous_topk.csv`): these confidence/margin review queues are a residual feature retained for future use with much larger, well-balanced datasets. At typical per-class spectra counts the predicted-probability and margin values are not well calibrated, so these queues are not informative. With the confidence threshold off (the default), `to_review_low_margin.csv` lists almost every spectrum. Tuning them meaningfully would require ground-truth labels and hyperparameter-style optimization that small datasets cannot support; we therefore do not recommend interpreting them for routine analysis.

**Known limitation.** Whatever the threshold setting, the run summary (`run_summary.json`, `run_summary.md` and the *Prediction summary* window) reports *Accepted* 0 and probability statistics of 0.000 (known issue 29). With the threshold off it also reports τ 0. The predicted labels and the per-class counts are correct (§7.2).

For composition-level profiling, read `class_summary.csv` together with the precursor-gate outputs; for any spectrum that matters individually (e.g., entries you intend to publish), trace it through `predictions_with_conf.csv` and inspect the spectrum itself.

That's the end-to-end CGA pipeline. You started with raw instrument files and ended with predictions on unlabeled spectra plus reviewer-friendly reports.

---

# Part 3 — Workflow guides

This Part covers the workflows that branch off from the Quickstart's CGA spine. Each subsection is independent — read the ones relevant to your data.

> **Shared across both routes (MAS and CGA).** Several steps are route-agnostic: they are configured once in the **Prepare Dataset** window and then applied to whichever labeling route you run. Specifically — **negative sampling** (§3.4, applied during *MAS Merge* and, for the ion-mining step, during *Start CGA analysis*), **ion mining** (§3.5, applied during *MAS Merge* and *Start CGA analysis*), the **selected-fragment-ion feature encoding** (§4.3), and the **ML training/prediction** stage (§2.7–2.9). The *CGA → Trainable* window has its own negative-sampling box (§2.6). In other words, MAS and CGA differ only in how labels are produced; everything downstream is the same pipeline.

## 3.1 The MAS pipeline with expert annotations

Use the MAS pipeline (Manual annotation support) when you have **expert-validated annotations** linking specific MS2 scans to specific glycan compositions. MAS skips the in-silico library and constraint scoring entirely; instead, it ingests your annotation Excel as ground truth, validates that the labels are consistent with the spectra, and passes everything through the same feature-extraction and ML stages as CGA.

**Inputs.** Beyond the converted scan CSV and metadata JSON (Step 2.2 of the Quickstart), MAS needs an annotation Excel with two sheets:

- **`MSlist` sheet** — one row per annotated scan. Required columns include the scan identifier and the assigned glycan composition label. Additional metadata columns (RT, intensity, manual notes) are preserved but not consumed.
- **`ionlist` sheet** — reference fragment ions for feature extraction. Must contain a numeric column named `mass`. **Link Sample** also accepts the header `fragmentation_mass` (used in Score B workbooks), but **Merge Sample** reads only a column named `mass`. Rename the column before you merge.

**Walkthrough.** From the Prepare Dataset window with your scan CSV(s) added per Quickstart Step 2.3:

1. Click **Add Excels (ion list/man annotation)** and select your annotation Excel. The application appends the file to the corresponding sample's `excel` slot.
2. Make sure each sample has its metadata attached (visible in the tree as a *Metadata: …* sub-row). If it has none, **Link Sample** asks "No metadata found. Would you like to create it now?" in a dialog titled *Metadata Missing*, once the other checks have passed. Answer *Yes* and GlycoMSP asks for an optional RAW file and an output folder, then opens the Metadata Editor without a conversion step. When the metadata is saved, GlycoMSP validates the sample again.
3. With the sample selected, click the **Link Sample** button in the action grid. The application runs six fail-closed validation stages: sample lookup, file presence, CSV structure, Excel sheet existence (the `MSlist` sheet), Excel content (validates the `MSlist` sheet structure and ensures the `ionlist` sheet has a numeric mass column), and metadata gate. On success, the sample's tree status changes to ✅ and a per-sample method JSON of family `MAS` is written. If a stage fails, an error dialog names the problem and the sample shows ⛔ (validation failed). A sample that has files but has not been validated shows ⚠️ (unvalidated), and a sample without metadata shows ❌ (metadata missing). A missing or non-numeric mass column in the `ionlist` sheet is reported in the log panel only, and the sample is not marked as validated.
4. With the sample now validated, click **Merge Sample**. The button is greyed out until validation succeeds. The application asks you for an output folder, reads `Derivatization Type` from the metadata (and prompts you for a custom mass shift if the value is `Others`), then builds a merged dataset by joining the converted CSV's spectra with the annotation Excel's labels.
5. If *Add Non-glycan entries* is enabled in **Negative Sampling Setting** (Section 3.4), low-evidence scans not annotated in the Excel are sampled as Non-glycan negatives at this stage.
6. The merged CSV is written to your chosen folder as `{sample_name}_merged_{timestamp}.csv`. The path is recorded in the sample's `files["trainable_csv"]` slot, and a method JSON of family `MAS` is refreshed.

**Output schema.** The merged CSV has `protonatedmass`, one feature column per ion-list mass, `Structure` (the expert-assigned label, or `Non-glycan` for sampled negatives), `IUPACname(optional)`, `Glycanannotation2`, `GlyToucan ID`, `WURCS` and `unique_ID` (the scan number). The CGA trainable CSV carries the scan number as `MS2scan_no` instead. Train and predict from this point onward exactly as in Quickstart Steps 2.7–2.8.

**Masses and the reference database.** The MAS merge computes precursor masses from permethylated residue masses, also when the metadata says `Native` or `Native(Reduced)` (§3.8). The KDN residue mass is 320.1471, the same value the CGA route uses (changed in GlycoMSP 1.12). Scan rows whose precursor charge is unknown keep charge 0 and mass 0 in the converted CSV and are not excluded downstream (§4.9). When *Fill GlyToucan ID / WURCS from reference DB (see Manage References)* is ticked in Prepare Dataset, the merge fills empty `GlyToucan ID` and `WURCS` cells from the local reference database, and values already in your sheet always win (added in GlycoMSP 1.11, §3.7).

**Differences worth noting compared with CGA.** MAS labels are not scored — the application trusts the expert. The application *does* validate the structure of the `MSlist` sheet and that the ion-list mass column is numeric. At merge, the annotated scans are joined to the converted CSV by `MS2scan_no`, so annotated scans that are not in the converted CSV are dropped without a message. The application does not double-check the assigned compositions. Score B is not applicable to MAS rows because Score B operates on candidate compositions; MAS rows have a single expert-assigned composition with no candidate set.

## 3.2 Combining trainable datasets (and when not to)

GlycoMSP can combine multiple trainable CSVs into one training pool via **Combine Datasets for Training** in the ML Analysis window's Train tab. Click **Select Datasets**, pick the CSVs (Cmd-/Ctrl-click for multi-select), and the application concatenates them and keeps the union of their columns. It keeps an existing `UID` column and otherwise builds one as `{file name without extension}#scan{scan number}`. It also adds `Origin_Basename` and `Origin_File` columns for traceability. A save dialog then proposes the name `combined_trainable_{timestamp}.csv`, and a `_combine_log.txt` is written next to the result. Load the combined file with *Select CSV File* in *Step 1*.

**The one hard requirement: the same fragment-ion list.** Trainable datasets are only safely mergeable when they were built from the **same fragment-ion (feature) list**, so that their feature columns mean the same thing. Combining datasets built from different ion lists produces a misaligned feature space and the model training will not identify features correctly (or will fail). This is the only compatibility GlycoMSP can guarantee at the file level.

**Scientific compatibility GlycoMSP does *not* enforce — your responsibility.** Even with a shared ion list, only merge datasets acquired under compatible conditions:

- **Same ionization mode** (positive vs. negative) and **same derivatization** (do not mix native and permethylated glycans — they fragment differently). Mixing these is improper and, because the feature ions differ, will generally fail at training unless you deliberately prepared a single ion list valid for both and built every dataset with it.
- **Different biological treatments** of the same system are fine to group.
- For building a comprehensive glycan database, even **different cell types or species** may be pooled — *provided the sample-preparation protocol and the spectral-acquisition workflow are identical*.

**Mixing CGA-derived with MAS-derived rows is generally not recommended.** CGA labels are simulated/constraint-scored pseudolabels and are less discriminative than expert MAS labels; mixing the two typically *reduces* performance rather than helping. In addition, if the MAS data uses customized (non-composition-string) labels, combining it with CGA output — which uses composition strings — leads to inconsistent label spaces and training errors. Mix the two routes only if you have a specific reason and have confirmed label and feature compatibility. The **CGA dataset?** checkbox in *Step 1* has no effect in this version, so it does not change training behavior or reconcile any of the above.

## 3.3 Score B enrichment — motif-aware re-ranking

Score B is a second-stage scoring pass that re-ranks CGA candidates using user-defined glycotopic motif rules. The manuscript Methods describes Score B as a *motif-aware re-ranking step that refines glycan composition assignments generated by CGA*. Use it when your samples include compositions that are isobaric (identical in mass, so Score A cannot discriminate them) but distinguishable in fragmentation pattern. For example, in permethylated glycans a Neu5Gc→Neu5Ac change (−CH₂O) combined with a Fuc→Hex change (+CH₂O) cancels out, so a `…G…F…`-containing composition is mass-degenerate with a `…S…H…` one — same precursor mass, different composition.

By default Score B is **inactive**. Enabling it requires a Score B workbook authored externally. The manuscript (Supplementary Data 2 in `manuscript/GlycoMSPdraft/`) provides example workbooks and authoring guidance; this manual covers only how to attach and consume the workbook.

**Workbook structure.** The Score B workbook is an Excel file (`.xlsx`, `.xlsm`, `.xltx`, or `.xltm`) that must contain six sheets. `ionlist` defines the diagnostic fragment ions. `interpretation_v2` maps observed ions to glyco-feature targets with a logic, a minimum number of hits and a weight. `TargetHierarchy` holds the parent-to-child relations between targets. `MotifPolicy` selects which targets are relevant for each candidate flag. `CompositionConsistency` holds penalties for evidence that is inconsistent with the candidate composition. `UnexpectedEvidencePenalty` holds penalties for evidence that does not fit an expected motif. Other sheets in the file are ignored. The exact column layouts are documented in `score_b_config_instructions.md` and in `docs/score_b_design.md`.

**Attaching the workbook.** From the CGA Setup window (Quickstart Step 2.4):

1. Click **Add Score B Workbook** and pick the Excel file. The path appears next to *Score B reference:*.
2. Click **Validate Score B Workbook**. The application runs a load test against `msp_CGA_structscore.ScoreBLoader` and pops up either a confirmation dialog titled *Score B Workbook* ("Workbook validated successfully.") or an error dialog with the reason, such as the names of missing required sheets. Validation is non-destructive — it only reads.
3. Optionally, click **Clear Score B Workbook** to detach. Detaching is the same as having never attached.

Then proceed with **Start CGA analysis** as in Quickstart Step 2.5. **Start CGA analysis** repeats the validation when a workbook is attached and stops with an error dialog if the workbook is not valid.

**Output columns when Score B is active.** The pseudolabel TSV gains nine columns:

- **`score_b`** — float in `[0, 1]`, the final Score B value for this candidate.
- **`score_b_rank`** — per-scan rank by Score B (1 = highest).
- **`score_b_selected`** — boolean, `true` only for the top-ranked candidate per scan that has a non-empty `selected_composition`.
- **`score_b_motif_summary`** — short text summary of which motifs contributed and which penalties applied.
- **`score_b_support_score`**, **`score_b_composition_penalty`** and **`score_b_unexpected_penalty`** — the three terms of the formula below for this candidate.
- **`selection_source`** — `score_b` on every candidate row that Score B scored.
- **`selected_composition`** — the composition Score B chose for the scan, repeated on every scored candidate row of that scan.

Plus the formula behind `score_b`:

`b = max(0, min(1, support − penalty_comp − penalty_unexpected))`

where *support* is motif-level evidence aggregated across the workbook's glycotopes, *penalty_comp* penalises evidence inconsistent with the candidate's composition, and *penalty_unexpected* penalises evidence not accounted for by any expected motif. The clamp ensures all candidates remain on a comparable `[0, 1]` scale.

**How composition labels are read (changed in GlycoMSP 1.12).** Score B reads the `composition` label case-sensitively and scores the residues `F`, `H`, `N`, `S`, `G` and `KDN` only. Lowercase `s` (sulfate), lowercase `p` (phosphate) and `A` (HexA) are not Score B residues. They are left out of scoring, and the labels that carry them are reported once at the end of the run in one warning line on the console. Score B matches fragments with a fixed 20 ppm window, which is independent of the *Score A fragment match* field in the *Mass tolerances (CGA matching)* frame of CGA Setup (§4.10).

**Selection rule in CGA → Trainable.** When you build the trainable CSV (Quickstart Step 2.6), the application detects the presence of `score_b_selected` and switches selection logic: per scan, the row with `score_b_selected = true` is taken as the assigned label. For scans where no candidate was Score B-selected (typically because no motif evidence was strong enough), the workflow reverts to Score A selection (§4.1). The manuscript names this *"Score B preferred, Score A fallback"*. Candidates must first pass the *min ion score* and *max |ppm_error|* filters of box *4) Thresholds / selection* (§2.6). If no row passes, the build stops with a message.

**Failure mode.** Score B execution is non-fatal. If `enrich_cga_tsv_with_score_b` raises any exception (workbook-load error, malformed rule, runtime issue), the application logs a warning, shows a warning dialog titled *Score B*, and returns the original Score-A-only TSV unchanged. Subsequent CGA → Trainable runs against that TSV simply use Score A selection. You will not lose your CGA work to a Score B failure.

**When motif-aware re-ranking helps.** Score B improves discrimination for isobaric same-mass candidates, particularly compositions that share monosaccharide counts but differ in glyco-feature expression (e.g., Lewis antigens versus their sialylated counterparts). It does not help when no motif evidence exists (the workbook rules don't match any peaks for the spectra in question), in which case the Score A fallback applies and Score B is effectively a no-op.

## 3.4 Negative sampling options

The MAS merge can sample low-evidence scans as Non-glycan negatives, providing the ML model with examples of the "not a glycan" class. The configuration is set in one dialog: **Negative Sampling Setting** in the Prepare Dataset window's action grid. The *CGA → Trainable* window has its own box *5) Negative sampling* (§2.6) and does not use this dialog.

> **[SCREENSHOT: `GlycoMSP_screenshots/Negative_Samplng_smoketest4_1.jpg` — Negative Sampling Options dialog, opened by the Negative Sampling Setting button.]**
> Modal dialog titled *Negative Sampling Options*, shown over the Prepare Dataset window. *Add Non-glycan entries (easy negatives)* checkbox (ticked in the screenshot and unticked by default). Below: *Max ratio (neg:pos):* Spinbox (range 0.0–10.0, step 0.5, **default 3.0**). *Keep scans with < min hits to ion list:* Spinbox (range 0–50, step 1, **default 3**). *Ion match tolerance (ppm):* Entry with input validation (**default 10**). *Load ion list now…* button (manually re-probes the feature count). Below: live hint label reading *"Gating non-glycan: ion matches < K hits (~P% ; K of N)"* where K is the gate value, P is the percentage of total features, and N is the auto-probed feature count from the current sample's ion list (the screenshot shows *ion matches < 3 hits (~1.1% ; 3 of 282)*). When N is unknown the hint ends with *(feature count unknown yet)*. *Close* button.

The dialog **does not run any workflow itself** — it only sets shared state. The next time you run **Merge Sample** (MAS), or **Start CGA analysis** with ion mining enabled (the negatives are then used for the ion suggestions only), the consumer reads these values.

**What the gate does.** A scan qualifies as a Non-glycan negative only if its ion-hit count against the reference ion list is **below** the gate (*Keep scans with < min hits to ion list:*). So if you set the gate to 3, scans matching 0–2 ion-list fragments within *Ion match tolerance (ppm):* are eligible negatives, and scans matching 3 or more are skipped. The *Max ratio (neg:pos):* field then caps the total negatives sampled at *ratio × n_positives*. Default settings (**ratio 3.0, gate 3, ppm 10**) cap negatives at three times the number of positives, sampled from non-positive scans matching fewer than 3 ion-list fragments.

**Recommended starting point.** Start with the defaults. If your model overfits to the positive class (test recall on Non-glycan is poor), raise *Max ratio (neg:pos):* above the default of 3.0, up to 10.0, to give the model more negative examples. If your negatives leak into the positive class (the model predicts Non-glycan for clear glycan spectra), lower *Keep scans with < min hits to ion list:* so only scans with very few ion-list matches qualify as negatives.

**Cross-pipeline state.** The dialog's settings affect the MAS merge and the negatives used for ion mining in CGA analysis. They do not affect *CGA → Trainable*. That window keeps its own box *5) Negative sampling* with its own defaults: *Add Non-glycan entries (easy negatives)* ticked, *Max ratio (neg:pos)* 3.0, *Keep scans with < min hits to ion list* 3, *Shuffle negatives (random)* ticked and *Seed:* 42. It matches negatives with its *ion list ppm* value (§2.6). Setting the dialog once and leaving it is fine.

## 3.5 Ion mining — discovering diagnostic fragments

Ion Mining inspects your data for fragment masses that distinguish glycan-class spectra from negative spectra and would not be in a generic ion list. The output is a ranked list of suggested m/z values you can review and optionally add to your ion list for future runs.

Open the configuration via **Ion Mining** in the Prepare Dataset action grid.

> **[SCREENSHOT: `GlycoMSP_screenshots/ion_mining.jpg` — Ion Mining dialog.]**
> Modal dialog titled *Ion Mining*, 5-row grid. *Generate ion suggestions on merge* checkbox (ticked in the screenshot and unticked by default, and it also enables mining in *Start CGA analysis*). Below: *ppm tolerance:* Spinbox (1.0–50.0, step 0.5, default 10.0). *DA floor (low-m/z merge):* Spinbox (0.00–0.10, step 0.005, format `%.3f`, default 0.030). *Min glycan support per ion:* Spinbox (1–100, default 5). *Top-K suggestions:* Spinbox (5–500, step 5, default 60). *Close* button.

Like the Negative Sampling dialog, this one only sets shared state. The actual mining runs during dataset building for **either** route when the enable checkbox is ticked — specifically during **MAS Merge** and during **Start CGA analysis** (it is not part of the CGA → Trainable step). The checkbox is unticked by default.

**How mining works at a high level.** The application takes your current spectra and assigns each one a binary label (`Glycan` for annotated/positive scans, `Non-glycan` for sampled negatives — the same negatives configured in Section 3.4). It then groups m/z peaks from the Glycan spectra at the configured PPM tolerance with the small-mass merge floor (so, e.g., peaks at 204.085 and 204.087 collapse to one cluster) and keeps clusters of at least three peaks. It counts how many Glycan vs Non-glycan spectra contain each cluster and computes statistics: support count per class, lift (ratio of class support to baseline support), odds ratio, and a `recommended` flag. The flag is `True` when the cluster is not already in your ion list. Clusters with fewer Glycan spectra than *Min glycan support per ion:* are dropped, and the top-K clusters by lift are exported. In *Start CGA analysis* the negatives for this step are sampled with the dialog of Section 3.4 when *Add Non-glycan entries (easy negatives)* is ticked. They are used for the suggestions only and are not written into the pseudolabel TSV. Without negatives, `support_non` is 0 for every cluster and the lift cannot separate class-specific ions from common ones.

**Reviewing results.** Click **Browse Ion Mining result** in the action grid. The application opens a modal viewer titled *Ion suggestions — {file name}* that shows the top 10 suggestions in a Treeview with columns `mz`, `support_glycan`, `support_non`, `lift`, `odds_ratio`, `already_in_list`, `recommended`. The rows are sorted by lift (descending), then by glycan support. Click a row and use **Copy selected m/z** to copy the value to clipboard for adding to your ion list. **Open folder** opens the OS file manager at the suggestions CSV's parent directory so you can audit the full list outside the GUI. **Close** closes the viewer.

**Output file.** Each enabled mining run writes `{sample_name}_ion_suggestions.csv` into the output folder of the merge (MAS) or into the folder of the converted scan CSV (CGA). The viewer shows the file of the most recent run in the current session. If no run has been made in the session, it asks you to pick a suggestions CSV.

**When to use ion mining.** Two natural cases. First, **early ion-list curation** — if you started with a generic ion list and want to expand it with class-specific diagnostic fragments observed in your own data. Second, **method troubleshooting** — if your model performs poorly, mining can reveal fragments your reference list missed but that consistently distinguish glycans in your dataset. In both cases, the suggestions are not added to the ion list automatically; you decide which to incorporate, edit the ion-list Excel, and re-run.

## 3.6 Building unlabeled datasets across multi-sample experiments

The Predict tab has two ways to load an unlabeled input. **Load unlabeled dataset** picks a single pre-built CSV and is what the Quickstart used. **Create unlabeled dataset of certain experiment** is the integrated multi-sample workflow.

**When to use it.** When you have an experiment containing several samples and want to predict on all of them in one go, with consistent feature spaces and traceable outputs. This is the typical setup for, e.g., a clinical cohort: the experiment is the cohort, the samples are individual patients, and you want one trained model applied uniformly.

**Walkthrough.** Click **Create unlabeled dataset of certain experiment** and pick the experiment JSON (`*.exp.json`) for the cohort. GlycoMSP then asks for the feature window in a dialog titled *PPM Tolerance*, which accepts 1 to 500 ppm with decimals. If you loaded a model first (Step 1) and it was trained with a non-default feature window, the dialog is prefilled with that value, and a notice appears if you enter a different one (§2.8). The application:

1. Reads the experiment JSON and traverses the `samples` array.
2. For each sample, finds the converted CSV and the ion list, either in the sample's own entry or in its first valid method JSON (listed in the sample's `methods` array). Samples with a missing or unreadable file are skipped and named in the final message.
3. Calls `extract_fragment_masses` to read the ion list and determine the column order, then `extract_ion_intensities` to compute the feature matrix for each spectrum (PPM-windowed, log10-transformed, identical to training-time feature extraction).
4. Adds a `UID` column of the form `{experiment hash}:{sample hash}:{scan number}`. Both hashes are the first eight characters of a SHA-1 hash of the experiment title and of the sample name, and the scan number has six digits, so each predicted row is traceable to source experiment, source sample, and source scan number.
5. Writes one `{csv_stem}_unlabeled_ppm{ppm}.csv` per sample, next to that sample's converted CSV.
6. Records the path of each unlabeled CSV as `unlabeled_csv` in the experiment JSON's per-sample blocks.

After this completes, click **Load unlabeled dataset**, pick one of the new files, and click **Run Prediction**. The prediction step reads one unlabeled CSV per run and writes its output folder next to the input. Repeat for each sample.

**A note on consistency.** The unlabeled feature space *must* align with the trained model's feature space. The application enforces this by reordering input columns to match the model's training feature list and filling any missing column with `0.0`, with extra columns dropped. A fragment that was searched for and not found in a spectrum has the value `1.0` in the feature matrix, so a missing column and an unmatched fragment differ. If your experiment uses a different ion list from the one used at training time, predictions will work but may be noisier — every "missing" reference fragment becomes an effective zero. For best results, use the same ion list (or a strict superset) at training and prediction time.

## 3.7 Reference database: build, select, inspect (Manage References)

Use the reference database when you want your annotation tables to carry a GlyTouCan accession and a WURCS string next to each glycan composition. It is a local SQLite file that maps a composition to its accession and WURCS string. It was added in GlycoMSP 1.11. The Manage References window that manages it was added in GlycoMSP 1.12.

The reference database is **off by default**. With the option off, every output is identical to v1.10. Only one component ever contacts the network, the *prebuild*, which fills the database from GlyCosmos. The fill that runs inside GlycoMSP reads the local database only and never calls the network. Once the database is built, your analysis runs work offline.

**Where it lives.** GlycoMSP keeps its per-user data in one folder. On Windows this is `%LOCALAPPDATA%\GlycoMSP`. On macOS and Linux it is `~/.glycomsp`.

| Item | Default location |
|---|---|
| Reference database | `<data folder>/references/glycomsp_reference.db` |
| Transaction log | `<data folder>/glytoucan_cache.jsonl` |
| Settings file | `<data folder>/settings.json` |

The transaction log is a JSONL file. Every successful fetch is appended to it. The database stays the source of truth and the log can be regenerated from it (see *Rebuild log* below).

Two environment variables override the defaults. Set `GLYCOMSP_REFERENCE_DB` to use another database file. Set `GLYCOMSP_GLYTOUCAN_LOG` to use another log file. While `GLYCOMSP_REFERENCE_DB` is set, a database path remembered in `settings.json` is ignored.

GlycoMSP does not create the database when it starts or when you open the Manage References window. The first prebuild creates it.

`settings.json` remembers the fill option and the chosen database between launches. GlycoMSP writes it only when you press *Set as default* or tick *Remember these settings*, and never by browsing. At launch GlycoMSP applies the file only if *Remember these settings* was ticked when it was saved. If the file cannot be read, GlycoMSP keeps the defaults and prints a warning in the console.

### Open the Manage References window

In the **Prepare Dataset** window, the last row of the action grid has a **Manage References…** button. Beside it, an *Active DB:* label shows which database is in use. The fill option sits in the row above as the checkbox *Fill GlyToucan ID / WURCS from reference DB (see Manage References)*. The checkbox is greyed out if the reference module could not be loaded.

The screenshot in §2.3 shows these controls.

The *Active DB:* label shows one of these texts:

- *Active DB: {name} {version} · {n} entries* when the database opens.
- *Active DB: not found ({path})* when the file does not exist.
- *Active DB: unreadable ({path})* when the file exists but cannot be opened as a reference database.
- *Active DB: (updating…)* while a prebuild or a log rebuild is running.
- *Active DB: resolver unavailable* when the reference module did not load.

The label changes when you select another database and when a prebuild or a log rebuild starts or ends.

Click **Manage References…**. The window opens once per session. Clicking the button again brings it to the front.

> **[SCREENSHOT: `GlycoMSP_screenshots/Manage_References.jpg` — Manage References window, no prebuild run yet.]**
> Title *Manage References — GlycoMSP v1.14 publication build*. Top frame *Reference database* with eight rows: *Path*, *Name* (here `glycomsp_reference    Version v0.01.20260921`), *Entries* (here `2,230   resolved 706 · wurcs_only 1,524 · unresolvable 0      by source api 2,230`), *Updated*, *Hash*, *Last prebuild* (here `never`), *Log* and *Resolver* (`available`). Three buttons sit at the right of the frame: *Select database…*, *Set as default* and *Copy hash*. Below, the frame *Options* holds two unticked checkboxes. The frame *Actions* holds the buttons *Prebuild / Update from source…*, *Rebuild log*, *Check again* (greyed out, its label shows a bracketed internal code), *Refresh* and *Cancel* (greyed out), then an empty log area and a status line reading *ready*. The attribution line and a *Close* button end the window.

**The *Reference database* frame.** The frame always shows the active database. Opening the window only reads it.

| Row | What it shows |
|---|---|
| *Path* | The file in use. |
| *Name* | The database name and its version. A red line appears under the rows if the database was written with a different schema than this GlycoMSP expects. |
| *Entries* | The number of entries, then the count per status (`resolved`, `wurcs_only`, `unresolvable`), then the count per data source after *by source* (for example `api`). |
| *Updated* | The time the database was last updated (UTC). |
| *Hash* | The first 12 characters of the file's SHA-256 hash, followed by an ellipsis. |
| *Last prebuild* | The source file name, the time (UTC) and the outcome of the most recent prebuild, or `never`. |
| *Log* | The log path, its number of lines and whether a `.bak` copy exists (`.bak present` or `no .bak`). It reads `not found` if there is no log. |
| *Resolver* | `available`, or `available; prebuild/rebuild disabled:` followed by the reason, or `unavailable:` followed by the error. |

If the database file does not exist, the rows *Name* to *Last prebuild* read *not found — run Prebuild to create it*.

The three buttons at the right of the frame:

- **Select database…** chooses another database file (see *Select a database* below).
- **Set as default** remembers the active database for the next launch (see *Select a database* below).
- **Copy hash** copies the full SHA-256 hash to the clipboard. It is active when a hash is shown.

**The *Options* frame.**

- *Fill GlyToucan ID / WURCS from reference DB on run* is the same setting as the Prepare Dataset checkbox *Fill GlyToucan ID / WURCS from reference DB (see Manage References)*. Ticking one ticks the other. It is off by default.
- *Remember these settings* writes the fill option and the active database to `settings.json` now and again each time either changes, and when you close the window. Unticking it leaves `settings.json` as it is. GlycoMSP still reads the file at the next launch.

**The *Actions* frame.** The buttons are *Prebuild / Update from source…*, *Rebuild log*, *Check again*, *Refresh* and *Cancel*. The label of *Check again* also shows an internal reference code in brackets. The button is greyed out (see *Not available yet*). *Refresh* reads the database and the log again and updates every row. Below the buttons, the log area shows the live output of a prebuild or a log rebuild. The status line under it reports what happened. It reads *ready* when the window opens. A coloured status line marks success, a warning or an error. Below the frame, the attribution line reads *GlyTouCan / GlyCosmos composition data, CC BY 4.0 (https://glytoucan.org, https://glycosmos.org)*.

While a prebuild or a log rebuild runs, the buttons *Select database…*, *Set as default*, *Prebuild / Update from source…*, *Rebuild log* and *Close* are greyed out. Only *Cancel* stays active. If you try to close the window by its title bar, the status line reads *a process is running — cancel or wait first*. Quitting GlycoMSP stops a running process.

### Select a database

Press **Select database…** and choose a `.db` file. GlycoMSP opens it read-only and checks that it is a GlycoMSP reference database. If the check fails, a dialog titled *Not a reference database* shows the path and the error. The status line reads *{file name} refused — previous database kept* and your previous selection stays active. If the check passes, the file becomes the active database and the status line reads *active database: {name} {version} ({path})*.

A selection applies to the current session. To use the database again at the next launch, press **Set as default** and tick **Remember these settings**. *Set as default* alone saves the path to `settings.json`. The next launch applies it only while *Remember these settings* is ticked. The status line says so after you press the button. *Set as default* is active only while the active database opens correctly.

A selected database that is not the default one keeps its own transaction log beside it, named `<database file name without .db>.jsonl`.

The `references/` folder of the GlycoMSP package is for shipped databases. It describes one showcase snapshot, built from a permethylated N-glycan CGA table. Select it with *Select database…* like any other database file. A shipped snapshot is a database that was in use at the time of publication. Later design or version changes may alter its structure and behaviour. There is **no compatibility promise** for shipped snapshots. A prebuild adds to the active database. Select your own database first if you want a shipped snapshot to stay unchanged.

### Build or update the database (prebuild)

Press **Prebuild / Update from source…** to add compositions to the active database. The dialog *Prebuild / Update from source* opens.

> **[SCREENSHOT: `GlycoMSP_screenshots/Prebuild_from_source.png` — *Prebuild / Update from source* dialog with a source file chosen.]**
> In this capture *Batch size (max 30)* was changed to 5. The default is 30. The two lines under the fields name the database and the log of the run.

| Field | What to enter |
|---|---|
| *Source file (CGA TSV / trainable CSV / MAS sheet / in-silico library)* | The file that lists your compositions. Use **Browse…** to choose it. The file must exist. |
| *Column (optional, auto-detect)* | The column that holds the compositions. Leave it empty to let GlycoMSP find it. |
| *Sheet (optional, Excel)* | The sheet of an Excel source. Empty means the sheet `MSlist` if the workbook has one, otherwise the first sheet. |
| *Diff only — no network (--no-network)* | Ticked: report what would be fetched and make no request. Unticked by default. |
| *Batch size (max 30)* | Compositions per request. A whole number from 1 to 30. Default 30. |
| *Pause between batches (s)* | Seconds between two requests. 0 or more. Default 10. |

Under the fields, the dialog shows the database and the log that the run will use. The buttons are **Start** and **Close**.

The source can be a CGA result TSV (pseudolabel TSV), a trainable CSV, a MAS workbook or an in-silico library. GlycoMSP looks for the columns `composition`, `Structure`, `pseudo compositions` and `Composition` in that order. For an in-silico library with residue count columns, it builds the compositions from the counts. The file browser offers `.tsv`, `.csv` and `.xlsx` files, and *All files* for anything else.

Press **Start**. The dialog closes and the log area shows the command that runs, then the output of the prebuild as it appears. The status line reads *running prebuild (pid {number}) — do not start a run until it finishes*. The prebuild runs the command-line tool `msp_CLI_glytoucan_prebuild.py` in a separate process, so the window stays responsive. §3.9 describes the same tool on the command line.

The prebuild compares the source with the database and fetches only the compositions that are missing. It sends at most 30 compositions per request and waits 10 s between requests. It stores the results batch by batch. If the network fails, it stops sending requests and keeps everything it already stored. Run it again later to resume from what is missing. A second run over an unchanged source makes zero requests.

The prebuild needs network access and the `certifi` Python package, which is listed in `src/requirements.txt`. On macOS, Python builds without a system certificate store fail with `CERTIFICATE_VERIFY_FAILED` until `certifi` is installed. Certificate verification is never switched off.

When the process ends, the log area shows `[exit {code}]`, the rows of the window are refreshed, and the status line reports the result:

| Status line | Meaning |
|---|---|
| *done: status={status} committed={n} api_calls={n}* | The prebuild finished. `committed` is the number of entries stored in this run. |
| *network fail-soft — rerun to resume (n_unsent={n})* | The network failed. Everything fetched so far is stored. Run the prebuild again. |
| *storage problem: {status}; see log pane* | The database or the log could not be written. Read the log area. |
| *cancelled — batches already committed are kept; rerun to resume* | You pressed **Cancel**. |

A diff-only run ends with `status=disabled` when compositions are missing and with `status=nothing_to_do` when none are. It still records itself under *Last prebuild*.

**Cancel** stops the process. Batches already committed stay in the database, and a rerun resumes from them.

**Rebuild log** regenerates the transaction log from the database. It asks you to confirm in a dialog titled *Rebuild log*. GlycoMSP then copies the current log to a `.bak` file next to it, replacing an older `.bak`, and replaces the live log as the last step. The database is opened read-only. Rebuild log never creates a database and never changes one. It is active only while the active database opens correctly. The status line then reads *rebuilt {n} lines → {log path} (previous kept as .bak)*. Use it after a prebuild reports that the log could not be written.

### Entry statuses

Every database entry has one of three statuses.

- `resolved`: the entry has a WURCS string and a GlyTouCan accession.
- `wurcs_only`: WURCS is available but the composition is not registered in GlyTouCan. This is seen for some KDN-, HexA- and sulfate-containing compositions. The accession stays empty until the composition is registered in GlyTouCan.
- `unresolvable`: the converter returned no WURCS string.

Values you entered yourself are never overwritten by the database.

### What the fill does at run time

With *Fill GlyToucan ID / WURCS from reference DB* ticked, GlycoMSP fills the cells `GlyToucan ID` and `WURCS` from the active database. The fill never uses the network. It applies in four places:

- the CGA annotation tables (the pseudolabel TSV),
- the MAS merge,
- the CGA → Trainable conversion,
- the prediction output, which gets two additional columns, `Predicted_GlyToucan_ID` and `Predicted_WURCS`. The column `Predicted_Label` is unchanged.

The fill follows these rules:

- It fills **empty cells only**. Your own values always win. In the MAS merge, values already in your sheet are kept.
- A row that already carries your own accession keeps its WURCS cell untouched.
- A row with your own WURCS keeps its accession cell untouched.

The database therefore never pairs a generic value with a structure-specific one. The composition column is never changed.

Each fill writes one line to the main-window log:

```
[glytoucan] fill ({column}): status={status} compositions={n} filled={n} missing={n} user_accession={n} db={version} entries={n} updated={time}
```

Here `compositions` is the number of distinct compositions in the table. `filled` is the number of them that the database has a WURCS string for. `missing` is the number it does not have. `user_accession` is the number of rows that keep an accession you entered. The `status` is `filled` when no composition is missing, `partial` when some are missing, `none` when the table has no composition, and `db_missing` when the database cannot be read. The prediction fill writes its own line starting with `[ML][Predict] GlyToucan/WURCS fill:`.

### Inspect what was filled

Right-click a sample in the Prepare Dataset tree and choose **Reference database info…**. The popup is titled *Reference database info — {sample name}*. For the CGA results and for the MAS results of that sample, it lists:

- the status of the fill,
- the database name, version and path,
- the number of entries in the database and its update time,
- the number of compositions, the number filled and the number not in the database,
- the number of user accessions kept and the number of user WURCS strings kept (marked *local annotation*),
- the time of the last fill.

If no fill is recorded for the sample, the popup reads *No reference fill recorded for this sample.* The popup reads the record of the fill only. It does not open the database.

> **[SCREENSHOT: `GlycoMSP_screenshots/Reference_database_info_popup.png` — *Reference database info* popup for a sample with a recorded fill.]**
> The popup is titled with the sample name. Its text starts with *MAS: reference fill filled*, then lists *database:*, *path:*, *entries:*, *compositions:* with *filled:* and *not in database:*, *user accessions kept:* with *user WURCS kept (local annotation):*, and *last resolved:*. An *OK* button closes it.

> **[SCREENSHOT: `GlycoMSP_screenshots/Sample_context_menu_refinfo.jpg` — right-click menu of a sample in the Prepare Dataset tree.]**
> The menu lists *Move to...*, *Link and Validate Sample*, *Reference database info…*, the disabled *Resolve missing references* and the *Remove* entry of the sample.

When a fill happened, the CGA or MAS method JSON also carries a block `artifacts.glytoucan_resolution`. §5.2 describes the method JSON.

### Annotation without registration

GlycoMSP does not ask you to register structures on GlyTouCan. A row with a blank accession and your own WURCS is a complete local annotation, not a gap.

### Not available yet

Three controls are visible but greyed out. Each label also shows an internal reference code in brackets.

- *Check again* in the *Actions* frame of the Manage References window, for unresolved rows.
- *Resolve missing references* on the right-click menu of a sample.
- *Update GlyToucan ID / WURCS in this file…* on the right-click menu of a CGA result TSV, a CGA-based trainable CSV or an Excel file.

A "check again" action for `wurcs_only` entries is planned. Merging a shipped snapshot into your own database is planned.

### Attribution

The composition data in the reference database come from GlyTouCan and GlyCosmos. They are licensed under CC BY 4.0 (https://glytoucan.org, https://glycosmos.org). The Manage References window shows this line at its foot. The `references/` folder carries it in `ATTRIBUTION.md`.

## 3.8 Native (underivatized) glycans and ion-trap data

Use this section when your glycans are native or native-reduced instead of permethylated. Such samples are often measured in negative mode and often on an ion-trap instrument. The two often come together but they are separate topics. Native masses decide which precursor masses the in-silico library holds. Ion-trap data decide how you convert the files and which tolerances you set. Everything here is additive: with every new option at its default, results equal those of v1.10, apart from the corrections listed in the changelog for GlycoMSP 1.13.

**Scope: the CGA route only.** Native masses apply to the CGA route only (added in GlycoMSP 1.13). The MAS route still uses permethylated residue masses for its theoretical masses, whatever *Derivatization Type* says, so a native project annotated through MAS gets wrong theoretical masses (known issue 22). Use the CGA route for native glycans.

**Choose Native or Native(Reduced).** You set the derivatization in two places.

- In the Metadata Editor (window title *Enter Experiment Metadata*, opened by **Convert Raw to CSV**, §2.2), the *Derivatization Type* list offers `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `Others`. The value is stored in the sample's metadata JSON.
- In CGA Setup, the *Metadata (can override here)* frame shows *Glycan Type*, *Mass Analyzer charge mode* and *Derivatization Type*, filled from the metadata. Its *Derivatization Type* list offers `PerMe`, `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `None`. The CGA Setup screenshot in §2.4 shows this list open. *Reset from file* restores the stored values.

**Generate In-Silico CSV** reads the values shown in CGA Setup. Only `Native` and `Native(Reduced)` give native masses. Any other value, including `None`, gives permethylated masses. `Native` is the free reducing end and `Native(Reduced)` is the reduced end (alditol). With either value the library holds underivatized precursor masses, always as the singly charged ion: [M+H]+ in positive mode and [M−H]− in negative mode (added in GlycoMSP 1.13). Before 1.13 the masses for native libraries were zero or silently permethylated. A library keeps the masses it was generated with, so generate a new one when you change the derivatization or the ionization mode (§2.4 (a)).

**Negative mode.** Set *Mass Analyzer charge mode* to `-` in the Metadata Editor or in CGA Setup. The library then holds the [M−H]− value of each composition. The readers take the polarity of each scan from the scan itself, not from this field. The measured precursor mass that they write to the scan CSV is the singly charged equivalent, so in negative mode it is the [M−H]− equivalent. Before GlycoMSP 1.13, negative-mode precursor masses used the positive-mode formula.

**Sulfate, phosphate and HexA in negative mode.** The *In-Silico Generation Flags* panel of CGA Setup holds these options. *SO3* and *PO3H* are checkboxes, with the ranges *SO3_range:* and *PO3H_range:*, in the frame *substituent(neg)*. HexA is the checkbox *allowHexA* in the frame *glyco-feature restraints* and the range *HexA_range:* in the frame *composition restraints*. For native glycans each group is a flat addition to the mass, and the ion stays singly charged.

The monosaccharide, reducing-end, sulfate, phosphate and HexA masses are documented in the module page `docs/modules/compnewv4.md`. To check the mass of a single composition, use `msp_CLI_precursor_mass.py` (§3.9). It prints the value exactly as the CGA library computes it.

**Convert ion-trap files.** Many ion-trap files do not report a precursor charge for MS2 spectra. GlycoMSP keeps these spectra (added in GlycoMSP 1.13). A spectrum without a charge is exported with its peaks and its precursor m/z, with charge `0` and mass `0`. The status file (`<name>.scan_status.tsv`) marks it as unknown charge (`charge_status` is `unknown`). While such spectra remain, the conversion ends with the outcome *complete with warnings*. The log panel shows the number of MS2 spectra without a charge. RAW files from ion-trap instruments are read from the vendor mass list. Versions before 1.13 exported empty peak lists for them. RAW conversion needs Windows (§2.2). mzML conversion works anywhere.

Before each GUI conversion, the *Conversion options* dialog opens. It has two options and both are off every time it opens.

The screenshot in §2.2 shows the dialog.

- **Isotope comb** (*Calculate the charge from the MS1 isotope pattern (charge 1–3)*). GlycoMSP reads the spacing of the precursor's isotope peaks in its parent MS1 scan and estimates the charge from it. The estimate is used only where the file gives no charge.
- **Dataset charge** (*Dataset charge for the rest (blank = off)*). You enter one whole number other than 0, for example `1`. GlycoMSP gives it to every MS2 spectrum that still has no charge after the comb. The polarity comes from each scan, so enter `1`, not `-1`.

A charge that the instrument reports is never overwritten. The order is: reported charge, then isotope comb, then dataset charge. Your choices are recorded in the project metadata (the keys `Conversion charge comb` and `Conversion dataset charge`). The same two options exist in `msp_extract_cli.py` (§3.9).

The isotope comb has limits.

- It tests charges 1 to 3 only and never calls z ≥ 4.
- It works on profile MS1 scans from non-FT acquisitions only. On mzML it needs an MS1 that was converted without peak picking.
- A mixture of charge states is reported as the stronger one without a flag. Of eight hand-checked z = 3 calls, three were mixtures.
- When the best charge does not beat the runner-up by the required margin, no charge is called. The status file records this as `below_margin`.
- It returns a charge and never a monoisotopic m/z. When the file has no monoisotopic m/z, the precursor mass of such a spectrum is calculated from the isolation m/z.

The comb is a feasibility-checked heuristic, not a validated charge assignment. The scoring steps downstream check its calls (known issue 20).

**How many spectra stay without a charge.** In the ion-trap native files that were checked, 20–29 % of the exported MS2 rows had an unknown charge (known issue 17). With the isotope comb on, 24–27 % of the MS2 scans in the files tried still had no charge called (`below_margin`, known issue 20). These are figures from a few files, so treat them as an order of magnitude. The counts for your own files are in the diagnostics report (`<name>.reader_diag.txt`).

Unknown-charge rows stay in the scan CSV with mass `0`. A mass of 0 never matches a library mass, so CGA analysis leaves such a scan unmatched. Downstream steps do not exclude these rows. A model trained with the mass feature reads mass 0 as the Non-glycan signature (known issue 17). The loss is bounded. It is not the reason why ion-trap data give few annotations with the default tolerances. The tolerances are the reason.

**Set tolerances for ion-trap data (added in GlycoMSP 1.14).** When the scan CSV has no monoisotopic precursor m/z, its precursor mass is calculated from the isolation m/z. On the ion-trap native files that were checked, the measured precursor errors were about −20 to −90 ppm at charge 1 and +150 to +250 ppm at charge 2. Fragment errors were 0.06–0.14 Da. With the defaults (*Precursor match* 20 ppm and *Score A fragment match* 10 ppm), CGA analysis then annotates almost nothing. In one check on such files, the 20 ppm window kept about 5 % of the real precursor matches, and Score A at 10 ppm removed nearly all of those (known issue 27). The error changes sign with charge, so a wider window is needed, not an offset.

The steps below are a worked example, not a recommendation for every dataset. In that check, about 300 ppm for the precursor and 0.3 Da for fragments gave tens of annotations per file against a low decoy rate (known issue 27). Your data can need other values. The decoy check in step 2 tells you.

1. In CGA Setup (§2.4 (d)), set *Precursor match* to `300` with the unit *ppm*, and *Score A fragment match* to `0.3` with the unit *Da*. Then click **Start CGA analysis**. The log panel marks both values *NON-DEFAULT* and the method JSON records them (§5.2).
2. Run the decoy check (§2.4 (d)) at these tolerances and at the defaults.
3. In CGA → Trainable (§2.6), widen the *precursor gate* and *max |ppm_error|* to at least 300 ppm. Choose an *ion list ppm* that covers your fragment errors, within the limit described under *The feature window* below.
4. On the Predict tab (§2.8), load the model. If the CGA → Trainable fields differed from their defaults, the *Model Loaded* dialog lists the values that were prefilled. Check *Apply precursor gate* and *Tolerance (ppm):*. The gate moves a prediction to the review queue when its precursor mass deviates from the predicted composition by more than this value, so keep it at least as wide as the precursor gate of step 3. When you click **Create unlabeled dataset of certain experiment**, the *PPM Tolerance* dialog asks for the feature window. Enter the *ion list ppm* of step 3.
5. Compare Score A only between runs with the same fragment window and the same ion list (§4.1, *How hits are counted*).

**Ion lists for native negative-mode N-glycans (added in GlycoMSP 1.14).** Two ion lists are in the `templates` folder.

| File | Ions | Role |
|---|---|---|
| `GlycoMSP_ionlist_native_neg_14.xlsx` | 14 | primary set |
| `GlycoMSP_ionlist_native_neg_candidate_97.xlsx` | 97 | extended set |

Use one of them with *Native(Reduced)* in negative mode. Attach it with **Add Ion List** in CGA Setup as the Score A ion list. Use the same file as the feature ion list in the field *2) Ion sheet (CSV/XLSX) or Manual Excel:* of CGA → Trainable. At prediction, use the same ion list or a superset of it (§3.6).

With a wide fragment window (0.3 Da or more) and a dense ion list, Score A can count one peak for two neighbouring ions (§4.1). With the 14-ion list no scan is affected. With the 97-ion list about one scan in eight is. The validation scripts report these counts (§8.5).

**The feature window.** The *ion list ppm* feature window in CGA → Trainable and in prediction is ppm only. The *PPM Tolerance* dialog of the prediction step accepts at most 500 ppm, so keep the CGA → Trainable value at 500 ppm or less. At 500 ppm the window is 0.13 Da at m/z 260 and 0.25 Da at m/z 500. A fragment error of 0.14 Da is inside it only above about m/z 280. Low-m/z ion-trap fragments therefore get a narrower window than a Da setting in CGA Setup gives them. A ppm-or-Da selector for feature extraction is planned. §4.10 and §7.2 treat this limit further.

**Always run the decoy check when you widen a window.** A wider window finds more real matches and also more chance matches. Make copies of the in-silico library and of the ion list with every mass increased by 3.37 Da, run CGA analysis against them with the same tolerances, and compare the counts with the real run (§2.4 (d)). Every decoy match is a chance match. The script `tools/validation/cga_tolerance_decoy_count.py` does the counting (§8.5). On a dense library with accurate-mass data, a wide window can give more decoy matches than real ones. In one check on a permethylated Orbitrap-class run at 300 ppm, the decoy count exceeded the real count. The choice of tolerance is therefore made per dataset.

**Other things you will meet.**

- **Size of the CGA result TSV.** The pseudolabel TSV carries every converted scan, so its size follows the raw file. A run of 100,000 MS2 scans gives about 1.7 GB. A matched-rows-only output is planned.
- **Profile MS2 in an mzML.** MS2 spectra stored as profile are written point by point, and the diagnostics warn with "profile data". This happens when MSConvert runs `zeroSamples` before `peakPicking vendor`. MSConvert then skips peak picking. Convert again with `peakPicking vendor msLevel=2-` as the first filter (known issue 19).
- **Score A comparisons.** Compare Score A only between runs at the same fragment window and with the same ion list (§4.1).

## 3.9 Command-line tools

GlycoMSP is a graphical application, and every step of the workflow is available in it. Three command-line tools in `src/` do work that is also useful without the GUI. Run them from the repository root with the Python environment you use for the GUI. Each prints its options with `--help`.

### `msp_extract_cli.py`: convert and inspect

This tool converts a RAW or mzML file to the scan CSVs without the GUI (`convert`) and probes the Thermo backend for selected scans (`inspect`). It calls the same readers as the GUI and offers the same two charge options (§3.8). General-purpose scripting support beyond these two commands is still in progress.

**convert.** RAW files need Windows and Thermo MSFileReader. mzML files convert on any system.

```bash
python src/msp_extract_cli.py convert --raw {input.raw} --out {output_folder} --alias {name} --charge-comb
```

| Option | Meaning |
|---|---|
| `--raw FILE` or `--mzml FILE` | The input. Give exactly one of the two. |
| `--out DIR` | Output folder (required). It is created if missing. |
| `--alias NAME` | Name for every output file instead of the input file stem. With an alias, the input name appears in no output. It has 1–64 characters (letters, digits, `_`, `-` and `.`) and starts with a letter or digit. |
| `--charge-comb` | The isotope comb (§3.8). Off by default. |
| `--dataset-charge N` | One charge for every MS2 spectrum that still has none. A whole number other than 0. Off by default. |

Outputs, in the `--out` folder, where `<name>` is the alias or the input file stem: `ms2_<name>.csv`, `ms3_<name>.csv`, `<name>.convert.json`, `<name>.scan_status.tsv`, `<name>.reader_diag.json` and `<name>.reader_diag.txt`. The two CSVs have the same 9-column MS2 and 7-column MS3 format as in the GUI. The terminal shows the outcome, the charge sources, any warnings and the names of the diagnostics files.

Exit codes of `convert` mirror the four outcomes of a conversion. `<name>.reader_diag.txt` states the same meanings.

| Exit code | Outcome | Meaning |
|---|---|---|
| 0 | complete | Converted, nothing to check. |
| 3 | complete with warnings | Converted. MS2 charges are unknown, no usable precursor mass exists, or the MS2 spectra are profile. |
| 4 | partial | Converted, but some spectra could not be decoded and were not exported. |
| 5 | failed | The input could not be opened, the conversion stopped, no MS2 spectrum with peaks was exported, or the diagnostics files could not be written. A CSV written before the failure is kept. |
| 2 | (usage error) | Bad arguments, an input file that is not found, or an invalid `--dataset-charge`. Nothing is converted. |

**inspect.** This is a read-only probe of the Thermo backend for selected scans of a RAW file. It needs Windows and Thermo MSFileReader. Use it to see which peak requests return data and what charge information the file holds.

```bash
python src/msp_extract_cli.py inspect {input.raw} --alias {name} --scans {scan_1} {scan_2} --out {output_folder}
```

| Option | Meaning |
|---|---|
| `input` | The Thermo `.raw` file (positional). |
| `--alias NAME` | Name used for every output file and inside the report instead of the input file name (required). |
| `--scans N [N ...]` | Native scan numbers to probe (required). |
| `--out DIR` | Output folder (required). It is created if missing. |
| `--window MZ` | m/z half-width kept around the precursor in the parent MS1 scan. Default 3.0. |
| `--full-points` | Also write every returned point of MS1 scans that are requested directly. |

Outputs, in the `--out` folder: `<alias>.inspect.json`, `<alias>.inspect.txt` and tables of the returned points. Exit codes: 0 when the probe ran. Individual backend requests may still have failed, so read the report. 2 for a usage error, and 5 when the input could not be opened or the backend is not available.

### `msp_CLI_precursor_mass.py`: theoretical precursor mass

This tool prints the theoretical precursor mass of a composition exactly as the CGA library computes it (added in GlycoMSP 1.13). It has no mass table of its own. Every value comes from the same function as the library, so it cannot differ from CGA. It does not start the GUI.

```bash
python src/msp_CLI_precursor_mass.py --composition H5N4F1S2 --deri native --reduced --mode -
```

A composition is written in the short CGA form: `F` (Fuc), `H` (Hex), `N` (HexNAc), `S` (NeuAc), `G` (NeuGc), `KDN` (or `K`), `A` (HexA), `s` (sulfate) and `p` (phosphate), each followed by a count, in any order. `H5N4F1S2` and `F1H5N4S2` are the same composition.

| Option | Meaning |
|---|---|
| `--composition LABEL` | The composition (see above). You can repeat it with `--table`. |
| `--deri perme\|native` | Derivatization. Default `perme`. |
| `--reduced` | Reduced (alditol) reducing end. Default is the free end. |
| `--mode +\|-` | Ionization mode. Default: print both ions. |
| `--hexa N`, `--so3 N`, `--po3h N` | HexA, sulfate and phosphate counts. Same as `A<N>`, `s<N>` and `p<N>` in the label. |
| `--json` | Print the result as JSON. |
| `--table FILE` | Write a table to a `.csv` or `.md` file. With `--composition`, the table lists those compositions. Without it, the table lists the 12 rows of the self-test. |
| `--selftest` | Check the built-in reference rows (see below). |

Outputs: a printed block with the composition, the derivatization, the neutral mass where one is defined, the singly charged ion or ions and the source of the value. With `--json` the same content is printed as JSON. With `--table` the tool writes the file you name.

`--selftest` computes 12 reference rows through the library function. The permethylated rows must equal the tool's stored reference values exactly. The native rows must be within 0.002 Da of their references, which were checked against GlycoMod and GlycoWorkBench. The tool prints one line per row and then a count. A shortened print looks like this:

```text
H5N4       PerMe  reduced=no  [M+H]+  2048.05329  ref 2048.05329  identical
H5N4       Native reduced=no  [M+H]+  1641.59942  ref 1641.59900  within 0.002 Da
...
selftest: 12 of 12 rows ok (compnewv4 1.08, 20261027)
```

Exit codes: 0 for success, 1 when a self-test row does not match, and 2 for a usage error or an input the tool does not support. It refuses permethylated compositions that carry both sulfate and phosphate, and permethylated compositions with more than two sulfates or more than two phosphates. The library has no mass rule for them (known issue 26). The tool prints a message and exits with 2. Native compositions have no such limit. The CGA library itself still computes the permethylated cases as in v1.10, and those values are not reliable.

### `msp_CLI_glytoucan_prebuild.py`: reference database prebuild

This tool builds or updates the local reference database of GlyTouCan accessions and WURCS strings (§3.7). It is the only component of GlycoMSP that contacts the network. The fill option in Prepare Dataset reads the database and never calls the network. The tool reads the compositions in a source file, fetches only those that the database does not hold yet from the GlyCosmos composition API (at most 30 compositions per request, with a 10 s pause between requests) and stores the results batch by batch. A network failure stops further requests and keeps everything already fetched. Run it again to resume. A second run over an unchanged source makes no API calls. The tool needs network access and the `certifi` package. On macOS Python builds without a system certificate store, a run fails with `CERTIFICATE_VERIFY_FAILED` until `certifi` is installed.

```bash
python src/msp_CLI_glytoucan_prebuild.py --source {sample_CGA_result.tsv}
```

| Option | Meaning |
|---|---|
| `--source FILE` | CGA TSV, trainable CSV, MAS workbook or in-silico library CSV. Required, unless `--summary` or `--rebuild-log` is given. |
| `--column NAME` | The composition column. It is detected when omitted. |
| `--sheet NAME` | Excel sheet name. Default is `MSlist` if present, else the first sheet. |
| `--db FILE` | The reference database. Default `glycomsp_reference.db` in the `references` folder of the GlycoMSP data folder. |
| `--log FILE` | The transaction log. Default `glytoucan_cache.jsonl` in the GlycoMSP data folder. |
| `--no-network` | Report what would be fetched. Makes no API call. |
| `--batch-size N` | Compositions per request. Default 30. This is also the maximum. |
| `--pause S` | Seconds between requests. Default 10. |
| `--timeout S` | Timeout of one request in seconds. Default 20. |
| `--summary` | Print the database status and stop. |
| `--rebuild-log` | Regenerate the transaction log from an existing database. The previous log is kept as `.bak`. It never creates or changes a database. |
| `--json` | Print the result as JSON. |

The GlycoMSP data folder is `~/.glycomsp/` on macOS and Linux and `%LOCALAPPDATA%\GlycoMSP\` on Windows. The environment variables `GLYCOMSP_REFERENCE_DB` and `GLYCOMSP_GLYTOUCAN_LOG` override the database and the log location.

Outputs: the database and the log, and lines starting with `[prebuild]` that report the status and the counts of committed, already present, unsent, unpersisted and unlogged compositions, and the number of API calls.

Exit codes: 0 when the run is done, when there is nothing to do, or when the network failed softly (status `offline` or `partial`, check the unsent count and run again). 2 for a storage problem (the database write failed, the transaction log is incomplete, or `--rebuild-log` was refused). A usage error such as a missing `--source` also exits with 2.

The same action is available in the GUI as *Prebuild / Update from source…* in the Manage References window (§3.7).

### Prediction from the command line

`src/msp_CLI_predict_skops.py` exists, but in GlycoMSP 1.14.1 it has known issues 30 to 32. It cannot load a model with the pinned skops version, its default `--drop-cols` removes the precursor-mass feature, and it writes integer class labels. Predict from the GUI (§2.8, and §7.1 for problems that can occur). The GUI is not affected by these issues.

---

# Part 4 — Concepts

This Part explains the design choices behind GlycoMSP. Reading is optional for a first run, but useful for thesis review and for understanding why a particular pipeline behaves the way it does. Formulas and naming follow the manuscript Methods.

## 4.1 Score A — ion-fitting score

Score A (manuscript variable *a*) is the **ion-fitting score**: the proportion of diagnostic ions detected in an MS2 spectrum, evaluated against the user-defined ion list. For each scan that has at least one candidate composition, the application:

1. Takes every reference ion of the ion list. The list does not depend on the candidate, so all candidate compositions of one scan receive the same Score A.
2. For each ion, searches the spectrum's peak list for peaks within the *Score A fragment match* tolerance (default 10 ppm). Since GlycoMSP 1.14 this tolerance is a field in CGA Setup, in ppm or Da. A peak with an intensity of 1 or less is not counted.
3. Computes the score as the number of hits divided by the number of ions in the list, limited to `[0, 1]`.

**How hits are counted.** Score A counts peak–ion matches, not distinct ions. For each ion in the list, every peak inside that ion's window counts as one hit. This has two consequences:

- **Several peaks in one ion's window count several times.** With a ±0.3 Da window, peaks at 299.80, 300.10 and 300.25 all lie in the window of ion 300.00, giving 3 hits for that one ion. `ion hits m/z` then lists 300.00 three times.
- **One peak can count for two ions.** If two ions are closer than twice the window (for ±0.3 Da: closer than 0.6 Da), one peak between them can fall in both windows. Ions 406.20 and 406.60 with one peak at 406.42 give 2 hits.

At the default 10 ppm (±0.003 Da at m/z 300), the first case almost never happens. The second depends on how closely your ion list is spaced. Both become common at wide windows: in the zebrafish brain test at 300 ppm / 0.3 Da, 32 % of matched rows had an ion counted more than once. GlycoMSP 1.14 keeps the v1.10 rule, so earlier results stay reproducible. What this means for you:

- Compare Score A only between runs with the same fragment window and the same ion list.
- Before widening the window, check your ion list for ions closer than twice the window.
- Check a widened window against chance with a decoy count (§2.4 (d)).

The manuscript frames Score A as serving "as a general filter to reject spectra having possible compositions without enough glycan-related fragments." It works at the **fragment-evidence level**: does the spectrum contain the glycan-related fragments of the ion list? It does not consider structural detail beyond fragment presence.

In practice, Score A decides which scans keep a label when Score B is inactive. A scan is kept when its Score A reaches the user-set minimum threshold. The candidates of that scan share one Score A, so Score A does not rank them against each other. CGA analysis writes the candidates of a scan in order of increasing precursor error, and CGA → Trainable takes the first candidate that passes its precursor filters (§2.6). Without Score B and with the default *Top N per scan* of 1, the assigned label is therefore the candidate with the smallest precursor error. The threshold is exposed in the **CGA → Trainable** modal (Quickstart Step 2.6). Score A is a ratio in `[0, 1]` (matched ions ÷ total ions in the list); GlycoMSP applies a default minimum of `0.07`. The manuscript does not prescribe a threshold value beyond this default — the choice is left to the user.

If the per-scan score attach fails, GlycoMSP falls back to a simpler scorer (`_fallback_simple_ion_scoring`) that matches the **same user-provided ion list** against each spectrum and computes the matched-fraction score. (Earlier builds used a small set of hardcoded anchor masses here; that was removed during the v14 scoring-logic fix, so the fallback now always scores against your configured ion list.) A properly configured ion list matching your derivatization chemistry is therefore required — there is no built-in default ion set.

## 4.2 Score B — motif-aware re-ranking

Score B (manuscript variable *b*) operates one level above Score A, on **structural motif evidence**. Where Score A asks *"do the expected fragments appear?"*, Score B asks *"do the appearing fragments imply the structural motifs that this candidate would predict?"*. The manuscript defines:

`b = max(0, min(1, support − penalty_comp − penalty_unexpected))`

where:

- **support** quantifies motif-level evidence assembled from user-defined glycotopic rules. Each rule contributes positively when its diagnostic fragments appear (subject to gate conditions if the rule is gated); contributions sum across all rules in the workbook.
- **penalty_comp** subtracts evidence that contradicts the candidate's composition — for example, a sialic-acid-specific fragment in a candidate that has no sialic acid.
- **penalty_unexpected** subtracts evidence that no expected motif accounts for — peaks that look diagnostic but match no rule applicable to this candidate.
- The outer `max(0, min(1, …))` clamps the result to the unit interval, so all candidates remain on a comparable `[0, 1]` scale regardless of motif-rule density.

The selection rule used downstream in CGA → Trainable is what the manuscript names *"Score B preferred, Score A fallback"*: the per-scan top-ranked Score B candidate is selected when it has a non-empty assigned composition; otherwise the workflow reverts to Score A's top-ranked candidate. In practice this means Score B refines decisions where motif evidence discriminates, and stays out of the way where it doesn't.

Score B is **inactive by default**. Activation requires a Score B Excel workbook authored by someone with glycotopic domain expertise; see Section 3.3 for attaching the workbook and the manuscript's Supplementary Data 2 for authoring instructions. The manuscript explicitly recommends that users without established glycotopic knowledge proceed with Score A alone (CGA-ML) and validate predictions manually rather than attempt Score B configuration without the prerequisite background.

## 4.3 Feature extraction

Both pipelines (CGA and MAS) produce feature matrices using the same per-cell formula. For each spectrum and each reference ion *n*:

`i_feature,n = log10(I_n) + 1`

where *I_n* is the maximum peak intensity within the configured ion-matching tolerance of fragment *n*'s reference mass. When no peak matches, the cell takes the default value `1.0`. In the MAS merge, a window that holds two or more peaks takes the last peak (the one with the highest m/z) instead of the most intense one (known issue 15). A matched peak with an intensity of 0 or below stops the feature build in CGA → Trainable with a math-domain error (known issue 33).

The window for this match is the **feature window**. For training sets it is **10 ppm** by default: fixed in the MAS merge, and the *ion list ppm* field in CGA → Trainable. For unlabeled **prediction** data it is the value you enter in the unlabeled-dataset dialog (§2.8). That dialog has no fixed value, and its prompt suggests 20 ppm, which is more permissive than the training default. At prediction a stricter precursor mass gate is also applied as a quality check, with a default of **10 ppm**. §4.10 lists every tolerance, its unit, and the mass it is measured against.

Two design choices worth understanding. First, **the log transform is applied before the +1, not after**: the formula is `log10(I) + 1`, not `log10(I + 1)`. The two are equivalent only at high intensity; at low intensity they differ. The application's choice means that a peak with raw intensity 1.0 produces a feature value of 1.0 (clean coincidence with the no-hit default), peaks below 1.0 produce values less than 1.0 (so weaker hits are distinguishable from no-hits but compressed), and peaks above 1.0 produce values greater than 1.0 (the dominant signal range). Second, **no-hit cells default to 1.0 rather than 0.0**: this lets the Random Forest treat "no peak in window" as a distinct, low-but-nonzero baseline rather than collapsing it into the same bin as a tiny hit, which would distort feature-importance estimates.

Feature extraction is identical for CGA-derived and MAS-derived rows. The model trained on a mix sees one feature space and treats all rows uniformly.

## 4.4 The ML model — Random Forest

GlycoMSP ships a single, well-tuned Random Forest classifier. The choice was deliberate: Random Forests handle mixed feature scales without normalisation (intensities span several orders of magnitude even after log-transform), are robust to small datasets and class imbalance, expose interpretable feature importances (useful when explaining model behavior to chemists), and are fast enough to train and predict on commodity laptops. The frozen v1.10 configuration is 400 decision trees, `max_depth=None` (unrestricted depth), `min_samples_split=2`, `min_samples_leaf=1`, `class_weight="balanced"`, and `random_state=42`; scikit-learn defaults apply for the rest (e.g. `max_features="sqrt"` and bootstrap sampling). These are the values reported in the manuscript.

**Threshold τ + margin.** The application supports a confidence step at prediction time, **disabled by default** (the *Enable confidence threshold → fallback to Majority label* checkbox in the Train/Test Parameters dialog). When you switch it on, it uses a threshold τ (*Threshold τ (0.00–0.99):*, default 0.65) and a margin (*Majority support margin δ (0.00–0.20):*, default 0.05). It acts in two places. In the predicted labels, a glycan prediction becomes the majority label (Non-glycan) when its top probability is below τ and the Non-glycan probability lies within the margin of the top probability. In the prediction report, rows with a top probability below τ are written to `to_review_low_conf.csv`, and rows whose margin (top minus second probability) is below the margin value are written to `to_review_low_margin.csv`. The two values work together. τ selects the predictions the model is unsure of. The margin limits the relabelling to those for which Non-glycan is nearly as probable as the top class.

The threshold step is disabled by default, which is the recommended setting for small datasets. With the checkbox off, no label is changed and the application passes τ = 0.0 and margin = 1.0 to the report, so the run summary shows τ 0. Whatever the threshold setting, the run summary shows *Accepted 0* and probability statistics of 0.000, which is a known limitation (known issue 29). The predicted labels and per-class counts are correct. Threshold demotion mostly helps once datasets exceed a few hundred spectra per class, where probability calibration becomes meaningful. The code holds four different τ defaults in total. The values above are the ones the GUI applies (§7.2).

**Class balancing.** Both Score-A-only and Score-B-active runs can produce class-imbalanced datasets (more candidates for common compositions, fewer for rare ones; lots of Non-glycan negatives sampled at default ratios). The application includes two mitigations: per-class minimum size (drop classes with fewer than 5 spectra; manuscript baseline) and majority-class capping (cap the dominant class — typically Non-glycan — at no more than 3× the largest glycan class; manuscript baseline). By default the cap applies to the training split only, so the held-out test split is not capped. Both are configurable in the Train/Test Parameters dialog. The defaults reflect what worked across the development datasets and are reasonable starting points.

## 4.5 Reproducibility — `training_run.json`

Each successful training run writes a JSON manifest next to the model file. The manifest captures everything needed to reproduce the run given the same input data. All fields sit under one top-level key, `training_run`:

- **`effective_params`** — the fully merged ML parameters that were actually used (UI vars + JSON editor + presets layered together).
- **`column_order`** — the exact ordered list of feature columns. Critical: at predict time the application reorders unlabeled-input columns to match this list; without it, Random Forest predictions become positional gibberish.
- **`classes`** — the ordered list of class labels seen during training.
- **`inputs`** — the path of the trainable CSV (`train_csv`) and of the linked experiment JSON (`linked_exp`).
- **`hashes`** — the SHA-256 hash of the trainable CSV (`train_csv`).
- **`versions`** — environment versions: Python, platform, sklearn, joblib, numpy, pandas, and the GlycoMSP build version (key `gms_gui`).
- **`seeds`** — the random seed of the train/test split.
- **`created_at`** — local-time ISO timestamp of when training completed.
- **`cga_tolerances`** (added in GlycoMSP 1.14, always written): the tolerances behind the trainable CSV. If GlycoMSP finds the CSV's CGA → Trainable record, this is `{"status": "found", "source": …, "trainable": {…}}`, with the fields listed in §5.2. Otherwise it is only a status: `not found` (MAS, pooled, moved or older-build trainables), `stale` (the CSV has changed since it was built), `ambiguous` (two different records for the same file), `trainable CSV changed during read`, or `lookup failed` followed by the error type in parentheses. GlycoMSP finds the record when the CSV was built in the same session, or when its CGA method JSON has been loaded or sits in the CSV's folder or the linked experiment JSON's folder, and only while the CSV still has the SHA-256 recorded when it was built. Prediction prefills from this block (§2.8).

This manifest accompanies the model file (`*.skops` plus `*_labelencoder.joblib` plus `*_features.json`). When you publish a model or share it with a collaborator, include all four files. The manifest is the document that lets reviewers verify a run was performed under the parameters reported in a paper.

## 4.6 Why two parallel labeling pipelines?

The architecture sits on a deliberate observation: the bottleneck for ML in glycomics is not the classifier — it's the labels. Most public glycomics datasets are not annotated in a way that supports cross-study ML, and even within a single lab the cost of expert annotation makes it impractical to manually label every spectrum.

GlycoMSP addresses this with two complementary entry points. **CGA** generates labels computationally from in-silico libraries and constraint scoring. The labels are imperfect by construction — Score A is a heuristic, Score B helps but requires expert workbook authoring — but they are *available*, in volume, for any dataset where you can specify the chemistry. **MAS** ingests expert-validated labels when you have them. The labels carry information that no rule-based system can match (subtle structural assignments, ambiguity resolution from human inspection), but they require expert time and tend to exist only for studies that have been carefully reviewed.

The design choice to share **feature extraction and ML stages** across both pipelines is what makes mixing possible. A researcher starting a new line of work might collect a small expert-annotated set (MAS rows) and combine it with a larger CGA-derived set on the same instrument, training a single classifier on the union. The expert labels anchor the model on confident assignments; the constraint labels expand its exposure. This is the *CGA-ML / MAS-ML* duality the manuscript describes — two ways into the same trained-model surface.

## 4.7 Ion mining (terminology)

The Prepare Dataset window has a button labelled **Ion Mining** (and a sister button **Browse Ion Mining result**). This manual uses "ion mining" throughout, matching the GUI label. The procedure is a frequency-and-lift analysis that ranks m/z peaks by how preferentially they appear in glycan vs Non-glycan spectra, with Laplace add-one smoothing to handle sparse counts.

Conceptually the same procedure could be called *feature mining*, since it surfaces candidate features for the ML model — but that name is not used in the manuscript and is not adopted as a formal term here. The GUI term **"ion mining"** is authoritative throughout this manual.

## 4.8 Evaluating MAS-ML: within-tissue, leave-one-tissue-out, and pooled CV

Because GlycoMSP trains *local* models, "performance" depends on which spectra the model has seen. The manuscript evaluates MAS-ML three complementary ways, and understanding the distinction matters when interpreting any GlycoMSP model:

- **Within-tissue held-out test** — split one tissue's MAS data into train/test (80/20, stratified, classes with ≥5 spectra). This asks: *can the model recover held-out scan-level labels from the same sample context?* It is the cleanest in-sample check and the strongest numbers (zebrafish glycan macro-F1 ≈ 0.93–0.95).
- **Leave-one-tissue-out (LOTO)** — train on two tissues, predict the third (held-out) tissue as an independent biological sample. Because tissues do not share identical glycan profiles, held-out compositions are split into **seen** classes (present in the training label space) and **out-of-label-space** classes (absent from training). Metrics are reported on seen classes only; out-of-label-space compositions are reported separately, **not** as classification errors (a supervised model cannot predict a label it never saw). LOTO performance is lower than within-tissue (≈ 0.60–0.83), which is expected and reflects label-space mismatch, not model failure.
- **Pooled cross-validation** — pool all tissues and run stratified 5-fold CV. Pooling raises the number of trainable compositions (≈14–25 per tissue → 43 pooled) and stabilises performance (macro-F1 ≈ 0.90 ± 0.04).

The practical lesson: a local RF model should be used **within its trained label space**, and accumulating compatible, manually reviewed spectra (pooling) is how you expand coverage. Primary metrics are always glycan-only — the operational Non-glycan class (an *estimated* negative, not ground truth) is reported separately. Part 8.2 gives the exact commands and expected numbers; the standalone reproduction script and its tech note document the metric definitions in full.

## 4.9 What the readers do and do not do

This section explains what the RAW and mzML converters write into the scan CSV and what they leave out. Converting files is described in §2.2, the command-line converter in §3.9, and the file layout in §5.1. The reader diagnostics and the two charge options were added in GlycoMSP 1.13.

**What a reader exports.** A reader writes one row per MS2 spectrum, with the peak list as the instrument or the converter stored it. MS3 spectra get their own rows, and each carries its own peaks. The 9-column MS2 and 7-column MS3 layouts are unchanged from v1.10. A spectrum with no peaks is not exported. A spectrum whose arrays cannot be decoded is not exported either, and the conversion then ends as *partial*. Neither case stops the conversion. The diagnostics report counts them.

**Vendor-centroided peaks, no peak picking.** GlycoMSP does not pick peaks and does not centroid. It writes the peaks that the source file holds.

- **RAW.** For FT scans, whose Thermo scan filter names FTMS, the reader takes the label data, which is the centroided peak list the instrument stores. Ion-trap scans have no label data. For them the reader takes the mass list from the Thermo library. If the scan is stored as profile, the Thermo library centroids it with its default peak width. If it is stored as centroid, the list is returned as stored. The result matches the vendor peak picking of MSConvert (known issue 16). Before 1.13, RAW conversion read label data only, so ion-trap scans were exported with empty peak lists.
- **mzML.** The peaks are whatever the mzML holds. If the MS2 spectra were stored as profile, every profile point is written as a peak, zeros included, and a row can hold hundreds of points. The diagnostics report then warns that exported MS2 spectra are profile data, and the outcome is *complete with warnings*. The remedy is to convert again with MSConvert and put `peakPicking vendor msLevel=2-` as the first filter. When `zeroSamples` runs before peak picking, MSConvert skips peak picking and prints "peakPicking is not the first filter" (known issue 19). With peak picking first, the mzML route reproduces the rows, peaks and charges of the RAW route.

**Where the precursor charge comes from.** For each MS2 spectrum the reader tries four sources in this order and takes the first that gives a charge.

1. **Reported by the file.** A charge of 0 in the file means "not determined" and counts as missing. A reported charge is never overwritten by any later source. The status file marks it `charge_source reported`.
2. **Isotope comb**, only if you switched it on. In the conversion options dialog this is the checkbox *Calculate the charge from the MS1 isotope pattern (charge 1–3)*. The status file marks the result `charge_source comb`.
3. **Dataset charge**, only if you entered one. In the conversion options dialog this is the field *Dataset charge for the rest (blank = off):*. One charge applies to every MS2 spectrum still without one. Enter 1, not −1, because the polarity comes from each scan. The status file marks it `charge_source dataset_value`.
4. **None.** The charge stays unknown and the status file marks `charge_source none`.

Both options are off by default. The conversion options dialog appears before each GUI conversion and applies to every file of the batch. Your choice is recorded in the project metadata. The reader never substitutes a charge silently, so no spectrum receives a charge of 1, −1 or the instrument method's default unless you asked for it through the dataset charge.

The status file also carries `charge_status` (`known` or `unknown`) and `mass_status` (`eligible` or `unavailable`) for each exported MS2 spectrum, and `comb_status` with the comb's result (for example `called` or `below_margin`, and `not_run` when the comb did not run). The diagnostics report counts the spectra for each charge source.

**What the isotope comb is.** The comb reads the spacing of the precursor's isotope peaks in its parent MS1 scan and estimates the charge from it. Its limits are part of its design (known issue 20):

- It tests z = 1, 2 and 3 only and never calls z of 4 or more.
- It works on profile MS1 scans that were not acquired by FT. In RAW files that means ion-trap profile MS1. In mzML files it means a profile MS1 spectrum, named in the MS2 spectrum, whose filter string does not state FTMS.
- A mixture of charge states is reported as the stronger one, without a flag. In a hand check, 3 of 8 z = 3 calls were mixtures.
- On the ion-trap files tried, 24–27 % of the MS2 scans stayed `below_margin`, so the comb gave them no charge.
- It returns a charge and never a monoisotopic m/z.

The comb is a heuristic that was checked for feasibility, not validated. Its calls are checked downstream by scoring.

**Rows with unknown charge: mass 0.** A spectrum whose charge stays unknown is still exported. The row keeps its peaks and its precursor m/z. It has `chargeState` 0 and `protonatedmass` 0. The status file marks it `charge_status unknown` and `mass_status unavailable`, the GUI log shows the count of MS2 spectra without a charge, and the outcome is *complete with warnings*. Downstream, mass 0 never matches a library mass in CGA, so such a scan stays unmatched and no reason is given. No downstream step excludes these rows. A model trained with the mass feature sees mass 0 as a Non-glycan signature (known issue 17). Known issue 17 records this behaviour as left unchanged on purpose, and measured 20–29 % of exported MS2 rows with unknown charge on the ion-trap native files it tested. The comb and the dataset charge reduce the number of affected rows. The status file lets you count them before you go on.

**Polarity and precursor fields in mzML.** The mzML reader takes the polarity from the polarity term of the spectrum. If the term is missing, it takes the polarity from the scan filter string. If neither states a polarity, the precursor mass uses the positive-mode formula and the status file leaves the polarity blank. The precursor m/z, its intensity and its charge always come from the first selected ion of the first precursor, all from the same ion. A charge is not taken from a later ion or precursor. GlycoMSP 1.13 also corrected the negative-mode precursor mass. Earlier versions used the positive-mode formula. In negative mode `protonatedmass` is now the [M−H]− equivalent, which is the form the in-silico library stores for negative mode.

**Outcome and exit code.** Every conversion ends with an outcome (*complete*, *complete with warnings*, *partial* or *failed*) and a matching exit code (0, 3, 4 and 5). The GUI shows the outcome in its log. A failed conversion is never saved as if it had succeeded, and a partial one is saved with a warning. §2.2 shows where the outcome appears, §3.9 lists the exit codes of the command-line converter, and §5.1 describes the status file and the diagnostics report.

**What the readers do not do.**

- They do not pick peaks or centroid. For profile ion-trap scans in RAW files the Thermo library does the centroiding.
- They do not deisotope. The peak list is written as stored.
- They do not export retention time. The scan CSV has no such column.
- They do not assign a charge beyond the two options above.

## 4.10 Precursor and fragment tolerances

GlycoMSP uses several tolerances. They are set in different windows and they do not follow each other. This section lists them, says what each one is measured against, and explains when to choose ppm and when Da. The steps are in §2.4 (d), §2.6 and §2.8. The two CGA Setup fields and the *precursor gate* field were added in GlycoMSP 1.14. With every field at its default, results equal those of v1.10, and CGA analysis writes the same TSV.

| Setting | Where | Unit | Default | Measured against |
|---|---|---|---|---|
| *Precursor match* | CGA Setup, *Mass tolerances (CGA matching)* | ppm or Da | 20 ppm | The observed precursor mass. A ppm window is the observed mass × ppm ÷ 10⁶ on both sides. A Da window is ± that many Da around it. |
| *Score A fragment match* | CGA Setup, *Mass tolerances (CGA matching)* | ppm or Da | 10 ppm | The ion-list mass. A Da window is symmetric. |
| *max \|ppm_error\|* | CGA → Trainable, box *4) Thresholds / selection* | ppm only | 20 | The theoretical mass. It tests the `ppm_error` column of the pseudolabel TSV. |
| *precursor gate* (checkbox and ppm entry) | CGA → Trainable, box *4) Thresholds / selection* | ppm only | on, 10 | The theoretical mass. It is recomputed from the scan CSV. |
| *ion list ppm* | CGA → Trainable, box *4) Thresholds / selection* | ppm only | 10 | The ion-list mass. This is the feature window. The Non-glycan sampling of this modal uses it too. |
| Feature window of the MAS merge | no field | ppm | 10, fixed | The ion-list mass. |
| Feature window at prediction | *PPM Tolerance* dialog of *Create unlabeled dataset of certain experiment* | ppm only, 1–500, decimals allowed | none, the prompt suggests 20 | The ion-list mass. |
| *Apply precursor gate* and *Tolerance (ppm):* | Predict tab, *Precursor (MS1) Gate for Predictions* | ppm only | on, 10 | The theoretical mass of the predicted composition. |
| Score B ion matching | no field | ppm | 20, fixed | The mass of the workbook ion. |

**ppm versus Da.** A ppm window grows with the mass. A Da window has the same width at every mass. At m/z 300, 10 ppm is ±0.003 Da, and at m/z 1,000 it is ±0.010 Da. A window of 0.3 Da is 1,000 ppm at m/z 300 and 300 ppm at m/z 1,000. Orbitrap-class data with a monoisotopic precursor m/z stay within the default windows, and mass errors as large as those of ion-trap data are not expected for them, so ppm suits them. For ion-trap fragments the measured errors are given in Da, 0.06 to 0.14 Da. In ppm that range is 300 to 700 ppm at m/z 200 and 60 to 140 ppm at m/z 1,000. No single ppm value fits the whole ion list, while one Da value does.

**Isolation m/z versus monoisotopic m/z.** When a RAW file has no monoisotopic precursor m/z, the precursor mass is calculated from the isolation m/z. This is common for ion-trap files. On the native ion-trap files measured, the precursor error was about −20 to −90 ppm at charge 1 and +150 to +250 ppm at charge 2, and the fragment error was 0.06 to 0.14 Da. The default 20 ppm kept about 5 % of the real precursor matches (known issue 27). The error changes sign with the charge, so one offset cannot correct it. A wider window is needed, in ppm or in Da. §3.8 gives a worked example.

**The fragment-match rule in ppm is slightly asymmetric.** For a peak above the ion mass, the error is taken relative to the peak mass instead of the ion mass. This moves the upper edge of the window by less than 0.3 ppm even at 500 ppm (known issue 10). The Da window is symmetric.

**The feature window is ppm only (known limitation).** The prediction dialog stops at 500 ppm, which is ±0.13 Da at m/z 260. A fragment window of 0.3 Da is about 1,150 ppm at m/z 260, so low-m/z ion-trap fragments get a narrower window at feature extraction than CGA gave them at scoring. A ppm-or-Da selector for feature extraction is planned. Until then, keep *ion list ppm* at 500 ppm or less, so that prediction features can be built with the same window (§2.6).

**Wider windows and chance matches.** A wider window finds more matches, and some of them are chance matches. Estimate the chance rate with the decoy check (§2.4 (d)). It repeats the run against a decoy library, which is a copy of the in-silico library and ion list with every mass shifted by 3.37 Da. The scripts that do this are described in §8.5. A wider fragment window can also raise Score A without added evidence, because one ion or one peak can be counted more than once (*How hits are counted*, §4.1).

**Score B enrichment.** Score B enrichment keeps a fixed 20 ppm. It has no field and does not follow *Score A fragment match*.

Non-default tolerances are recorded in the CGA method JSON and in `training_run.json` (§5.2, §4.5), and they prefill the prediction dialogs (§2.8).

---

# Part 5 — Reference

This Part is the reference. It is denser than Parts 1–4 and intended to be searched, not read linearly. Schemas, default values, and output column specifications are listed exhaustively where useful for reproducibility.

With every option added after v1.10 left at its default, the data files and the method JSON files that GlycoMSP writes equal those of v1.10, except where a section below says otherwise.

## 5.1 Input file formats and column requirements

### Scan CSV (output of Step 2.2, input to everything else)

Produced by GlycoMSP's converter (`mspext.convert_raw_to_csv` for RAW, `mspmzmlext.extract_mzML` for mzML). One row per exported MS2 scan. The file is tab-separated although its name ends in `.csv`. The conversion also writes a second file with one row per exported MS3 scan. The 9-column MS2 format and the 7-column MS3 format are unchanged since v1.10.

**MS2 file** (`ms2_{savename}.csv`), columns in file order:

| Column | Type | Notes |
|--------|------|-------|
| `entry_no` | int | Running number of the exported MS2 rows, starting at 1. |
| `MS1scan_no` | int | Parent MS1 scan number. `0` when the file does not give one. |
| `MS1_isolationmass` | float | MS1 isolation (precursor) m/z. `0` when unknown. |
| `MS1_monoisolationmass` | float | Monoisotopic isolation m/z reported by the file. For mzML the isolation m/z is used when the file has no monoisotopic term. `0` when unknown. |
| `chargeState` | int | Precursor charge state. `0` means unknown (see below). |
| `protonatedmass` | float | Precursor mass in monoisotopic protonated form, calculated from the precursor m/z and the charge. For RAW files the monoisotopic m/z is used when the file gives one, otherwise the isolation m/z. Used by the precursor MS1 gate. `0` when the charge is unknown. |
| `MS2scan_no` | int | Scan number from the instrument file. Primary key for downstream joins. |
| `peaklist` | str | Peak m/z values of the spectrum, written as a Python tuple in parentheses, for example `(94.03, 95.39, …)`. |
| `peakintensity` | str | Peak intensities in the same format and of the same length as `peaklist`. |

The MAS merge checks this exact header and these types when it loads a scan CSV: `entry_no`, `MS1scan_no`, `chargeState` and `MS2scan_no` as integers, the mass columns as floats, and `peaklist` and `peakintensity` as text.

**MS3 file** (`ms3_{savename}.csv`), columns in file order:

| Column | Type | Notes |
|--------|------|-------|
| `entry no` | int | Running number of the exported MS3 rows. The header contains a space, not an underscore. |
| `MS2scan_no` | int | Parent MS2 scan number. |
| `MS2_isolationmass` | float | Isolation m/z of the parent MS2 scan. |
| `MS2_monoisomass` | float | Monoisotopic isolation m/z of the parent MS2 scan. |
| `MS3scan_no` | int | Scan number of the MS3 spectrum. |
| `peaklist` | str | The MS3 spectrum's own peaks, same format as above. |
| `peakintensity` | str | Intensities, same format. |

A spectrum without peaks is not exported. It is counted in the status file and the diagnostics report (added in GlycoMSP 1.13).

Calculation of `protonatedmass`: for positive mode with charge 1 the value equals the precursor m/z. For positive mode with charge z of 2 or more it is `m/z · z − (z − 1) · 1.00784`. For negative mode (the scan polarity decides) it is `m/z · z + (z − 1) · 1.007276`, the [M−H]− equivalent (the positive formula was used in negative mode before GlycoMSP 1.13).

**Unknown charge (`chargeState 0` and `protonatedmass 0`, added in GlycoMSP 1.13).** When the instrument file reports no precursor charge, the spectrum is still exported with its peaks and its precursor m/z. Its `chargeState` is `0` and its `protonatedmass` is `0`. You recognise these rows by `chargeState` equal to `0`. The status file marks the same rows with `charge_status` `unknown`, and the diagnostics report counts them as *unknown charge*. A charge supplied by one of the conversion options (isotope comb or dataset charge) replaces the `0`, and the status file records which source gave it. Steps that need a precursor mass cannot use rows with `protonatedmass` equal to `0`. The downstream steps do not yet exclude rows with unknown charge.

Retention time (RT) and other LC/EIC-related parameters are **not** retained — GlycoMSP is built around MS2 scan-level traceability tied to the peak list and intensities, not chromatographic features.

### Status file `<name>.scan_status.tsv` (added in GlycoMSP 1.13)

Every RAW or mzML conversion writes a per-scan status table next to the two CSV files. It is tab-separated, has a header row, and has one row for every MS2 and MS3 spectrum the reader read, including the spectra that were not exported. MS1 spectra have no row. A blank cell means the value is unknown or does not apply. A true value is written as `yes` and a false value as `no`. After you save the conversion the file is named `{savename}.scan_status.tsv` (see §5.4).

Columns in file order:

| Column | Allowed values | Meaning |
|--------|----------------|---------|
| `ms_level` | `2`, `3` | MS level of the spectrum. |
| `scan` | number | Scan number as the reader reports it. |
| `native_id` | text, or blank | Full native identifier of the spectrum (mzML only). Blank for RAW. |
| `exported` | `yes`, `no` | Whether the spectrum is in the CSV file. |
| `parent_scan` | number, or blank | Parent scan number. |
| `precursor_mz` | number, or blank | Precursor (isolation) m/z. |
| `charge` | integer, or blank | Charge used for the mass. Blank when unknown (the CSV holds `0`). |
| `charge_source` | `reported`, `comb`, `dataset_value`, `none` | Where the charge came from. `reported` is the charge in the instrument file. `comb` is the isotope-pattern estimate. `dataset_value` is the dataset charge you entered. `none` means unknown. Filled for exported MS2 rows. |
| `charge_status` | `known`, `unknown` | Filled for exported MS2 rows only. |
| `mass_status` | `eligible`, `unavailable` | `eligible` when a known charge gave a finite, non-zero precursor mass. Filled for exported MS2 rows only. |
| `arrays_status` | `ok`, `empty`, `length_mismatch`, `decode_error` | State of the peak arrays. `empty` marks a readable spectrum without peaks. `decode_error` marks a spectrum whose arrays could not be decoded. |
| `n_peaks` | integer, or blank | Number of peaks. Blank for `decode_error`. |
| `polarity` | `+`, `-`, or blank | Scan polarity. Blank when the file states none. |
| `polarity_source` | `cv`, `filter`, or blank | mzML only. `cv` is the polarity term of the spectrum. `filter` is the polarity token of the Thermo filter string. |
| `representation` | `centroid`, `profile`, or blank | Peak representation. RAW rows are always `centroid`. |
| `peak_source` | `label_data`, `mass_list`, `mass_list_vendor_centroid`, or blank | RAW only. Where the peaks came from. `label_data` is the Thermo label stream of an FTMS scan. `mass_list` is the vendor mass list of another analyzer. `mass_list_vendor_centroid` is that list after the vendor library centroided a profile scan. |
| `comb_status` | `not_run`, `no_isolation_mz`, `no_parent_ms1`, `no_profile_ms1`, `ft_ms1_not_tried`, `ms1_not_decodable`, `called`, `below_margin`, `no_signal`, `invalid_profile` | Result of the isotope comb for the spectrum. `not_run` means the comb is off or the charge was reported. `ms1_not_decodable` occurs for mzML only. Filled for exported MS2 rows. |
| `comb_z` | `1`, `2`, `3`, or blank | Charge the comb called. Blank when the comb did not call a charge. |
| `comb_best_z` | `1`, `2`, `3`, or blank | Charge with the highest comb score. |
| `comb_score` | number, or blank | Score of the called charge, rounded to four decimals. |
| `comb_margin` | number, or blank | Margin between the best and the next charge, rounded to four decimals. |
| `comb_anchor` | number, or blank | m/z of the anchor that gave the best comb score, rounded to four decimals. |
| `mass_mz_source` | `monoisotopic`, `isolation`, or blank | RAW only. The m/z used for the mass. `isolation` is used when the file gives no monoisotopic m/z. |
| `issues` | comma-separated list, or blank | Issue flags of the spectrum: `charge_unknown`, `precursor_mz_missing`, `parent_ref_missing`, `precursor_metadata_error`, `derived_mass_unavailable`, `arrays_empty`, `arrays_length_mismatch`, `decode_error`, `not_written`. `not_written` marks spectra that were read before a failed conversion stopped and whose CSV was never written. |

The comb columns `comb_z`, `comb_best_z`, `comb_score`, `comb_margin` and `comb_anchor` are filled only for spectra whose charge was not reported and for which the comb computed a score.

### Diagnostics report `<name>.reader_diag.json` and `.txt` (added in GlycoMSP 1.13)

Every conversion also writes a diagnostics report, even when the conversion fails. The `.json` file is the complete record. The `.txt` file is a human-readable summary of the same data: the outcome and what its exit code means, the explanation, the warnings and notes, the MS2 and MS3 counts, the charge sources, the exported peak representation and polarity, the coverage, the instrument, the settings, the file names, the failure details when there are any, and the run time. The reports carry the input file's stem and base names only. The input path is not written into them.

**Outcome and exit code.** The `outcome` field of the `.json` file holds one of three values. The fourth outcome, *complete with warnings*, is `complete` with exit code `3`. The `.txt` file and the log panel of the main window print it as "complete with warnings".

| `outcome` in the file | Exit code | Label in `.txt` and log | Meaning |
|-----------------------|-----------|-------------------------|---------|
| `complete` | `0` | complete | Converted, nothing to check. |
| `complete` | `3` | complete with warnings | Converted. Read the warnings: MS2 charges unknown, no usable precursor mass, or profile MS2 data. |
| `partial` | `4` | partial | Converted, but some spectra could not be decoded and were not exported. |
| `failed` | `5` | failed | The input could not be read, the conversion stopped, or no MS2 spectrum with peaks was exported. |

Exit code `2` (usage error, nothing converted) exists in the command-line tool only. No report is written for it. The graphical interface never saves a failed conversion as if it had succeeded. It saves a partial conversion with a warning.

**Top-level fields of the `.json` file**, in file order:

| Field | Content |
|-------|---------|
| `json_type` | `glycomsp.reader_diagnostics` |
| `schema_version` | `1.1.0` |
| `outcome` | `complete`, `partial` or `failed` (see above). |
| `exit_code` | `0`, `3`, `4` or `5`. |
| `exit_code_meaning` | The meaning of the exit code as a sentence. |
| `warnings` | List of warning sentences. Empty for exit code `0`. |
| `notes` | List of notes that are not warnings: spectra without peaks that were not exported, MS2 spectra without polarity, and charges that came from the isotope comb or the dataset charge. |
| `explanation` | Reason for a failed outcome or for a missing precursor mass. `null` otherwise. |
| `input` | `alias` (the file stem), `format` (`raw` or `mzml`) and `size_bytes`. |
| `tool` | `reader` (`msprawextractor` or `mzmlreader`), `reader_version`, `reader_last_update`, `shared_rules` and `entry` (`direct` for a conversion from the graphical interface, `cli-convert` for the command-line tool). |
| `environment` | `python`, `platform` and `backend` (the reading library and its version). |
| `settings` | The options of this run: `charge_comb`, `dataset_charge` and the comb parameters. For mzML also `round`. |
| `instrument` | Instrument information read from the file. |
| `representation` | Counts of exported spectra by `centroid`, `profile` and `unknown`, for MS2 and MS3. |
| `polarity` | Counts of exported spectra by `+`, `-` and `unknown`, for MS2 and MS3. |
| `counts` | `spectra_by_ms_level`, and for `ms2` and `ms3` the counts `seen`, `exported`, `empty`, `failed` and `not_written`. For `ms2` also `charge_known`, `charge_unknown`, `mass_eligible` and `charge_sources`. |
| `charge` | `sources`, `comb_status`, `comb_vs_dataset_disagree` and `mass_mz_source`. |
| `affected_scans` | For MS2 and MS3, each issue flag with the list of scans it affects. |
| `coverage` | `spectra_declared`, `spectra_processed`, `last_scan_processed` and `aborted`. |
| `failure` | `null`, or `stage`, `exception`, `message`, `last_scan_started`, `last_ms_level` and `traceback` of the error that stopped the conversion. |
| `outputs` | File names `ms2_csv`, `ms3_csv` and `status_tsv`, and `partial` (true when the outputs are incomplete). After you save the conversion the names are the final names. |
| `timing` | `started_utc`, `finished_utc` and `duration_s`. |

An exit code `3` warning is raised when at least one exported MS2 spectrum has no charge, when no exported MS2 spectrum has a usable precursor mass, or when exported MS2 spectra are profile data. Profile MS2 rows hold every profile point instead of peaks. Re-convert such a file with peak picking on MS2.

### Metadata sidecar

The metadata sidecar `{savename}.json` is written next to the scan CSV files when you save a conversion in the *Enter Experiment Metadata* window. It is a JSON object with these keys in file order:

| Key | Value |
|-----|-------|
| `Experiment Title` | Text you entered. |
| `Experiment Description` | Text you entered. |
| `Author running this analysis` | Text you entered. |
| `Raw data acquired date` | Text you entered. |
| `Glycan Type` | `N`, `O` or `N+O` in the conversion window. The CGA Setup window offers `N` and `O`. |
| `Mass Analyzer charge mode` | `+` or `-`. |
| `Derivatization Type` | One of the values below. |
| `Parameters when GlycoMSP launched` | List holding `Autofilled by GUI`. |
| `Date of file extracted from raw file` | Time stamp `YYYYMMDD_HHMMSS`. |
| `Conversion charge comb` | `true` or `false`: whether the isotope-comb option was on (added in GlycoMSP 1.13). |
| `Conversion dataset charge` | The dataset charge as an integer, or `null` when it was off (added in GlycoMSP 1.13). |
| `Original raw file path` | Path of the input file, or `not linked`. |
| `Raw filename` | File name without extension. |
| `json_type` | `glycomsp.metadata` |
| `schema_version` | `1.0.0` |

The two conversion entries are written only by a conversion and are absent from files that were produced without one. A conversion made with both options off records `false` and `null`.

*Derivatization Type* values: the conversion window offers `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `Others`. The CGA Setup window offers `PerMe`, `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `None`. When the value is `Others`, the MAS merge asks for a custom derivatization mass.

### Annotation Excel (MAS pipeline)

The expert-prepared annotation file consumed by the MAS pipeline. GlycoMSP reads two sheets from it.

- **`MSlist`** — manually assigned glycan compositions paired with MS2 scan numbers. The header row must equal one of four layouts, written here in column order. The newest layout has 12 columns: `Structure`, `MS2scan_no`, `peak`, `charge`, `mass shift`, `adduct`, `diff profile (integer for groups if different from same annotation)`, `IUPACname(optional)`, `Glycanannotation2`, `note`, `GlyToucan Access Number`, `WURCS`. The 11-column layout is the same up to `note` and then has `GlyToucan ID`. The 10-column layout stops at `note`. The basic 6-column layout is `Structure`, `MS2scan_no`, `charge`, `mass shift`, `adduct`, `note`. `Structure`, `MS2scan_no`, `charge`, `mass shift` and `adduct` are type-checked (text, integer, integer, integer, text). A column that is partly empty is rejected, and a column that is entirely empty is accepted. `GlyToucan Access Number` is read as `GlyToucan ID`. Consumed by `directassign_files`.
- **`ionlist`** — fragmentation ions used by GlycoMSP for feature extraction. Required column: `mass` as a floating-point number. Other columns are ignored.

The template workbooks also contain sheets named `Glycanlist` and `information_sheet`. GlycoMSP does not read them. A template is included in the GlycoMSP project folder.

### Ion list (CGA pipeline, also acceptable as MAS `ionlist` sheet)

CSV or Excel. One row per reference fragment. For Excel files GlycoMSP reads the sheet named `ionlist`, `ions` or `ion_list` (not case-sensitive), or else the first sheet that has a mass column.

| Column | Required | Recognised header variants | Purpose |
|--------|----------|---------------------------|---------|
| Mass | Yes | `mass`, `mz`, `ion_mz`, `m/z` (not case-sensitive) | Reference m/z for PPM or Da matching |

Other columns are ignored.

Two ion lists for negative-mode native N-glycans are included in the `templates/` folder (added in GlycoMSP 1.14): `GlycoMSP_ionlist_native_neg_14.xlsx` (the primary set) and `GlycoMSP_ionlist_native_neg_candidate_97.xlsx` (the extended set). Use them as the Score A ion list and as the feature ion list with *Native(Reduced)* in negative mode.

### In-silico glycan library (CGA pipeline)

Generated by **Generate In-Silico CSV** in the CGA Setup window, or supplied externally. The generated file is named `{sample_name}_insilico_{YYYYMMDD_HHMMSS}.csv`. It is a comma-separated file with the columns `Hex`, `HexNAc`, `NeuAc`, `NeuGc`, `KDN`, `Fuc`, `HexA`, `Mode`, `SO3`, `PO3H`, `Mass` in this order. The same columns are written for positive and negative mode. `Mode` holds `+` or `-`. `SO3` and `PO3H` hold the sulfate and phosphate counts and are `0` when the run has none. `Mass` is the theoretical precursor mass with four decimals. One row per theoretical composition consistent with the user-supplied flag set. Reusable across samples sharing derivatization and ionization mode.

### Score B workbook (optional, CGA pipeline)

Six-sheet Excel file (`.xlsx`, `.xlsm`, `.xltx`, or `.xltm`). See Section 5.6 for the high-level sheet structure; full authoring guidance is in `score_b_config_instructions.md` (manuscript supplementary).

## 5.2 Method JSON v1 schema

Produced by `build_method_v1_from_tree`; written by `save_method_v1_for_sample` and `_export_method_v1`. One file per (sample, family) pair.

```jsonc
{
  "json_type": "glycomsp.method",
  "schema_version": "1.0.0",
  "uid": "<uuid4>",
  "created_utc": "<iso8601>",
  "updated_utc": "<iso8601>",

  "method": {
    "family": "MAS" | "CGA",
    "name": "<exp_name>:<sample_name>:<family>",
    "description": "",
    "tags": []
  },

  "sample": {
    "sample_name": "<sample_name>",
    "experiment_title": "<exp_title>"
  },

  "inputs": {
    "converted_csv":         {"path": "<native-normpath>"},
    "metadata_json":         {"path": "<native-normpath>"},
    // MAS only:
    "annotation_excel":      {"path": "<native-normpath>"},
    // CGA only:
    "ion_list":              {"path": "<native-normpath>"},
    "insilico_glycan_list":  {"path": "<native-normpath>"},
    // CGA optional:
    "score_b_workbook":      {"path": "<native-normpath>"}
  },

  "parameters": {
    "mas": {},
    "cga": {},          // GlycoMSP 1.14, CGA family: see "CGA tolerances" below; {} from older builds
    "scoring": {},
    "ml": {}
  },

  "artifacts": {
    "reports": [],
    // CGA only — populated when corresponding sample fields are set:
    "CGAresult_tsv": {"path": "<native-normpath>"},
    "trainable_csv": {"path": "<native-normpath>"},
    "unlabeled_csv": {"path": "<native-normpath>"},
    // GlycoMSP 1.11, CGA or MAS — only when the reference-database fill ran: see "Reference-database provenance" below
    "glytoucan_resolution": {…}
  },

  "validation": {
    "status": "unknown",
    "checked_utc": null,
    "items": []
  },

  "software": {
    "glycomsp": {"version": "", "commit": ""},
    "extractor": {"name": "", "version": ""}
  },

  "operator": {
    "name": "",
    "note": ""
  }
}
```

**Key contracts:**

- `method.family` determines which `inputs` keys are populated. MAS reads `annotation_excel`; CGA reads `ion_list` + `insilico_glycan_list`. Both read `converted_csv` + `metadata_json`. `score_b_workbook` is optional and CGA-only.
- Path values are stored as **native** OS-specific paths via `os.path.normpath`. Windows produces `C:\path\file`; macOS / Linux produce `/path/file`. Cross-platform load is handled at read-time by `_resolve_path` (Windows-drive regex detection).
- Hash fields are **not currently written** by the writer; the schema reserves them for a future addition.
- `artifacts.CGAresult_tsv` was previously named `pseudolabels_tsv`; reads accept either key for backward compatibility, but writes always use the new name.
- The `parameters.mas`, `parameters.scoring` and `parameters.ml` sub-objects are placeholder shells in the current schema. They are reserved for future per-method parameter snapshots. GlycoMSP 1.14 fills `parameters.cga` (see below).

**CGA tolerances in `parameters.cga` (added in GlycoMSP 1.14).** A CGA method JSON written by GlycoMSP 1.14 records the tolerances used, at the defaults too, as soon as a CGA run or a *CGA → Trainable* build has been recorded for the sample. Its content therefore differs from a v1.10 method JSON for the same sample: the record is added provenance, and the data files (pseudolabel TSV, trainable CSV) are unchanged at the defaults. A method JSON from an older build keeps `"cga": {}` and loads as before.

```jsonc
"cga": {
  "run": {                                   // written by Start CGA analysis
    "tsv": "<pseudolabel TSV path>",
    "precursor_tolerance": {"value": 20.0, "unit": "ppm"},   // or "Da"
    "fragment_tolerance":  {"value": 10.0, "unit": "ppm"},   // Score A
    "tolerance_non_default": false
  },
  "trainable": {                             // written by Build Trainable CSV
    "csv": "<trainable CSV path>",
    "csv_sha256": "<hex>",                   // the record is used only while the CSV has this hash
    "cga_tsv": "<TSV the CSV was built from>",
    "cga_run": {"precursor_tolerance": {…}, "fragment_tolerance": {…}, "tolerance_non_default": false},
                                             // null when cga_tsv is not the TSV of the recorded run
    "precursor_gate_ppm": 10.0,              // null = gate unticked
    "max_abs_ppm": 20.0,
    "min_ion_score": 0.07,
    "feature_ion_ppm": 10.0,                 // "ion list ppm"
    "negatives": {"enabled": true, "ion_ppm": 10.0},
    "tolerance_non_default": false           // gate, max |ppm_error|, feature window or the CGA run non-default
  }
}
```

On load, GlycoMSP keeps only the well-formed parts of this block. A `run` part needs two valid tolerances. A `trainable` part needs `csv`, `csv_sha256`, a positive `max_abs_ppm`, a positive `feature_ion_ppm` and a `precursor_gate_ppm` that is `null` or positive. An unusable part is dropped and its values are treated as unknown.

The run log of the sample (`pseudolabel_run` entries) also carries `precursor_tolerance`, `fragment_tolerance` and `tolerance_non_default`. Its old `ppm_value` / `ion_ppm` keys are `null` when that window was set in Da.

**Reference-database provenance in `artifacts.glytoucan_resolution` (added in GlycoMSP 1.11).** When the option *Fill GlyToucan ID / WURCS from reference DB* is on and a fill ran, the CGA or MAS method JSON of the sample gains this block. The fill runs in the CGA analysis, in the MAS merge, and in *CGA → Trainable*. The block describes the most recent fill for the family. A new CGA analysis or MAS merge that is produced with the option off removes the block. A *CGA → Trainable* build with the option off leaves it in place. `schema_version` stays `1.0.0`, and a method JSON produced with the option off is identical to one from v1.10.

```jsonc
"glytoucan_resolution": {
  "status": "filled",               // "filled" | "partial" | "none" | "db_missing"
  "route": "database",              // the fill reads the local database only, never the network
  "db": {                           // the database that filled the cells
    "name": "<database name>",
    "version": "<database version>",
    "path": "<database path>",
    "sha256": "<hex>",              // hash of the database file at the time of the fill
    "n_entries": <integer>,
    "updated_utc": "<iso8601>"
  },
  "n_compositions": <integer>,      // distinct compositions found in the composition column
  "n_filled": <integer>,            // of these, compositions that have a database row with a WURCS
  "n_missing": <integer>,           // compositions without such a row
  "n_user_accession": <integer>,    // rows whose own accession differs from the database record: their WURCS is left as it is
  "n_user_wurcs": <integer>,        // rows with their own WURCS and a blank accession: the accession stays blank
  "last_resolved_utc": "<iso8601>"
}
```

`status` is `filled` when every composition has a database row, `partial` when some have none, `none` when the table has no composition to look up, and `db_missing` when the database could not be read. Cells that already hold a value are never overwritten. The prediction fill (§5.4) does not write this block.

A loaded method whose family is recognised as UNKNOWN hides its `artifacts` until *Export Method v1 (MAS)* is used once. This is a known issue that predates GlycoMSP 1.11.

## 5.3 Experiment JSON schema

Produced by `export_experiment_json`; one file per experiment. Holds the experiment's organisational state — sample list, per-sample method references, ML defaults, and policies. GlycoMSP 1.11 to 1.14.1 added no keys to this file.

```jsonc
{
  "json_type": "glycomsp.experiment",
  "schema_version": "1.0.0",
  "created_utc": "<iso8601>",
  "updated_utc": "<iso8601>",

  "experiment": {
    "title": "<exp_title>",
    "description": "",
    "tags": []
  },

  "workspace": {
    "preferred_method_folder": "",
    "preferred_exp_path": "",
    "ui_order": ["<sample_name>", ...]
  },

  "samples": [
    {
      "sample_name": "<sample_name>",
      "validated": true | false,
      "methods": [
        {
          "method_id": "<uuid4>",
          "method_name": "<basename>",
          "family": "MAS" | "CGA" | "UNKNOWN",
          "method_ref": {
            "path": "<posix-str>",
            "hash": {"alg": "sha256", "value": "<hex>"},
            "last_seen_utc": "<iso8601>"
          },
          "validation": {
            "status": "ok" | "warn" | "missing" | "incomplete" | "unknown",
            "checked_utc": "<iso8601>",
            "notes": "<str>"
          }
        },
        ...
      ]
    },
    ...
  ],

  "validation_policy": {
    "hash_alg": "sha256",
    "allow_missing_method": true,
    "allow_hash_mismatch": true
  },

  "ml_defaults": {
    // optional; written/read by the ML params editor's
    // "Apply to exp" / "Pull from exp" buttons
  },

  "train_parameters": {
    // optional; written by train_model after a successful training run
  }
}
```

**Key contracts:**

- Path values inside `method_ref.path` are stored as **POSIX** strings via `pathcanon.to_posix_str`, regardless of host OS. This keeps experiment JSON files diff-clean across hosts.
- `hash` is omitted when the method file does not exist at save time, and empty values are not written.
- `validation_policy.allow_hash_mismatch: true` is the default — a hash mismatch on import sets the method's status to `"warn"` but loads the method anyway. Set to `false` to require exact hash match.
- `ml_defaults` and `train_parameters` are optional top-level keys populated by the ML Analysis window. Their internal structure follows the merged ML-params dictionary schema and is not validated by GlycoMSP beyond JSON parseability. `train_parameters` holds `split_ratio_test`, `split_ratio_val`, `min_samples`, `balancing_enabled`, `majority_label`, `majority_factor`, `stratify`, `rf_n_estimators`, `rf_class_weight`, `real_world_test_cap_training_only`, `threshold_enabled`, `threshold_tau`, `threshold_margin`, `model_path`, `labelencoder_path` and `report_path`.
- `linked_validated_samples` (in-memory) round-trips via the per-sample `validated: true/false` flag. Renaming a sample post-validation without re-saving the exp JSON loses the flag.

## 5.4 Output file specifications

### Conversion outputs (Step 2.2)

A conversion that you save in the *Enter Experiment Metadata* window gives these files in the output folder. `{savename}` is `{project name}_{YYYYMMDD_HHMMSS}_{input file stem}`, where the project name is the one you type when the window asks for it.

| File | Content |
|------|---------|
| `ms2_{savename}.csv` | MS2 scan CSV (§5.1). |
| `ms3_{savename}.csv` | MS3 scan CSV (§5.1). |
| `{savename}.json` | Metadata sidecar (§5.1). |
| `{savename}.scan_status.tsv` | Per-scan status file (§5.1). Added in GlycoMSP 1.13. |
| `{savename}.reader_diag.json` | Diagnostics report, complete record (§5.1). Added in GlycoMSP 1.13. |
| `{savename}.reader_diag.txt` | Diagnostics report, text summary (§5.1). Added in GlycoMSP 1.13. |
| `{savename}.log` | Log of the conversion. |

While the conversion runs, the readers write `ms2tmp_{stem}.csv` and `ms3tmp_{stem}.csv` and the diagnostics files under the input file's stem. Saving moves them to the names above, and the `outputs` entries of the `.json` report are rewritten to the final names. A failed conversion is not saved. A partial conversion is saved with its warning.

### Pseudolabel TSV (CGA, output of `run_pseudolabeling`)

Tab-separated. The file carries every converted scan, so its size follows the converted file (about 1.7 GB for a run of 100,000 MS2 scans). A matched-rows-only output is planned. A scan with precursor matches has one row per candidate composition, at most the top 3 per scan. A scan without a match keeps one row with empty annotation cells.

| Column | When present | Notes |
|--------|--------------|-------|
| `entry_no`, `MS1scan_no`, `MS1_isolationmass`, `MS1_monoisolationmass`, `chargeState`, `protonatedmass`, `MS2scan_no`, `peaklist`, `peakintensity` | always | The columns of the converted scan CSV, unchanged. |
| `composition` | always | Compact composition label of the candidate, e.g. `F1H5N4`. Empty for a scan without a match. |
| `theoretical_mass` | always | Theoretical precursor mass of the candidate. |
| `ppm_error` | always | Precursor error against this candidate's theoretical mass, in ppm. |
| `ion score` | an ion list is set | Score A: ion-fitting fraction in `[0, 1]` |
| `ion hit count` | an ion list is set | Integer count of ion-list fragments hit within tolerance (counted as peak–ion matches, see §4.1, *How hits are counted*) |
| `ion hits m/z` | an ion list is set | Semicolon-separated m/z list of hits |
| `Glycanannotation2` | always | Empty unless the input table already carried a value. |
| `GlyToucan ID` | always | Empty unless filled. With the option *Fill GlyToucan ID / WURCS from reference DB* on, empty cells are filled from the reference database. Values already present are kept. Added in GlycoMSP 1.11. |
| `WURCS` | always | As `GlyToucan ID`. |
| `score_b` | Score B active | Float in `[0, 1]` |
| `score_b_rank` | Score B active | Per-scan rank by `score_b` (1 = highest) |
| `score_b_selected` | Score B active | Boolean. `true` for the per-scan top-ranked candidate |
| `score_b_support_score` | Score B active | Weighted motif support of the candidate. |
| `score_b_composition_penalty` | Score B active | Composition-consistency penalty of the candidate. |
| `score_b_unexpected_penalty` | Score B active | Unexpected-evidence penalty of the candidate. |
| `score_b_motif_summary` | Score B active | Short text rationale |
| `selection_source` | Score B active | `score_b` for scored rows. |
| `selected_composition` | Score B active | Composition of the top-ranked candidate of the scan, repeated on every scored row of that scan. |

Score B adds its columns to the same file. The columns appear in the order listed. Rows that Score B does not score keep `score_b` empty and `score_b_selected` `false`.

Filename: `{sample_name}_CGA_{YYYYMMDD_HHMMSS}.tsv` next to the converted CSV.

### Trainable CSV (output of `build_trainable_from_CGA` or MAS merge)

Comma-separated. One row per spectrum.

**CGA path** (*CGA → Trainable*), columns in file order:

| Column | Notes |
|--------|-------|
| `MS2scan_no` | From scan CSV |
| `Structure` | Composition label (e.g., `F1H5N4`) for positive rows; literal `Non-glycan` for sampled negatives. |
| `Glycanannotation2` | Empty unless the CGA TSV carried a value. |
| `GlyToucan ID` | Empty unless the CGA TSV carried a value or the fill option is on. With the option on, empty cells are filled from the reference database (added in GlycoMSP 1.11). |
| `WURCS` | As `GlyToucan ID`. |
| `protonatedmass` | Present only when *Include precursor mass as feature* is enabled, and then first among the feature columns. Blank for Non-glycan rows. The Train step reads blank and Non-glycan values as `0.0` when *Use protonated mass as a model feature (Non-glycan = 0)* is on. |
| `<ion_mz_1>`, `<ion_mz_2>`, ..., `<ion_mz_N>` | Per-fragment feature columns. Cell value is `log10(I) + 1` where I is the max intensity in the PPM window, or `1.0` if no peak matched. Column names are the reference ion m/z values formatted as floats. |

Filename: `{sample_name}_trainable_fromCGA_{YYYYMMDD_HHMMSS}{_withMass|_noMass}.csv`.

**MAS path** (MAS merge), columns in file order: `protonatedmass`, one column per ion-list mass (`<ion_mz_1>` … `<ion_mz_N>`, value as above), `Structure`, `IUPACname(optional)`, `Glycanannotation2`, `GlyToucan ID`, `WURCS`, `unique_ID`. `unique_ID` holds the MS2 scan number. Values already in your annotation sheet always win over the reference database. With the fill option on, only empty `GlyToucan ID` and `WURCS` cells are filled. Filename: `{sample_name}_merged_{YYYYMMDD_HHMMSS}.csv`.

**Combined datasets** (*Combine Datasets for Training*): the combined file is the union of the columns of its sources. It starts with the columns `UID`, `Origin_Basename` and `Origin_File`. `UID` has the format `{origin_basename}#scan{scan_no}`, unless the source already had a `uid` column. A short text log `{combined name}_combine_log.txt` is written beside it.

### Predictions CSV and supporting outputs

**Unlabeled prediction input (the input of a prediction run).** *Create unlabeled dataset of certain experiment* writes `{converted CSV stem}_unlabeled_ppm{ppm}.csv` next to the converted CSV. Its columns are `protonatedmass`, `MS2scan_no`, one feature column per ion (named by the ion m/z, as in the trainable CSV), and `UID`.

Each successful prediction run creates a folder next to the unlabeled input CSV. Its name is `{method_base}{ms1_suffix}{gate_suffix}`. `{method_base}` is the name of the sample's `.method.json` file without `.json`, or the stem of the unlabeled CSV when no method file is found. `{ms1_suffix}` is `_MS1feat` when the model uses the mass feature and `_noMS1` otherwise. `{gate_suffix}` is `_PG{N}ppm` when the precursor gate is on, where `N` is the gate tolerance as a whole number with the decimals dropped, and `_noPG` when it is off. An example is `{method_base}_MS1feat_PG10ppm/`. In the file names below, `{input_stem}` is the stem of the unlabeled CSV.

Folder contents:

- **`{input_stem}_predicted{ms1_suffix}{gate_suffix}.csv`** — the complete prediction table with every row of the input. It holds the input columns plus `Predicted_Label`, `ppm_precursor` and `pred_ok`. With the option *Fill GlyToucan ID / WURCS from reference DB* on, it also holds `Predicted_GlyToucan_ID` and `Predicted_WURCS` (added in GlycoMSP 1.11). The fill uses the local reference database only and leaves `Predicted_Label` unchanged. `ppm_precursor` is the precursor error in ppm against the closest library mass of the predicted composition, and is empty when no gate ran. With the gate on, `pred_ok` is `True` for glycan rows that pass the gate and `False` for Non-glycan rows and for glycan rows that fail it. With the gate off, `pred_ok` is `True` for every row.
- **`predictions_with_conf.csv`** — main prediction report table. When the precursor gate ran it holds only the glycan rows that passed the gate, and the counts in `class_summary.csv` and `run_summary.json` refer to these rows. Columns: the columns of the full table with `Predicted_Label` renamed `pred_label`, then `proba_vector`, `max_proba` (top-class probability), `second_label`, `second_proba`, `margin` (`max_proba` minus `second_proba`), `entropy`, `accepted`, `reject_reason` (`low_conf`, `low_margin`, `low_conf_and_margin`, or empty) and, as two comma-separated strings, `topk_labels` and `topk_probas`. A `UID` column from the input is passed through.
- **`class_summary.csv`** — one row per predicted class with `pred_label`, `pred_count`, `accepted_count`, `mean_max_proba`, `median_max_proba`, `mean_margin`, `median_margin`, `mean_entropy`, `median_entropy` and `accepted_rate_within_class`.
- **`by_sample_summary.csv`** — written only when the input has an `experiment_title` or `sample_name` column. One row per experiment and sample with counts, accepted rate, mean `max_proba`, margin and entropy, and the three most frequent classes.
- **`run_summary.json`** — summary of the run: `n_total`, `n_accepted`, `n_rejected`, `tau`, `margin`, `max_proba_stats`, `margin_stats` and `entropy_stats` (each with `mean`, `median`, `p25` and `p75`), `top_classes_by_count`, `top_classes_by_accepted`, `ambiguous_count`, `low_conf_count`, `class_names`, `class_summary`, `method_suffix` and `model_tolerances`. GlycoMSP 1.14 adds `model_tolerances` with these keys: `prefilled` (the values prefilled from the model's record, non-default values only), `gate_requested`, `gate_ppm_requested`, `gate_applied`, `gate_ppm_used`, `gate_request_differs_from_model`, `gate_outcome_matches_model` (`null` when the model has no non-default gate record) and `feature_window`. The feature window is not checked here. It was set when the unlabeled CSV was built.
- **`run_summary.md`** — short human-readable summary: totals, τ and margin, statistics of `max_proba`, margin and entropy, and the eight most frequent classes.
- **`to_review_low_conf.csv`** — rows whose `max_proba` is below τ.
- **`to_review_low_margin.csv`** — rows whose `margin` is below the configured margin.
- **`to_review_ambiguous_topk.csv`** — accepted rows whose runner-up probability lies within the margin of the top probability.

When thresholding is off (the default), the summary reports τ 0 and margin 1. Whatever the threshold setting, it counts no row as accepted (*Accepted 0*) and shows probability statistics of 0.000 (known issue 29). The predicted labels and per-class counts are correct. The summary will be corrected (changelog 1.14, known limitations).

**Prediction gate files.** These are written only when the precursor gate is on and ran. They are listed here by exact name pattern.

| File | Content |
|------|---------|
| `{input_stem}_predicted{ms1_suffix}_PG{N}ppm_only.csv` | The rows of the complete prediction table with `pred_ok` equal to `True`. |
| `{input_stem}_gate_passed.csv` | The glycan rows that passed the gate, with the columns `MS2scan_no`, `pred_label`, `canon_pred`, `protonatedmass`, `ppm_precursor`, `masses_found` and `kept`. `canon_pred` is the canonical composition label and `masses_found` the number of library masses for it. |
| `{input_stem}_gate_debug.csv` | The glycan rows that were dropped by the gate, with the same columns. Use it to see why spectra were demoted. |
| `{input_stem}_gate_nolut_bycounts_debug.csv` | Written only when some predicted compositions have no entry in the library lookup. Lists up to the first 30 such rows with `MS2scan_no`, `pred_label`, `canon_pred`, `protonatedmass` and `Mass_candidates` (library masses found by residue counts). |
| `{input_stem}_PG_counts.csv` | One row per class among the rows that passed the gate, with `pred_label` and `count`, sorted by count. |

### Trained model and sidecars (output of `train_model`)

`<dataset_stem>` is the path of the trainable CSV without `.csv`, followed by `_withMass` when the model uses the mass feature.

| File | Contents |
|------|----------|
| `<dataset_stem>_rf_model.skops` | The trained `RandomForestClassifier` serialized via skops. If skops is not installed, the model is saved as `<dataset_stem>_rf_model.joblib`. |
| `<dataset_stem>_labelencoder.joblib` | Fitted `LabelEncoder` mapping class strings to integer labels. |
| `<dataset_stem>_features.json` | Ordered list of feature column names. Critical for predict-time alignment. |
| `<dataset_stem>_rf_performance.txt` | Performance report on the held-out test set: the per-class classification report, followed by the glycan-only report that excludes the majority label. |
| `training_run.json` | Reproducibility manifest, written into the same folder under this fixed name. A second training run in the same folder overwrites it. See Section 4.5. Its `training_run` object holds `effective_params`, `column_order`, `classes`, `versions`, `seeds`, `inputs`, `hashes`, `created_at` and `cga_tolerances`. |

Models saved by older builds may use `.joblib` or `.pkl` extensions instead of `.skops`. The application's loader handles all three; sidecar paths are derived from the model path with the appropriate extension.

### Random Forest timing report `rf_timing_<alias>.txt`

This report is written by the validation script `tools/validation/rf_timing.py`, not by the application. The script trains a Random Forest on a trainable CSV exactly as the *Train* button does for a trainable with the mass feature, with all other Train settings at their defaults, and times it.

```bash
python tools/validation/rf_timing.py --trainable <trainable CSV> --alias <name> [--repeats 3] [--trees 400] [--out <folder>]
```

The file name is `rf_timing_{alias}.txt`. It is written to the folder given by `--out`, or to the folder of the trainable CSV when `--out` is not given. The alias is 1 to 64 characters of letters, digits, `_`, `.` and `-`, and starts with a letter or digit. The report carries the alias only and never the input path. It has these blocks:

| Block | Content |
|-------|---------|
| `[data]` | Rows in the file, rows after the class-size filter, number of features, whether the mass feature is included, classes kept and dropped (below 5 spectra), training rows before and after the cap on the majority class, and test rows. |
| `[random forest, as the product sets it]` | `n_estimators` (400 unless `--trees` is given), `class_weight`, `random_state`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `n_jobs` and the split. |
| `[environment]` | Python, scikit-learn, numpy and pandas versions, platform, processor string, CPU model and number of logical cores. |
| `[timing, N repeat(s), seconds, time.perf_counter()]` | Seconds per repeat for the fit alone and for the whole step from CSV load to fitted model, with the mean and maximum of both. |
| `[held-out accuracy of the last fitted model, threshold off]` | Multiclass accuracy and, when a Non-glycan class is kept, the glycan versus Non-glycan binary accuracy. |

### Reference database files

`<data dir>` is `%LOCALAPPDATA%\GlycoMSP` on Windows and `~/.glycomsp` on macOS and Linux, where `~` is your home folder.

| File | Location | Content |
|------|----------|---------|
| `glycomsp_reference.db` | `<data dir>/references/` (default). Override with the environment variable `GLYCOMSP_REFERENCE_DB`. | The local reference database. A SQLite file with a `meta` table and an `entries` table. An entry holds the canonical composition key, WURCS, accession, status (`resolved`, `wurcs_only` or `unresolvable`), source and the times it was added and last checked. The *Manage References window* can select another database file for the session. |
| `glytoucan_cache.jsonl` | `<data dir>/` (default). Override with the environment variable `GLYCOMSP_GLYTOUCAN_LOG`. For a selected database that is not the default one, the log is `<database folder>/<database name>.jsonl`. | JSONL transaction log. Every successful fetch of the prebuild is appended as one line with `ts`, `op`, `composition_key`, `query`, `accession`, `wurcs`, `status`, `endpoint`, `db_name` and `db_version`. |
| `glytoucan_cache.jsonl.bak` | Next to the log. | Copy of the previous log, kept when *Rebuild log* (or the `--rebuild-log` option of the prebuild tool) regenerates the log from the database. An older `.bak` is overwritten. |
| `settings.json` | `<data dir>/` | Written only when you press *Set as default* or tick *Remember these settings* in the Manage References window, never by browsing. Keys: `schema_version`, `updated_utc` and a `glytoucan` object with `remember`, `fill_on_run` and `db_path`. At start-up the file is applied only while `remember` is `true`. A remembered `db_path` is ignored while `GLYCOMSP_REFERENCE_DB` is set. |

## 5.5 Default parameters

A single reference table for every default value the application uses out-of-the-box. Source: the application code and screenshots. Manuscript-published values are flagged where they diverge. Every row of the table was re-verified against the GlycoMSP 1.14.1 code (the application code of 1.14.1 is identical to 1.14).

| Parameter | Default | Where set | Notes |
|-----------|---------|-----------|-------|
| Precursor match tolerance (CGA matching) | 20 ppm | CGA Setup *Precursor match* → `run_pseudolabeling` `ppm_value` / `tol_da` | Relative to the observed precursor mass. GlycoMSP 1.14: field, ppm or Da. v1.10: fixed, no field |
| Top-N candidates per scan (CGA) | 3 | `run_pseudolabeling` `keep_top_n_per_scan` | Caps candidates by smallest `\|ppm_error\|` |
| Score A fragment match tolerance | 10 ppm | CGA Setup *Score A fragment match* → `run_pseudolabeling` `ion_ppm` / `ion_tol_da` | GlycoMSP 1.14: field, ppm or Da. v1.10: fixed, no field. Hits are peak–ion matches (§4.1) |
| Score A anchors required | not used | `run_pseudolabeling` `anchors_required` (value 2) | The parameter exists but is not read. Score A has no anchor gate (§4.1) |
| Score B fragment match tolerance | 20 ppm | fixed in the call of `run_score_b_enrichment` | No field. Not changed by the *Score A fragment match* field |
| CGA → Trainable precursor gate | on, 10 ppm | *precursor gate* checkbox and entry | Relative to the theoretical mass; ppm only. GlycoMSP 1.14: field. v1.10: applied with no field (§2.6) |
| CGA → Trainable max \|ppm_error\| | 20 ppm | *max \|ppm_error\|* field | Filters the TSV `ppm_error`; stays active when the gate is unticked; ppm only |
| CGA → Trainable min ion score | 0.07 | *min ion score* field | |
| CGA → Trainable Top N per scan | 1 | *Top N per scan* field | |
| CGA → Trainable feature window (*ion list ppm*) | 10 ppm | *ion list ppm* field → `build_features_from_peaks_log10_plus1` `ppm` | MS2 features and this modal's Non-glycan sampling; ppm only |
| CGA → Trainable negatives | on, ratio 3.0, < 3 ion hits | *5) Negative sampling* box | Its own settings, not the shared dialog below |
| CGA → Trainable mass feature | off | *Include precursor mass as feature* | |
| `build_features_from_peaks_log10_plus1` function default | 20 ppm | `ml_ng_utils` | Not used by the GUI: CGA → Trainable always passes its *ion list ppm* |
| No-peak default cell value | 1.0 | feature extraction | `log10(1) + 0 = 0` would lose info; +1 keeps cells distinguishable |
| Negative sampling: enabled | False | Negative Sampling Setting (`add_negatives_var`) | Shared dialog, used by MAS merge and the ion-mining pre-negatives of CGA analysis; not by CGA → Trainable |
| Negative sampling: ratio | 3.0 | `neg_ratio_var` | neg:pos cap |
| Negative sampling: min ion-hit gate | 3 | `min_hits_var` | Scans below this become negatives |
| Negative sampling: PPM tolerance | 10 ppm | `ppm_tol_var` | Shared dialog only; ppm only |
| Ion mining: enabled | False | `ion_suggest_enable_var` | Opt-in |
| Ion mining: PPM tolerance | 10.0 ppm | `ion_suggest_ppm_var` | |
| Ion mining: Da floor | 0.030 Da | `ion_suggest_dafloor_var` | Low-m/z merge floor |
| Ion mining: min support | 5 | `ion_suggest_minsupp_var` | Per glycan class |
| Ion mining: top-K | 60 | `ion_suggest_topk_var` | Suggestions exported |
| Min samples per class | 5 | `train_model` | Drop classes smaller than this |
| Max negative cap | 3× largest glycan class | `train_model` | Manuscript baseline; configurable |
| Test split | 0.20 | `test_split_var` | |
| Validation split | 0.00 | `val_split_var` | No validation fold by default |
| Stratified split | True | `stratify_var` | |
| RF n_estimators | 400 | `n_estimators_var` | |
| RF max_depth | None | RF defaults | Unrestricted depth |
| RF min_samples_split | 2 | RF defaults | |
| RF min_samples_leaf | 1 | RF defaults | |
| RF max_features | `sqrt` | RF defaults | |
| RF class_weight | `balanced` | `class_weight_var` | |
| RF random_state | 42 | fixed | |
| Training mass feature | off | *Use protonated mass as a model feature (Non-glycan = 0)* (`include_mass_train_var`) | |
| Confidence threshold τ | 0.65 | `thresh_val_var` | Only when threshold enabled. Four τ defaults exist in the code. The table gives the GUI value (§7.2) |
| Confidence margin | 0.05 | `margin_val_var` | Only when threshold enabled |
| Threshold enabled | False | `enable_thresh_var` | Disabled by default |
| Feature encoding ppm (training) | 10 ppm | CGA → Trainable *ion list ppm*; MAS merge fixed (no field) | ppm only |
| Feature encoding ppm (prediction) | none: you type it (the dialog suggests 20) | Unlabeled-dataset dialog → `extract_ion_intensities` (function default 20) | GlycoMSP 1.14: 1–500 ppm with decimals, prefilled from a model's non-default record. v1.10: whole numbers 1–200 |
| Precursor gate (predict) | True | `apply_pred_precursor_gate_var` | Opt-out at prediction time. GlycoMSP 1.14: unticked on load when the model's record has the gate off |
| Precursor gate PPM (predict) | 10 ppm | `pred_precursor_ppm_var` | Relative to the theoretical mass of the predicted composition. GlycoMSP 1.14: prefilled from the model's non-default gate record |
| Validation policy: allow missing method | True | exp JSON | |
| Validation policy: allow hash mismatch | True | exp JSON | |
| Conversion option: *Calculate the charge from the MS1 isotope pattern (charge 1–3)* | off | Conversion options dialog | Added in GlycoMSP 1.13. The dialog starts with both options off every time. |
| Conversion option: *Dataset charge for the rest (blank = off)* | blank (off) | Conversion options dialog | Added in GlycoMSP 1.13. A whole number other than 0. |
| *Fill GlyToucan ID / WURCS from reference DB* | off | Prepare Dataset checkbox and Manage References window | Added in GlycoMSP 1.11. Remembered only through *Remember these settings*. |
| Reference database file | `glycomsp_reference.db` in `<data dir>/references/` | `GLYCOMSP_REFERENCE_DB`, or *Select database…* and *Set as default* in Manage References | See §5.4. |
| Prebuild batch size | 30 (dialog and command line) | *Batch size (max 30)* in the *Prebuild / Update from source* dialog, `--batch-size` | The dialog accepts 1 to 30. |
| Prebuild pause between batches | 10 s (dialog and command line) | *Pause between batches (s)* in the dialog, `--pause` | The dialog accepts 0 or more. |
| Prebuild request timeout | 20 s | `--timeout` of `msp_CLI_glytoucan_prebuild.py` | The dialog has no field for it and uses the command-line default. |
| Prebuild *Diff only — no network* | off | *Diff only — no network (--no-network)* in the dialog, `--no-network` | |
| Decoy offset (validation script only) | +3.37 Da | `--shift` of `tools/validation/cga_tolerance_decoy_count.py` | Not part of the application. |

Known issue 34 concerns the four different τ defaults in the code. The table states the GUI values. §7.2 describes the issue.

## 5.6 Score B workbook structure (overview)

The Score B Excel workbook contains six required sheets. The schema is sufficient detail for **identifying** that you have a workbook in the right shape; **authoring** a workbook from scratch requires reading `score_b_config_instructions.md` (manuscript supplementary).

| Sheet | Purpose |
|-------|---------|
| `ionlist` | Defines the named glycotopes (motifs) and their reference fragment masses used by all other sheets. Required columns: `glycotope_name`, `ion_type`, `fragmentation_mass`, `structural_identifier`, `charge_mode`, `derivatization`. |
| `interpretation_v2` | The rules. Direct evidence rules score positively when their listed fragments are observed. Gated rules apply only when a specified precondition has fired, which encodes evidence chains. Required columns: `display_name`, `target_id`, `target_level`, `target_group`, `parent_target_id`, `ion_struct_id`, `ion_struct_target_id`, `ion_role`, `logic`, `min_hits`, `weight`, `charge_mode_required`, `derivatization_required`. |
| `TargetHierarchy` | Parent and child relations between targets. Required columns: `parent_target_id`, `child_target_id`. |
| `MotifPolicy` | Which targets a candidate's flags activate. Required columns: `flag_id`, `target_id`, `activation_mode`, `include_descendants`. |
| `CompositionConsistency` | Penalties subtracted when observed evidence contradicts the candidate's composition. Required columns: `rule_id`, `composition_field`, `when_value`, `evidence_target_id`, `penalty`. |
| `UnexpectedEvidencePenalty` | Penalties subtracted when observed evidence is not accounted for by any expected motif rule applicable to the candidate. Required columns: `flag_id`, `target_id`, `include_descendants`, `penalty`. |

Workbook-level tags for the derivatization and the charge mode are columns of `ionlist` and `interpretation_v2`. The application validates the sheets and columns on workbook load.

Validation runs via `msp_CGA_structscore.ScoreBLoader.load`. Click **Validate Score B Workbook** in the CGA Setup window to surface load errors before running CGA analysis.

Two changes were made in GlycoMSP 1.12. Composition parsing is case-sensitive, so a lowercase `s` or `p` in a composition label is sulfate or phosphate and is not read as NeuAc or another residue. Tokens that Score B does not know are left out of scoring. The labels that carry them are reported once at the end of the run, in one warning line on the console, and the run does not stop (§3.3).

---

# Part 6 — Configuration & advanced

For most users, accepting the defaults from Part 2 is enough. This Part covers the controls you reach for when the defaults don't match your data — typically because your dataset is unusually large, unusually small, drawn from a non-typical instrument configuration, or being prepared for a comparison study where parameter choices matter.

## 6.1 The ML parameters editor

The Train/Test Parameters dialog (Quickstart Step 2.7) covers the most common adjustments. For deeper configuration — anything beyond what the dialog exposes, or any parameter you'd like to version-control as a JSON preset — the application provides a richer editor.

From inside the Train/Test Parameters dialog, click the *Edit ML Parameters…* button. The Train tab's *Set Parameters / Train the Model* button opens that dialog. This opens a 1200×600 modal window titled *ML Parameters* with two side-by-side panels.

> **[SCREENSHOT: `GlycoMSP_screenshots/ML_Params_Editor.jpg` — ML Parameters window.]**
> Window titled *ML Parameters*, two panels. Left: editable JSON text area in a frame titled *Editor (JSON)*. Right: read-only JSON preview in a frame titled *Effective (Builtins ← Files ← Editor)*, recomputed live as you type, with the line *live: test=0.20, val=0.00, trees=400* below it. Buttons at the top right in three rows: *Load from method…*, *Save to method…*, *Close*, then *Validate & Use*, *Load preset…*, *Save preset as…*, then *Pull from experiment*, *Apply to experiment*, *Link exp manually*.

**The merge model.** The "effective" parameters used at training time are computed from up to three layers:

1. The current values of the Train/Test Parameters dialog.
2. The experiment JSON's `ml` block (its `ml.parameters` entry if present), if an experiment JSON is linked.
3. The editor pane's JSON, if non-empty.

Each layer overrides the previous. The right-hand pane shows the final merged result so you can see exactly what the model will train with before you commit. Training uses the result of *Validate & Use* if you clicked it, and otherwise the current values of the Train/Test Parameters dialog.

**Buttons.**

- **Pull from experiment** — reads `ml_defaults` from the linked experiment JSON into the editor pane. Lets you start from a previously-saved exp-level baseline.
- **Apply to experiment** — writes the editor JSON back into `ml_defaults` of the linked experiment JSON. Useful when an experiment-wide standard parameter set is what you want. The merge layer listed above reads the `ml` block and not `ml_defaults`, so load such a set back with *Pull from experiment*.
- **Load preset…** — picks any JSON file via dialog, loads its content into the editor (accepts both whole-exp payloads — `ml_defaults` is extracted — and plain ML-params files).
- **Save preset as…** — writes the current editor JSON to a file. Use for sharing parameters across labs or version-controlling them in git.
- **Link exp manually** — re-binds the editor to a different `.exp.json` so subsequent Pull/Apply operate on it. Selecting a training CSV in the Train tab links an `.exp.json` from the same folder automatically.
- **Load from method…** — reads `ml.parameters` from a method JSON.
- **Save to method…** — writes the merged effective parameters to a method JSON's `ml.parameters` block (with a local-time `ml.updated` timestamp).
- **Validate & Use** — re-runs the merge with the editor as the rightmost layer, rebinds `effective_ml_params` (the closure variable read by the trainer), and persists into the in-memory ML state. The window does NOT auto-close after — you can verify, adjust, and re-validate as needed.
- **Close** — closes the editor without persisting (your changes since the last *Validate & Use* are lost).

**Recommended workflow for parameter studies.** Author a JSON preset with the variation under test, save it, load it with *Load preset…*, click *Validate & Use*, then train. Repeat with the next preset. *Apply to experiment* stores a parameter set in the experiment JSON as `ml_defaults` so that you can load it again with *Pull from experiment*. The trained model's `training_run.json` records the parameters that were used, and it also holds a `cga_tolerances` block (§4.5).

## 6.2 Precursor (MS1) gate

The **Precursor (MS1) Gate for Predictions** panel on the Predict tab (Quickstart Step 2.8) compares each spectrum's measured precursor mass against the theoretical mass of the predicted composition. A glycan prediction whose measured-vs-theoretical PPM error exceeds the configured tolerance fails the gate. It keeps `pred_ok` = False and a `ppm_precursor` value in the full prediction table, and it is left out of the report files (§2.9).

**When the gate helps.** The gate is most useful for CGA-derived models where the labels were assigned at a known precursor PPM tolerance. If a prediction's composition doesn't agree with the precursor mass within a similar tolerance, something is off — either the model is wrong, or the spectrum's precursor was assigned to the wrong scan, or the in-silico CSV is incomplete. The gate surfaces these cases in `*_gate_debug.csv` rather than mixing them into the main predictions.

**When to disable.** Disable the gate when working with high-mass-error datasets (some instrument configurations don't reach 10 ppm reliably), when the in-silico CSV is known incomplete, or when running an exploratory prediction where you want the model's raw output regardless of mass agreement. Untick the *Apply precursor gate* checkbox; the application then sets `pred_ok = True` for all rows and skips the gate step.

**Not the CGA → Trainable gate (GlycoMSP 1.14).** This prediction gate and the CGA → Trainable precursor gate (§2.6) are separate settings with the same 10 ppm default. If you widened the CGA tolerances and the CGA → Trainable gate, a 10 ppm prediction gate demotes nearly every prediction. Loading a model whose training record has a non-default CGA → Trainable gate (or the gate off) therefore prefills this panel with that setting (§2.8). You can still change it; a notice and `run_summary.json` record the difference.

**Required input when enabled.** The gate needs an in-silico CSV mapping composition → theoretical_mass, set via the *Select in-silico CSV file* button. The same in-silico CSV used in CGA Setup (Step 2.4) works directly; if you trained on multiple CGA samples with different libraries, pick the library that covers the broadest composition range. The library needs a column named `Mass` or `theoretical_mass`, and either a composition column (`composition`, `structure`, `label` or `comp_str`) or the monomer count columns.

**Gate vs confidence threshold.** The two checks are independent and they act differently.
- Confidence threshold (when enabled, §6.3): a prediction of a glycan class is relabelled as the majority label (*Non-glycan* by default) when its top probability is below τ and the majority-class probability is within the majority support margin of the top probability.
- Precursor gate (when enabled): a glycan prediction passes when |measured_mass − theoretical_mass| / theoretical_mass × 1e6 ≤ `gate_ppm`. If the library lists several masses for a composition, the closest one is used.

Predictions that fail the gate keep `pred_ok` = False and a `ppm_precursor` value in the full table `{input stem}_predicted…csv`, and they are left out of the report files. The report's own review files are `to_review_low_conf.csv` (top probability below τ), `to_review_low_margin.csv` (top probability minus second probability below the margin) and `to_review_ambiguous_topk.csv` (accepted predictions whose runner-up is within the margin). The output folder also contains `*_gate_debug.csv` (the glycan predictions that failed the gate), `*_gate_passed.csv`, `*_PG{n}ppm_only.csv`, `*_PG_counts.csv` and `*_gate_nolut_bycounts_debug.csv`. The gate files show the per-row gate decision and are useful when investigating cases that surprised you.

**Counts-based fallback for composition lookup.** When the in-silico CSV's composition strings do not match the predicted labels (for example because the tokens are in a different order), the gate's string lookup fails. Both sides are first converted to one canonical form, and the legacy token `K` is read as `KDN`. If the string lookup still fails, the application matches the predicted label by monomer counts against the in-silico CSV's monomer columns (`Hex`, `HexNAc`, `Fuc`, `NeuAc`, `NeuGc`, `KDN`, `HexA`, `SO3`, `PO3H`). When the CSV has no composition column, the strings are built from these columns. The file `*_gate_nolut_bycounts_debug.csv` lists up to 30 predictions that had no string match, with the candidate masses found by counts. A prediction that matches by neither path fails the gate. The log panel reports the number of string matches, counts matches and misses.

With the reference-database fill ticked in Prepare Dataset, the full prediction table also carries the columns `Predicted_GlyToucan_ID` and `Predicted_WURCS` (added in GlycoMSP 1.11, §3.7).

## 6.3 Threshold sweeps and τ selection

The GUI does not produce a τ-sweep summary. To choose τ for your data, train with different values of *Threshold τ (0.00–0.99)* in the Train/Test Parameters dialog and compare the reports. Each training run writes `{trainable CSV name}_rf_performance.txt` next to the trainable CSV. It holds the per-class precision, recall, F1-score and support on the held-out test fold with the accuracy, the macro average and the weighted average, followed by a glycan-only report that leaves out the majority label. When the confidence threshold is enabled, the report is computed after the threshold has relabelled the test predictions. A new run on the same trainable CSV overwrites the file, so copy it before you train again.

The general pattern: as τ rises, more low-confidence glycan predictions are relabelled as the majority label, so fewer glycan predictions remain (coverage drops) but those that remain have higher confidence (precision rises). Macro-F1 typically peaks at a moderate τ — too low and confused predictions degrade per-class recall. Too high and rare classes get relabelled as the majority label entirely.

**Practical rules of thumb.**

- **Small datasets (< 100 spectra/class):** disable threshold demotion entirely. Set the *Enable confidence threshold → fallback to Majority label* checkbox off. The model's probability calibration is too coarse for thresholds to be meaningful at this scale.
- **Mid-size datasets (100–1000 spectra/class):** enable the threshold and start near τ = 0.65 (the GUI default value), then adjust based on the macro average F1-score in the `_rf_performance.txt` reports of your trial runs — use whichever τ maximises it.
- **Large datasets (> 1000 spectra/class):** thresholds become meaningful, and aggressive values (τ ≥ 0.70) start to make sense for production prediction runs where you'd rather relabel uncertain glycan predictions as Non-glycan than keep them in the final output.

**Margin selection.** The *Majority support margin δ (0.00–0.20)* has a default of 0.05. It sets how close the majority-class probability must be to the top probability before a low-confidence glycan prediction is relabelled as the majority label. Raise δ if you want more low-confidence glycan predictions relabelled as Non-glycan, up to the maximum of 0.20. Lower δ if you want fewer. The margin compares the majority-class probability with the top probability only. It does not look at ties between two glycan classes. The prediction report also lists a `margin` column, which is the top probability minus the second probability.

The τ and margin you settle on can be saved into the experiment JSON via the *Apply to experiment* button and loaded later with *Pull from experiment*.

**What the summaries show.** Whatever the threshold setting, the run summary reports *Accepted 0* and probability statistics of 0.000 (known issue 29). With the threshold off, which is the default, it also reports τ 0. The predicted labels and the per-class counts are correct. The code contains four different default values for τ (known issue 34). The rule described here, set in the Train/Test Parameters dialog, is the one the GUI applies at training and at prediction.

---

# Part 7 — Troubleshooting & known limitations

This Part covers the most common issues users encounter, plus the limitations baked into the current version that are good to know upfront so they don't surprise you mid-analysis. Entries marked "known issue NN" refer to the numbered list in `docs/known-issues.md`.

## 7.1 Common errors and their fixes

### Conversion

**RAW conversion fails with a "Raw Conversion Not Available" dialog.** GlycoMSP needs the Thermo RawFileReader library plus the `pymsfilereader` Python package to read RAW files directly. This works only on Windows. If you're on macOS or Linux, convert your RAW files to mzML externally (e.g., MSConvert) and use **Add mzML File** instead.

**The "Convert Raw to CSV" button shows a "pymzML Module Missing" dialog when only mzML files are selected.** The mzML reader needs the `pymzml` library. The dialog opens with "Could not import mzmlreader." followed by the import error. Install it with `pip install pymzml` (it is listed in `src/requirements.txt`) and start the conversion again. A RAW conversion does not need `pymzml`.

**The "Convert Raw to CSV" button shows an "Ambiguous Input" dialog.** RAW and mzML files are both selected. Convert one format at a time: clear one selection and retry.

**Conversion fails with "could not write to output directory".** The conversion writes a temporary CSV before promoting it to the final filename. The first file of a batch is converted into the folder of the input file, and the finished CSVs are moved to the folder you choose when you save. Later files of the batch are converted into the folder you chose. Check that the folders exist, that your user account has write permission, and that any antivirus / sync client (OneDrive, Dropbox) isn't holding the destination file open. Watch the application's log panel for the specific path that failed. If only the move fails, the log shows a `[WARNING] MS2 finalize failed` or `[WARNING] MS3 finalize failed` line (see the entry on "All metadata have been completed" below).

**Conversion ends "failed".** The reader could not produce a usable result, for example because the file could not be opened, a decoding stage stopped, or no MS2 spectrum with peaks was exported. A failed conversion is never saved as if it had succeeded. You see:

- A dialog titled *Conversion failed*. It names the file, gives the exit code (5), the reason, and the sentence "Nothing is promoted to the output folder; the reader's files stay where it wrote them." The last line gives the path of the diagnostics report. If the error is not one the reader reported itself, the dialog is titled *Thermo Library Error* instead and reads "Raw file could not be processed." followed by the details.
- In the log panel, `[reader] outcome: failed (exit code 5)` and a `[reader] reason:` line, followed by the path of the diagnostics report.
- In the editor, the title changes to *Error:* and the file name, and the status line turns red. The editor then closes by itself and the main window status reads *Idle (error)*. The rest of the batch is not converted. Leftover temporary files are removed (`[Cleanup]` lines in the log).

The diagnostics report (`<name>.reader_diag.txt` and `.json`, §5.1) stays next to the input file, or in the output folder if you had already chosen one. It states the stage that stopped and the exception. Click **Save Log** on the main window and read the log for the underlying error message. Most often the cause is a corrupt RAW file, malformed metadata, or a temp-path collision.

**Conversion ends "partial".** Some spectra could not be decoded and were not exported. The CSVs hold all the other spectra and are saved with a warning. You see:

- A warning dialog titled *Partial conversion*. It names the file, gives the exit code (4), says that some spectra could not be decoded and were not exported, lists the warnings (the numbers of MS2 and MS3 spectra), and states that the CSVs hold the rest and will be saved with this warning.
- In the log panel, `[reader] outcome: partial (exit code 4)` and the warning lines. When you save, the log adds a `[finalize] partial conversion promoted with its warning` line that names the diagnostics report.
- In the editor, the status line is orange: "Conversion partial (see log). Ready for metadata."

Open the diagnostics report and the status file (`<name>.scan_status.tsv`, §5.1). Spectra that were not decoded carry `arrays_status decode_error` there. Decide whether the loss matters for your analysis. If it does, re-convert the file.

**Conversion ends "complete with warnings".** The conversion finished and saved every spectrum it could export, but the reader has something to report. The usual cause is that MS2 spectra come without a precursor charge, which is common for ion-trap data. These spectra are kept with their peaks and precursor m/z, charge `0` and mass `0`, and are marked as unknown charge in the status file. You see no dialog. You see:

- In the log panel, `[reader] outcome: complete with warnings (exit code 3)`, then `[reader] warning:` lines. The charge warning reads "… exported MS2 spectra have no charge (written as chargeState 0, protonatedmass 0); charge-dependent steps cannot use them."
- In the editor, the status line is orange: "Conversion complete with warnings (see log). Ready for metadata."

Look at the warnings in the log, then at the diagnostics report and the status file (§5.1). The reader also gives exit code 3 when exported MS2 spectra are profile data (see the next entry). If many spectra have no charge, you can ask the conversion to supply one. In the conversion options dialog, tick *Calculate the charge from the MS1 isotope pattern (charge 1–3)* or enter a *Dataset charge for the rest*. Both options are off by default. A charge that the file reports is never changed. §3.8 describes when to use them and §4.9 describes what they do and do not do.

**Converted mzML rows hold hundreds of points per scan, or the diagnostics warn that exported MS2 spectra are "profile data" (known issue 19).** The mzML file holds profile MS2 spectra, and GlycoMSP does not centroid them. Each row then holds every profile point, zeros included, not peaks. This happens when MSConvert runs `zeroSamples` before `peakPicking vendor`. MSConvert then skips peak picking altogether. Re-convert the mzML file with peak picking on MS2 as the first filter: `peakPicking vendor msLevel=2-`. With that, the mzML route gives the same rows, peaks and charges as the RAW route.

**RAW files from an ion-trap instrument gave empty peak lists in older versions.** Before GlycoMSP 1.13, the RAW converter read only the label data of each MS2 scan, which exist only for Fourier-transform scans. Ion-trap scans therefore gave rows with empty peak lists. Since 1.13 the converter reads the vendor mass list when no label data exist. Re-convert such files with 1.13 or later. MS3 rows also carried the parent's peaks in older versions and now carry their own.

**CGA annotates almost nothing on ion-trap data.** The CGA tolerances at their defaults (20 ppm for the precursor, 10 ppm for Score A fragments) are too narrow when the scan CSV has no monoisotopic precursor m/z and the precursor mass is the isolation m/z. Widen the tolerances in the *Mass tolerances (CGA matching)* frame of CGA Setup and check the choice with the decoy check. §3.8 gives a worked example, §4.10 explains the tolerances, and §2.4 (d) describes the frame.

### Reference database

**Prebuild fails with `CERTIFICATE_VERIFY_FAILED`.** The prebuild step contacts the network and verifies the server certificate. Some macOS Python builds have no system certificate store. Install the `certifi` package (`pip install certifi`) and run the prebuild again. `certifi` is listed in `src/requirements.txt`. GlycoMSP never switches certificate verification off.

**Selecting a file in the Manage References window shows "Not a reference database".** The file you chose could not be opened as a GlycoMSP reference database. The dialog gives the path and the error. The status line of the window adds that the file was refused and that the previous database was kept. Choose the correct `.db` file with *Select database…* (§3.7).

**The *Active DB:* label in Prepare Dataset reads "unreadable (…)".** The path in brackets exists but could not be opened as a reference database. Open **Manage References…** and select a valid database file, or run *Prebuild / Update from source…* to create one. The label reads "not found (…)" when the file does not exist at that path.

### Metadata and batches

**A batch ends with "All metadata have been completed" although one file failed to finalize (known issue 23).** When the move of a file's CSVs to the output folder fails, the application still advances the batch. After the last file, the dialog titled *Metadata* reads "All metadata have been completed." and does not mention the failure. Look in the log for a `[WARNING] MS2 finalize failed` or `[WARNING] MS3 finalize failed` line and check the output folder for the `ms2_…` and `ms3_…` CSVs of each file. Click **Save Log** before you close the application.

### Manual annotation support

**Right-click "Link and Validate Sample" stops with a "Missing File" dialog.** The Sample needs both a CSV and an Excel attached before validation can run. Check that both files are attached (visible as `csv` and `excel` slots on the sample node). Metadata is validated separately at stage 6 of the link-and-validate workflow. A missing metadata file prompts you ("No metadata found. Would you like to create it now?") rather than blocking.

**MAS validation fails on the Excel file with "Excel Validation Failed".** The dialog names the file. If the workbook cannot be read or has no `MSlist` sheet, the dialog also shows the error. If the `MSlist` structure is invalid, the dialog says "See log for details", and the log panel shows a `[Validation Failed] MAS Annotation failed` line with the path. The specific reason (a header mismatch, a missing column, an empty or non-string column) is printed to the console window from which you started GlycoMSP, not to the log panel. The most common reasons: missing `MSlist` sheet, an `MSlist` header that matches none of the accepted header sets, missing or non-numeric mass column in the `ionlist` sheet.

**Drag-drop in the Dataset Explorer doesn't work on rows above the sample level.** Drag-drop only reassigns Sample nodes between Experiment nodes. Dragging an Experiment node, a file slot, or a method reference is silently ignored. To rename or reorganise at the Experiment level, edit the experiment JSON directly or use the right-click menu on Sample nodes.

### Machine learning and prediction

**ML training fails with an "Invalid Label" message.** The Label Column setting in the Train tab holds a value that is not a supported label column. Set the *Step 2: Select Label Column* combobox to one of `Structure` (the default), `IUPACname(optional)`, `Glycanannotation2` or `GlyToucan ID`, and make sure that column exists in your CSV and holds labels.

**Predictions look random — every spectrum gets the same class.** Most likely cause: feature columns at predict time don't align with the model's training column order. Check that you're using the same ion list at predict time as at train time. The application reorders columns to match `train_feats` from the model's `_features.json`. If your unlabeled CSV is missing too many columns, the log reports "N features missing from input, filled with 0". The model then receives 0 for those features, which is not the value of 1.0 that a scan with no matching peak receives, and it can fall back to predicting the majority class.

**Predict tab gate writes lots of `*_gate_debug.csv` rows but few `*_gate_passed.csv` rows.** The model is predicting compositions whose theoretical masses don't match the measured precursor masses. Check: (a) is your in-silico CSV the right one for this dataset? (b) is the *Tolerance (ppm)* of the precursor gate reasonable for your instrument (the default is 10)? (c) is the dataset actually a different glycan class from what the model was trained on? Disable the gate temporarily to see what the model would predict without filtering.

**Feature building stops with a math-domain error (known issue 33).** A converted CSV holds a peak with zero or negative intensity that matches an ion of your list, and the run stops with `ValueError: math domain error`. Whether real data contain such peaks has not been checked. Remove the peaks with intensity 0 or below from the converted CSV before training or prediction. A fix is planned for 1.15.

### Command-line prediction tool

**The command-line prediction tool cannot load a model (known issue 30).** `msp_CLI_predict_skops.py` stops with `TypeError: trusted must be a list of strings`. The tool loads the model with `trusted=True`, which skops 0.10 and later reject, and `src/requirements.txt` pins a newer skops. Two workarounds exist. Predict from the GUI (Predict tab), or install an older library for the command line: `pip install "skops<0.10"`. The GUI is unaffected and was used for every result reported in the manuscript. The compatibility converter `msp_CLI_joblibtoskops_compat.py` fails in the same way. A fix is planned for 1.15.

**The command-line prediction tool drops the precursor-mass feature by default (known issue 31).** The default `--drop-cols` list (`MS2scan_no protonatedmass`) removes the precursor-mass column before prediction. A model trained with the mass feature then finds the column missing and receives 0 for it, which is also the mass value of Non-glycan rows. The predictions are made on a changed input and no message is shown. Predict in the GUI, or pass `--drop-cols` with an explicit list that does not name the mass column the model was trained with. A fix is planned for 1.15.

**The command-line prediction tool writes numbers instead of glycan compositions (known issue 32).** The GUI trains the Random Forest on integer class codes and saves the label encoder next to the model. The command-line tool does not read that file, so its output column holds the integer codes. Predict in the GUI, or map the integers back with the saved `_labelencoder.joblib` file. A fix is planned for 1.15.

### Platform

**On Windows, the application's console flashes briefly after an error then closes before you can read it.** If an unhandled error occurs while GlycoMSP starts or imports its modules, the application prints the traceback and a short environment check, saves them to a file named `startup_error_<date>_<time>.log` in `%LOCALAPPDATA%\GlycoMSP` (in `~/.glycomsp` on macOS and Linux), and exits. On Windows it also waits for you to press Enter before it closes, but only when its standard input is not an interactive terminal. To keep the console open in every case, launch GlycoMSP from a persistent terminal (Command Prompt or PowerShell). The `startup_error` file holds the same text. To check your environment without opening the GUI, run `python src/mspfileloaderv14.py --diag`. For errors during a session, use **Save Log** before the application closes.

**Lower CGA Setup flags are cut off.** On a small or remote display, the CGA Setup window can be shorter than its content. The file links, the *Mass tolerances* frame and both button rows always stay visible; the lowest rows of the *In-Silico Generation Flags* panel are cut off instead. Their values are still used. Maximise the window to see or edit them.

**macOS file dialogs occasionally crash the application during multi-file selection.** A platform issue with tk on certain macOS versions. The application catches the underlying exception and degrades to the **Add MS2 CSV(s)** / **Add Excels** path. If repeated, prefer single-file selections.

## 7.2 Known limitations

Entries follow the workflow: conversion, library and masses, CGA and tolerances, training and prediction, reference database, command-line tools, general. With every option added since v1.10 left at its default, results equal those of v1.10, except for the corrections listed under *Fixed* in the changelog for 1.13.

### Conversion

**GlycoMSP does not centroid.** If an mzML file holds profile MS2 spectra, each row holds every profile point. GlycoMSP warns (outcome "complete with warnings", exit code 3) and does not change the data. Re-convert the file with peak picking on MS2 (§7.1, Known issue 19).

**Rows with unknown charge are kept with mass 0 and are not excluded downstream (known issue 17).** When the instrument does not assign a charge, the spectrum is exported with charge `0` and mass `0`. The CGA matching, the CGA → Trainable precursor gate, the mass feature, training and prediction do not exclude these rows. CGA keeps such a scan as unmatched without giving a reason, and a model trained with the mass feature sees mass 0 as a Non-glycan signature. On the ion-trap native files measured, 20 to 29 % of the exported MS2 rows have unknown charge. This loss is bounded and is not the cause of the near-zero CGA result on those files (see the tolerance entry below). This is left unchanged on purpose. The *Calculate the charge from the MS1 isotope pattern* and *Dataset charge* options reduce the number of affected rows (§3.8).

**The isotope comb has limits (known issue 20).** The isotope comb tests charges 1 to 3 on profile MS1 scans that are not Fourier-transform scans. It does not call charge 4 or higher. A mixture of charge states is reported as the stronger one without a flag (3 of 8 hand-checked charge 3 calls were mixtures). On the ion-trap files tried, 24 to 27 % of MS2 scans stayed without a call (status `below_margin`). The comb is a heuristic that was checked for feasibility, not a validation. Scoring checks its calls downstream.

**Array integrity is not checked by the steps that read a converted CSV (known issue 21).** GlycoMSP's own converters did not write a malformed row in any file checked. A scan CSV that you edited by hand, or that another tool wrote, can hold a row whose m/z list and intensity list differ in length. The status file marks such a row `arrays_status length_mismatch`, and a step that scores it can stop with an `IndexError`. Check such a file before you use it.

**Two simultaneous conversions of files with the same name into one folder can mix their diagnostics (known issue 24).** This needs two GUI windows, or the GUI and the command-line converter, converting files with the same name into the same folder at the same moment. Only the explanation shown and the diagnostics files moved can be affected. Whether the CSVs are saved follows the reader's own result. Convert such files one after the other.

### Library and masses

**N-glycans and O-GalNAc only.** The in-silico library generators (`compnewv4.NGlaunch` and `compnewv4.OGlaunch`) cover N-glycans and O-GalNAc glycans. Other glycan types (O-mannose, glycolipids, GAGs) are not currently supported. Adding new glycan classes requires extending the composition generator's flag schema and partition sets.

**Native masses apply to the CGA route only (added in GlycoMSP 1.13, Known issue 22).** With derivatization *Native* or *Native(Reduced)*, the CGA library gives underivatized precursor masses (added in GlycoMSP 1.13). Manual annotation (MAS) still sums permethylated residue masses, so a native project annotated through MAS gets permethylated theoretical masses. Use the CGA route for native glycans.

**Permethylated masses are not modelled for mixed sulfate and phosphate, or for three or more of either (known issue 26).** The sulfate and phosphate corrections were fitted to observed species, one group type at a time. Compositions with both sulfate and phosphate, or with three or more sulfate or phosphate groups, get a wrong theoretical mass. These compositions were not observed in the analyses. The command-line mass tool refuses them (exit code 2), and the CGA library still computes them as in v1.10. Native masses are not affected.

### CGA and tolerances

**CGA tolerances are fixed at 20 ppm and 10 ppm unless you change them (known issue 27).** In v1.10 the precursor match (20 ppm) and the Score A fragment match (10 ppm) could not be changed, and ion-trap data whose precursor mass is the isolation m/z received almost no annotations. Since GlycoMSP 1.14 both are fields in CGA Setup, in ppm or Da. They default to the v1.10 values. §3.8 gives the worked example for ion-trap data and §4.10 explains the concept. Score B enrichment keeps its fixed 20 ppm.

**Score A counts peak–ion pairs (known issue 28).** With a wide fragment window (0.3 Da or more) and a dense ion list, one peak can be counted for two neighbouring ions. With the 14-ion native list no scan is affected. With a 97-ion list about one scan in eight is. At the default 10 ppm window the effect is low. The v1.10 counting rule is kept on purpose so that earlier results stay reproducible. Compare Score A only between runs with the same fragment window and ion list (§4.1, *How hits are counted*).

**The ppm window is slightly asymmetric (known issue 10).** The ppm error is divided by the larger of the peak m/z and the ion m/z, so the window above an ion is not exactly the size of the window below it. The difference is very small and applies uniformly to every comparison. The window of a Da tolerance is symmetric.

**The feature window is ppm only.** The *ion list ppm* field in *CGA → Trainable* and the *PPM Tolerance* dialog of prediction take ppm, not Da. The prediction dialog accepts at most 500 ppm (0.13 Da at m/z 260). The *CGA → Trainable* field has no limit, so keep it at 500 ppm or less so that prediction can use the same window. Low-m/z ion-trap fragments therefore get a narrower window at feature extraction than the CGA step's Da setting allows. A ppm-or-Da selector for feature extraction is planned.

**The CGA result TSV carries every converted scan.** Its size follows the raw file. A run of 100,000 MS2 scans gives about 1.7 GB. A matched-rows-only output is planned.

**Score B is independent of MAS.** The motif-aware re-ranking step operates on candidate compositions, of which MAS rows have only one (the expert label). Score B is therefore CGA-only. There is no analogue for refining expert labels with motif evidence.

**Saved CGA settings do not include the mass tolerances.** The *Save CGA settings* button writes the flag panel state to JSON, including the O-GalNAc core types. The *Mass tolerances (CGA matching)* fields are not part of that file. Enter them again after you load a saved preset. A run with non-default tolerances is recorded in the CGA method JSON (§5.2).

### Training and prediction

**Single classifier in v1.** The Train tab's *Choose Classifier* combobox lists Random Forest as the only option. The infrastructure (label encoding, balancing, splitting, threshold demotion) is classifier-agnostic. Adding alternative classifiers (gradient-boosted trees, support vector machines, neural networks) is feasible but not in this build.

**Ion list must align across mixed-source training.** Combining CGA-derived and MAS-derived rows requires matching feature spaces. *Combine Datasets for Training* joins the files by column name and leaves empty cells where a file lacks a column. Training fills empty cells with 0, which is not the value of 1.0 that a scan with no matching peak receives. Mismatched ion lists therefore work but degrade signal. Use the same ion list across samples you intend to combine.

**Negative sampling: the shared dialog does not apply to CGA → Trainable.** The Negative Sampling Setting dialog (§3.4) is shared by the MAS merge and by the negatives added during CGA ion mining. *CGA → Trainable* has its own *5) Negative sampling* box (§2.6) and ignores the dialog. Configure the dialog before each MAS merge or CGA analysis with ion mining, and set the box of *CGA → Trainable* separately.

**The MAS feature writer keeps the last in-window peak, not the most intense (known issue 15).** When two or more peaks fall inside the ppm window of one ion, the MAS-ML feature is taken from the last peak (the one with the highest m/z). CGA-ML training and prediction take the most intense peak. The severity is very low. No such case occurred on the published zebrafish MAS-ML data.

**Prediction run summaries show *Accepted 0* and probability statistics of 0.000 whatever the threshold setting (known issue 29).** With thresholding off, which is the default, they also show τ 0. The per-row columns `max_proba`, `second_proba`, `margin` and `entropy` of `predictions_with_conf.csv` are 0.000 in every row for the same reason. The predicted labels and the per-class counts are correct. The summary will be corrected.

**Four different default confidence thresholds exist in the code (known issue 34).** The GUI uses τ = 0.65 and margin = 0.05 as its defaults, and its margin-aware rule is the one applied. The Predict tab falls back to 0.60 if the τ entry cannot be read, the prediction report object defaults to 0.60, and a helper function used only for the τ sweep defaults to 0.55. A script that uses one of those library defaults gets a different cut-off than the GUI. Set τ and margin explicitly in the GUI and record them (they are written to the run summary). One documented default is planned for 1.15.

**Feature building stops on a zero or negative peak intensity (known issue 33).** See §7.1. A fix is planned for 1.15.

**`Create unlabeled dataset of certain experiment` is per-experiment.** The current build runs this workflow against one experiment JSON at a time. Predicting across multiple experiments simultaneously requires running it once per experiment.

**Older models (`.joblib`, `.pkl`) load but cannot be re-saved with full provenance.** The application loads `.joblib` and `.pkl` models for predict-time use, but the `training_run.json` reproducibility manifest is only generated during *training*. Models trained in older builds may lack the manifest. Treat them as legacy artifacts and re-train when reproducibility matters.

### Reference database

**Shipped snapshots carry no compatibility promise.** This applies to the showcase snapshot in the `references/` folder. Merging a shipped snapshot into your own database is planned for a later preview.

**`wurcs_only` rows keep an empty accession.** A composition with a WURCS string but no registration in GlyTouCan has no accession. It stays empty. Such rows were seen for some KDN-, HexA- and sulfate-containing compositions.

**Three actions of the Manage References window and the sample menu are disabled placeholders.** They are *Check again*, *Resolve missing references* and *Update GlyToucan ID / WURCS in this file…*. Each label also shows an internal reference code in brackets. They do nothing in this version.

**Prebuild needs network access.** It is the only part of GlycoMSP that contacts the network. Filling from the reference database during a run never does. The option is off by default, and with it off all outputs equal those of v1.10 (§3.7).

### Command-line tools

**The command-line prediction tool (known issues 30–32).** It cannot load a model with the pinned skops version, drops the precursor-mass feature by default, and writes integer class labels. §7.1 gives the error text and the workarounds. Fixes are planned for 1.15.

**Command-line scripting support beyond `convert` and `inspect` is still in progress.** Command-line conversion (`convert`) and a read-only probe of the Thermo backend (`inspect`) exist (§3.9). General-purpose scripting support beyond these two commands is still in progress.

### General

**Multi-method tracking is partial.** Each sample's `_methods` array can hold multiple method-JSON references (one per family — typically one MAS and one CGA), but the active-method slot `files["json"]` holds only the most recently saved one. UI surfaces (Export Method v1 buttons, exp-JSON writebacks) act on the active slot. If your workflow requires switching between MAS and CGA methods on the same sample, save each method explicitly and re-link before running the relevant pipeline.

**A loaded method of unrecognised family hides its `artifacts`.** The `artifacts` block of the method JSON stays hidden until *Export Method v1 (MAS)* is used once. This predates the reference-database work.

**Method JSON paths are stored OS-natively.** Method JSON files written on Windows have backslash separators. macOS / Linux files use forward slashes. The application's read-time path resolver handles cross-platform load (Windows-drive regex), but diff-ing a method JSON across hosts shows path-separator changes that aren't semantic. Experiment JSON files use POSIX paths uniformly and diff cleanly across hosts.

---

# Part 8 — Reproducing the manuscript results

This Part is for reviewers and readers who want to reproduce the results reported in the GlycoMSP manuscript. Each subsection maps to a package in **Supplementary Table S2c** (Zenodo), so you can find the exact input/output files. The reported analyses used a pre-release build of **GlycoMSP v1.10**. Sections 8.0 to 8.4 refer to GlycoMSP v1.10, which stays the release the manuscript cites. Later versions keep the v1.10 defaults, so these results are unchanged with every new option at its default. Sections 8.5 and 8.6 describe what was added for the revision.

> **Note.** Package ID numbers below follow Supplementary Table S2c. The MAS-ML LOTO and pooled cross-validation analyses are bundled into a single package (ID 7); downstream IDs shift accordingly (GlycoGenius → ID 8, CandyCrunch → ID 9).

## 8.0 Quick start for reviewers

If you only do one thing, reproduce the **MAS-ML evaluation (Result 2)** — it is a single standalone script with bundled data and fixed seeds:

```bash
python mas_ml_loto_manuscript_analysis.py \
  --brain zf_sPerMeNG_brain_merged.csv \
  --intestine zf_sPerMeNG_intestine_merged.csv \
  --ovary zf_sPerMeNG_ovary_merged.csv \
  --features training_run.json \
  --gms-src /path/to/GlycoMSP/src \
  --outdir MAS_ML_LOTO_results
```

Expected (glycan-only macro-F1): within-tissue **0.931 / 0.935 / 0.952** (brain / intestine / ovary); leave-one-tissue-out **0.635 / 0.599 / 0.832**; pooled 5-fold CV **0.904 ± 0.038**. The brain within-tissue run reproduces the GUI `rf_performance.txt` exactly (test support 219, accuracy 0.99). The other results are GUI-driven workflows described below.

## 8.1 Setup and environment

- **Source code:** GitHub `https://github.com/henry4204aaa/GlycoMSParser` (v1.10, MIT License).
- **Install:** Python 3.10+ with `pandas`, `numpy`, `scikit-learn`, `openpyxl`, `pymzml`, `joblib`/`skops` (and `pymsfilereader` on Windows for RAW). See the repository `requirements` files.
- **Reanalysis data:** the Zenodo archive (DOI in the manuscript Availability section), organised into the packages listed in Supplementary Table S2c.
- **Reference environment:** the reported metrics were generated with scikit-learn 1.2.0; outputs were verified byte-identical on scikit-learn 1.7.2, so a current install reproduces the same numbers.
- **Paths caveat:** GlycoMSP method/experiment JSON files may store absolute paths. When reproducing on another machine, update the saved paths to your local working directory before loading them (see the Quick Start Guide bundled in the Zenodo Misc package).

## 8.2 MAS-ML evaluation — Result 2 (Supp S2c ID 7)

The MAS-ML held-out, leave-one-tissue-out (LOTO), and pooled cross-validation results are reproduced by the standalone script in `manuscript_reproduction/` (also archived as Zenodo package ID 7 — the MAS-ML standalone, bundling held-out, LOTO, and pooled cross-validation). Run the command in §8.0 against the three zebrafish MAS trainable CSVs and the `training_run.json` from the package. The script reads GlycoMSP-exported trainable CSVs and reproduces the GUI pipeline without modifying core code; full design, metric definitions, and the GUI-equivalence evidence are in the package's `MAS_ML_LOTO_technote.md` and `MAS_ML_LOTO_validation_report.md`. Primary metrics are glycan-only (the operational Non-glycan class and out-of-label-space compositions are reported separately).

## 8.3 CGA / CGA-ML — Results 3–4 (Supp S2c ID 4c, 1c)

CGA/CGA-ML results are reproduced through the GUI using the provided in-silico glycan library, fragment-ion list, and CGA settings:

1. Convert the provided mzML/RAW to an MS2-indexed CSV (or use the bundled converted CSV).
2. In **CGA manager**, link the bundled in-silico CSV and ion list, then **Start CGA analysis**.
3. **CGA → Trainable** to build the trainable CSV, then **Train** and **Predict** as in Part 2.

For **U-937 O-GalNAc (Result 3, package 4c)** this recovers the nine MALDI-profile compositions plus the four CGA-ML-supported minor compositions (Fig. 3). For the **zebrafish brain F2H4N5 rescue (Result 4, package 1c)** the CGA-ML output recovers the omitted composition traceable to its MS2 scans (Fig. 4). CGA-only and negative-mode contexts correspond to Supplementary Figs. 2–3.

## 8.4 External-tool comparison — Result 5 (Supp S2c ID 8–9, Supp Table S3)

The GlycoGenius and CandyCrunch comparison (Fig. 5) is reproduced from the archived tool settings and reports: **GlycoGenius** (ID 8 — `.ini`/`.gg` project files + results; Supp S3a settings) and **CandyCrunch** (ID 9 — exact run commands + reports; Supp S3b). Reported values: MAS-ML recovered 19 of 20 eligible compositions; GlycoGenius re-identified all 38 imported MAS compositions; CandyCrunch overlapped 8 and produced 11 additional; in the larger search-space analysis GlycoMSP CGA assigned 573 (88 retained by CGA-ML), GlycoGenius reported 1042 (230 good-quality), and CandyCrunch reported 19.

> The archived packages also include runtime/error observations and the minimal workarounds we applied, recorded as *issues encountered in our tested environment* (with OS and version details). These are provided for reproducibility and transparency, not as a judgement of the other tools.

## 8.5 Automated tests and validation tools

GlycoMSP 1.14.1 is a test release. It adds an automated test suite and a folder of validation tools to the public repository. It makes no change to the application, so `src/` is identical to 1.14. This section tells you how to run the tests, what they cover, which tests are withheld and how to regenerate what they need, and how to use the validation tools.

### Run the test suite

The suite runs from a clean checkout of the repository. It needs no data files. Run these commands from the repository root:

```bash
pip install -r src/requirements.txt      # Python 3.10-3.12 (3.13: src/requirementspy313.txt)
pip install pytest
python -m pytest tests
```

The repository READMEs print the last command as `python3 -m pytest tests -q -rs`. The option `-q` shortens the output and `-rs` prints the reason for every skipped test. Add them if you want those reasons on screen.

**Result to expect.** From a clean checkout the run gives 626 tests with 2 skipped, 1 expected failure and 0 failed. A test that cannot run on your machine skips with a reason and does not fail. The table lists the skip reasons.

| Skip | Reason |
|---|---|
| A few graphical-interface tests that need a real Tk window | The machine has no display (a headless server). The same logic is tested without Tk elsewhere. |
| One test that rebuilds an excerpt of a real data file | The excerpt is optional, local to the development machine, and not distributed. |
| One test that compares the MAS writer with the v1.10 development commit | That commit exists only in the development repository. |
| One test of process termination in the Manage References window | It runs on Windows only. |

**The expected failure.** One test is a strict expected failure (`xfail`). It records known issue 30: the command-line prediction tool cannot load a model with the pinned skops version. The test is kept on purpose. It fails loudly once the defect is fixed, so that its marker is removed deliberately and the fix is not missed. Two further tests in the same file are labelled CHARACTERISATION. They record what the command-line prediction tool does today and are not a claim that this behaviour is wanted. Known issues 30 to 34 were found by these tests. Fixes for them are planned for 1.15.

### What the tests build and what they cover

Every test builds its own input. The inputs are synthetic and seeded: a synthetic mzML file, a fake Thermo reader, and small tables. Each test file finishes in seconds. Nothing is downloaded, and files are written to temporary folders only. The suite was developed and run on Python 3.12 with pytest 9 on macOS. Runs on Windows, Linux or a continuous-integration service have not been recorded.

| Workflow stage | What is tested |
|---|---|
| RAW and mzML conversion | Charge and polarity handling, the diagnostics files, exit codes, and the conversion options dialog and reader result in the GUI |
| Composition library | Enumeration of N-glycan constraints and O-glycan cores |
| Machine learning | Features, Random Forest training, the `.skops` round trip, command-line prediction, the confidence rule and the prediction report |
| Precursor mass calculation | Native, PerMe, and sulfate and phosphate (SO3/PO3H) compositions, and the KDN residue mass |
| CGA | Score A, tolerances in ppm or Da, and the recording of the settings |
| Score B | The composition tokenizer |
| Trainable CSV command line | Import and code-hygiene guards |
| Reference database | The GlyTouCan / WURCS reference database and the Manage References window, with a stub command line in place of the network |
| Validation scripts | The Random Forest timing script, the decoy count and the shared-peak count (31 tests for the decoy script) |

Three test files carry the work added for the revision. `tests/test_ml_features_train_predict.py` runs the machine-learning path end to end on known inputs. `tests/test_compnewv4_enumeration.py` checks the in-silico library enumeration: N- and O-glycan composition sets, ranges, switches and modes. `tests/test_cga_tolerance_decoy_count.py` holds the 31 tests of the decoy-count and Score A peak-counting scripts. For the machine-learning tests, nine deliberate breakages of the product code were made and all nine are detected.

Several test files read the source text of the main GUI file and look for fixed anchors, because the GUI cannot be imported without starting Tk. If a code change moves an anchor, the test skips with the message "coverage LOST" and does not pass silently.

### The test register

`docs/TEST_REGISTER.md` is the register of the tests. It lists every test file with its purpose, the modules of `src/` it exercises, and whether its last run was confirmed by hand. The column headed *Confirmed* is filled only after a manual run on a clean checkout. An entry gives the date and the setup of that run. An empty cell means the file has not been confirmed by a manual run. The script-style checks described below are in this state because they were not part of that run. The register also records the earlier counts of each development session. Quote the totals of the previous subsection and not the totals of those sessions.

### Script-style checks

Five older checks run as scripts. pytest does not collect them. They need no data files and each finishes in under a second:

```bash
python tests/test_b01_kdn_smoke.py
python tests/test_b03_apply_to_exp_smoke.py
python tests/test_b04_excepthook_smoke.py
python tests/test_b06_status_label_smoke.py      # needs a display
python tests/test_b17_glytoucan_columns_smoke.py
```

Some of them mirror small GUI helpers instead of importing them, so they can drift from the product. Treat them as regression guards for the bug each one fixed. They are not proof of current behaviour.

### Withheld tests and how to regenerate their reference outputs

Eight test files exist in the development repository and are not part of the public set. They are listed in `tests/public_exclude.txt`.

| Withheld file | What it does |
|---|---|
| `test_smoke_02_raw_to_csv.py` | Compares the stored RAW-to-CSV conversion output of the manuscript run with a reference (stored outputs only) |
| `test_smoke_06_ml_training.py` | Compares the stored outputs of the GUI training runs (stored outputs only) |
| `test_smoke_08_cga_insilico.py` | Compares the stored in-silico libraries (stored outputs only) |
| `test_smoke_09_cga_composition_fitting.py` | Compares the stored CGA result tables (stored outputs only) |
| `test_smoke_02a_mzml_to_csv.py` | Re-runs the mzML converter on the reference inputs |
| `test_smoke_03_metadata_json.py` | Re-runs the metadata step on the reference inputs |
| `test_smoke_04_mas_to_trainable.py` | Re-runs the validator and the MAS merge on the reference inputs |
| `test_json_upgrade.py` | Tests the JSON format upgrade for versions before v1.10. It does not apply to the released code. |

The seven smoke checks need the reference outputs of the manuscript runs, about 8 GB, which are not distributed. Four of them (02, 06, 08 and 09) compare stored reference outputs and do not execute GlycoMSP code. The other three (02a, 03 and 04) re-run the converter or the validator on the reference inputs.

The reference outputs were derived from the public raw-data deposit named in the manuscript's data-availability statement. The smoke scripts compare a regenerated output with the stored one. They do not create it. To regenerate:

1. Convert the deposited RAW files with the current converter. RAW conversion needs Windows and Thermo MSFileReader. An mzML file converts on any system with `--mzml` in place of `--raw`.

   ```bash
   python src/msp_extract_cli.py convert --raw <file.raw> --out <folder> --alias <name>
   ```

2. Run the later stages (annotation workbook, trainable CSV, Random Forest, CGA) in the GUI on those files. `docs/TEST_REGISTER.md`, section "Regenerating the withheld reference outputs", gives the recipe and the expected values for each script. It also marks the settings that are not recorded in the repository yet, and it states that the recipes have not been executed end to end.

A conversion with the current converters can differ from the stored reference, which was made at the v1.10 freeze. MS3 peaks and charge handling changed, so a count that differs in the conversion check needs interpretation before you treat it as a failure.

One reference cannot be regenerated. The O-glycan run used by the CGA smoke check comes from an unpublished dataset and is not in the deposit.

### Validation tools

`tools/validation/` holds scripts that you or an assistant run on inputs you choose. They are validation tools and not application code. Nothing in `src/` imports them and they change no GlycoMSP output. They print counts or timings and answer one question each. They need Python 3.9 or later with numpy and pandas (`openpyxl` only for an `.xlsx` ion list). No script prints an input file name. `rf_timing.py` and the decoy-count script report under an `--alias` that you choose, and `B39_dup_hits.py` takes an `alias=file` argument. Where a script has tests, they are part of the normal suite. The scripts stay in this folder because they are tools and not tests.

| Script | Question it answers |
|---|---|
| `rf_timing.py` | How long does the GUI's Random Forest training step take on a given trainable CSV? It times the fit alone and the whole step from CSV load to fitted model. |
| `cga_tolerance_decoy_count.py` | When a CGA mass tolerance is widened, how many of the extra matches are chance? |
| `B39_dup_hits.py` | How many Score A rows list the same ion more than once, because several peaks fall inside one ion's window? |
| `B39_shared_peaks.py` | Which ion pairs are close enough to share one peak, and how many peaks does GlycoMSP count as a hit for two or more ions? |

`B39_dup_hits.py` is run as `python tools/validation/B39_dup_hits.py <alias>=<CGA TSV>`. `B39_shared_peaks.py` is run as `python tools/validation/B39_shared_peaks.py <ion list> [<converted MS2 CSV>]`. Both measure the counting rule described in §4.1 under *How hits are counted*.

**Decoy count.** The script matches the converted MS2 CSV against the in-silico library and against a decoy library with the same tolerance. The decoy library is the real library with every mass shifted by +3.37 Da. No decoy mass is a plausible glycan of the sample, so every decoy match is a chance match. A run with several tolerances reports every precursor and fragment pair:

```bash
python tools/validation/cga_tolerance_decoy_count.py --alias run1 \
    --converted <converted MS2 CSV> --library <in silico CSV> \
    --ions <ion list CSV/XLSX> \
    --precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da --out run1_decoy.csv
```

Per tolerance the script counts the real scans (at least one real candidate), the decoy scans (at least one decoy candidate) and the candidate pairs, capped at the CGA top 3 per scan. With `--ions` it also computes Score A against the real ion list and against a +3.37 Da shifted copy. A scan passes at Score A of 0.07 or higher (`--min-score`), or at a hit count of N or more (`--min-hits`).

**Reading the result.** Use the decoy count as the chance estimate at the same tolerance. Real minus decoy is roughly the number of matches above chance. The ratio decoy / real is roughly the fraction expected to be chance. Prefer the tolerance where real minus decoy stops growing while decoy keeps growing. The script reports counts only. It does not choose a tolerance and it does not touch any GlycoMSP output. The decision is made per dataset, as described in §2.4 (d) and §3.8.

The script was checked against GlycoMSP itself. The check is described in §8.6.

**Random Forest timing.** `rf_timing.py` builds the training matrix as the GUI's *Train* button does for a trainable CSV that carries the precursor-mass feature. Every other Train-panel setting stays at its default. It writes `rf_timing_<alias>.txt` into `--out`. The default for `--out` is the folder that holds the CSV.

```bash
python tools/validation/rf_timing.py --trainable <trainable CSV> --alias {name} \
    --repeats 3 --trees 400 --out <folder>
```

The report lists the data (rows, features, classes kept, the training-fold size before and after the majority-class cap, and the test rows), the Random Forest parameters as the product sets them, and the environment (Python, scikit-learn, numpy and pandas versions, platform, processor and logical cores). It then lists the fit time and the full-step time for each repeat with the mean and maximum, and the held-out accuracy of the last fitted model. The report carries the alias and never the input path. An alias is 1 to 64 characters of letters, digits, `_`, `.` and `-`.

## 8.6 External-dataset validation runs

The two runs in this section used data that the manuscript's own analyses did not cover. A further permethylated Orbitrap run checked the validation tools. They are described by chemistry and instrument only. They are not new manuscript results. This manual does not list their files. Their identifiers and files are listed with the reproducibility archive and not here.

### Ion-trap native negative-mode files

These files hold native (underivatized) glycans recorded in negative mode on a linear ion-trap instrument. They exposed four properties of low-resolution data. Each led to a change in GlycoMSP 1.13 or 1.14.

**Peak lists.** Before the fix, RAW conversion of ion-trap scans gave empty peak lists. The RAW reader took peaks from label data, which exist only for Fourier-transform scans, so for ion-trap scans the call succeeded and returned no points. A read-only probe showed that the mass list holds the points and that the vendor-centroided mass list is bit-identical to the vendor peak picking of MSConvert. Since GlycoMSP 1.13 the reader takes peaks from the vendor mass list when no label data exist. Two later re-runs gave peaks equal to the vendor peak picking of MSConvert (known issue 16).

**Unknown precursor charge.** Ion-trap data often carry no precursor charge. On these files 20 to 29 % of the exported MS2 rows have unknown charge. These rows keep their peaks with charge `0` and mass `0` (known issue 17). The loss is bounded and is not the cause of the near-zero CGA result described below. The isotope comb (the MS1 isotope-pattern charge estimate) can supply a charge. It tests charges 1 to 3 on profile, non-Fourier-transform MS1 scans only and never calls a charge of 4 or more. A mixture of charge states is reported as the stronger one without a flag, and 3 of 8 hand-checked calls of charge 3 were mixtures. 24 to 27 % of the MS2 scans in the files tried stayed below the comb's margin. The comb is a heuristic and not a validation (known issue 20).

**Precursor and fragment errors.** The mono-isotopic m/z of these files is 0, so the precursor mass is built from the isolation m/z. The precursor error against the library was −20 to −90 ppm at charge 1 and +150 to +250 ppm at charge 2. The fragment errors were 0.06 to 0.14 Da. The error changes sign with the charge, so a wider window is needed and an offset does not help (known issue 27).

**Effect of the tolerances.** With the fixed 20 ppm precursor and 10 ppm Score A windows of GlycoMSP 1.13, these files received almost no annotations (changelog 1.13). The 20 ppm window kept about 5 % of the real precursor matches, and Score A at 10 ppm removed nearly all of those, so almost no trainable CSV could be built. The estimate recorded in known issue 27 is that about 300 ppm and 0.3 Da give tens of annotations per file against a low decoy rate. GlycoMSP 1.14 makes both windows configurable in ppm or Da (§4.10). The decoy-count script of §8.5 is the check for a chosen window.

**The 14-ion and 97-ion lists.** For these files GlycoMSP 1.14 ships two negative-mode native N-glycan ion lists in `templates/` for use with *Native(Reduced)*: a primary 14-ion set and an extended 97-ion set. A wide fragment window of 0.3 Da or more lets one peak count for two neighbouring ions, because Score A counts peak and ion pairs. With the 14-ion list no scan is affected. With the 97-ion list about one scan in eight is (changelog 1.14, known issue 28).

### Permethylated positive-mode run with ion-trap MS2 spectra

The second run is a permethylated, positive-mode run from a public deposit, acquired on an Orbitrap instrument with ion-trap MS2 spectra, run with the Score A fragment tolerance at 0.3 Da. This manual reports no counts for it. §4.10 describes the Da fragment tolerance, and the decoy-count script of §8.5 is the check for a chosen window.

### Check of the validation tools on a permethylated positive-mode Orbitrap run

The validation tools were checked on a separate permethylated positive-mode run with Orbitrap MS2 spectra, with its in-silico library and a 282-ion list. The decoy-count script gave real counts equal to the CGA result tables of GlycoMSP on the same inputs at two tolerance pairs, 20 ppm / 10 ppm and 300 ppm / 0.3 Da. At 300 ppm the decoy count exceeded the real count. A wide window is therefore a choice made per dataset, and the decoy check is the means to make it.

A permethylated Orbitrap run with a 282-ion list gave the Score A counting numbers of known issue 28. Rows that list one ion more than once were 0 % of the matched rows at 20 ppm / 10 ppm and 31.8 % at 300 ppm / 0.3 Da. Scans with one peak inside the windows of two ions were 323 at 10 ppm, 1,285 at 0.3 Da and 8,844 at 0.5 Da, out of 23,193 scans.

### What these runs are for

They show how the reader diagnostics and the tolerance options behave on low-resolution ion-trap data, and how the validation tools behave on accurate-mass Orbitrap data. Known issue 27 rates the fixed-tolerance problem as none for Orbitrap-class data with a monoisotopic m/z and as high for ion-trap data whose mass is the isolation m/z. The results above are checks of behaviour. They do not add a result to the manuscript.

---

# Appendix A — Glossary

Terms specific to GlycoMSP, plus selected machine-learning and mass-spectrometry terms whose meaning may not be obvious from context. The manual assumes mass-spectrometry literacy and defines machine-learning terms in-prose at first use; this glossary is a backstop for forgetting.

**Anchor (in Score A).** An earlier design of Score A required a minimum number of anchor ions before the score was computed. The current Score A does not apply an anchor gate and treats every ion of the ion list in the same way (§4.1).

**CGA — Constraint-based Glycan Annotation.** GlycoMSP's computed-labeling pipeline. Generates candidate glycan compositions from an in-silico library, filters by precursor mass agreement, scores by fragment-ion support (Score A), optionally re-ranks by motif evidence (Score B), and produces per-spectrum labels. Manuscript-canonical name. The end-to-end pipeline including ML training is named **CGA-ML**.

**Class weight.** A per-class multiplier applied during Random Forest training that compensates for class imbalance. The `balanced` setting computes weights inversely proportional to class frequency, so rare classes contribute proportionally more to the loss function during training.

**Confidence.** The probability of the top class assigned by the Random Forest's `predict_proba`. Used by the threshold-demotion step at prediction time.

**Constraint flags.** User-supplied configuration to the in-silico library generator (`compnewv4.NGlaunch` / `compnewv4.OGlaunch`), specifying which glycan structural features to include or exclude in enumeration. Examples: hybrid-type N-glycans (`hybrid`), core fucose (`corefuc`), fucose (`allowfuc`), sialic-acid limits (`allow5ac`, `Neu5Ac_range`).

**Dataset charge.** One precursor charge that you enter in the conversion options dialog (or with `--dataset-charge` on the command line) for every MS2 spectrum that still has no charge after the instrument's own value and the isotope comb. Blank by default, which means off. A charge that the file reports is never replaced (§4.9). Added in GlycoMSP 1.13.

**Decoy library.** A copy of the in-silico library, and of the ion list, with every mass shifted by 3.37 Da. No decoy mass is a glycan or a fragment of your sample, so every match against the decoy copies is a chance match. Comparing the real and the decoy counts at the same tolerance estimates how many matches are chance (§2.4 (d), §8.5).

**Derivatization.** The chemical modification applied to glycans during sample preparation prior to mass-spec analysis (e.g., permethylation). The `Derivatization Type` metadata field tells GlycoMSP how to compute precursor mass shifts and which fragment-ion patterns to expect. The Metadata Editor offers `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `Others`. The CGA Setup window offers `PerMe`, `PerMe(Freeend)`, `PerMe(Reduced)`, `Native`, `Native(Reduced)` and `None`. With `Others`, the MAS merge asks for a custom mass shift.

**Diagnostics report.** The files `<name>.reader_diag.json` and `<name>.reader_diag.txt` that every conversion writes next to the scan CSVs. They hold the outcome, the counts of exported and skipped spectra, the counts of known and unknown precursor charges, the warnings, and the reason for any failure (§5.1). Added in GlycoMSP 1.13.

**Exit code.** The number a conversion ends with: 0 (complete), 3 (complete with warnings), 4 (partial) or 5 (failed). The command-line converter returns it, the GUI prints it in the log, and the diagnostics report records it. The command-line tools use 2 for a usage error (§3.9).

**Feature matrix.** The wide-format CSV consumed by the ML model. Rows are spectra; columns are reference-ion features (one per fragment in the ion list); cells are `log10(I) + 1` where I is the max peak intensity within the PPM window for that fragment, or `1.0` when no peak matched.

**Feature window.** The tolerance within which a peak gives the feature value of an ion-list mass. It is the *ion list ppm* field in CGA → Trainable, a fixed 10 ppm in the MAS merge, and the value you enter in the *PPM Tolerance* dialog when you build an unlabeled dataset. It is in ppm only (§4.10).

**Gate.** A filter applied to a prediction. The application has two: (1) the **confidence threshold**, which relabels a glycan prediction as the majority label (*Non-glycan* by default) when its top probability is below τ and the majority-class probability is within the majority support margin of the top probability, and (2) the **precursor MS1 gate**, which fails a glycan prediction whose predicted-vs-measured precursor mass error exceeds `gate_ppm`. Predictions that fail the precursor gate are left out of the report files.

**Glycotope.** A structural sub-element of a glycan recognisable by characteristic fragmentation. Score B operates on glycotope-level evidence rather than whole-composition evidence, which is what makes motif-aware re-ranking possible.

**In-silico library.** A computationally generated catalog of theoretical glycan compositions consistent with user-specified constraints. Each entry has a name, a monosaccharide-count tuple, and a theoretical mass. Reusable across samples sharing derivatization and ionization mode.

**Ion mining.** The GUI procedure (button "Ion Mining") that identifies m/z peaks appearing preferentially in glycan vs Non-glycan spectra, ranks them by lift, and outputs candidate diagnostic ions. (Informally a "feature-mining" step since it surfaces candidate ML features; the manuscript does not use that term.)

**Isotope comb.** The optional charge estimate of the readers. It reads the spacing of the precursor's isotope peaks in the parent MS1 scan and calls a charge of 1, 2 or 3. It is used only where the file reports no charge, and it is off by default (§4.9). Added in GlycoMSP 1.13.

**Label encoder.** A deterministic mapping of class strings (e.g., `F1H5N4`, `Non-glycan`) to integer labels (0, 1, 2, ...). Random Forest works internally with integer labels; the encoder is fitted at training time and saved alongside the model so prediction-time outputs can be decoded back to strings.

**Manual Annotation (MA) / MAS.** "MA" is the input format — an Excel workbook with the sheets `MSlist`, `ionlist`, `Glycanlist` and `information_sheet` (as in `templates/GlycoMSP_MAS_template.xlsx`). The pipeline reads the `MSlist` and `ionlist` sheets. **MAS** is "MA-Supported" and refers to the pipeline that consumes MA inputs. The end-to-end pipeline including ML training is named **MAS-ML**.

**Margin.** In the Train/Test Parameters dialog the margin is the *Majority support margin δ*. A low-confidence glycan prediction is relabelled as the majority label when the majority-class probability is within δ of the top probability. Default 0.05. The prediction report also lists a `margin` column, which is the top-class probability minus the second-class probability of one prediction.

**Motif.** A structural pattern matched by Score B rules. Each glycotope defined in the Score B workbook corresponds to a motif; rules fire when the spectrum's fragment evidence matches the motif's diagnostic ions.

**Native / Native(Reduced).** *Derivatization Type* values for underivatized glycans with a free or a reduced reducing end. With either value the CGA in-silico library holds underivatized precursor masses. The MAS route still uses permethylated residue masses (§3.8). Added in GlycoMSP 1.13.

**Outcome.** How a conversion ended: *complete*, *complete with warnings*, *partial* or *failed*. The GUI prints the outcome in the log panel, and the diagnostics report records it together with the exit code (§2.2, §5.1). Added in GlycoMSP 1.13.

**Prebuild.** The step that fills the reference database. It reads the compositions of a source file, fetches the ones the database does not hold yet from the GlyCosmos composition API, and stores them. It is the only part of GlycoMSP that contacts the network. You start it from the Manage References window or with `msp_CLI_glytoucan_prebuild.py` (§3.7, §3.9).

**Precursor.** The intact glycan ion that was selected and fragmented to produce an MS2 spectrum. Its measured mass (`protonatedmass`) is the basis for precursor-PPM matching against the in-silico library.

**Pseudolabel.** Original term for CGA-derived computed labels. Persists in the source code (file naming, function names like `run_pseudolabeling`) but is being phased out in favor of "constraint-based annotation" / "CGA result". The manuscript and this manual use the new naming.

**Random Forest (RF).** An ensemble classifier consisting of many decision trees trained on bootstrap samples of the data with random feature subsets at each split. Predictions are the majority vote across trees. Probabilities (`predict_proba`) are the per-class fractions of trees voting for that class.

**Reference database.** A local SQLite file that maps a glycan composition to its GlyTouCan accession and WURCS string. With the fill option ticked, GlycoMSP fills empty `GlyToucan ID` and `WURCS` cells from it. The option is off by default (§3.7). Added in GlycoMSP 1.11.

**Score A.** The ion-fitting score: the number of peak–ion matches divided by the number of ions in the ion list, computed against the whole ion list and limited to `[0, 1]`. A match is a peak within the Score A fragment match tolerance of an ion (default 10 ppm, in ppm or Da since GlycoMSP 1.14). Manuscript variable *a*. In CGA it decides which scans keep a label when Score B is inactive (§4.1).

**Score B.** The motif-aware score: a [0, 1]-clamped combination of motif-level support minus composition-consistency and unexpected-evidence penalties. Manuscript variable *b*. Optional; enabled by attaching a Score B workbook in the CGA Setup window.

**Status file.** The file `<name>.scan_status.tsv` that every conversion writes next to the scan CSVs. It has one row for every MS2 and MS3 spectrum the reader read and records, among other things, whether the spectrum was exported, where its charge came from and whether its precursor mass is usable (§5.1). Added in GlycoMSP 1.13.

**Stratification.** Maintaining class proportions across train / validation / test splits. With `n` total spectra and a class with `k` examples, stratification ensures each split contains approximately `(k/n) × split_fraction` examples of that class. Avoids "test fold has no rare-class examples" pathologies.

**τ (tau).** The confidence threshold for glycan predictions. Default value 0.65, with the threshold **disabled by default** in v1.10 (the checkbox *Enable confidence threshold → fallback to Majority label* is unticked). When it is enabled, a glycan prediction with a top probability below τ is relabelled as the majority label if the majority-class probability is within the majority support margin of the top probability (§6.3).

**Trainable CSV.** The wide-format CSV that ML training consumes — schema in Section 5.4. Both CGA and MAS pipelines produce this same shape so models can train on a mix.

**Training / validation / test split.** The standard ML data partition. Training data fits model parameters; an optional validation fold can tune hyperparameters and select τ during the τ-sweep; test data is held out entirely until final reporting. The GlycoMSP v1.10 default is an **80 / 20 train / test split** (`test_split = 0.20`, `val_split = 0.00` — no validation fold by default), matching the manuscript; raise `val_split` only if you want a validation fold for τ tuning.

**Unknown-charge spectra.** MS2 spectra for which neither the file nor a conversion option gave a precursor charge. They are exported with `chargeState` 0 and `protonatedmass` 0 and are marked `charge_status unknown` in the status file. They never match a library mass, and later steps do not exclude them (§4.9, §7.2).

**Unlabeled dataset.** A feature-extracted CSV without `Structure` labels. The input to prediction. *Load unlabeled dataset* picks a pre-built file, and *Create unlabeled dataset of certain experiment* traverses an experiment JSON and produces one unlabeled CSV per sample.

**WYSIWYG (in CGA design).** The manuscript's framing for CGA's enumeration philosophy: the in-silico library encodes only structural motifs the user expects to observe in their samples ("what you see is what you get"), avoiding reliance on glycogene databases that may be incomplete or organism-specific.

---

# Appendix B — Keyboard and right-click reference

GlycoMSP is primarily mouse-driven, but a handful of keyboard and right-click conventions speed up routine work.

## Dataset Explorer (Prepare Dataset window)

**Right-click on a Sample node** opens a context menu. Its entries are, in order:

- *Move to...* — a submenu that lists every experiment and sample as `{experiment} → {sample}`.
- *Link and Validate Sample* — to link and validate a sample, use the **Link Sample** button in the action grid (Section 3.1).
- *Reference database info…* — shows which reference database filled the sample's CGA and MAS results, how many compositions were filled and when (added in GlycoMSP 1.12, Section 3.7).
- *Resolve missing references* — disabled in this version. The label also shows an internal reference code in brackets.
- A *Remove …* item. The text after *Remove* names the node type.

**Right-click on a file inside a Sample** (the CSV, Excel and Metadata rows) offers *Move to...* and a *Remove* action that clears just that slot (*Remove CSV*, *Remove EXCEL*, *Remove METADATA*, for example to remove the attached Excel without removing the CSV). **Bookkeeping only — does not delete files on disk.** The CSV / Excel / metadata files remain in their folders, and you can re-add them via the *Add* buttons. On Excel rows, *CGA result tsv* rows and *CGA-based trainable csv* rows the menu also shows *Update GlyToucan ID / WURCS in this file…*, which is disabled in this version. Its label also shows an internal reference code in brackets.

**Right-click on an Experiment node** opens a menu with *Move to...* and a *Remove …* item. Experiment-level operations go through the bottom toolbar (*Set Experiment Method File*, *Save .exp.json*, etc.).

**Drag-and-drop** moves a file between Samples. Click and hold on a file row (CSV, Excel or Metadata) under a Sample, drag onto a different Sample row, release. The file moves to that Sample, within the same Experiment or into another one. Dropping anywhere that is not a Sample row is silently ignored.

The action grid of Prepare Dataset also holds **Manage References…** with the label *Active DB:* beside it (Section 3.7).

The *macOS variant* binds right-click to both `Button-3` and the macOS-conventional `Button-2` plus `Control-Button-1`, so the right-click context menu works regardless of mouse hardware.

## Modal dialog conventions

These dialogs grab focus: *Conversion options*, *Negative Sampling Options*, *Ion Mining*, *Ion suggestions*, *CGA → Trainable (one-pass)*, *ML Parameters* and *Prebuild / Update from source*. While one is open, clicks on the parent window are ignored until the dialog closes. If you find a window unresponsive, check whether a modal dialog is open behind it.

The CGA Setup window and the *Train/Test Parameters* dialog do not grab focus and stay open beside the main window. Close any of these windows with its *Close* button (or *Cancel*), or with the close control of the window frame. GlycoMSP binds no keyboard shortcuts in these dialogs.

## Main window

GlycoMSP binds no keyboard shortcuts. Buttons (*Add Raw File*, *Add mzML File*, *Clear*, *About*, *Save Log*, *Convert Raw to CSV*, *Prepare Dataset*, *Run ML Analysis*) are click-only. Closing the main window prompts a confirmation dialog titled *Quit* ("Do you really want to quit?") before exiting.

## Ion Suggestions Viewer

In the *Browse Ion Mining result* viewer:

- **Click a row** to select it.
- **Copy selected m/z** copies the selected row's first column (m/z value) to the clipboard. Useful for pasting into your ion-list Excel.
- **Open folder** opens the OS file manager at the suggestions CSV's parent directory, so you can audit the full CSV outside the viewer.

---

# Appendix C — Pointers to related documents

- **Manuscript (v1.0)** — the canonical source for naming (Score A / Score B, CGA / MAS, CGA-ML / MAS-ML), the Score B formula `b = max(0, min(1, support − penalty_comp − penalty_unexpected))`, the "Score B preferred, Score A fallback" rule, and published parameters (RF 400 trees, `max_depth=None`, `class_weight="balanced"`, 80/20 train/test). Note: the manuscript does **not** define "feature mining/FM" — see §4.7.
- **Supplementary information (v0.99)** — datasets and file mapping (Supp Table S2), external-tool settings and reports (Supp Table S3), and the short Score B overview (Supplementary Data 2).
- **Release notes** — `CHANGELOG.md` at the repository root lists every user-visible change from v1.10 to 1.14.1, with the known limitations of each version. This manual claims nothing about a version that the changelog does not.
- **Known issues** — `docs/known-issues.md` is the numbered list of known issues. This manual refers to its entries as "known issue NN" (§7.1, §7.2).
- **Tests and validation tools** — `docs/TEST_REGISTER.md` (the register of every test file), `tests/README.md` (how to run the suite and what is withheld) and `tools/validation/README.md` (the validation scripts and the decoy count method). §8.5 summarises them.
- **Reference database** — `references/README.md` (shipped database snapshots and how to select one) and `references/ATTRIBUTION.md` (data source and licence). §3.7 describes the workflow.
- **Score B detailed design** — `docs/score_b_design.md` describes the current Score B implementation. The example workbook is `templates/GlycoMSP_scoring_update_example_v7.xlsx`.
- **Developer reference (repository `docs/`)** — `docs/index.md` (overview), `docs/indexes/mspfileloaderv14.md` (function and variable maps), `docs/indexes/mspfileloaderv14_docs.md` (function reference for the main GUI), `docs/modules/<module>.md` (module pages) and `docs/routes.md` (data flow for both pipelines and scoring). Line numbers in these pages may lag the source, so confirm a location in the code before you cite it.
- **MAS-ML LOTO reproduction** — `manuscript_reproduction/`: the standalone script, `MAS_ML_LOTO_technote.md`, and the validation report; also archived in the Zenodo MAS-ML package (ID 7).
- **Reproducibility archive (Zenodo)** — all versions: <https://doi.org/10.5281/zenodo.20823041>. The v1.10 archive: <https://doi.org/10.5281/zenodo.20823042>.

---

# Appendix D — Manual changelog (v1.0 → v1.1)

Manual v1.0 described GlycoMSP v1.10, with a few passages on the CGA mass tolerances of 1.14. Manual v1.1 describes GlycoMSP 1.14.1. It adds the features of 1.11 to 1.14.1 and it corrects statements of v1.0 that did not match the code. The corrections describe behaviour that v1.10 already had, unless a version is named. They change the description and not the software.

| Section | What changed in manual v1.1 |
|---|---|
| Front matter | New version table (v1.10 to 1.14.1), which version the manuscript cites, where to get each version. |
| Part 1 | Title bar and launch command updated. Metadata fields corrected. One-sentence pointers to the reference database, the conversion diagnostics, native and ion-trap data, and the command-line tools. |
| §2.1 | Title bar of 1.14, new main-window screenshot. |
| §2.2 | New: conversion options dialog, the four outcomes and the `[reader]` log lines, the status file and the diagnostics report, unknown-charge spectra. Corrected: order of the steps (the conversion starts when the Metadata Editor opens), the Metadata Editor fields and buttons, the save steps and file names, the scan CSV columns (nine, tab-separated) and the MS3 file, the error state. |
| §2.3 | New: *Manage References…*, the *Active DB:* label and the fill checkbox. New screenshot. |
| §2.4 | CGA Setup re-described from the 1.14 window (title, frames, buttons, file rows). Ion-list header rules corrected. *Native* and *Native(Reduced)* mentioned. |
| §2.5 | Corrected: Score A has no anchor gate and is the same for all candidates of a scan. The real file name and column names of the pseudolabel TSV. New: the TSV carries every converted scan. |
| §2.6 | CGA → Trainable re-described from the 1.14 modal (boxes 1 to 7, real labels). Corrected: selection rule, negatives, file name and column order of the trainable CSV. New: the fill option. |
| §2.7 | Corrected: the *CGA dataset?* checkbox has no effect, the *Trainable File Info* panel, the labels of the Train/Test Parameters dialog, the training steps, the sidecar files (`_rf_performance.txt`), the *Training Complete* dialog. The τ-sweep statement is removed. |
| §2.8 | Corrected: the feature window at prediction is the value you enter (no fixed 20 ppm), the confidence threshold relabels to *Non-glycan*, what the precursor gate does with a prediction that fails, the list of output files. |
| §2.9 | Corrected: folder name, file names and columns of the prediction outputs, the gate files. New: the reporter files hold only gate-passed glycan predictions when the gate ran, the fill columns, the run-summary limitation. |
| §3.1–§3.6 | Corrected against the code: Score B workbook sheets and output columns, MAS merged CSV columns, ion-list header rules, validation symbols, *Combine Datasets for Training*, which steps read the shared Negative Sampling dialog, the advice on ratio and gate, ion-mining labels and details, the unlabeled-dataset walkthrough, the fill value for missing model features. |
| §3.7 | New section: reference database and the Manage References window. |
| §3.8 | New section: native glycans and ion-trap data. The worked example of v1.0 §7.2 moved here. |
| §3.9 | New section: command-line tools. |
| §4.1 | Corrected: no anchor gate, Score A is computed against the whole ion list, how the candidate of a scan is chosen without Score B. |
| §4.3 | Corrected: feature window at prediction, the MAS merge takes the last in-window peak (known issue 15), the zero-intensity error (known issue 33). |
| §4.4 | Corrected: what the confidence threshold does and what the run summary shows (known issue 29). |
| §4.5 | Corrected: the fields of `training_run.json`. |
| §4.9 | New section: what the readers do and do not do. |
| §4.10 | New section: precursor and fragment tolerances. |
| §5.1 | Scan CSV and MS3 tables rebuilt. New: unknown-charge convention, status file, diagnostics report, conversion entries of the metadata sidecar. Corrected: metadata keys, annotation sheets, ion-list headers, in-silico library columns. |
| §5.2 | New: `artifacts.glytoucan_resolution`. The `parameters.cga` block re-verified. |
| §5.3 | Corrected: header values of the experiment JSON, the keys of `train_parameters`. |
| §5.4 | Rebuilt from the code: conversion outputs, pseudolabel TSV, trainable CSV, prediction folder and files, gate files, training sidecars. New: Random Forest timing report, reference database files. |
| §5.5 | Every row re-verified. New rows for the options of 1.11 to 1.14.1. |
| §5.6 | Corrected: the six required sheets and their columns. New: the two Score B parsing changes of 1.12. |
| Part 6 | Corrected: the ML parameters editor (button names, merge order), the gate formula, the gate and review files, the counts-based fallback, the τ-sweep statement, the threshold label. |
| §7.1 | New: conversion outcomes, profile mzML, ion-trap peak lists, reference-database errors, the batch message after a failed finalize, the command-line prediction tool (known issues 30 to 32), the zero-intensity error (known issue 33). Several v1.0 entries corrected against the dialogs the code shows. Sub-headings added. |
| §7.2 | Brought to 1.14.1 and ordered by workflow stage. The ion-trap entry is shortened and points to §3.8. Three v1.0 entries corrected (negative sampling, mixed-source training, saved CGA settings). New entries include the ppm-only feature window and the run-summary limitation (known issue 29). |
| §8.0–§8.4 | Unchanged, apart from one sentence that these sections refer to v1.10 and the source-code line of §8.1, which no longer names a tag (the repository has no release tags). |
| §8.5 | New section: automated tests and validation tools. |
| §8.6 | New section: external-dataset validation runs. |
| Appendix A | New terms for 1.11 to 1.14.1. Corrected entries: Anchor, Derivatization, Gate, Manual Annotation, Margin, Score A, τ. |
| Appendix B | Right-click menus, drag-and-drop and dialog conventions corrected against the code. |
| Appendix C | New pointers: changelog, known issues, test register, reference database, Score B design page. The Zenodo entry now gives the DOI for all versions and the DOI of the v1.10 archive. |
| Appendix D | New. |
| Screenshots | Main window, Prepare Dataset, CGA Setup and CGA → Trainable replaced by 1.14 captures. New: conversion options dialog, Manage References window, *Prebuild / Update from source* dialog, sample right-click menu, reference-database info popup. |

---

*This PDF was exported from the source Markdown, with revisions requested and confirmed by the author. It is provided as a reference copy in both the Zenodo archive and the GitHub repository. The canonical source Markdown is maintained on GitHub and will continue to be updated with future maintenance and version releases.*
