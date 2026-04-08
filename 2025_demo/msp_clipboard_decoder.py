version = "1.00"
last_update = 20260409
# I decided the spec and needed components from GlycoMSP, then let the Claude Code (opus4.6) do the vibe coding
# QA and user tests are performed manually beside cli tests, and the code review has been done by author to confirm the behaviors as expected

"""Decode mass spectrum data pasted from vendor software clipboard exports.

Currently supports:
- Thermo Xcalibur "Exact Mass" clipboard export (SPECTRUM - MS format)

Other vendor formats (Waters MassLynx, Bruker DataAnalysis, Agilent MassHunter)
are not yet supported and will produce a warning.

Usage::

    from msp_clipboard_decoder import parse_clipboard

    result = parse_clipboard(clipboard_text)
    if result["format"] == "unsupported":
        print(result["message"])
    else:
        entry_dict = result["entry"]
"""

from __future__ import annotations

import re


def detect_format(text: str) -> str:
    """Sniff clipboard text and return the detected vendor format.

    Returns:
        ``"thermo"`` — Thermo Xcalibur Exact Mass export
        ``"unsupported"`` — unrecognised format
    """
    lines = text.strip().splitlines()
    if not lines:
        return "unsupported"

    # Thermo: first line is "SPECTRUM - MS" (possibly with trailing whitespace)
    if lines[0].strip().upper().startswith("SPECTRUM"):
        return "thermo"

    # Thermo variant: sometimes starts with the filter string directly
    # e.g. "FTMS + p NSI ..."
    if lines[0].strip().upper().startswith("FTMS") or lines[0].strip().upper().startswith("ITMS"):
        return "thermo"

    return "unsupported"


def parse_thermo_clipboard(text: str) -> dict:
    """Parse Thermo Xcalibur 'Exact Mass' clipboard export.

    Expected format::

        SPECTRUM - MS
        sample_file.raw
        FTMS + p NSI d Full ms2 1041.2157@hcd15.00 [90.0000-2000.0000]
        Scan #: 5467-32694
        RT: 24.26-61.27
        AV: 12
        m/z\tIntensity\tRelative\tCharge

        111.04\t146633.0\t12.26\t1.00
        ...

    Returns:
        Dict with keys matching curated library entry fields, plus extra
        metadata under ``"_meta"``.
    """
    lines = text.strip().splitlines()

    # ---- Parse header lines ----
    sample_id = ""
    filter_line = ""
    scan_range = ""
    rt_range = ""
    avg_count = ""
    header_end = 0  # line index where peak data starts

    i = 0
    # Skip the "SPECTRUM - MS" header if present
    if lines[i].strip().upper().startswith("SPECTRUM"):
        i += 1

    # Next non-empty line: RAW filename (or filter if no filename)
    while i < len(lines) and not lines[i].strip():
        i += 1

    if i < len(lines):
        candidate = lines[i].strip()
        # If it looks like a filename (has a dot extension), it's the sample ID
        if re.search(r"\.\w{2,5}$", candidate) and not candidate.upper().startswith(("FTMS", "ITMS")):
            sample_id = candidate
            i += 1

    # Filter string line (FTMS/ITMS ...)
    while i < len(lines) and not lines[i].strip():
        i += 1

    if i < len(lines) and re.match(r"(FTMS|ITMS)", lines[i].strip(), re.IGNORECASE):
        filter_line = lines[i].strip()
        i += 1
    elif i < len(lines) and not filter_line:
        # Could still be a filter line with different format
        filter_line = lines[i].strip()
        i += 1

    # Optional metadata lines: Scan #, RT, AV
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.lower().startswith("scan"):
            scan_range = re.sub(r"^scan\s*#?\s*:\s*", "", line, flags=re.IGNORECASE).strip()
            i += 1
        elif line.lower().startswith("rt"):
            rt_range = re.sub(r"^rt\s*:\s*", "", line, flags=re.IGNORECASE).strip()
            i += 1
        elif line.lower().startswith("av"):
            avg_count = re.sub(r"^av\s*:\s*", "", line, flags=re.IGNORECASE).strip()
            i += 1
        elif line.lower().startswith("m/z") or "\t" in line:
            # Column header or data line — we've reached the table
            if line.lower().startswith("m/z"):
                i += 1  # skip column header
            break
        else:
            i += 1

    # ---- Parse filter string for instrument metadata ----
    precursor_mz = 0.0
    ion_mode = ""
    collision_energy = ""
    ms_level = ""
    scan_range_mz = ""

    if filter_line:
        # Polarity: "+" or "-" after FTMS/ITMS
        if re.search(r"\s\+\s", filter_line):
            ion_mode = "positive"
        elif re.search(r"\s-\s", filter_line):
            ion_mode = "negative"

        # MS level
        if "ms2" in filter_line.lower() or "ms²" in filter_line.lower():
            ms_level = "MS2"
        elif "ms3" in filter_line.lower():
            ms_level = "MS3"
        else:
            ms_level = "MS1"

        # Precursor m/z: number before @ symbol
        m_prec = re.search(r"(\d+\.?\d*)\s*@", filter_line)
        if m_prec:
            precursor_mz = float(m_prec.group(1))

        # Collision energy: hcdNN.NN or cidNN.NN
        m_ce = re.search(r"(hcd|cid|etd)(\d+\.?\d*)", filter_line, re.IGNORECASE)
        if m_ce:
            collision_energy = f"{m_ce.group(1).lower()}{m_ce.group(2)}"

        # Scan range: [low-high]
        m_range = re.search(r"\[(\d+\.?\d*)\s*-\s*(\d+\.?\d*)\]", filter_line)
        if m_range:
            scan_range_mz = f"{m_range.group(1)}-{m_range.group(2)}"

    # ---- Parse peak table ----
    mzs: list[float] = []
    intensities: list[float] = []
    charges: list[float] = []

    # Skip blank lines before data
    while i < len(lines) and not lines[i].strip():
        i += 1

    for line_idx in range(i, len(lines)):
        line = lines[line_idx].strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            # Try whitespace split as fallback
            parts = line.split()
        if len(parts) >= 2:
            try:
                mz_val = float(parts[0])
                int_val = float(parts[1])
                mzs.append(mz_val)
                intensities.append(int_val)
                if len(parts) >= 4:
                    try:
                        charges.append(float(parts[3]))
                    except ValueError:
                        charges.append(0.0)
            except ValueError:
                continue  # skip unparseable lines

    # ---- Build scan identifier ----
    # Composite spectra don't have a single scan number.
    # Use "AVnn_range" format to mark it as averaged.
    if avg_count:
        scan_id = f"AV{avg_count}_{scan_range}" if scan_range else f"AV{avg_count}"
    elif scan_range and "-" in scan_range:
        scan_id = f"range_{scan_range}"
    else:
        scan_id = scan_range or "0"

    # ---- Build notes ----
    notes_parts = []
    if avg_count:
        notes_parts.append(f"Averaged spectrum (AV:{avg_count})")
    if rt_range:
        notes_parts.append(f"RT {rt_range}")
    if scan_range:
        notes_parts.append(f"Scans {scan_range}")
    if collision_energy:
        notes_parts.append(f"CE: {collision_energy}")
    if scan_range_mz:
        notes_parts.append(f"Scan range m/z: {scan_range_mz}")

    # Convert peaks to semicolon-separated strings for DB storage
    peaklist_str = ";".join(f"{mz:.4f}" for mz in mzs)
    peakintensity_str = ";".join(f"{inten:.1f}" for inten in intensities)

    # ---- Build entry dict ----
    entry = {
        "sample_id": sample_id,
        "ms2_scan_no": scan_id,
        "precursor_mz": precursor_mz,
        "ion_mode": ion_mode,
        "collision_energy": collision_energy,
        "peaklist": peaklist_str,
        "peakintensity": peakintensity_str,
        "peak_count": len(mzs),
        "annotation_source": "manual_clipboard",
        "confidence": "tentative",
        "diagnostic_notes": "; ".join(notes_parts) if notes_parts else "",
        # Left empty for user to fill:
        "composition": "",
        "structural_notes": "",
        "precursor_adduct": "",
        "precursor_charge": "",
        "derivatization": "",
    }

    # Set default adduct based on polarity
    if ion_mode == "negative":
        entry["precursor_adduct"] = "[M-H]-"
        entry["precursor_charge"] = "-1"
    elif ion_mode == "positive":
        entry["precursor_adduct"] = "[M+Na]+"
        entry["precursor_charge"] = "1"

    meta = {
        "format": "thermo",
        "filter_string": filter_line,
        "scan_range": scan_range,
        "rt_range": rt_range,
        "avg_count": avg_count,
        "ms_level": ms_level,
        "scan_range_mz": scan_range_mz,
        "charges": charges,
    }

    return {"entry": entry, "_meta": meta}


def parse_clipboard(text: str) -> dict:
    """Parse clipboard text from any supported vendor format.

    Returns:
        Dict with:
        - ``"format"``: ``"thermo"`` or ``"unsupported"``
        - ``"entry"``: Dict of curated library fields (if supported)
        - ``"_meta"``: Extra metadata not stored in DB (if supported)
        - ``"message"``: Human-readable warning (if unsupported)
    """
    if not text or not text.strip():
        return {
            "format": "unsupported",
            "message": "Clipboard is empty.",
        }

    fmt = detect_format(text)

    if fmt == "thermo":
        result = parse_thermo_clipboard(text)
        result["format"] = "thermo"
        return result

    return {
        "format": "unsupported",
        "message": (
            "Unrecognised clipboard format. Currently only Thermo Xcalibur "
            "'Exact Mass' clipboard export is supported.\n\n"
            "Tip: In Xcalibur, right-click the spectrum → Export → Clipboard (Exact Mass)."
        ),
    }
