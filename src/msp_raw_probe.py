"""Read-only probe of the Thermo RAW backend for selected scans (READER-DIAG CP-C).

version 0.3, 20261023 (v0.2 READER-DIAG CP-D.1: the alias replacement runs in one pass; v0.3 R2-F5: the shared
`msp_spectrum_record.Redactor` (stem as a whole word, escaped paths, case) and dict keys are no longer redacted, so a
short stem cannot corrupt the trailer labels "Charge State" / "Master Scan Number" before they are read)

Purpose: find out, for a few native scans, which backend requests return peaks and what charge
information the file holds, before any production change. Nothing here is used by conversion.

Per requested scan it records:
- MS order, scan filter, profile/centroid flags, scan header and trailer-extra values
  (the trailer carries the instrument's "Charge State" and "Monoisotopic M/Z"),
- precursor info as returned by GetPrecursorInfoFromScanNum (the call conversion uses),
- three independent peak requests: GetLabelData (what conversion uses today),
  GetMassListFromScanNum with centroiding off, and with centroiding on.
For an MS2 scan it also probes the parent MS1 scan and keeps the points within +/- window m/z of
the precursor (profile and centroid), which is the input any isotope-spacing charge check needs.

Every backend call is attempted on its own; an exception in one is recorded and does not stop
the others. The wrapper's real call signatures are recorded too, since this probe was written
without the Windows backend at hand.

Privacy: outputs are named by a caller-supplied alias. The input path, file name and stem are
replaced by the alias in every recorded string; traceback frames keep only module:line. Creator,
serial number, acquisition date and file name accessors are not called.
"""
import importlib
import inspect as _inspect
import json
import os
import platform
import re
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msp_spectrum_record import Redactor as _SharedRedactor  # noqa: E402

version = "0.3"
last_update = 20261023

REQUESTS = ("label_data", "mass_list_centroid_off", "mass_list_centroid_on")
_SAFE_FILE_CALLS = ("GetNumSpectra", "GetFirstSpectrumNumber", "GetLastSpectrumNumber",
                    "GetVersionNumber", "GetInstModel", "GetInstName", "GetInstSoftwareVersion",
                    "GetInstHardwareVersion", "GetNumberOfControllers", "GetMassResolution")
_CALLS_RECORDED = ("GetLabelData", "GetMassListFromScanNum", "GetPrecursorInfoFromScanNum",
                   "GetTrailerExtraForScanNum", "GetScanHeaderInfoForScanNum", "GetFilterForScanNum",
                   "GetMSOrderForScanNum", "IsProfileScanForScanNum", "IsCentroidScanForScanNum")
_MASS_LIST_KW = ("scanNumber", "scanFilter", "intensityCutoffType", "intensityCutoffValue",
                 "nMaxNumberOfPeaks", "centroidResult", "centroidPeakWidth")


class Redactor(_SharedRedactor):
    """The shared redactor (msp_spectrum_record); values that are not strings pass through unchanged."""

    def __init__(self, input_path, alias):
        super().__init__(input_path, alias)
        self.alias = alias

    def __call__(self, value):
        return super().__call__(value) if isinstance(value, str) else value


def _jsonable(value, redact, depth=0):
    """Convert backend return values (tuples, COM/ctypes objects) into JSON-safe data."""
    if depth > 6:
        return redact(repr(value))
    if value is None or isinstance(value, (bool, int)):
        return value
    if isinstance(value, float):
        return value if value == value and value not in (float("inf"), float("-inf")) else repr(value)
    if isinstance(value, str):
        return redact(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v, redact, depth + 1) for k, v in value.items()}   # keys: backend labels
    if isinstance(value, (list, tuple)):
        return [_jsonable(v, redact, depth + 1) for v in value]
    try:
        return _jsonable(float(value), redact, depth + 1)
    except (TypeError, ValueError):
        return redact(repr(value))


def _short_traceback(exc):
    frames = traceback.extract_tb(exc.__traceback__)
    return [f"{os.path.basename(f.filename)}:{f.lineno} in {f.name}" for f in frames]


def _attempt(reader, name, args, kwargs, redact):
    """Call reader.<name>(*args, **kwargs); never raises. Returns (record, raw_value)."""
    rec = {"api": name, "args": _jsonable(list(args), redact), "kwargs": _jsonable(kwargs, redact)}
    fn = getattr(reader, name, None)
    if fn is None:
        rec.update(ok=False, exception="AttributeError", message=f"backend has no {name}")
        return rec, None
    started = time.perf_counter()
    try:
        value = fn(*args, **kwargs)
    except Exception as exc:  # backend (COM) errors are Exception subclasses
        rec.update(ok=False, exception=type(exc).__name__, message=redact(str(exc)),
                   traceback=_short_traceback(exc))
        value = None
    else:
        rec["ok"] = True
    rec["seconds"] = round(time.perf_counter() - started, 4)
    return rec, value


def _shape(value, depth=0):
    """Describe the nesting and lengths of a return value, e.g. 'tuple[2](tuple[6](tuple[15]..'."""
    if isinstance(value, (list, tuple)) and depth < 3:
        inner = _shape(value[0], depth + 1) if value else ""
        return f"{type(value).__name__}[{len(value)}]" + (f"({inner})" if inner else "")
    return type(value).__name__


def _two_rows(value):
    """Pull (masses, intensities) out of a label-data or mass-list return.

    Both wrappers return (array_2xN, flags); conversion reads value[0][0] and value[0][1].
    """
    try:
        rows = value[0]
        masses = [float(x) for x in rows[0]]
        intensities = [float(x) for x in rows[1]]
        return masses, intensities, rows
    except Exception:
        return None, None, None


def _mass_list_call(reader, scan, centroid):
    """Arguments for GetMassListFromScanNum: no filter, no intensity cutoff, all points."""
    values = (scan, "", 0, 0, 0, bool(centroid), 0.0)
    fn = getattr(reader, "GetMassListFromScanNum", None)
    try:
        params = list(_inspect.signature(fn).parameters) if fn is not None else []
    except (TypeError, ValueError):
        params = []
    if all(p in params for p in _MASS_LIST_KW):
        return (), dict(zip(_MASS_LIST_KW, values))
    return values, {}


def _peak_request(reader, scan, request, redact):
    if request == "label_data":
        rec, value = _attempt(reader, "GetLabelData", (scan,), {}, redact)
    else:
        args, kwargs = _mass_list_call(reader, scan, request == "mass_list_centroid_on")
        rec, value = _attempt(reader, "GetMassListFromScanNum", args, kwargs, redact)
    rec["requested_scan"] = scan
    rec["returned_scan"] = None   # neither wrapper call returns a scan number
    masses = intensities = None
    if rec["ok"]:
        rec["return_shape"] = _shape(value)
        masses, intensities, rows = _two_rows(value)
        if masses is None:
            rec["parse_error"] = "could not read value[0][0] / value[0][1]"
        else:
            rec["n_mz"] = len(masses)
            rec["n_intensity"] = len(intensities)
            if request == "label_data" and rows is not None and len(rows) >= 6:
                # label rows: mass, intensity, resolution, baseline, noise, charge
                charges = sorted({int(round(float(c))) for c in rows[5]}) if len(rows[5]) else []
                rec["label_charges_present"] = charges
            try:
                flags = value[1]
                rec["n_flags"] = len(flags) if flags is not None else None
            except Exception:
                rec["n_flags"] = None
    return rec, masses, intensities


def _scan_info(reader, scan, redact):
    info = {"scan": scan, "calls": []}
    for name, key in (("GetMSOrderForScanNum", "ms_order"), ("GetFilterForScanNum", "filter"),
                      ("IsProfileScanForScanNum", "is_profile"), ("IsCentroidScanForScanNum", "is_centroid"),
                      ("GetPrecursorInfoFromScanNum", "precursor_info"),
                      ("GetTrailerExtraForScanNum", "trailer_extra"),
                      ("GetScanHeaderInfoForScanNum", "scan_header")):
        rec, value = _attempt(reader, name, (scan,), {}, redact)
        info["calls"].append(rec)
        info[key] = _jsonable(value, redact) if rec["ok"] else None
    return info


def _trailer_value(trailer, *labels):
    if not isinstance(trailer, dict):
        return None
    norm = {re.sub(r"[\s:]+$", "", str(k)).strip().lower(): v for k, v in trailer.items()}
    for label in labels:
        if label.lower() in norm:
            return norm[label.lower()]
    return None


def _number(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None


def _precursor_summary(info):
    """Collect every charge / precursor m/z / parent scan the file offers for one MS2 scan."""
    pinfo = info.get("precursor_info")
    flat = pinfo
    if isinstance(pinfo, list) and len(pinfo) == 1 and isinstance(pinfo[0], list):
        flat = pinfo[0]   # some wrappers return a list of precursor records
    out = {"from_precursor_info": None, "from_trailer": None, "from_filter": None}
    if isinstance(flat, list) and len(flat) >= 4:
        # order as used by convert_raw_to_csv: [0] isolation m/z, [1] mono m/z, [2] charge, [3] parent scan
        out["from_precursor_info"] = {"isolation_mz": _number(flat[0]), "mono_mz": _number(flat[1]),
                                      "charge": _number(flat[2]), "parent_scan": _number(flat[3])}
    trailer = info.get("trailer_extra")
    out["from_trailer"] = {"charge_state": _trailer_value(trailer, "charge state"),
                           "monoisotopic_mz": _trailer_value(trailer, "monoisotopic m/z"),
                           "master_scan": _trailer_value(trailer, "master scan number")}
    match = re.search(r"([0-9]+\.[0-9]+)@", info.get("filter") or "")
    out["from_filter"] = {"precursor_mz": float(match.group(1)) if match else None}
    return out


def _parent_scan(reader, scan, summary, redact, max_back=50):
    """Parent MS1 scan: precursor info, else trailer master scan, else nearest earlier MS1."""
    pi = summary.get("from_precursor_info") or {}
    if pi.get("parent_scan"):
        return int(pi["parent_scan"]), "precursor_info"
    master = _number((summary.get("from_trailer") or {}).get("master_scan"))
    if master:
        return int(master), "trailer_master_scan"
    for back in range(scan - 1, max(0, scan - max_back), -1):
        rec, order = _attempt(reader, "GetMSOrderForScanNum", (back,), {}, redact)
        if rec["ok"] and order == 1:
            return back, "nearest_earlier_ms1"
    return None, "not_found"


def _precursor_mz(summary):
    pi = summary.get("from_precursor_info") or {}
    for value in (pi.get("mono_mz"), pi.get("isolation_mz"),
                  _number((summary.get("from_trailer") or {}).get("monoisotopic_mz")),
                  (summary.get("from_filter") or {}).get("precursor_mz")):
        if value:
            return value
    return None


def _write_points(path, masses, intensities):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("mz\tintensity\n")
        for m, i in zip(masses, intensities):
            fh.write(f"{m!r}\t{i!r}\n")


def _environment(reader, redact):
    env = {"python": sys.version.split()[0], "platform": platform.platform(),
           "probe_version": version, "packages": {}}
    for pkg in ("pymsfilereader", "comtypes"):
        try:
            from importlib import metadata
            env["packages"][pkg] = metadata.version(pkg)
        except Exception:
            env["packages"][pkg] = None
    rec, value = _attempt(reader, "Version", (), {}, redact) if reader is not None else ({"ok": False}, None)
    env["backend_library_version"] = _jsonable(value, redact) if rec.get("ok") else None
    env["signatures"] = {}
    for name in _CALLS_RECORDED + ("Version",):
        fn = getattr(reader, name, None) if reader is not None else None
        try:
            env["signatures"][name] = str(_inspect.signature(fn)) if fn is not None else None
        except (TypeError, ValueError):
            env["signatures"][name] = "unavailable"
    return env


def inspect_raw_scans(raw_path, scans, alias, outdir, window=3.0, full_points=False, reader_factory=None):
    """Probe selected native scans of a RAW file. Writes <alias>.inspect.json/.txt and point TSVs.

    Returns the report dict. Never raises for backend errors; an open failure gives outcome 'failed'.
    `reader_factory` defaults to the MSFileReader class conversion uses (msprawextractor).
    """
    redact = Redactor(raw_path, alias)
    os.makedirs(outdir, exist_ok=True)
    report = {"json_type": "glycomsp.raw_probe", "schema_version": "1.0.0", "alias": alias,
              "outcome": None, "settings": {"scans": list(scans), "window_mz": window,
                                            "full_points": bool(full_points)},
              "input": {"size_bytes": os.path.getsize(raw_path) if os.path.exists(raw_path) else None},
              "file": {}, "scans": [], "files_written": [],
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if reader_factory is None:
        rawext = importlib.import_module("msprawextractor")
        reader_factory = getattr(rawext, "MSFileReader", None) if getattr(rawext, "pymsreader", False) else None
    reader = None
    if reader_factory is None:
        report["open"] = {"ok": False, "exception": "ImportError",
                          "message": "pymsfilereader / Thermo MSFileReader not available on this machine"}
    else:
        try:
            reader = reader_factory(raw_path)
            report["open"] = {"ok": True}
        except Exception as exc:
            report["open"] = {"ok": False, "exception": type(exc).__name__, "message": redact(str(exc)),
                              "traceback": _short_traceback(exc)}
    report["environment"] = _environment(reader, redact)

    if reader is not None:
        for name in _SAFE_FILE_CALLS:
            rec, value = _attempt(reader, name, (), {}, redact)
            report["file"][name] = _jsonable(value, redact) if rec["ok"] else {"error": rec.get("exception")}
        for scan in scans:
            report["scans"].append(_probe_one(reader, int(scan), alias, outdir, window, full_points,
                                              redact, report["files_written"]))
        rec, _ = _attempt(reader, "Close", (), {}, redact)
        report["close"] = rec
    report["outcome"] = "complete" if reader is not None else "failed"
    report["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    json_path = os.path.join(outdir, f"{alias}.inspect.json")
    with open(json_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(report, fh, indent=2)
    txt_path = os.path.join(outdir, f"{alias}.inspect.txt")
    with open(txt_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(summary_text(report))
    return report


def _probe_one(reader, scan, alias, outdir, window, full_points, redact, written):
    info = _scan_info(reader, scan, redact)
    entry = {"scan": scan, "ms_order": info["ms_order"], "filter": info["filter"],
             "is_profile": info["is_profile"], "is_centroid": info["is_centroid"],
             "precursor_info_raw": info["precursor_info"], "trailer_extra": info["trailer_extra"],
             "scan_header": info["scan_header"], "info_calls": info["calls"], "requests": {}}
    for request in REQUESTS:
        rec, masses, intensities = _peak_request(reader, scan, request, redact)
        if masses is not None and (full_points or entry["ms_order"] != 1):
            name = f"{alias}_scan{scan}_{request}.tsv"
            _write_points(os.path.join(outdir, name), masses, intensities)
            rec["points_file"] = name
            written.append(name)
        entry["requests"][request] = rec

    if entry["ms_order"] is not None and entry["ms_order"] >= 2:
        summary = _precursor_summary(info)
        entry["charge_sources"] = summary
        parent, how = _parent_scan(reader, scan, summary, redact)
        target = _precursor_mz(summary)
        entry["parent"] = {"scan": parent, "found_by": how, "precursor_mz_used": target,
                           "window_mz": window, "requests": {}}
        if parent is not None:
            pinfo = _scan_info(reader, parent, redact)
            entry["parent"].update(ms_order=pinfo["ms_order"], filter=pinfo["filter"],
                                   is_profile=pinfo["is_profile"], is_centroid=pinfo["is_centroid"])
            for request in REQUESTS:
                rec, masses, intensities = _peak_request(reader, parent, request, redact)
                if masses is not None and target is not None:
                    keep = [(m, i) for m, i in zip(masses, intensities) if abs(m - target) <= window]
                    rec["n_in_window"] = len(keep)
                    name = f"{alias}_scan{scan}_parent{parent}_{request}_window.tsv"
                    _write_points(os.path.join(outdir, name), [m for m, _ in keep], [i for _, i in keep])
                    rec["window_points_file"] = name
                    written.append(name)
                entry["parent"]["requests"][request] = rec
    return entry


def summary_text(report):
    lines = [f"GlycoMSP RAW probe v{version} - alias {report['alias']} - outcome {report['outcome']}",
             f"open: {report.get('open')}",
             f"backend library: {report['environment'].get('backend_library_version')}; "
             f"packages: {report['environment'].get('packages')}",
             f"file: {report.get('file')}", ""]
    for entry in report["scans"]:
        lines.append(f"scan {entry['scan']}: MS{entry['ms_order']} profile={entry['is_profile']} "
                     f"centroid={entry['is_centroid']} filter={entry['filter']!r}")
        for request, rec in entry["requests"].items():
            lines.append(_request_line(request, rec))
        if "charge_sources" in entry:
            lines.append(f"  charge sources: {entry['charge_sources']}")
            par = entry.get("parent") or {}
            lines.append(f"  parent MS1 {par.get('scan')} (by {par.get('found_by')}), precursor m/z used "
                         f"{par.get('precursor_mz_used')}, profile={par.get('is_profile')}")
            for request, rec in (par.get("requests") or {}).items():
                lines.append(_request_line("parent " + request, rec) + f", in +/-{par.get('window_mz')} m/z: "
                             f"{rec.get('n_in_window')}")
        lines.append("")
    return "\n".join(lines) + "\n"


def _request_line(label, rec):
    if not rec.get("ok"):
        return f"  {label}: FAILED {rec.get('exception')}: {rec.get('message')}"
    return (f"  {label}: ok, shape {rec.get('return_shape')}, m/z {rec.get('n_mz')}, "
            f"intensity {rec.get('n_intensity')}" + (f", parse_error {rec['parse_error']}" if rec.get("parse_error") else ""))
