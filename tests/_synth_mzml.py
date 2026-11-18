"""Synthetic mzML builder for reader tests (READER-DIAG, 2026-10-13).

All spectra are made up here; no third-party data. Produces a small non-indexed mzML 1.1 file
that pymzml can iterate.
"""
import base64
import struct
import zlib

_NS = "http://psi.hupo.org/ms/mzml"
_NATIVE = "controllerType=0 controllerNumber=1 scan={}"


def encode_array(values, zlib_compress=True):
    raw = struct.pack("<%df" % len(values), *values)
    if zlib_compress:
        raw = zlib.compress(raw)
    return base64.b64encode(raw).decode("ascii")


def _cv(acc, name, value="", extra=""):
    return f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value="{value}"{extra}/>'


def _array_xml(kind, values, binary=None):
    acc, name = ("MS:1000514", "m/z array") if kind == "mz" else ("MS:1000515", "intensity array")
    payload = encode_array(values) if binary is None else binary
    return (f'<binaryDataArray encodedLength="{len(payload)}">'
            f'{_cv("MS:1000521", "32-bit float")}{_cv("MS:1000574", "zlib compression")}'
            f'{_cv(acc, name)}<binary>{payload}</binary></binaryDataArray>')


def spectrum_xml(index, scan, ms_level, mz, intensity, *, polarity="-", centroid=True,
                 precursor=None, mz_binary=None, filter_string=None, representation=True):
    """One <spectrum>. `precursor` is a dict with optional keys parent, mz, charge, intensity, and more_ions
    (a list of dicts with mz, charge, intensity: further <selectedIon>s after the first; R3-F3, 2026-10-26).

    A key that is absent (or None) is left out of the XML, which is how a converter writes
    missing metadata. `mz_binary` overrides the encoded m/z payload (for decode-error cases).
    `filter_string` adds a Thermo filter string to the scan; representation=False omits the
    centroid/profile term.
    """
    parts = [
        f'<spectrum index="{index}" id="{_NATIVE.format(scan)}" defaultArrayLength="{len(mz)}">',
        _cv("MS:1000511", "ms level", ms_level),
        _cv("MS:1000130", "positive scan") if polarity == "+" else _cv("MS:1000129", "negative scan"),
        ("" if not representation else
         _cv("MS:1000127", "centroid spectrum") if centroid else _cv("MS:1000128", "profile spectrum")),
        '<scanList count="1">', _cv("MS:1795069", "no combination"),
        f'<scan>{_cv("MS:1000016", "scan start time", 1.0 + index / 100.0)}'
        f'{_cv("MS:1000512", "filter string", filter_string) if filter_string else ""}</scan></scanList>',
    ]
    if precursor is not None:
        ref = precursor.get("parent")
        ref_attr = f' spectrumRef="{_NATIVE.format(ref)}"' if ref is not None else ""
        ions = []
        for spec in [precursor] + list(precursor.get("more_ions") or ()):
            ion = []
            if spec.get("mz") is not None:
                ion.append(_cv("MS:1000744", "selected ion m/z", spec["mz"]))
            if spec.get("charge") is not None:
                ion.append(_cv("MS:1000041", "charge state", spec["charge"]))
            if spec.get("intensity") is not None:
                ion.append(_cv("MS:1000042", "peak intensity", spec["intensity"]))
            ions.append(f'<selectedIon>{"".join(ion)}</selectedIon>')
        parts.append(
            f'<precursorList count="1"><precursor{ref_attr}>'
            f'<selectedIonList count="{len(ions)}">{"".join(ions)}</selectedIonList>'
            f'<activation>{_cv("MS:1000133", "collision-induced dissociation")}</activation>'
            f'</precursor></precursorList>')
    parts.append('<binaryDataArrayList count="2">')
    parts.append(_array_xml("mz", mz, mz_binary))
    parts.append(_array_xml("i", intensity))
    parts.append('</binaryDataArrayList></spectrum>')
    return "".join(parts)


def build_mzml(spectra_xml, run_id="synthetic_run"):
    body = "\n".join(spectra_xml)
    return (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        f'<mzML xmlns="{_NS}" version="1.1.0">\n'
        '<cvList count="1"><cv id="MS" fullName="PSI-MS" URI="https://example.invalid/psi-ms.obo"/></cvList>\n'
        f'<run id="{run_id}">\n'
        f'<spectrumList count="{len(spectra_xml)}">\n{body}\n</spectrumList>\n'
        '</run>\n</mzML>\n'
    )


def write_mzml(path, spectra_xml):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(build_mzml(spectra_xml))
    return str(path)
