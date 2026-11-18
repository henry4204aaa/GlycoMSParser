
We've tested some datasets on native glycans using LTQ, and noticed the MS1 precursor error might be larger than suggested in user manual.

> 100 - 200 ppm

Things to know before you review:
- The mass now uses the isolation m/z when the file has no monoisotopic m/z, with or without the options. It was part of the plan; I'm repeating it because it's the only change that applies with the options off.
  - Before, such rows got mass 0, or a negative mass if a charge was reported.
  - Rows that have a monoisotopic m/z are unchanged.
- Deploying the RAW converter to Drive now takes three files: msprawextractor.py, plus msp_charge_comb.py and msp_spectrum_record.py. Both are new to Drive, and the converter won't import without them. v14 itself is unaffected.
- Expect many scans the comb leaves unknown on LTQ data. On synthetic isotope patterns with peaks as wide as on your LTQ scans:
  - z = 3 is never called; it is left unknown, never given a wrong charge;
  - z = 2 only just clears the margin (0.12);
  - your dataset value fills those scans.