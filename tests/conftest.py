"""Session-wide bookkeeping for the GlycoMSP test suite (added 2026-09-23, Codex S1 R3-F3).

Records the value of GLYCOMSP_REFERENCE_DB at session start so `tests/test_s1_zz_env_isolation.py`
(collected last by name) can assert that no earlier test module leaked a database override.
"""
import os

GLYCOMSP_ENV_AT_SESSION_START = {"GLYCOMSP_REFERENCE_DB": os.environ.get("GLYCOMSP_REFERENCE_DB"),
                                 "GLYCOMSP_GLYTOUCAN_LOG": os.environ.get("GLYCOMSP_GLYTOUCAN_LOG")}
