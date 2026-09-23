"""
msp_glytoucan_resolver.py — GlyTouCan accession / WURCS resolution for GlycoMSP (v1.20, S1 / T1).

Design spec: handoff/WURCS_designspec20260921.md (read it before changing behaviour here).

Three layers, stdlib only (sqlite3, urllib, json):

  1. Canonicalisation   canonical_composition / composition_to_api_query
  2. Store              ReferenceDB  (SQLite; the ONLY lookup source at run time)
                        append_log   (JSONL transaction log; audit / backup, never queried)
  3. Client + prebuild  fetch_compositions / prebuild   (the ONLY code that touches the network)
     Run-time fill      annotate_dataframe             (database only; never network)

Rules of record (Henry, 2026-09-20/21):
  - One composition -> one composition-level (linkage-undefined) GlyTouCan record in 2026.
  - A WURCS with no accession is a legitimate row (status 'wurcs_only'); never auto-register.
  - Run-time fill never overwrites a non-empty cell (MAS sheet value always wins).
  - WURCS / accession are annotations only. Nothing here may feed scoring or ML.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

__version__ = "0.1.0"
SCHEMA_VERSION = "1.0.0"
API_ENDPOINT = "https://api.glycosmos.org/glycancompositionconverter/1.0.0/composition2wurcs"
ATTRIBUTION = "GlyTouCan / GlyCosmos composition data, CC BY 4.0 (https://glytoucan.org, https://glycosmos.org)"
DEFAULT_DB_NAME = "glycomsp_reference"
DEFAULT_DB_VERSION = "v0.01.20260921"

ENV_DB = "GLYCOMSP_REFERENCE_DB"
ENV_LOG = "GLYCOMSP_GLYTOUCAN_LOG"

# ---------------------------------------------------------------------------
# 1. Canonicalisation and symbol mapping
# ---------------------------------------------------------------------------

# Same alternation and canonical order as the B-01 helper in mspfileloaderv14.run_prediction
# (KDN-first greedy; bare 'K' aliased to KDN). Case-sensitive: 's'/'p' are sulfate/phosphate.
_COMP_PAT = re.compile(r"(KDN|K|F|H|N|S|G|A|s|p)\s*([0-9]+)")
_CANON_ORDER = ("F", "H", "N", "S", "G", "KDN", "A", "s", "p")

# GlycoMSP symbol -> GlyCosmos composition2wurcs body key. Note the 'S' collision:
# GlycoMSP S = NeuAc (-> neu5ac); API "S" = sulfate (<- GlycoMSP 's').
SYMBOL_TO_API = {
    "F": "dhex",
    "H": "hex",
    "N": "hexnac",
    "S": "neu5ac",
    "G": "neu5gc",
    "KDN": "kdn",
    "A": "hexa",
    "s": "S",
    "p": "P",
}
API_QUERY_KEYS = ("hex", "hexnac", "dhex", "neu5ac", "neu5gc", "kdn", "hexa", "P", "S", "Ac")

_NON_COMPOSITION = {"", "non-glycan", "nonglycan", "non_glycan", "nan", "none"}


def canonical_composition(s) -> str:
    """Canonical GlycoMSP composition key ('F1H5N4S2'), or '' if `s` is not a composition.

    Returns '' for empty/NaN/'Non-glycan'. Strings whose characters are not fully consumed by
    the composition grammar (e.g. IUPAC names) also return ''.
    """
    if s is None:
        return ""
    if not isinstance(s, str):
        try:
            if s != s:  # NaN
                return ""
        except Exception:
            pass
        s = str(s)
    t = s.strip()
    if t.lower() in _NON_COMPOSITION:
        return ""
    # the whole string must be composition tokens (ignoring whitespace)
    if re.sub(_COMP_PAT, "", t).strip():
        return ""
    counts = {k: 0 for k in _CANON_ORDER}
    for k, v in _COMP_PAT.findall(t):
        if k == "K":
            k = "KDN"
        counts[k] += int(v)
    parts = [f"{k}{counts[k]}" for k in _CANON_ORDER if counts[k] > 0]
    return "".join(parts)


def composition_to_api_query(key: str) -> Dict[str, int]:
    """Canonical key -> composition2wurcs body object (ints; zero keys omitted)."""
    key = canonical_composition(key)
    if not key:
        raise ValueError(f"not a composition: {key!r}")
    q: Dict[str, int] = {}
    for k, v in _COMP_PAT.findall(key):
        if k == "K":
            k = "KDN"
        q[SYMBOL_TO_API[k]] = q.get(SYMBOL_TO_API[k], 0) + int(v)
    return q


def _api_body_item(q: Dict[str, int]) -> Dict[str, str]:
    """Full body item as the API documents it (all keys present, as strings)."""
    return {k: str(int(q.get(k, 0))) for k in API_QUERY_KEYS}


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# 2. Store
# ---------------------------------------------------------------------------

STATUS_RESOLVED = "resolved"       # wurcs + accession
STATUS_WURCS_ONLY = "wurcs_only"   # wurcs, accession not registered in GlyTouCan
STATUS_UNRESOLVABLE = "unresolvable"  # converter returned no wurcs

SOURCE_API = "api"
SOURCE_USER = "user"
SOURCE_SHIPPED = "shipped"


@dataclass
class Entry:
    composition_key: str
    query_json: str
    wurcs: Optional[str]
    accession: Optional[str]
    status: str
    source: str = SOURCE_API
    source_version: Optional[str] = None
    level: str = "composition"
    added_utc: str = ""
    last_checked_utc: str = ""
    note: Optional[str] = None

    @staticmethod
    def from_api(key: str, query: Dict[str, int], wurcs: Optional[str], accession: Optional[str],
                 *, now: Optional[str] = None) -> "Entry":
        """Build an Entry from one API item. `wurcs`/`accession` must be str or None; anything else
        (number, dict, list) raises TypeError so the client can mark the item malformed (Codex S1 F4)."""
        now = now or _utc_now()
        for name, val in (("wurcs", wurcs), ("id", accession)):
            if val is not None and not isinstance(val, str):
                raise TypeError(f"{name}: expected str or null, got {type(val).__name__}")
        wurcs = (wurcs or "").strip() or None
        accession = (accession or "").strip() or None
        if wurcs and not wurcs.startswith("WURCS="):
            raise TypeError(f"wurcs: unexpected value {wurcs[:20]!r}")
        if wurcs and accession:
            status = STATUS_RESOLVED
        elif wurcs:
            status = STATUS_WURCS_ONLY
        else:
            status = STATUS_UNRESOLVABLE
        return Entry(composition_key=key, query_json=json.dumps(query, sort_keys=True),
                     wurcs=wurcs, accession=accession, status=status, source=SOURCE_API,
                     added_utc=now, last_checked_utc=now)


def default_data_dir() -> Path:
    """Per-user data directory; same rule as the startup-error logger in mspfileloaderv14."""
    if os.name == "nt":
        return Path(os.environ.get("LOCALAPPDATA", os.path.expanduser("~"))) / "GlycoMSP"
    return Path(os.path.expanduser("~")) / ".glycomsp"


def default_db_path() -> Path:
    env = os.environ.get(ENV_DB)
    if env:
        return Path(env)
    return default_data_dir() / "references" / f"{DEFAULT_DB_NAME}.db"


def default_log_path() -> Path:
    env = os.environ.get(ENV_LOG)
    if env:
        return Path(env)
    return default_data_dir() / "glytoucan_cache.jsonl"


_DDL = f"""
CREATE TABLE IF NOT EXISTS meta (
  key   TEXT PRIMARY KEY,
  value TEXT
);
CREATE TABLE IF NOT EXISTS entries (
  composition_key  TEXT PRIMARY KEY,
  query_json       TEXT NOT NULL,
  wurcs            TEXT,
  accession        TEXT,
  level            TEXT NOT NULL DEFAULT 'composition',
  status           TEXT NOT NULL,
  source           TEXT NOT NULL,
  source_version   TEXT,
  added_utc        TEXT NOT NULL,
  last_checked_utc TEXT NOT NULL,
  note             TEXT
);
CREATE INDEX IF NOT EXISTS idx_entries_status ON entries(status);
"""


class ReferenceDB:
    """SQLite reference database. Opened lazily; created with schema on first use.

    Use as a context manager or call close(). Only successful conversions are stored
    (a network failure stores nothing, so the composition is retried next prebuild).
    """

    def __init__(self, path: Optional[os.PathLike | str] = None, *,
                 db_name: str = DEFAULT_DB_NAME, db_version: str = DEFAULT_DB_VERSION,
                 create: bool = True, readonly: bool = False):
        """create=True  : writable; file + schema + meta are created on first use (prebuild / manager).
        readonly=True   : SQLite `mode=ro` URI, no DDL, no meta writes, no mkdir — the run-time fill path.
                          A missing, empty or foreign file surfaces as an error at query time (callers map it
                          to status 'db_missing'); the file's bytes are never touched. (Codex S1 F5)
        create=False, readonly=False: writable but refuses to create a missing file."""
        self.path = Path(path) if path else default_db_path()
        self._conn: Optional[sqlite3.Connection] = None
        self._db_name = db_name
        self._db_version = db_version
        self._create = bool(create) and not readonly
        self._readonly = bool(readonly)

    # -- lifecycle --
    def _connect(self) -> sqlite3.Connection:
        if self._conn is None:
            if self._readonly:
                if not self.path.is_file():
                    raise FileNotFoundError(str(self.path))
                # Codex R2-F4: percent-encode via Path.as_uri() so '%', '?', '#', spaces and Unicode in the
                # filesystem path all address the SAME file that exists()/sha256() see.
                uri = self.path.resolve().as_uri() + "?mode=ro"
                conn = sqlite3.connect(uri, uri=True)
                conn.row_factory = sqlite3.Row
                # touch nothing; fail early with a clear error if this is not a reference DB
                if conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='entries'").fetchone() is None:
                    conn.close()
                    raise sqlite3.OperationalError(f"not a GlycoMSP reference database (no 'entries' table): {self.path}")
                self._conn = conn
                return conn
            if not self.path.exists() and not self._create:
                raise FileNotFoundError(str(self.path))
            self.path.parent.mkdir(parents=True, exist_ok=True)
            conn = sqlite3.connect(str(self.path))
            conn.row_factory = sqlite3.Row
            conn.executescript(_DDL)
            self._ensure_meta(conn)
            self._conn = conn
        return self._conn

    def _ensure_meta(self, conn: sqlite3.Connection) -> None:
        defaults = {
            "schema_version": SCHEMA_VERSION,
            "db_name": self._db_name,
            "db_version": self._db_version,
            "created_utc": _utc_now(),
            "attribution": ATTRIBUTION,
            "endpoint": API_ENDPOINT,
        }
        for k, v in defaults.items():
            conn.execute("INSERT OR IGNORE INTO meta(key, value) VALUES (?, ?)", (k, v))
        conn.commit()

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> "ReferenceDB":
        self._connect()
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # -- meta --
    def meta(self) -> Dict[str, str]:
        conn = self._connect()
        return {r["key"]: r["value"] for r in conn.execute("SELECT key, value FROM meta")}

    def set_meta(self, key: str, value: str) -> None:
        if self._readonly:
            raise sqlite3.OperationalError("read-only reference database")
        conn = self._connect()
        conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)", (key, value))
        conn.commit()

    def sha256(self) -> str:
        """Content hash of the database file (after flushing). '' if the file does not exist."""
        # Codex R2-F1: never commit here — an inspection must not make an unfinished write durable.
        if not self.path.exists():
            return ""
        h = hashlib.sha256()
        with open(self.path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    def identity(self) -> Dict[str, object]:
        m = self.meta()
        return {
            "name": m.get("db_name", ""),
            "version": m.get("db_version", ""),
            "schema_version": m.get("schema_version", ""),
            "path": str(self.path),
            "sha256": self.sha256(),
            "n_entries": self.count(),
            "updated_utc": m.get("updated_utc", m.get("created_utc", "")),   # Henry obs. c: same-day updates distinguishable
        }

    # -- queries --
    def count(self) -> int:
        conn = self._connect()
        return int(conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0])

    def summary(self) -> Dict[str, object]:
        conn = self._connect()
        by_status = {r[0]: int(r[1]) for r in conn.execute(
            "SELECT status, COUNT(*) FROM entries GROUP BY status")}
        by_source = {r[0]: int(r[1]) for r in conn.execute(
            "SELECT source, COUNT(*) FROM entries GROUP BY source")}
        return {"n_entries": self.count(), "by_status": by_status, "by_source": by_source,
                **{k: v for k, v in self.meta().items() if k in ("db_name", "db_version", "schema_version")}}

    def lookup(self, keys: Iterable[str]) -> Dict[str, Entry]:
        """Canonical keys -> Entry for rows that exist. Non-composition keys are skipped."""
        want = sorted({canonical_composition(k) for k in keys} - {""})
        if not want:
            return {}
        conn = self._connect()
        out: Dict[str, Entry] = {}
        # chunk to stay under SQLite's variable limit
        for i in range(0, len(want), 500):
            chunk = want[i:i + 500]
            q = f"SELECT * FROM entries WHERE composition_key IN ({','.join('?' * len(chunk))})"
            for r in conn.execute(q, chunk):
                out[r["composition_key"]] = Entry(**{k: r[k] for k in r.keys()})
        return out

    def missing(self, keys: Iterable[str]) -> List[str]:
        """Canonical keys with no row at all (any status counts as present)."""
        want = sorted({canonical_composition(k) for k in keys} - {""})
        have = self.lookup(want)
        return [k for k in want if k not in have]

    def get(self, key: str) -> Optional[Entry]:
        return self.lookup([key]).get(canonical_composition(key))

    # -- writes --
    def upsert(self, entries: Iterable[Entry]) -> int:
        """Write entries as ONE transaction. On any failure the transaction is rolled back and the
        exception propagates: either every row (plus `updated_utc`) is durable, or none is (Codex R2-F1)."""
        if self._readonly:
            raise sqlite3.OperationalError("read-only reference database")
        conn = self._connect()
        n = 0
        try:
            if not conn.in_transaction:
                conn.execute("BEGIN")
            for e in entries:
                if not e.composition_key:
                    continue
                conn.execute(
                    """INSERT INTO entries(composition_key, query_json, wurcs, accession, level, status,
                                           source, source_version, added_utc, last_checked_utc, note)
                       VALUES (?,?,?,?,?,?,?,?,?,?,?)
                       ON CONFLICT(composition_key) DO UPDATE SET
                           query_json=excluded.query_json, wurcs=excluded.wurcs,
                           accession=excluded.accession, level=excluded.level, status=excluded.status,
                           source=excluded.source, source_version=excluded.source_version,
                           last_checked_utc=excluded.last_checked_utc, note=excluded.note""",
                    (e.composition_key, e.query_json, e.wurcs, e.accession, e.level, e.status,
                     e.source, e.source_version, e.added_utc or _utc_now(),
                     e.last_checked_utc or _utc_now(), e.note))
                n += 1
            if n:
                conn.execute("INSERT OR REPLACE INTO meta(key, value) VALUES ('updated_utc', ?)", (_utc_now(),))
            conn.commit()
        except Exception:
            try:
                conn.rollback()
            except Exception:
                pass
            raise
        return n

    # -- log reconciliation (Codex R2-F2, light option) --
    def log_records(self) -> List[Dict[str, object]]:
        """Every entry rendered as a transaction-log record (the database is a superset of the log)."""
        conn = self._connect()
        m = self.meta()
        out: List[Dict[str, object]] = []
        for r in conn.execute("SELECT * FROM entries ORDER BY added_utc, composition_key"):
            src = r["source"] or SOURCE_API
            try:
                q = json.loads(r["query_json"] or "{}")
            except Exception:
                q = {}
            out.append({
                "ts": r["last_checked_utc"] or r["added_utc"], "op": "api_resolve" if src == SOURCE_API else f"{src}_entry",
                "composition_key": r["composition_key"], "query": q,
                "accession": r["accession"] or "", "wurcs": r["wurcs"] or "", "status": r["status"],
                "endpoint": m.get("endpoint", API_ENDPOINT),
                "db_name": m.get("db_name", ""), "db_version": m.get("db_version", ""),
                "rebuilt": True,
            })
        return out


class LogWriteError(OSError):
    """Raised by append_log when the JSONL could not be (fully) written. `.written` = lines that made it."""
    def __init__(self, msg: str, written: int):
        super().__init__(msg)
        self.written = written


def append_log(log_path: Optional[os.PathLike | str], records: Iterable[Dict[str, object]]) -> int:
    """Append JSON lines to the transaction log. Returns lines written.
    Codex R2-F2: raises LogWriteError (carrying the partial count) instead of swallowing I/O errors, so callers
    can report a degraded status and the CLI can rebuild the log from the database (`rebuild_log`)."""
    p = Path(log_path) if log_path else default_log_path()
    n = 0
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "a", encoding="utf-8") as f:
            for rec in records:
                f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True) + "\n")
                f.flush()
                n += 1
    except Exception as ex:
        raise LogWriteError(f"could not append to log {p}: {type(ex).__name__}: {ex}", n) from ex
    return n


class RebuildLogError(OSError):
    """rebuild_log refused or failed; the live log and its backup are untouched unless stated."""


def _same_file(a: Path, b: Path) -> bool:
    """True if the two paths name the same existing filesystem object (follows symlinks, catches hard links)."""
    try:
        return os.path.samefile(a, b)
    except (FileNotFoundError, OSError):
        return False


def rebuild_log(db: "ReferenceDB", log_path: Optional[os.PathLike | str] = None) -> int:
    """Regenerate the whole JSONL from the database (recovery for a lost/partial log).

    Safety contract (Codex R3-F1/F2 + R4-F1/F2):
      1. The source is ALWAYS read through a fresh read-only connection on `db.path`, whatever handle the caller
         passed (writable / unconnected handles are never used). A missing, empty or foreign file raises
         RebuildLogError before anything is written; nothing is created or initialised by recovery.
      2. Destination aliasing is refused before any mutation: the database must not be the same file as the live
         log, its backup or the temp file (real-file identity, not string equality); a live log or backup that is a
         symlink is refused (a copy would write through it).
      3. The projection is fully written to an exclusively created unique temp file next to the log.
      4. If a live log exists it is COPIED to `<log>.bak` (the live file is never moved away); a backup failure
         raises RebuildLogError and leaves the live log untouched.
      5. The live target is then replaced atomically by `os.replace(tmp, log)`; that is the only mutation of it.
    Returns the number of lines written."""
    import shutil
    import tempfile
    p = Path(log_path) if log_path else default_log_path()
    src_path = Path(db.path)
    bak = p.with_suffix(p.suffix + ".bak")
    # 1. read the source through our OWN read-only handle (R4-F1): never the caller's, never writable
    if not src_path.is_file():
        raise RebuildLogError(f"reference database does not exist: {src_path} (recovery never creates one)")
    try:
        ro = ReferenceDB(src_path, readonly=True)
        try:
            recs = ro.log_records()
        finally:
            ro.close()
    except Exception as ex:
        raise RebuildLogError(f"cannot read reference database {src_path}: {type(ex).__name__}: {ex}") from ex
    # 2. refuse destination aliasing (R4-F2) before any mutation
    for label, dest in (("live log", p), ("backup", bak)):
        if _same_file(src_path, dest):
            raise RebuildLogError(f"refused: the {label} {dest} is the same file as the reference database {src_path}")
        if dest.is_symlink():
            raise RebuildLogError(f"refused: the {label} {dest} is a symlink (a copy/replace would write through it)")
    if p.exists() and bak.exists() and _same_file(p, bak):
        raise RebuildLogError(f"refused: live log {p} and backup {bak} are the same file")
    try:
        if p.parent.exists() and _same_file(src_path.parent, p.parent) and src_path.name.startswith(p.name + "."):
            # a database named like one of our sibling temp/backup files would be a foot-gun; refuse
            if src_path.suffix in (".tmp", ".bak"):
                raise RebuildLogError(f"refused: reference database {src_path} collides with the log's sibling files")
    except RebuildLogError:
        raise
    except Exception:
        pass
    # 3. complete projection into an EXCLUSIVELY created unique temp file
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd, tmp_name = tempfile.mkstemp(prefix=p.name + ".", suffix=".tmp", dir=str(p.parent))
        tmp = Path(tmp_name)
    except Exception as ex:
        raise RebuildLogError(f"could not create a temp file next to {p}: {type(ex).__name__}: {ex}") from ex
    try:
        if _same_file(src_path, tmp):
            raise RebuildLogError(f"refused: temp file {tmp} aliases the reference database")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            for rec in recs:
                f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True) + "\n")
            f.flush()
            os.fsync(f.fileno())
    except Exception as ex:
        try:
            tmp.unlink()
        except Exception:
            pass
        if isinstance(ex, RebuildLogError):
            raise
        raise RebuildLogError(f"could not write {tmp}: {type(ex).__name__}: {ex}") from ex
    # 4. backup by COPY (live file stays in place); abort if the backup cannot be established
    if p.exists():
        try:
            if bak.is_dir():
                raise IsADirectoryError(str(bak))
            shutil.copy2(p, bak)
            if _same_file(bak, src_path):     # belt and braces: the copy must not have landed on the source
                raise RebuildLogError(f"backup {bak} aliases the reference database")
        except Exception as ex:
            try:
                tmp.unlink()
            except Exception:
                pass
            raise RebuildLogError(f"could not back up {p} to {bak}: {type(ex).__name__}: {ex} — live log untouched") from ex
    # 5. atomic replace of the live target
    try:
        os.replace(tmp, p)
    except Exception as ex:
        try:
            tmp.unlink()
        except Exception:
            pass
        raise RebuildLogError(f"could not replace {p}: {type(ex).__name__}: {ex} — live log untouched, backup kept") from ex
    return len(recs)


# ---------------------------------------------------------------------------
# 3a. Client (network) — used ONLY by prebuild / update actions
# ---------------------------------------------------------------------------

HttpPost = Callable[[str, bytes, float], bytes]


def _ssl_context():
    """Default verifying context; prefers certifi's bundle when installed (python.org builds on
    macOS ship without system CA access until 'Install Certificates.command' is run).
    Verification is never disabled."""
    import ssl
    try:
        import certifi  # optional; not in requirements.txt today (see handoff note)
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


def _default_http_post(url: str, body: bytes, timeout_s: float) -> bytes:
    req = urllib.request.Request(
        url, data=body, method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": f"GlycoMSP-glytoucan-resolver/{__version__}"})
    with urllib.request.urlopen(req, timeout=timeout_s, context=_ssl_context()) as resp:
        return resp.read()


@dataclass
class FetchResult:
    entries: List[Entry] = field(default_factory=list)   # parsed OK (may include rows the callback failed to persist)
    status: str = "ok"            # 'ok' | 'partial' | 'offline' | 'disabled' | 'persist_failed'
    api_calls: int = 0
    n_requested: int = 0
    errors: List[str] = field(default_factory=list)
    unsent: List[str] = field(default_factory=list)       # keys never sent (after failure / invalid)
    unpersisted: List[str] = field(default_factory=list)  # parsed but on_batch failed (Codex R2-F1); retry next run


def fetch_compositions(keys: Iterable[str], *, batch_size: int = 30, pause_s: float = 10.0,
                       timeout_s: float = 20.0, http_post: Optional[HttpPost] = None,
                       endpoint: str = API_ENDPOINT, progress: Optional[Callable[[str], None]] = None,
                       sleep: Callable[[float], None] = time.sleep,
                       on_batch: Optional[Callable[[List["Entry"]], None]] = None) -> FetchResult:
    """POST canonical compositions to the converter in batches; fail soft.

    - Batches of `batch_size` (<=30 per the courtesy rule) with `pause_s` between batches.
    - Any transport/HTTP/parse failure ends the run: entries already obtained are kept,
      remaining keys are reported in `unsent`, status becomes 'partial' (or 'offline' if the
      very first batch failed).
    - Response items are matched to inputs by position (the API echoes counts; we verify the
      echoed counts agree with what we sent and skip mismatches).
    - `on_batch(entries)` is called after each successfully parsed batch so callers can persist
      incrementally; a later batch failure then cannot lose earlier results (Codex S1 F4).
    """
    post = http_post or _default_http_post
    say = progress or (lambda s: None)
    res = FetchResult()
    want: List[str] = []
    for k in keys:
        ck = canonical_composition(k)
        if ck and ck not in want:
            want.append(ck)
    res.n_requested = len(want)
    if not want:
        return res
    batch_size = max(1, min(int(batch_size), 30))
    batches = [want[i:i + batch_size] for i in range(0, len(want), batch_size)]
    for bi, batch in enumerate(batches):
        if bi > 0 and pause_s > 0:
            say(f"pausing {pause_s:.0f}s before batch {bi + 1}/{len(batches)}")
            sleep(pause_s)
        queries = [composition_to_api_query(k) for k in batch]
        body = json.dumps([_api_body_item(q) for q in queries]).encode("utf-8")
        say(f"batch {bi + 1}/{len(batches)}: {len(batch)} compositions")
        try:
            raw = post(endpoint, body, timeout_s)
            res.api_calls += 1
            payload = json.loads(raw.decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError) as ex:
            res.errors.append(f"batch {bi + 1}: {type(ex).__name__}: {ex}")
            res.unsent.extend(k for b in batches[bi:] for k in b)
            res.status = "offline" if not res.entries else "partial"
            say(f"batch {bi + 1} failed: {ex} — stopping (fail soft)")
            return res
        data = payload.get("data") if isinstance(payload, dict) else payload
        if not isinstance(data, list) or len(data) != len(batch):
            res.errors.append(f"batch {bi + 1}: unexpected response shape "
                              f"(got {type(data).__name__} len={len(data) if isinstance(data, list) else 'n/a'}, want {len(batch)})")
            res.unsent.extend(k for b in batches[bi:] for k in b)
            res.status = "offline" if not res.entries else "partial"
            return res
        now = _utc_now()
        batch_entries: List[Entry] = []
        for key, q, item in zip(batch, queries, data):
            # Codex S1 F4: every per-item failure (shape, echo mismatch, bad field type) marks THAT item
            # unsent + an error line and never escapes; earlier items and batches are kept.
            try:
                if not isinstance(item, dict):
                    raise TypeError(f"item is {type(item).__name__}, expected object")
                # verify the echoed counts match what we sent (zero keys are dropped in the echo)
                echoed = {k: int(v) for k, v in item.items()
                          if k in API_QUERY_KEYS and str(v).strip().lstrip("-").isdigit()}
                if {k: v for k, v in echoed.items() if v} != {k: v for k, v in q.items() if v}:
                    raise ValueError(f"echoed composition {echoed} != sent {q}")
                batch_entries.append(Entry.from_api(key, q, item.get("wurcs"), item.get("id"), now=now))
            except Exception as ex:
                res.errors.append(f"{key}: malformed response item ({type(ex).__name__}: {ex}); skipped")
                res.unsent.append(key)
        res.entries.extend(batch_entries)
        if on_batch is not None and batch_entries:
            try:
                on_batch(batch_entries)      # prebuild persists per batch (Codex S1 F4)
            except Exception as ex:
                # Codex R2-F1: storage failure is NOT fail-soft-and-continue — record the keys, stop requesting.
                res.errors.append(f"batch {bi + 1}: persistence failed: {type(ex).__name__}: {ex}")
                res.unpersisted.extend(e.composition_key for e in batch_entries)
                res.unsent.extend(k for b in batches[bi + 1:] for k in b)
                res.status = "persist_failed"
                say(f"batch {bi + 1}: persistence failed: {ex} — stopping")
                return res
    if res.unsent and res.status == "ok":
        res.status = "partial"
    return res


# ---------------------------------------------------------------------------
# 3b. Prebuild: diff against the database -> fetch misses -> upsert -> log
# ---------------------------------------------------------------------------

@dataclass
class PrebuildResult:
    n_source: int = 0            # distinct valid compositions in the source
    n_already: int = 0           # present in the database before this run
    n_fetched: int = 0           # rows COMMITTED this run (not merely parsed)
    api_calls: int = 0
    status: str = "ok"           # 'ok' | 'partial' | 'offline' | 'disabled' | 'nothing_to_do'
                                 # | 'persist_failed' (DB write failed; Codex R2-F1)
                                 # | 'log_failed' (DB ok, JSONL not fully written; Codex R2-F2 → run --rebuild-log)
    errors: List[str] = field(default_factory=list)
    unsent: List[str] = field(default_factory=list)
    unpersisted: List[str] = field(default_factory=list)   # parsed, not committed → retried next run
    unlogged: List[str] = field(default_factory=list)      # committed, not logged → rebuild_log recovers
    db: Dict[str, object] = field(default_factory=dict)

    @property
    def storage_ok(self) -> bool:
        return self.status not in ("persist_failed", "log_failed")


def prebuild(keys: Iterable[str], db: ReferenceDB, log_path: Optional[os.PathLike | str] = None,
             *, allow_network: bool = True, progress: Optional[Callable[[str], None]] = None,
             **fetch_kw) -> PrebuildResult:
    say = progress or (lambda s: None)
    res = PrebuildResult()
    want = sorted({canonical_composition(k) for k in keys} - {""})
    res.n_source = len(want)
    miss = db.missing(want)
    res.n_already = len(want) - len(miss)
    say(f"{len(want)} compositions in source; {res.n_already} already in database; {len(miss)} to fetch")
    def _identity_into(res_: PrebuildResult) -> None:
        # Codex R2-F1 / R3-F5: inspection (hash/meta/count) must never turn a finished run into a crash — on ANY path.
        try:
            res_.db = db.identity()
        except Exception as ex:
            res_.errors.append(f"identity() failed: {type(ex).__name__}: {ex}")
            res_.db = {"path": str(db.path)}

    if not miss:
        res.status = "nothing_to_do"
        _identity_into(res)
        return res
    if not allow_network:
        res.status = "disabled"
        res.unsent = miss
        _identity_into(res)
        return res
    m = db.meta()
    endpoint = fetch_kw.get("endpoint", API_ENDPOINT)

    def _persist(entries: List[Entry]) -> None:
        # per-batch persistence (Codex S1 F4): upsert (one transaction) then log, so a later failure keeps
        # earlier rows. A DB failure propagates (R2-F1); a log failure is recorded per key and the run
        # continues because the database is the source of truth and the log is rebuildable (R2-F2).
        res.n_fetched += db.upsert(entries)
        recs = [{
            "ts": e.last_checked_utc, "op": "api_resolve", "composition_key": e.composition_key,
            "query": json.loads(e.query_json), "accession": e.accession or "", "wurcs": e.wurcs or "",
            "status": e.status, "endpoint": endpoint,
            "db_name": m.get("db_name", ""), "db_version": m.get("db_version", ""),
        } for e in entries]
        try:
            append_log(log_path, recs)
        except LogWriteError as lex:
            res.unlogged.extend(e.composition_key for e in entries[lex.written:])
            res.errors.append(str(lex))
            say(f"log write failed ({lex.written}/{len(entries)} lines) — database committed; run --rebuild-log later")

    fr = fetch_compositions(miss, progress=say, on_batch=_persist, **fetch_kw)
    res.api_calls = fr.api_calls
    res.errors = fr.errors + [e for e in res.errors if e not in fr.errors]
    res.unsent = fr.unsent
    res.unpersisted = fr.unpersisted
    res.status = fr.status
    if res.unlogged and res.status in ("ok", "partial", "offline"):
        res.status = "log_failed"
    _identity_into(res)
    say(f"done: {res.n_fetched} rows committed, {res.api_calls} API calls, status={res.status}")
    return res


# ---------------------------------------------------------------------------
# 3c. Run-time fill — database only, never network
# ---------------------------------------------------------------------------

@dataclass
class FillResult:
    status: str = "db_missing"   # 'filled' | 'partial' | 'none' | 'db_missing'
    route: str = "database"
    db: Dict[str, object] = field(default_factory=dict)
    n_compositions: int = 0      # distinct valid compositions seen in comp_col
    n_filled: int = 0            # distinct compositions that had a row (any wurcs)
    n_missing: int = 0           # distinct compositions with no row
    n_cells_written: int = 0
    n_user_accession: int = 0    # rows whose GlyToucan ID differs from the composition-level record (WURCS left as-is)
    n_user_wurcs: int = 0        # rows whose WURCS differs from the record while the ID is blank (ID left blank; D1)
    last_resolved_utc: str = ""

    def to_artifact_block(self) -> Dict[str, object]:
        """Shape written to method JSON artifacts.glytoucan_resolution (see spec §4)."""
        return {
            "status": self.status,
            "route": self.route,
            "db": {k: self.db.get(k, "") for k in ("name", "version", "path", "sha256", "n_entries", "updated_utc")},
            "n_compositions": self.n_compositions,
            "n_filled": self.n_filled,
            "n_missing": self.n_missing,
            "n_user_accession": self.n_user_accession,
            "n_user_wurcs": self.n_user_wurcs,
            "last_resolved_utc": self.last_resolved_utc,
        }


def _is_empty_cell(v) -> bool:
    if v is None:
        return True
    try:
        if v != v:  # NaN
            return True
    except Exception:
        pass
    return str(v).strip() == ""


def annotate_dataframe(df, comp_col: str, db: Optional[ReferenceDB], *,
                       id_col: str = "GlyToucan ID", wurcs_col: str = "WURCS") -> FillResult:
    """Fill EMPTY `id_col` / `wurcs_col` cells from the reference database, in place.

    - Never overwrites a non-empty cell (a MAS sheet value always wins).
    - A row whose `id_col` already holds an accession DIFFERENT from the composition-level record is a
      structure-specific user assignment: its WURCS cell is left exactly as-is (empty or not), because the
      composition-level WURCS would not describe that structure. Counted in `n_user_accession`.
      (Henry manual test 2026-09-22, obs. d / run 9.) Resolving those needs the accession→WURCS client (later).
    - Mirror rule (Codex R2 D1, Henry approved 2026-09-23): a row whose `wurcs_col` already holds a WURCS DIFFERENT
      from the composition record while `id_col` is blank is a structure-specific user WURCS: the accession is left
      blank (the generic accession would not name that structure). Counted in `n_user_wurcs`; T1-h resolves it later.
    - Never touches `comp_col`.
    - Creates `id_col` / `wurcs_col` if absent (empty strings), matching the B-17 seeding.
    - Database only. If `db` is None or unreadable, returns status 'db_missing' with the
      frame unchanged apart from column creation.
    """
    res = FillResult(last_resolved_utc=_utc_now())
    for c in (id_col, wurcs_col):
        if c not in df.columns:
            df[c] = ""
    if comp_col not in df.columns:
        res.status = "none"
        return res
    keys = {canonical_composition(v) for v in df[comp_col].tolist()} - {""}
    res.n_compositions = len(keys)
    if db is None:
        return res
    try:
        found = db.lookup(keys)
        res.db = db.identity()
    except Exception as ex:
        print(f"[glytoucan][WARN] reference database unreadable: {ex}", file=sys.stderr)
        return res
    res.n_filled = sum(1 for k in keys if k in found and found[k].wurcs)
    res.n_missing = res.n_compositions - res.n_filled
    if res.n_compositions == 0:
        res.status = "none"
        return res
    ids = df[id_col].tolist()
    wur = df[wurcs_col].tolist()
    comps = df[comp_col].tolist()
    n_written = 0
    n_user = 0
    n_user_w = 0
    for i, c in enumerate(comps):
        e = found.get(canonical_composition(c))
        if e is None:
            continue
        id_empty = _is_empty_cell(ids[i])
        if not id_empty and str(ids[i]).strip() != (e.accession or ""):
            n_user += 1              # user-assigned accession: keep the row's WURCS untouched
            continue
        wur_empty = _is_empty_cell(wur[i])
        if id_empty and not wur_empty and str(wur[i]).strip() != (e.wurcs or ""):
            n_user_w += 1            # D1: user-assigned WURCS: leave the accession blank for the reverse resolver
            continue
        if e.accession and id_empty:
            ids[i] = e.accession
            n_written += 1
        if e.wurcs and _is_empty_cell(wur[i]):
            wur[i] = e.wurcs
            n_written += 1
    res.n_user_accession = n_user
    res.n_user_wurcs = n_user_w
    if n_written:
        df[id_col] = ids
        df[wurcs_col] = wur
    res.n_cells_written = n_written
    res.status = "filled" if res.n_missing == 0 else "partial"
    return res


# ---------------------------------------------------------------------------
# Convenience: composition sources (used by the CLI and, later, the manager panel)
# ---------------------------------------------------------------------------

_SOURCE_COLUMNS = ("composition", "Structure", "pseudo compositions", "Composition")

# Native compnewv4 library schema (save_glycan_pseudocomp_to_csv): residue COUNT columns, no composition string.
# Column -> GlycoMSP symbol; canonical order is applied by canonical_composition(). (Codex S1 F3)
_COUNT_COLUMNS = (("Fuc", "F"), ("Hex", "H"), ("HexNAc", "N"), ("NeuAc", "S"), ("NeuGc", "G"),
                  ("KDN", "KDN"), ("HexA", "A"), ("SO3", "s"), ("PO3H", "p"))
_COUNT_REQUIRED = ("Hex", "HexNAc")


def compositions_from_count_frame(df) -> List[str]:
    """Distinct canonical keys from a residue-count table (compnewv4 in-silico library). Columns other
    than Hex/HexNAc are optional and treated as 0 when absent. Non-numeric cells are ignored."""
    cols = [(c, sym) for c, sym in _COUNT_COLUMNS if c in df.columns]
    if not all(c in df.columns for c in _COUNT_REQUIRED):
        raise ValueError("not a residue-count table (needs at least Hex and HexNAc columns)")
    keys = set()
    for row in df[[c for c, _ in cols]].itertuples(index=False, name=None):
        parts = []
        for (c, sym), v in zip(cols, row):
            try:
                n = int(float(v))
            except (TypeError, ValueError):
                n = 0
            if n > 0:
                parts.append(f"{sym}{n}")
        k = canonical_composition("".join(parts))
        if k:
            keys.add(k)
    return sorted(keys)


def compositions_from_file(path: os.PathLike | str, column: Optional[str] = None, *, sheet: Optional[str] = None) -> List[str]:
    """Distinct canonical compositions from a CGA TSV / trainable CSV / MAS sheet /
    in-silico library. Column auto-detected from `_SOURCE_COLUMNS` unless given."""
    import pandas as pd  # local import: keep module importable without pandas for the DB layer
    p = Path(path)
    suf = p.suffix.lower()
    if suf in (".xlsx", ".xlsm", ".xls"):
        # MAS workbooks keep annotations on the 'MSlist' sheet (mspvalidator_merger.directassign_files);
        # fall back to the first sheet for other workbooks. `sheet` overrides.
        xl = pd.ExcelFile(p)
        sheet_name = sheet if sheet is not None else ("MSlist" if "MSlist" in xl.sheet_names else xl.sheet_names[0])
        df = xl.parse(sheet_name, dtype=str)
    elif suf == ".tsv":
        df = pd.read_csv(p, sep="\t", dtype=str)
    else:
        # sniff: CGA TSVs may carry .csv; try tab if a comma read yields one column
        df = pd.read_csv(p, dtype=str)
        if df.shape[1] == 1:
            df = pd.read_csv(p, sep="\t", dtype=str)
    col = column
    if col is None:
        for c in _SOURCE_COLUMNS:
            if c in df.columns:
                col = c
                break
    if col is None and all(c in df.columns for c in _COUNT_REQUIRED):
        return compositions_from_count_frame(df)          # native in-silico library (Codex S1 F3)
    if col is None or col not in df.columns:
        raise ValueError(f"no composition column found in {p.name}; columns: {list(df.columns)[:15]}")
    keys = sorted({canonical_composition(v) for v in df[col].tolist()} - {""})
    return keys
