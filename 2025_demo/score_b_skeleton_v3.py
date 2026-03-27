from __future__ import annotations
#202603231826(GMT+9) - Latest timestamp
#202603231800(GMT+9) 
#20260322 validation finished
#20260322 biological intents check 
#20260323 adjust the score - start from v3
#20260323-24 binding current workflow
#20260324-26 add into GUI
#will move the code to msp_CGA_structscore.py when we start implementing CGA score B

"""
Score B skeleton for GlycoMSP CGA motif-based scoring.

Purpose
-------
This file is a structured starting point for implementing the Excel-driven
Score B engine discussed in the design session.

Current scope
-------------
Phase 1
    - workbook loading
    - sheet/column validation
    - dataclass construction
    - TargetHierarchy compilation

Phase 2
    - Pass-1 evaluation (direct evidence)

Phase 3
    - Pass-2 evaluation (gated / derived evidence)

Placeholders are intentionally kept for:
    - MotifPolicy handling
    - CompositionConsistency penalties
    - UnexpectedEvidencePenalty
    - candidate-level final Score B aggregation

Notes
-----
1. parent_target_id is treated as a context gate, not hierarchy inheritance.
2. TargetHierarchy is only used for group expansion / gate resolution / penalties.
3. interpretation_v2 contains positive support rules only.
4. Non-MS3 targets should be explicitly defined in Excel; no implicit inheritance.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class ScoreBError(Exception):
    """Base exception for Score B loader / compiler / engine errors."""


class WorkbookValidationError(ScoreBError):
    """Raised when the workbook structure or values are invalid."""


class HierarchyValidationError(ScoreBError):
    """Raised when TargetHierarchy contains invalid references or cycles."""


# ---------------------------------------------------------------------------
# Dataclasses: raw workbook objects
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IonDef:
    row_idx: int
    glycotope_name: str
    ion_type: str
    fragmentation_mass: float
    structural_identifier: str
    charge_mode: str
    derivatization: str
    notes: Optional[str] = None


@dataclass(frozen=True)
class InterpretationRule:
    row_idx: int
    display_name: str
    target_id: str
    target_level: str
    target_group: Optional[str]
    parent_target_id: Optional[str]
    ion_struct_id: str
    ion_struct_target_id: Optional[str]
    ion_role: Optional[str]
    logic: str
    min_hits: int
    weight: float
    charge_mode_required: str
    derivatization_required: str
    notes: Optional[str] = None


@dataclass(frozen=True)
class HierarchyEdge:
    row_idx: int
    parent_target_id: str
    child_target_id: str
    notes: Optional[str] = None


@dataclass(frozen=True)
class MotifPolicyRule:
    row_idx: int
    flag_id: str
    target_id: str
    activation_mode: str
    include_descendants: bool
    notes: Optional[str] = None


@dataclass(frozen=True)
class CompositionConsistencyRule:
    row_idx: int
    rule_id: str
    composition_field: str
    when_value: str
    evidence_target_id: str
    penalty: float
    mode: str = "auto"
    notes: Optional[str] = None


@dataclass(frozen=True)
class UnexpectedEvidencePenaltyRule:
    row_idx: int
    flag_id: str
    target_id: str
    include_descendants: bool
    penalty: float
    notes: Optional[str] = None


# ---------------------------------------------------------------------------
# Dataclasses: compiled / runtime objects
# ---------------------------------------------------------------------------


@dataclass
class HierarchyIndex:
    children_by_parent: dict[str, set[str]] = field(default_factory=dict)
    descendants_by_target: dict[str, set[str]] = field(default_factory=dict)
    concrete_targets: set[str] = field(default_factory=set)
    group_nodes: set[str] = field(default_factory=set)

    def expand_target(self, target_id: str) -> set[str]:
        """
        Expand a target/group to concrete target_ids.

        If target_id is already a concrete target, it returns {target_id}.
        If target_id is a group node, it returns all concrete descendants.
        """
        if target_id in self.concrete_targets:
            return {target_id}
        return set(self.descendants_by_target.get(target_id, set()))


@dataclass
class ScoreBConfig:
    workbook_path: Path
    ions_by_struct_id: dict[str, IonDef]
    interpretation_rules: list[InterpretationRule]
    pass1_rules: list[InterpretationRule]
    pass2_rules: list[InterpretationRule]
    hierarchy_edges: list[HierarchyEdge]
    hierarchy_index: HierarchyIndex
    motif_policy_rules: list[MotifPolicyRule]
    composition_rules: list[CompositionConsistencyRule]
    unexpected_penalty_rules: list[UnexpectedEvidencePenaltyRule]
    rules_by_target_id: dict[str, list[InterpretationRule]] = field(default_factory=dict)
    rules_by_ion_struct_id: dict[str, list[InterpretationRule]] = field(default_factory=dict)


@dataclass(frozen=True)
class ObservedIonHit:
    ion_struct_id: str
    mz: float
    intensity: float
    ppm_error: float
    scan_level: Optional[str] = None
    matched: bool = True


@dataclass
class SpectrumEvidence:
    observed_hits: list[ObservedIonHit]
    charge_mode: str
    derivatization: str

    def hit_count_by_ion_struct_id(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for hit in self.observed_hits:
            if not hit.matched:
                continue
            counts[hit.ion_struct_id] = counts.get(hit.ion_struct_id, 0) + 1
        return counts


@dataclass
class RuleEvaluation:
    rule: InterpretationRule
    gate_open: Optional[bool]
    charge_ok: bool
    derivatization_ok: bool
    hit_count: int
    supported: bool
    row_score: float
    status: str


@dataclass
class TargetEvidence:
    target_id: str
    pass1_score: float = 0.0
    pass2_score: float = 0.0
    pass1_max: float = 0.0
    pass2_max: float = 0.0
    supported: bool = False
    contributing_rule_rows: list[int] = field(default_factory=list)

    # Runtime rows that were actually evaluated in this run.
    # Pass 1 rows are always active; Pass 2 rows are active only if gate opens.
    active_rule_rows: list[int] = field(default_factory=list)

    # Runtime rows that positively contributed score in this run.
    supported_rule_rows: list[int] = field(default_factory=list)

    # Phase 5A runtime summary fields
    logic_mode: str = "OR"
    partial_score: float = 0.0
    logic_satisfied: bool = False
    active_rule_count: int = 0
    supported_rule_count: int = 0

    @property
    def total_raw_score(self) -> float:
        return self.pass1_score + self.pass2_score

    @property
    def total_max_score(self) -> float:
        return self.pass1_max + self.pass2_max

    @property
    def normalized_score(self) -> float:
        if self.total_max_score <= 0:
            return 0.0
        return self.total_raw_score / self.total_max_score


@dataclass
class GlobalEvidenceResult:
    target_evidence: dict[str, TargetEvidence]
    pass1_rule_evaluations: list[RuleEvaluation]
    pass2_rule_evaluations: list[RuleEvaluation]


@dataclass
class ScoreBResult:
    candidate_composition: tuple[int, int, int, int, int, int]
    motif_support_score: float
    composition_penalty: float
    unexpected_penalty: float
    final_score_b: float
    target_evidence: dict[str, TargetEvidence]
    applied_comp_rules: list[dict[str, Any]] = field(default_factory=list)
    applied_unexpected_rules: list[dict[str, Any]] = field(default_factory=list)

    # Phase 5A debug payloads
    candidate_mask_trace: Optional[CandidateTargetMaskTrace] = None
    scoring_units: list[ScoringUnit] = field(default_factory=list)
    aggregation_debug: dict[str, float] = field(default_factory=dict)

@dataclass(frozen=True)
class MotifPolicyMatch:
    row_idx: int
    flag_id: str
    target_id: str
    activation_mode: str
    include_descendants: bool
    flag_value: Optional[bool]
    is_active: bool
    expanded_targets: list[str]


@dataclass
class CandidateTargetMaskTrace:
    input_flags: dict[str, bool]
    matched_policy_rows: list[MotifPolicyMatch] = field(default_factory=list)
    selected_root_targets: list[str] = field(default_factory=list)
    selected_concrete_targets: list[str] = field(default_factory=list)
    unmatched_flags: list[str] = field(default_factory=list)


@dataclass
class ScoringUnit:
    unit_id: str
    source_target_ids: list[str]
    score: float
    logic_satisfied: bool
    weight: float = 1.0
    unit_kind: str = "leaf"


# ---------------------------------------------------------------------------
# Sheet schema definitions
# ---------------------------------------------------------------------------


REQUIRED_SHEETS: tuple[str, ...] = (
    "ionlist",
    "interpretation_v2",
    "TargetHierarchy",
    "MotifPolicy",
    "CompositionConsistency",
    "UnexpectedEvidencePenalty",
)

# NOTE:
# These are normalized internal names. Adjust aliases below if workbook column
# headers differ slightly from this idealized schema.
REQUIRED_COLUMNS: dict[str, set[str]] = {
    "ionlist": {
        "glycotope_name",
        "ion_type",
        "fragmentation_mass",
        "structural_identifier",
        "charge_mode",
        "derivatization",
    },
    "interpretation_v2": {
        "display_name",
        "target_id",
        "target_level",
        "target_group",
        "parent_target_id",
        "ion_struct_id",
        "ion_struct_target_id",
        "ion_role",
        "logic",
        "min_hits",
        "weight",
        "charge_mode_required",
        "derivatization_required",
    },
    "TargetHierarchy": {
        "parent_target_id",
        "child_target_id",
    },
    "MotifPolicy": {
        "flag_id",
        "target_id",
        "activation_mode",
        "include_descendants",
    },
    "CompositionConsistency": {
        "rule_id",
        "composition_field",
        "when_value",
        "evidence_target_id",
        "penalty",
    },
    "UnexpectedEvidencePenalty": {
        "flag_id",
        "target_id",
        "include_descendants",
        "penalty",
    },
}

# Optional alias mapping: Excel header -> internal normalized header
COLUMN_ALIASES: dict[str, str] = {
    "Structural identifiers": "structural_identifier",
    "Structural identifier": "structural_identifier",
    "structural identifiers": "structural_identifier",
    "Glycotope name": "glycotope_name",
    "Fragmentation mass": "fragmentation_mass",
    "Charge mode": "charge_mode",
    "Display name": "display_name",
    "Target ID": "target_id",
    "Target level": "target_level",
    "Target group": "target_group",
    "Parent target ID": "parent_target_id",
    "Ion struct ID": "ion_struct_id",
    "Ion struct target ID": "ion_struct_target_id",
    "Ion role": "ion_role",
    "Min hits": "min_hits",
    "Charge mode required": "charge_mode_required",
    "Derivatization required": "derivatization_required",
    "Parent target": "parent_target_id",
    "Child target": "child_target_id",
    "Flag ID": "flag_id",
    "Activation mode": "activation_mode",
    "Include descendants": "include_descendants",
    "Rule ID": "rule_id",
    "Composition field": "composition_field",
    "When value": "when_value",
    "Evidence target ID": "evidence_target_id",
}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def _normalize_column_name(name: str) -> str:
    stripped = str(name).strip()
    if stripped in COLUMN_ALIASES:
        return COLUMN_ALIASES[stripped]
    return stripped



def _normalize_blank(value: Any) -> Any:
    if pd.isna(value):
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    return value



def _normalize_bool(value: Any, *, field_name: str) -> bool:
    v = _normalize_blank(value)
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)) and v in (0, 1):
        return bool(v)
    if isinstance(v, str):
        text = v.strip().lower()
        if text in {"true", "yes", "y", "1"}:
            return True
        if text in {"false", "no", "n", "0"}:
            return False
    raise WorkbookValidationError(f"Invalid boolean for {field_name}: {value!r}")



def _coerce_float(value: Any, *, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise WorkbookValidationError(f"Invalid float for {field_name}: {value!r}") from exc
    return out



def _coerce_int(value: Any, *, field_name: str) -> int:
    try:
        out = int(value)
    except (TypeError, ValueError) as exc:
        raise WorkbookValidationError(f"Invalid int for {field_name}: {value!r}") from exc
    return out



def _require_columns(sheet_name: str, df: pd.DataFrame) -> None:
    required = REQUIRED_COLUMNS[sheet_name]
    missing = required - set(df.columns)
    if missing:
        raise WorkbookValidationError(
            f"Sheet '{sheet_name}' is missing required columns: {sorted(missing)}"
        )



def _assert_nonempty(value: Any, *, field_name: str, row_idx: int, sheet_name: str) -> str:
    v = _normalize_blank(value)
    if v is None:
        raise WorkbookValidationError(
            f"Sheet '{sheet_name}' row {row_idx}: field '{field_name}' must not be empty"
        )
    return str(v).strip()


def _rule_is_active(ev: RuleEvaluation) -> bool:
    """
    A row is considered active if it was actually evaluated at runtime.

    Pass 1 rows have gate_open=None and are always active.
    Pass 2 rows are active only when gate_open is True.
    """
    if ev.gate_open is None:
        return True
    return bool(ev.gate_open)

# ---------------------------------------------------------------------------
# Workbook loading and parsing
# ---------------------------------------------------------------------------


class ScoreBLoader:
    def __init__(self, workbook_path: str | Path):
        self.workbook_path = Path(workbook_path)
        if not self.workbook_path.exists():
            raise FileNotFoundError(f"Workbook not found: {self.workbook_path}")

    def load(self) -> ScoreBConfig:
        """
        Phase 1 entry point.

        Steps:
            1. Read workbook sheets.
            2. Normalize column names.
            3. Validate required sheets/columns.
            4. Parse dataclasses.
            5. Validate cross-sheet references.
            6. Compile hierarchy and rule indexes.
        """
        sheets = self._read_required_sheets()

        ions = self._parse_ionlist(sheets["ionlist"])
        rules = self._parse_interpretation_rules(sheets["interpretation_v2"])
        hierarchy_edges = self._parse_hierarchy_edges(sheets["TargetHierarchy"])
        motif_policy_rules = self._parse_motif_policy_rules(sheets["MotifPolicy"])
        composition_rules = self._parse_composition_rules(sheets["CompositionConsistency"])
        unexpected_rules = self._parse_unexpected_penalty_rules(
            sheets["UnexpectedEvidencePenalty"]
        )

        hierarchy_index = self._compile_hierarchy(rules, hierarchy_edges)

        self._validate_cross_references(
            ions=ions,
            rules=rules,
            hierarchy_index=hierarchy_index,
            composition_rules=composition_rules,
            unexpected_rules=unexpected_rules,
            motif_policy_rules=motif_policy_rules,
        )

        pass1_rules = [r for r in rules if not r.parent_target_id]
        pass2_rules = [r for r in rules if r.parent_target_id]

        rules_by_target_id: dict[str, list[InterpretationRule]] = {}
        rules_by_ion_struct_id: dict[str, list[InterpretationRule]] = {}
        for rule in rules:
            rules_by_target_id.setdefault(rule.target_id, []).append(rule)
            rules_by_ion_struct_id.setdefault(rule.ion_struct_id, []).append(rule)

        self._validate_per_target_logic(rules_by_target_id)

        return ScoreBConfig(
            workbook_path=self.workbook_path,
            ions_by_struct_id={ion.structural_identifier: ion for ion in ions},
            interpretation_rules=rules,
            pass1_rules=pass1_rules,
            pass2_rules=pass2_rules,
            hierarchy_edges=hierarchy_edges,
            hierarchy_index=hierarchy_index,
            motif_policy_rules=motif_policy_rules,
            composition_rules=composition_rules,
            unexpected_penalty_rules=unexpected_rules,
            rules_by_target_id=rules_by_target_id,
            rules_by_ion_struct_id=rules_by_ion_struct_id,
        )

    def _read_required_sheets(self) -> dict[str, pd.DataFrame]:
        xls = pd.ExcelFile(self.workbook_path)
        available = set(xls.sheet_names)
        missing = set(REQUIRED_SHEETS) - available
        if missing:
            raise WorkbookValidationError(
                f"Workbook missing required sheets: {sorted(missing)}"
            )

        out: dict[str, pd.DataFrame] = {}
        for sheet_name in REQUIRED_SHEETS:
            df = pd.read_excel(self.workbook_path, sheet_name=sheet_name)
            df = df.rename(columns={c: _normalize_column_name(c) for c in df.columns})
            _require_columns(sheet_name, df)
            out[sheet_name] = df
        return out

    def _parse_ionlist(self, df: pd.DataFrame) -> list[IonDef]:
        ions: list[IonDef] = []
        seen_ids: set[str] = set()
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            struct_id = _assert_nonempty(
                row["structural_identifier"],
                field_name="structural_identifier",
                row_idx=actual_row,
                sheet_name="ionlist",
            )
            if struct_id in seen_ids:
                raise WorkbookValidationError(
                    f"Duplicate ion structural_identifier in ionlist: {struct_id!r}"
                )
            seen_ids.add(struct_id)
            ions.append(
                IonDef(
                    row_idx=actual_row,
                    glycotope_name=_assert_nonempty(
                        row["glycotope_name"],
                        field_name="glycotope_name",
                        row_idx=actual_row,
                        sheet_name="ionlist",
                    ),
                    ion_type=_assert_nonempty(
                        row["ion_type"],
                        field_name="ion_type",
                        row_idx=actual_row,
                        sheet_name="ionlist",
                    ),
                    fragmentation_mass=_coerce_float(
                        row["fragmentation_mass"],
                        field_name="fragmentation_mass",
                    ),
                    structural_identifier=struct_id,
                    charge_mode=_assert_nonempty(
                        row["charge_mode"],
                        field_name="charge_mode",
                        row_idx=actual_row,
                        sheet_name="ionlist",
                    ),
                    derivatization=_assert_nonempty(
                        row["derivatization"],
                        field_name="derivatization",
                        row_idx=actual_row,
                        sheet_name="ionlist",
                    ),
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return ions

    def _parse_interpretation_rules(self, df: pd.DataFrame) -> list[InterpretationRule]:
        rules: list[InterpretationRule] = []
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            min_hits = _coerce_int(row["min_hits"], field_name="min_hits")
            if min_hits < 1:
                raise WorkbookValidationError(
                    f"interpretation_v2 row {actual_row}: min_hits must be >= 1"
                )
            weight = _coerce_float(row["weight"], field_name="weight")
            if weight < 0:
                raise WorkbookValidationError(
                    f"interpretation_v2 row {actual_row}: weight must be >= 0"
                )

            rules.append(
                InterpretationRule(
                    row_idx=actual_row,
                    display_name=_assert_nonempty(
                        row["display_name"],
                        field_name="display_name",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    target_id=_assert_nonempty(
                        row["target_id"],
                        field_name="target_id",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    target_level=_assert_nonempty(
                        row["target_level"],
                        field_name="target_level",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    target_group=_normalize_blank(row["target_group"]),
                    parent_target_id=_normalize_blank(row["parent_target_id"]),
                    ion_struct_id=_assert_nonempty(
                        row["ion_struct_id"],
                        field_name="ion_struct_id",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    ion_struct_target_id=_normalize_blank(row["ion_struct_target_id"]),
                    ion_role=_normalize_blank(row["ion_role"]),
                    logic=_assert_nonempty(
                        row["logic"],
                        field_name="logic",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ).upper(),
                    min_hits=min_hits,
                    weight=weight,
                    charge_mode_required=_assert_nonempty(
                        row["charge_mode_required"],
                        field_name="charge_mode_required",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    derivatization_required=_assert_nonempty(
                        row["derivatization_required"],
                        field_name="derivatization_required",
                        row_idx=actual_row,
                        sheet_name="interpretation_v2",
                    ),
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return rules

    def _parse_hierarchy_edges(self, df: pd.DataFrame) -> list[HierarchyEdge]:
        edges: list[HierarchyEdge] = []
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            parent = _normalize_blank(row["parent_target_id"])
            child = _normalize_blank(row["child_target_id"])
            if parent is None or child is None:
                continue
            edges.append(
                HierarchyEdge(
                    row_idx=actual_row,
                    parent_target_id=str(parent).strip(),
                    child_target_id=str(child).strip(),
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return edges

    def _parse_motif_policy_rules(self, df: pd.DataFrame) -> list[MotifPolicyRule]:
        rules: list[MotifPolicyRule] = []
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            flag_id = _normalize_blank(row["flag_id"])
            target_id = _normalize_blank(row["target_id"])
            if flag_id is None or target_id is None:
                continue
            rules.append(
                MotifPolicyRule(
                    row_idx=actual_row,
                    flag_id=str(flag_id).strip(),
                    target_id=str(target_id).strip(),
                    activation_mode=_assert_nonempty(
                        row["activation_mode"],
                        field_name="activation_mode",
                        row_idx=actual_row,
                        sheet_name="MotifPolicy",
                    ),
                    include_descendants=_normalize_bool(
                        row["include_descendants"],
                        field_name="include_descendants",
                    ),
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return rules

    def _parse_composition_rules(
        self, df: pd.DataFrame
    ) -> list[CompositionConsistencyRule]:
        rules: list[CompositionConsistencyRule] = []
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            rule_id = _normalize_blank(row["rule_id"])
            evidence_target = _normalize_blank(row["evidence_target_id"])
            if rule_id is None or evidence_target is None:
                continue
            penalty = _coerce_float(row["penalty"], field_name="penalty")
            if not (0 <= penalty <= 1):
                raise WorkbookValidationError(
                    f"CompositionConsistency row {actual_row}: penalty must be between 0 and 1"
                )
            rules.append(
                CompositionConsistencyRule(
                    row_idx=actual_row,
                    rule_id=str(rule_id).strip(),
                    composition_field=_assert_nonempty(
                        row["composition_field"],
                        field_name="composition_field",
                        row_idx=actual_row,
                        sheet_name="CompositionConsistency",
                    ),
                    when_value=_assert_nonempty(
                        row["when_value"],
                        field_name="when_value",
                        row_idx=actual_row,
                        sheet_name="CompositionConsistency",
                    ),
                    evidence_target_id=str(evidence_target).strip(),
                    penalty=penalty,
                    mode=str(_normalize_blank(row.get("mode")) or "auto"),
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return rules

    def _parse_unexpected_penalty_rules(
        self, df: pd.DataFrame
    ) -> list[UnexpectedEvidencePenaltyRule]:
        rules: list[UnexpectedEvidencePenaltyRule] = []
        for row_idx, row in df.iterrows():
            actual_row = row_idx + 2
            flag_id = _normalize_blank(row["flag_id"])
            target_id = _normalize_blank(row["target_id"])
            if flag_id is None or target_id is None:
                continue
            penalty = _coerce_float(row["penalty"], field_name="penalty")
            if not (0 <= penalty <= 1):
                raise WorkbookValidationError(
                    f"UnexpectedEvidencePenalty row {actual_row}: penalty must be between 0 and 1"
                )
            rules.append(
                UnexpectedEvidencePenaltyRule(
                    row_idx=actual_row,
                    flag_id=str(flag_id).strip(),
                    target_id=str(target_id).strip(),
                    include_descendants=_normalize_bool(
                        row["include_descendants"],
                        field_name="include_descendants",
                    ),
                    penalty=penalty,
                    notes=_normalize_blank(row.get("notes")),
                )
            )
        return rules

    def _compile_hierarchy(
        self,
        rules: list[InterpretationRule],
        edges: list[HierarchyEdge],
    ) -> HierarchyIndex:
        concrete_targets = {r.target_id for r in rules}
        children_by_parent: dict[str, set[str]] = {}
        group_nodes: set[str] = set()

        for edge in edges:
            children_by_parent.setdefault(edge.parent_target_id, set()).add(edge.child_target_id)
            group_nodes.add(edge.parent_target_id)

        descendants_by_target: dict[str, set[str]] = {}

        def dfs(root: str, node: str, path: list[str], seen: set[str]) -> set[str]:
            if node in path:
                cycle = " -> ".join(path + [node])
                raise HierarchyValidationError(f"Cycle detected in TargetHierarchy: {cycle}")
            out: set[str] = set()
            if node in concrete_targets:
                out.add(node)
            for child in children_by_parent.get(node, set()):
                if child in seen:
                    continue
                out.update(dfs(root, child, path + [node], seen | {child}))
            return out

        all_nodes = group_nodes | concrete_targets | {e.child_target_id for e in edges}
        for node in all_nodes:
            descendants_by_target[node] = dfs(node, node, [], {node})

        return HierarchyIndex(
            children_by_parent=children_by_parent,
            descendants_by_target=descendants_by_target,
            concrete_targets=concrete_targets,
            group_nodes=group_nodes,
        )

    def _validate_cross_references(
        self,
        *,
        ions: list[IonDef],
        rules: list[InterpretationRule],
        hierarchy_index: HierarchyIndex,
        composition_rules: list[CompositionConsistencyRule],
        unexpected_rules: list[UnexpectedEvidencePenaltyRule],
        motif_policy_rules: list[MotifPolicyRule],
    ) -> None:
        ion_ids = {ion.structural_identifier for ion in ions}
        concrete_targets = hierarchy_index.concrete_targets
        group_nodes = hierarchy_index.group_nodes
        known_targets = concrete_targets | group_nodes

        for rule in rules:
            if rule.ion_struct_id not in ion_ids:
                raise WorkbookValidationError(
                    f"interpretation_v2 row {rule.row_idx}: ion_struct_id {rule.ion_struct_id!r} not found in ionlist"
                )
            if rule.parent_target_id and rule.parent_target_id not in known_targets:
                raise WorkbookValidationError(
                    f"interpretation_v2 row {rule.row_idx}: parent_target_id {rule.parent_target_id!r} is unknown"
                )

        valid_comp_fields = {"H", "N", "S", "G", "K", "F"}
        for rule in composition_rules:
            if rule.composition_field not in valid_comp_fields:
                raise WorkbookValidationError(
                    f"CompositionConsistency row {rule.row_idx}: invalid composition_field {rule.composition_field!r}"
                )
            if rule.evidence_target_id not in known_targets:
                raise WorkbookValidationError(
                    f"CompositionConsistency row {rule.row_idx}: unknown evidence_target_id {rule.evidence_target_id!r}"
                )

        for rule in unexpected_rules:
            if rule.target_id not in known_targets:
                raise WorkbookValidationError(
                    f"UnexpectedEvidencePenalty row {rule.row_idx}: unknown target_id {rule.target_id!r}"
                )

        for rule in motif_policy_rules:
            if rule.target_id not in known_targets:
                raise WorkbookValidationError(
                    f"MotifPolicy row {rule.row_idx}: unknown target_id {rule.target_id!r}"
                )

    def _validate_per_target_logic(
        self,
        rules_by_target_id: dict[str, list[InterpretationRule]],
    ) -> None:
        for target_id, target_rules in rules_by_target_id.items():
            logics = {r.logic for r in target_rules}
            if len(logics) > 1:
                raise WorkbookValidationError(
                    f"target_id {target_id!r} mixes multiple logic modes: {sorted(logics)}"
                )


# ---------------------------------------------------------------------------
# Pass-1 / Pass-2 evaluation engine
# ---------------------------------------------------------------------------


class ScoreBEngine:
    def __init__(self, config: ScoreBConfig, *, gate_threshold: float = 0.0):
        self.config = config
        self.gate_threshold = gate_threshold

    def evaluate_global_motif_evidence(
        self,
        spectrum: SpectrumEvidence,
    ) -> GlobalEvidenceResult:
        """
        Evaluate motif evidence in two passes.

        Pass 1
            - rows with empty parent_target_id
            - direct evidence only

        Pass 2
            - rows with parent_target_id set
            - only evaluated when gate is open based on Pass-1 evidence
        """
        hit_counts = spectrum.hit_count_by_ion_struct_id()
        target_evidence = self._initialize_target_evidence()

        pass1_evals = self._evaluate_rules(
            rules=self.config.pass1_rules,
            hit_counts=hit_counts,
            spectrum=spectrum,
            target_evidence=target_evidence,
            pass_name="pass1",
        )

        pass2_evals = self._evaluate_rules(
            rules=self.config.pass2_rules,
            hit_counts=hit_counts,
            spectrum=spectrum,
            target_evidence=target_evidence,
            pass_name="pass2",
        )
        self._finalize_target_evidence_logic(target_evidence)
        #self._finalize_support_flags(target_evidence)

        return GlobalEvidenceResult(
            target_evidence=target_evidence,
            pass1_rule_evaluations=pass1_evals,
            pass2_rule_evaluations=pass2_evals,
        )

    def _initialize_target_evidence(self) -> dict[str, TargetEvidence]:
        targets = {r.target_id for r in self.config.interpretation_rules}
        return {target_id: TargetEvidence(target_id=target_id) for target_id in targets}

    def _evaluate_rules(
        self,
        *,
        rules: list[InterpretationRule],
        hit_counts: dict[str, int],
        spectrum: SpectrumEvidence,
        target_evidence: dict[str, TargetEvidence],
        pass_name: str,
    ) -> list[RuleEvaluation]:
        evaluations: list[RuleEvaluation] = []
        for rule in rules:
            gate_open: Optional[bool] = None
            if pass_name == "pass2":
                assert rule.parent_target_id is not None
                gate_open = self._gate_is_open(rule.parent_target_id, target_evidence)
                if not gate_open:
                    eval_result = RuleEvaluation(
                        rule=rule,
                        gate_open=False,
                        charge_ok=False,
                        derivatization_ok=False,
                        hit_count=0,
                        supported=False,
                        row_score=0.0,
                        status="blocked_by_gate",
                    )
                    evaluations.append(eval_result)
                    self._accumulate_target_evidence(
                        rule=rule,
                        target_evidence=target_evidence,
                        row_score=0.0,
                        pass_name=pass_name,
                        counted_as_max=True,
                        is_active=False,
                        is_supported=False,
                    )
                    continue

            charge_ok = self._mode_matches(
                required=rule.charge_mode_required,
                observed=spectrum.charge_mode,
            )
            deriv_ok = self._mode_matches(
                required=rule.derivatization_required,
                observed=spectrum.derivatization,
            )
            hit_count = hit_counts.get(rule.ion_struct_id, 0)
            supported = charge_ok and deriv_ok and (hit_count >= rule.min_hits)
            row_score = rule.weight if supported else 0.0

            status = "supported" if supported else "not_supported"
            if not charge_ok:
                status = "charge_mode_mismatch"
            elif not deriv_ok:
                status = "derivatization_mismatch"
            elif hit_count < rule.min_hits:
                status = "insufficient_hits"

            eval_result = RuleEvaluation(
                rule=rule,
                gate_open=gate_open,
                charge_ok=charge_ok,
                derivatization_ok=deriv_ok,
                hit_count=hit_count,
                supported=supported,
                row_score=row_score,
                status=status,
            )
            evaluations.append(eval_result)

            self._accumulate_target_evidence(
                rule=rule,
                target_evidence=target_evidence,
                row_score=row_score,
                pass_name=pass_name,
                counted_as_max=True,
                is_active=_rule_is_active(eval_result),
                is_supported=eval_result.supported and eval_result.row_score > 0,
            )

        return evaluations

    def _mode_matches(self, *, required: str, observed: str) -> bool:
        """
        Minimal matcher.

        Extend later if workbook uses wildcard / ALL / comma-separated values.
        """
        req = required.strip().upper()
        obs = observed.strip().upper()
        if req in {"ALL", "ANY", "*"}:
            return True
        return req == obs

    def _gate_is_open(
        self,
        parent_target_id: str,
        target_evidence: dict[str, TargetEvidence],
    ) -> bool:
        """
        parent_target_id is a context gate, not hierarchy inheritance.

        Gate source is Pass-1 evidence only.
        If parent_target_id is a concrete target:
            gate opens when that target has pass1 evidence > threshold.
        If parent_target_id is a group node:
            gate opens when any expanded concrete descendant has pass1 evidence > threshold.
        """
        expanded_targets = self.config.hierarchy_index.expand_target(parent_target_id)
        if not expanded_targets:
            return False
        for target_id in expanded_targets:
            te = target_evidence[target_id]
            if te.pass1_max > 0:
                norm = te.pass1_score / te.pass1_max
            else:
                norm = 0.0
            if norm > self.gate_threshold:
                return True
        return False

    def _accumulate_target_evidence(
        self,
        *,
        rule: InterpretationRule,
        target_evidence: dict[str, TargetEvidence],
        row_score: float,
        pass_name: str,
        counted_as_max: bool,
        is_active: bool,
        is_supported: bool,
    ) -> None:
        te = target_evidence[rule.target_id]

        if pass_name == "pass1":
            te.pass1_score += row_score
            if counted_as_max:
                te.pass1_max += rule.weight
        elif pass_name == "pass2":
            te.pass2_score += row_score
            if counted_as_max:
                te.pass2_max += rule.weight
        else:
            raise ValueError(f"Unknown pass_name: {pass_name}")

        te.contributing_rule_rows.append(rule.row_idx)

        if is_active:
            te.active_rule_rows.append(rule.row_idx)

        if is_supported:
            te.supported_rule_rows.append(rule.row_idx)

    def _resolve_target_logic_mode(self, target_id: str) -> str:
        rules = self.config.rules_by_target_id.get(target_id, [])
        if not rules:
            return "OR"

        logics = {str(r.logic).strip().upper() for r in rules}
        if len(logics) > 1:
            raise WorkbookValidationError(
                f"target_id {target_id!r} mixes multiple logic modes at runtime: {sorted(logics)}"
            )
        return next(iter(logics)) if logics else "OR"

    def _finalize_target_evidence_logic(
        self,
        target_evidence: dict[str, TargetEvidence],
    ) -> None:
        """
        Phase 5A target-level finalization.

        - partial_score keeps the continuous evidence strength
        - logic_satisfied enforces target-level AND / OR semantics
        - supported is aligned to logic_satisfied
        """
        for target_id, te in target_evidence.items():
            te.contributing_rule_rows = sorted(set(te.contributing_rule_rows))
            te.active_rule_rows = sorted(set(te.active_rule_rows))
            te.supported_rule_rows = sorted(set(te.supported_rule_rows))

            te.logic_mode = self._resolve_target_logic_mode(target_id)
            te.partial_score = te.normalized_score
            te.active_rule_count = len(te.active_rule_rows)
            te.supported_rule_count = len(te.supported_rule_rows)

            if te.logic_mode == "AND":
                te.logic_satisfied = (
                    te.active_rule_count > 0
                    and te.supported_rule_count == te.active_rule_count
                )
            else:
                te.logic_satisfied = te.supported_rule_count > 0

            te.supported = te.logic_satisfied


# ---------------------------------------------------------------------------
# Candidate-level penalties and final Score B (placeholders for later phases)
# ---------------------------------------------------------------------------


COMPOSITION_INDEX = {
    "H": 0,
    "N": 1,
    "S": 2,
    "G": 3,
    "K": 4,
    "F": 5,
}


class CandidateScoreBScorer:
    def __init__(self, config: ScoreBConfig):
        self.config = config

    def evaluate_candidate(
        self,
        *,
        global_evidence: GlobalEvidenceResult,
        candidate_composition: tuple[int, int, int, int, int, int],
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> ScoreBResult:
        """
        Phase 5A flow:
            1. Resolve MotifPolicy candidate mask trace.
            2. Build policy-selected scoring units.
            3. Compute default motif support.
            4. Apply CompositionConsistency penalties.
            5. Apply UnexpectedEvidencePenalty.
            6. Clamp final score to [0, 1].
        """
        candidate_mask_trace = self.build_candidate_target_mask_trace(candidate_flags)
        scoring_units = self._build_scoring_units_from_mask(
            global_evidence.target_evidence,
            candidate_mask_trace,
        )

        if scoring_units:
            motif_support = self._compute_motif_support_score_default(scoring_units)
        else:
            motif_support = self._compute_motif_support_score(
                global_evidence.target_evidence,
                candidate_flags=candidate_flags,
            )

        aggregation_debug = self._compute_strength_coverage_summary(scoring_units)

        composition_penalty, applied_comp_rules = self._compute_composition_penalty(
            target_evidence=global_evidence.target_evidence,
            candidate_composition=candidate_composition,
        )
        unexpected_penalty, applied_unexpected_rules = self._compute_unexpected_penalty(
            target_evidence=global_evidence.target_evidence,
            candidate_flags=candidate_flags,
        )

        final_score = max(
            0.0,
            min(1.0, motif_support - composition_penalty - unexpected_penalty),
        )
        return ScoreBResult(
            candidate_composition=candidate_composition,
            motif_support_score=motif_support,
            composition_penalty=composition_penalty,
            unexpected_penalty=unexpected_penalty,
            final_score_b=final_score,
            target_evidence=global_evidence.target_evidence,
            applied_comp_rules=applied_comp_rules,
            applied_unexpected_rules=applied_unexpected_rules,
            candidate_mask_trace=candidate_mask_trace,
            scoring_units=scoring_units,
            aggregation_debug=aggregation_debug,
        )

    def _normalize_candidate_flags(
        self,
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> dict[str, bool]:
        if not candidate_flags:
            return {}
        return {str(k): bool(v) for k, v in candidate_flags.items()}

    def _is_policy_rule_active(
        self,
        rule: MotifPolicyRule,
        candidate_flags: dict[str, bool],
    ) -> tuple[Optional[bool], bool]:
        mode = str(rule.activation_mode).strip().lower()
        flag_value = candidate_flags.get(rule.flag_id)

        if mode in {"always"}:
            return flag_value, True
        if mode in {"when_true", "true", "if_true", "enable"}: #currently in sheet, it is enable when flags are selected
            return flag_value, bool(flag_value) is True
        if mode in {"when_false", "false", "if_false", "disable"}:
            return flag_value, (flag_value is False)

        raise WorkbookValidationError(
            f"Unsupported MotifPolicy activation_mode {rule.activation_mode!r} "
            f"at row {rule.row_idx}"
        )

    def build_candidate_target_mask_trace(
        self,
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> CandidateTargetMaskTrace:
        flags = self._normalize_candidate_flags(candidate_flags)

        matched_policy_rows: list[MotifPolicyMatch] = []
        selected_root_targets: set[str] = set()
        selected_concrete_targets: set[str] = set()

        for rule in sorted(self.config.motif_policy_rules, key=lambda r: r.row_idx):
            flag_value, is_active = self._is_policy_rule_active(rule, flags)

            expanded_targets: list[str] = []
            if is_active:
                selected_root_targets.add(rule.target_id)
                if rule.include_descendants:
                    expanded = self.config.hierarchy_index.expand_target(rule.target_id)
                    expanded_targets = sorted(expanded) if expanded else [rule.target_id]
                else:
                    expanded_targets = [rule.target_id]

                for target_id in expanded_targets:
                    selected_concrete_targets.add(target_id)

            matched_policy_rows.append(
                MotifPolicyMatch(
                    row_idx=rule.row_idx,
                    flag_id=rule.flag_id,
                    target_id=rule.target_id,
                    activation_mode=rule.activation_mode,
                    include_descendants=rule.include_descendants,
                    flag_value=flag_value,
                    is_active=is_active,
                    expanded_targets=expanded_targets,
                )
            )

        known_flag_ids = {r.flag_id for r in self.config.motif_policy_rules}
        unmatched_flags = sorted(set(flags) - known_flag_ids)

        return CandidateTargetMaskTrace(
            input_flags=flags,
            matched_policy_rows=matched_policy_rows,
            selected_root_targets=sorted(selected_root_targets),
            selected_concrete_targets=sorted(selected_concrete_targets),
            unmatched_flags=unmatched_flags,
        )

    def _get_selected_targets_from_trace(
        self,
        trace: CandidateTargetMaskTrace,
    ) -> list[str]:
        return list(trace.selected_concrete_targets)

    def _build_scoring_units_from_mask(
        self,
        target_evidence: dict[str, TargetEvidence],
        trace: CandidateTargetMaskTrace,
    ) -> list[ScoringUnit]:
        units: list[ScoringUnit] = []

        for match in trace.matched_policy_rows:
            if not match.is_active:
                continue

            source_target_ids = list(match.expanded_targets)
            if not source_target_ids:
                continue

            if len(source_target_ids) == 1:
                target_id = source_target_ids[0]
                te = target_evidence.get(target_id)
                score = te.partial_score if te is not None else 0.0
                logic_satisfied = te.logic_satisfied if te is not None else False
                unit_kind = "leaf"
            else:
                tes = [target_evidence.get(tid) for tid in source_target_ids]
                valid_tes = [te for te in tes if te is not None]

                if valid_tes:
                    score = max(te.partial_score for te in valid_tes)
                    logic_satisfied = any(te.logic_satisfied for te in valid_tes)
                else:
                    score = 0.0
                    logic_satisfied = False
                unit_kind = "group"

            units.append(
                ScoringUnit(
                    unit_id=f"policy_row_{match.row_idx}:{match.target_id}",
                    source_target_ids=source_target_ids,
                    score=score,
                    logic_satisfied=logic_satisfied,
                    weight=1.0,
                    unit_kind=unit_kind,
                )
            )

        return units

    def _compute_motif_support_score_default(
        self,
        scoring_units: list[ScoringUnit],
    ) -> float:
        if not scoring_units:
            return 0.0

        total_weight = sum(unit.weight for unit in scoring_units)
        if total_weight <= 0:
            return 0.0

        return sum(unit.score * unit.weight for unit in scoring_units) / total_weight

    def _compute_motif_support_score_experimental(
        self,
        scoring_units: list[ScoringUnit],
    ) -> float:
        # Placeholder hook for future comparison designs.
        return self._compute_motif_support_score_default(scoring_units)

    def _compute_strength_coverage_summary(
        self,
        scoring_units: list[ScoringUnit],
    ) -> dict[str, float]:
        if not scoring_units:
            return {
                "strength_mean": 0.0,
                "coverage_mean": 0.0,
                "hybrid_score": 0.0,
            }

        total_weight = sum(unit.weight for unit in scoring_units)
        if total_weight <= 0:
            return {
                "strength_mean": 0.0,
                "coverage_mean": 0.0,
                "hybrid_score": 0.0,
            }

        strength_mean = (
            sum(unit.score * unit.weight for unit in scoring_units) / total_weight
        )
        coverage_mean = (
            sum((1.0 if unit.logic_satisfied else 0.0) * unit.weight for unit in scoring_units)
            / total_weight
        )
        hybrid_score = 0.7 * strength_mean + 0.3 * coverage_mean

        return {
            "strength_mean": strength_mean,
            "coverage_mean": coverage_mean,
            "hybrid_score": hybrid_score,
        }

    def _compute_motif_support_score(
        self,
        target_evidence: dict[str, TargetEvidence],
        *,
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> float:
        """
        Phase 5A default implementation.

        Behavior:
            - resolve candidate-relevant policy rows
            - expand to concrete targets
            - collapse descendants per active MotifPolicy row
            - average policy-selected scoring units
        """
        if not target_evidence:
            return 0.0

        trace = self.build_candidate_target_mask_trace(candidate_flags)
        scoring_units = self._build_scoring_units_from_mask(target_evidence, trace)

        if scoring_units:
            return self._compute_motif_support_score_default(scoring_units)

        # Fallback for cases where no MotifPolicy row is active:
        vals = [te.partial_score for te in target_evidence.values()]
        return sum(vals) / len(vals) if vals else 0.0

    def _compute_composition_penalty(
        self,
        *,
        target_evidence: dict[str, TargetEvidence],
        candidate_composition: tuple[int, int, int, int, int, int],
    ) -> tuple[float, list[dict[str, Any]]]:
        """
        Phase 4B implementation.

        Rule behavior:
            - composition_field selects one value from candidate_composition
            - when_value defines the activation condition
            - if rule mode is "auto":
                * when_value == 0 / ==0  -> penalize_present
                * when_value like >=1, >0 -> penalize_absent
            - evidence_target_id may be a concrete target or a hierarchy group
            - present penalty  = penalty * max_evidence
            - absent penalty   = penalty * (1 - max_evidence)
            - total penalty is summed and clamped to 1.0
        """
        total_penalty = 0.0
        applied_rules: list[dict[str, Any]] = []

        for rule in sorted(self.config.composition_rules, key=lambda r: r.row_idx):
            candidate_value = self._candidate_comp_value(
                candidate_composition,
                rule.composition_field,
            )

            condition_matched = self._matches_when_value(
                candidate_value,
                rule.when_value,
            )
            if not condition_matched:
                continue

            resolved_targets = self._resolve_comp_targets(rule.evidence_target_id)

            target_scores: list[tuple[str, float]] = []
            for target_id in sorted(resolved_targets):
                te = target_evidence.get(target_id)
                score = te.normalized_score if te is not None else 0.0
                target_scores.append((target_id, score))

            if target_scores:
                trigger_target_id, max_evidence = max(target_scores, key=lambda x: x[1])
            else:
                trigger_target_id, max_evidence = None, 0.0

            mode = self._infer_comp_rule_mode(rule)

            if mode == "penalize_present":
                evidence_term = max_evidence
            elif mode == "penalize_absent":
                evidence_term = 1.0 - max_evidence
            else:
                raise ValueError(
                    f"Unsupported CompositionConsistency mode '{mode}' "
                    f"for rule row {rule.row_idx}"
                )

            # defensive clamp
            evidence_term = max(0.0, min(1.0, evidence_term))
            rule_penalty = rule.penalty * evidence_term
            total_penalty += rule_penalty

            applied_rules.append(
                {
                    "row_idx": rule.row_idx,
                    "rule_id": rule.rule_id,
                    "composition_field": rule.composition_field,
                    "candidate_value": candidate_value,
                    "when_value": rule.when_value,
                    "mode": mode,
                    "evidence_target_id": rule.evidence_target_id,
                    "resolved_targets": sorted(resolved_targets),
                    "trigger_target_id": trigger_target_id,
                    "max_evidence": max_evidence,
                    "evidence_term": evidence_term,
                    "base_penalty": rule.penalty,
                    "applied_penalty": rule_penalty,
                }
            )

        total_penalty = min(1.0, total_penalty)
        return total_penalty, applied_rules

    def _candidate_comp_value(
        self,
        candidate_composition: tuple[int, int, int, int, int, int],
        composition_field: str,
    ) -> int:
        key = str(composition_field).strip().upper()
        if key not in COMPOSITION_INDEX:
            raise ValueError(f"Unknown composition_field: {composition_field!r}")
        return int(candidate_composition[COMPOSITION_INDEX[key]])

    def _matches_when_value(
        self,
        candidate_value: int,
        when_value: str,
    ) -> bool:
        text = str(when_value).strip().replace(" ", "")

        if text.startswith(">="):
            return candidate_value >= int(text[2:])
        if text.startswith("<="):
            return candidate_value <= int(text[2:])
        if text.startswith("=="):
            return candidate_value == int(text[2:])
        if text.startswith("!="):
            return candidate_value != int(text[2:])
        if text.startswith(">"):
            return candidate_value > int(text[1:])
        if text.startswith("<"):
            return candidate_value < int(text[1:])

        # plain integer like "0", "1", "2"
        return candidate_value == int(text)

    def _infer_comp_rule_mode(
        self,
        rule: CompositionConsistencyRule,
    ) -> str:
        if rule.mode and rule.mode != "auto":
            return rule.mode

        text = str(rule.when_value).strip().replace(" ", "")

        # For current workbook semantics:
        #   0 / ==0       -> penalize present evidence
        #   >=1 / >0      -> penalize missing evidence (1 - E)
        if text in {"0", "==0"}:
            return "penalize_present"
        if text.startswith(">=") or text.startswith(">"):
            return "penalize_absent"

        # conservative default
        return "penalize_present"

    def _resolve_comp_targets(
        self,
        evidence_target_id: str,
    ) -> set[str]:
        expanded = self.config.hierarchy_index.expand_target(evidence_target_id)
        if expanded:
            return expanded
        return {evidence_target_id}

    def _compute_unexpected_penalty(
        self,
        *,
        target_evidence: dict[str, TargetEvidence],
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> tuple[float, list[dict[str, Any]]]:
        """
        Phase 4A implementation.

        Behavior:
            - active rules are selected by candidate_flags[flag_id] == True
            - forbidden target/group is resolved to concrete targets
            - penalty is scaled by the strongest forbidden evidence
            - multiple penalties are summed and clamped to 1.0
        """
        if not candidate_flags:
            return 0.0, []

        active_flag_ids = {
            flag_id
            for flag_id, enabled in candidate_flags.items()
            if bool(enabled)
        }
        if not active_flag_ids:
            return 0.0, []

        total_penalty = 0.0
        applied_rules: list[dict[str, Any]] = []

        for rule in sorted(self.config.unexpected_penalty_rules, key=lambda r: r.row_idx):
            if rule.flag_id not in active_flag_ids:
                continue

            resolved_targets = self._resolve_penalty_targets(
                rule.target_id,
                include_descendants=rule.include_descendants,
            )

            target_scores: list[tuple[str, float]] = []
            for target_id in sorted(resolved_targets):
                te = target_evidence.get(target_id)
                score = te.normalized_score if te is not None else 0.0
                target_scores.append((target_id, score))

            if target_scores:
                trigger_target_id, max_evidence = max(target_scores, key=lambda x: x[1])
            else:
                trigger_target_id, max_evidence = None, 0.0

            rule_penalty = rule.penalty * max_evidence
            total_penalty += rule_penalty

            applied_rules.append(
                {
                    "row_idx": rule.row_idx,
                    "flag_id": rule.flag_id,
                    "target_id": rule.target_id,
                    "include_descendants": rule.include_descendants,
                    "resolved_targets": sorted(resolved_targets),
                    "trigger_target_id": trigger_target_id,
                    "max_evidence": max_evidence,
                    "base_penalty": rule.penalty,
                    "applied_penalty": rule_penalty,
                }
            )

        total_penalty = min(1.0, total_penalty)
        return total_penalty, applied_rules

    def _resolve_penalty_targets(
        self,
        target_id: str,
        *,
        include_descendants: bool,
    ) -> set[str]:
        if include_descendants:
            expanded = self.config.hierarchy_index.expand_target(target_id)
            if expanded:
                return expanded
        return {target_id}


# ---------------------------------------------------------------------------
# Debug helpers
# ---------------------------------------------------------------------------


def summarize_target_evidence(evidence: GlobalEvidenceResult) -> pd.DataFrame:
    rows = []

    for target_id in sorted(evidence.target_evidence):
        s = evidence.target_evidence[target_id]
        rows.append(
            {
                "target_id": s.target_id,
                "pass1_score": s.pass1_score,
                "pass2_score": s.pass2_score,
                "total_raw_score": s.total_raw_score,
                "total_max_score": s.total_max_score,
                "normalized_score": s.normalized_score,
                "partial_score": s.partial_score,
                "logic_mode": s.logic_mode,
                "active_rule_count": s.active_rule_count,
                "supported_rule_count": s.supported_rule_count,
                "logic_satisfied": s.logic_satisfied,
                "supported": s.supported,
                "contributing_rule_rows": ",".join(map(str, s.contributing_rule_rows)),
                "active_rule_rows": ",".join(map(str, s.active_rule_rows)),
                "supported_rule_rows": ",".join(map(str, s.supported_rule_rows)),
            }
        )

    return pd.DataFrame(rows)



def summarize_rule_evaluations(evaluations: Iterable[RuleEvaluation]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for ev in evaluations:
        rows.append(
            {
                "row_idx": ev.rule.row_idx,
                "display_name": ev.rule.display_name,
                "target_id": ev.rule.target_id,
                "parent_target_id": ev.rule.parent_target_id,
                "ion_struct_id": ev.rule.ion_struct_id,
                "min_hits": ev.rule.min_hits,
                "weight": ev.rule.weight,
                "gate_open": ev.gate_open,
                "charge_ok": ev.charge_ok,
                "derivatization_ok": ev.derivatization_ok,
                "hit_count": ev.hit_count,
                "supported": ev.supported,
                "row_score": ev.row_score,
                "status": ev.status,
            }
        )
    return pd.DataFrame(rows)

def summarize_candidate_mask_trace(trace: CandidateTargetMaskTrace) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for match in trace.matched_policy_rows:
        rows.append(
            {
                "row_idx": match.row_idx,
                "flag_id": match.flag_id,
                "flag_value": match.flag_value,
                "activation_mode": match.activation_mode,
                "is_active": match.is_active,
                "target_id": match.target_id,
                "include_descendants": match.include_descendants,
                "expanded_targets": ",".join(match.expanded_targets),
            }
        )
    return pd.DataFrame(rows)


def summarize_scoring_units(scoring_units: list[ScoringUnit]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for unit in scoring_units:
        rows.append(
            {
                "unit_id": unit.unit_id,
                "unit_kind": unit.unit_kind,
                "source_target_ids": ",".join(unit.source_target_ids),
                "score": unit.score,
                "logic_satisfied": unit.logic_satisfied,
                "weight": unit.weight,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# High-level convenience API
# ---------------------------------------------------------------------------



def load_scoreb_workbook(path: str | Path) -> ScoreBConfig:
    loader = ScoreBLoader(path)
    return loader.load()



def evaluate_global_motif_evidence(
    config: ScoreBConfig,
    spectrum: SpectrumEvidence,
    *,
    gate_threshold: float = 0.0,
) -> GlobalEvidenceResult:
    engine = ScoreBEngine(config, gate_threshold=gate_threshold)
    return engine.evaluate_global_motif_evidence(spectrum)



def evaluate_candidate_score_b(
    config: ScoreBConfig,
    global_evidence: GlobalEvidenceResult,
    candidate_composition: tuple[int, int, int, int, int, int],
    candidate_flags: Optional[dict[str, bool]] = None,
) -> ScoreBResult:
    scorer = CandidateScoreBScorer(config)
    return scorer.evaluate_candidate(
        global_evidence=global_evidence,
        candidate_composition=candidate_composition,
        candidate_flags=candidate_flags,
    )

# ---------------------------------------------------------------------------
# Smoke Test/ QA function
# ---------------------------------------------------------------------------
#def run_gate_case(cfg: ScoreBConfig, case_name: str, observed_hits: list[ObservedIonHit]) -> None:
def run_gate_case(
    cfg: ScoreBConfig,
    case_name: str,
    observed_hits: list[ObservedIonHit],
    focus_target_ids: Optional[list[str]] = None,
    focus_row_indices: Optional[list[int]] = None,
) -> GlobalEvidenceResult:
    spectrum = SpectrumEvidence(
        observed_hits=observed_hits,
        charge_mode="any",
        derivatization="any",
    )
    evidence = evaluate_global_motif_evidence(cfg, spectrum)

    focus_target_ids = focus_target_ids or []
    focus_row_indices = focus_row_indices or []

    print(f"\n=== {case_name} : target evidence ===")
    target_df = summarize_target_evidence(evidence)
    if focus_target_ids:
        target_df = target_df[target_df["target_id"].isin(focus_target_ids)]
    print(target_df.to_string(index=False) if not target_df.empty else "(no matching target rows)")

    print(f"\n=== {case_name} : pass1 evals ===")
    pass1_df = summarize_rule_evaluations(evidence.pass1_rule_evaluations)
    if focus_row_indices:
        pass1_df = pass1_df[pass1_df["row_idx"].isin(focus_row_indices)]
    print(pass1_df.to_string(index=False) if not pass1_df.empty else "(no matching pass1 rows)")

    print(f"\n=== {case_name} : pass2 evals ===")
    pass2_df = summarize_rule_evaluations(evidence.pass2_rule_evaluations)
    if focus_row_indices:
        pass2_df = pass2_df[pass2_df["row_idx"].isin(focus_row_indices)]
    print(pass2_df.to_string(index=False) if not pass2_df.empty else "(no matching pass2 rows)")

    return evidence

#def run_candidate_case(cfg: ScoreBConfig, case_name: str, observed_hits: list[ObservedIonHit], candidate_composition: tuple[int, int, int, int, int, int], candidate_flags: dict[str, bool]) -> None:
def run_candidate_case(
    cfg: ScoreBConfig,
    case_name: str,
    observed_hits: list[ObservedIonHit],
    candidate_composition: tuple[int, int, int, int, int, int],
    candidate_flags: dict[str, bool],
    focus_target_ids: Optional[list[str]] = None,
    focus_row_indices: Optional[list[int]] = None,
) -> ScoreBResult:
    spectrum = SpectrumEvidence(
        observed_hits=observed_hits,
        charge_mode="any",
        derivatization="any",
    )
    evidence = evaluate_global_motif_evidence(cfg, spectrum)
    result = evaluate_candidate_score_b(
        cfg,
        evidence,
        candidate_composition=candidate_composition,
        candidate_flags=candidate_flags,
    )

    print(f"\n=== {case_name} : pass1 evals ===")
    print(summarize_rule_evaluations(evidence.pass1_rule_evaluations).to_string(index=False))

    print(f"\n=== {case_name} : pass2 evals ===")
    print(summarize_rule_evaluations(evidence.pass2_rule_evaluations).to_string(index=False))

    print(f"\n=== {case_name} : target evidence ===")
    print(summarize_target_evidence(evidence).to_string(index=False))

    print_candidate_debug(result)


def print_candidate_debug(result: ScoreBResult) -> None:
    print("\n=== Candidate target mask trace ===")
    if result.candidate_mask_trace is None:
        print("(no trace)")
    else:
        print("input_flags:", result.candidate_mask_trace.input_flags)
        print("selected_root_targets:", result.candidate_mask_trace.selected_root_targets)
        print("selected_concrete_targets:", result.candidate_mask_trace.selected_concrete_targets)
        print("unmatched_flags:", result.candidate_mask_trace.unmatched_flags)
        trace_df = summarize_candidate_mask_trace(result.candidate_mask_trace)
        if not trace_df.empty:
            print(trace_df.to_string(index=False))

    print("\n=== Scoring units ===")
    units_df = summarize_scoring_units(result.scoring_units)
    if units_df.empty:
        print("(no scoring units)")
    else:
        print(units_df.to_string(index=False))

    print("\n=== Aggregation debug ===")
    print(result.aggregation_debug)

    print("\n=== Candidate score B ===")
    print("motif_support_score:", result.motif_support_score)
    print("composition_penalty:", result.composition_penalty)
    print("unexpected_penalty:", result.unexpected_penalty)
    print("final_score_b:", result.final_score_b)

    if result.applied_comp_rules:
        print("\n=== Applied CompositionConsistency rules ===")
        print(pd.DataFrame(result.applied_comp_rules).to_string(index=False))

    if result.applied_unexpected_rules:
        print("\n=== Applied UnexpectedEvidencePenalty rules ===")
        print(pd.DataFrame(result.applied_unexpected_rules).to_string(index=False))

# ---------------------------------------------------------------------------
# Minimal smoke-test block
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    # Example manual smoke test.
    # Replace the path and sample hits with your real development case.
    #workbook = Path(r"G:\其他電腦\My Computer\GlycoMSParser\2025_demo\GlycoMSP_scoring_update_example_v5.xlsx")
    workbook = Path("GlycoMSP_scoring_update_example_v6.xlsx")
    if workbook.exists():
        cfg = load_scoreb_workbook(workbook)
        OPEN_GATE_ION = "Neu5Ac(a2-6)GalNAc(b1-4)GlcNAc-"
        CHILD_ION = "Neu5Ac(a2-6)GalNAc(b1-4)-"
        #Test T07
        run_candidate_case(
    cfg,
    "T07_integrated_candidate_score",
    [
        ObservedIonHit("Neu5Ac-", mz=0.0, intensity=700000.0, ppm_error=0.0),
        ObservedIonHit("Neu5Ac[dCH3OH]-", mz=0.0, intensity=650000.0, ppm_error=0.0),
        ObservedIonHit("Neu5AcGal(b1-4)GlcNAc-", mz=0.0, intensity=800000.0, ppm_error=0.0),
    ],
    candidate_composition=(3, 3, 1, 0, 0, 0),
    candidate_flags={
        "Allow5Ac": True,
        "AllowLewis": False,
        "Allow5Gc": False,
        "AllowLDNC": False,
        "AllowBG": False,
        "AllowFuc": False,
        "Forbid5Gc": False,
        "Forbid5Ac": False,
        "LeA_penalize_LeX": False,
        "LeX_penalize_LeA": False,
        "LeY_penalize_LeB": False,
        "LeB_penalize_LeY": False,
    },
    focus_target_ids=[
        "Neu5Ac",
        "SialylAcLacNAc",
        "SialylAcLDNC",
        "LeA",
        "LeX",
        "Neu5Gc",
    ],
)
    else:
        print("Workbook not found for smoke test. Adjust __main__ path before testing.")
"""
        #Test T01
        CHILD_ION = "Neu5Ac(a2-6)GalNAc(b1-4)-"
        run_gate_case(
            cfg,
            "T01_gate_child_only",
            [
                ObservedIonHit("Neu5Ac(a2-6)GalNAc(b1-4)-", mz=0.0, intensity=900000.0, ppm_error=0.0),
            ],
            focus_target_ids=["SialylAcLDNC", "SialylAcGalNAc"],
            focus_row_indices=[23, 26],
        )
        #Test T02
        run_gate_case(
    cfg,
    "T02_gate_parent_only",
    [
        ObservedIonHit("Neu5Ac(a2-6)GalNAc(b1-4)GlcNAc-", mz=0.0, intensity=900000.0, ppm_error=0.0),
    ],
    focus_target_ids=["SialylAcLDNC", "SialylAcGalNAc"],
    focus_row_indices=[23, 26],
)
        #Test T03  
        run_gate_case(
            cfg,
            "T03_gate_parent_plus_child",
            [
                ObservedIonHit("Neu5Ac(a2-6)GalNAc(b1-4)GlcNAc-", mz=0.0, intensity=900000.0, ppm_error=0.0),
                ObservedIonHit("Neu5Ac(a2-6)GalNAc(b1-4)-", mz=0.0, intensity=700000.0, ppm_error=0.0),
            ],
            focus_target_ids=["SialylAcLDNC", "SialylAcGalNAc"],
            focus_row_indices=[23, 26],
        )

        #Test T04
        run_gate_case(
            cfg,
            "T04_and_partial_neu5ac",
            [
                ObservedIonHit("Neu5Ac-", mz=0.0, intensity=500000.0, ppm_error=0.0),
            ],
            focus_target_ids=["Neu5Ac"],
            focus_row_indices=[19, 20],
        )


        #Test T05
        run_candidate_case(
            cfg,
            "T05_policy_trace_lewis",
            [],
            candidate_composition=(3, 3, 0, 0, 0, 1),
            candidate_flags={
                "AllowLewis": True,
            },
        )

        #Test T06
        run_candidate_case(
    cfg,
    "T06_policy_trace_lewis_5ac",
    [],
    candidate_composition=(3, 3, 1, 0, 0, 1),
    candidate_flags={
        "AllowLewis": True,
        "Allow5Ac": True,
    },
)

"""


"""

        OPEN_GATE_ION = "Neu5Ac(a2-6)GalNAc(b1-4)GlcNAc-"
        CHILD_ION = "Neu5Ac(a2-6)GalNAc(b1-4)-"   # row 23 ion
        # Case A: child ion only -> expected blocked_by_gate
        run_gate_case(cfg,
            "CASE_A_child_only",
            [
                ObservedIonHit(
                    CHILD_ION,
                    mz=0.0,
                    intensity=900000.0,
                    ppm_error=0.0,
                ),
            ],
        )
        spectrum = SpectrumEvidence(
            observed_hits=[
                ObservedIonHit(
                    "Neu5AcGal-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
            ],
            charge_mode="any",
            derivatization="any",
        )

        evidence = evaluate_global_motif_evidence(cfg, spectrum)

        print("=== Pass 1 ===")
        print(summarize_rule_evaluations(evidence.pass1_rule_evaluations).to_string(index=False))
        print("=== Pass 2 ===")
        print(summarize_rule_evaluations(evidence.pass2_rule_evaluations).to_string(index=False))
        print("=== Target evidence ===")
        print(summarize_target_evidence(evidence).to_string(index=False))

        print("=== Target evidence ===")
        print(summarize_target_evidence(evidence).to_string(index=False))

        result = evaluate_candidate_score_b(
            cfg,
            evidence,
            candidate_composition=(3, 3, 1, 0, 0, 0),   # S>=1
            candidate_flags={
                "Forbid5Gc": False,
                "Forbid5Ac": False,
                "LeA_penalize_LeX": False,
                "LeX_penalize_LeA": False,
                "LeY_penalize_LeB": False,
                "LeB_penalize_LeY": False,
            },
        )

        print(summarize_target_evidence(evidence).to_string(index=False))
        print_candidate_debug(result)
"""