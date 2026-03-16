from __future__ import annotations

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

    # Optional future debug/count fields
    # rule_count_total: int = 0
    # rule_count_active: int = 0
    # rule_count_supported: int = 0

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

        self._finalize_support_flags(target_evidence)

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

    def _finalize_support_flags(self, target_evidence: dict[str, TargetEvidence]) -> None:
        """
        Current simple support definition:
            supported if total_raw_score > 0

        Later, if target-level AND / OR logic requires stricter aggregation,
        this should be replaced by target-wise logic evaluation using the stored
        per-rule results.
        """
        for te in target_evidence.values():
            te.supported = te.total_raw_score > 0

            te.contributing_rule_rows = sorted(set(te.contributing_rule_rows))
            te.active_rule_rows = sorted(set(te.active_rule_rows))
            te.supported_rule_rows = sorted(set(te.supported_rule_rows))


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
        Placeholder for later phases.

        Planned flow:
            1. Apply MotifPolicy / candidate target mask.
            2. Compute motif support score from target evidence.
            3. Apply CompositionConsistency penalties.
            4. Apply UnexpectedEvidencePenalty.
            5. Clamp final score to [0, 1].
        """
        motif_support = self._compute_motif_support_score(
            global_evidence.target_evidence,
            candidate_flags=candidate_flags,
        )
        composition_penalty, applied_comp_rules = self._compute_composition_penalty(
            target_evidence=global_evidence.target_evidence,
            candidate_composition=candidate_composition,
        )
        unexpected_penalty, applied_unexpected_rules = self._compute_unexpected_penalty(
            target_evidence=global_evidence.target_evidence,
            candidate_flags=candidate_flags,
        )

        final_score = max(0.0, min(1.0, motif_support - composition_penalty - unexpected_penalty))
        return ScoreBResult(
            candidate_composition=candidate_composition,
            motif_support_score=motif_support,
            composition_penalty=composition_penalty,
            unexpected_penalty=unexpected_penalty,
            final_score_b=final_score,
            target_evidence=global_evidence.target_evidence,
            applied_comp_rules=applied_comp_rules,
            applied_unexpected_rules=applied_unexpected_rules,
        )

    def _compute_motif_support_score(
        self,
        target_evidence: dict[str, TargetEvidence],
        *,
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> float:
        """
        Temporary default implementation.

        Current behavior:
            - average normalized score across all concrete targets

        Replace later with policy-aware target masking / weighting.
        """
        if not target_evidence:
            return 0.0
        vals = [te.normalized_score for te in target_evidence.values()]
        return sum(vals) / len(vals)

    def _compute_composition_penalty(
        self,
        *,
        target_evidence: dict[str, TargetEvidence],
        candidate_composition: tuple[int, int, int, int, int, int],
    ) -> tuple[float, list[dict[str, Any]]]:
        """
        Placeholder.

        Planned behavior:
            - compile when_value predicates
            - evaluate candidate composition field
            - expand evidence_target_id if group
            - scale penalty by evidence strength or missing evidence (1-E)
        """
        _ = target_evidence, candidate_composition
        return 0.0, []

    def _compute_unexpected_penalty(
        self,
        *,
        target_evidence: dict[str, TargetEvidence],
        candidate_flags: Optional[dict[str, bool]] = None,
    ) -> tuple[float, list[dict[str, Any]]]:
        """
        Placeholder.

        Planned behavior:
            - resolve active forbid flags
            - expand forbidden groups if include_descendants
            - compute penalty from forbidden target evidence
        """
        _ = target_evidence, candidate_flags
        return 0.0, []


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
                "supported": s.supported,
                "contributing_rule_rows": ",".join(map(str, s.contributing_rule_rows)),
                "active_rule_rows": ",".join(map(str, s.active_rule_rows)),
                "supported_rule_rows": ",".join(map(str, s.supported_rule_rows)),
                # Optional future count fields
                # "rule_count_total": s.rule_count_total,
                # "rule_count_active": s.rule_count_active,
                # "rule_count_supported": s.rule_count_supported,
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
# Minimal smoke-test block
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    # Example manual smoke test.
    # Replace the path and sample hits with your real development case.
    workbook = Path(r"G:\其他電腦\My Computer\GlycoMSParser\2025_demo\GlycoMSP_scoring_update_example_v5.xlsx")

    if workbook.exists():
        cfg = load_scoreb_workbook(workbook)

        #test 2G  -  LeY parent evidence + one LeX BY ion
        spectrum = SpectrumEvidence(
            observed_hits=[

                ObservedIonHit(
                    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
                ObservedIonHit(
                    "-Gal(b1-4)[Fuc(a1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),

            ],
            charge_mode="any",
            derivatization="any",
        )

        """
        #test 1
        spectrum = SpectrumEvidence(
            observed_hits=[
                ObservedIonHit(
                    "GalNAc(b1-4)GlcNAc-",
                    mz=0.0,              # placeholder for now
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
            ],
            charge_mode="any",
            derivatization="any",
        )
        #test for any entry violating the criteria
        
        spectrum = SpectrumEvidence(
            observed_hits=[
                ObservedIonHit("EXAMPLE_ION_A", mz=366.14, intensity=1000, ppm_error=2.1),
                ObservedIonHit("EXAMPLE_ION_B", mz=528.19, intensity=800, ppm_error=1.5),
            ],
            charge_mode="POS",
            derivatization="PERMETHYL",
        )
        #test 2A — child hit only, gate should stay closed
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
        #test 2B — parent + child hit, gate should open
        spectrum = SpectrumEvidence(
            observed_hits=[
                ObservedIonHit(
                    "Neu5AcGal(b1-4)GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
            ],
            charge_mode="any",
            derivatization="any",
        )
        #test 2C — parent + child hit, gate should open
                        ObservedIonHit(
                    "Neu5Ac-",
                    mz=0.0,
                    intensity=5000.0,
                    ppm_error=0.0,
                ),
        #test 2D - LeX only
        spectrum = SpectrumEvidence(
            observed_hits=[

                ObservedIonHit(
                    "Gal(b1-4)[Fuc(a1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),

            ],
            charge_mode="any",
            derivatization="any",
        )
        #test 2E  -  sLeX parent evidence + one LeX BY ion
        Test 2-F Direct LeX + sLeX + LeX BY ion
        spectrum = SpectrumEvidence(
            observed_hits=[

                ObservedIonHit(
                    "Gal(b1-4)[Fuc(a1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
                ObservedIonHit(
                    "Neu5Ac(a2-3)Gal(b1-4)[Fuca(1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
                ObservedIonHit(
                    "-Gal(b1-4)[Fuc(a1-3)]GlcNAc-",
                    mz=0.0,
                    intensity=1000.0,
                    ppm_error=0.0,
                ),
        """
        evidence = evaluate_global_motif_evidence(cfg, spectrum)
        print("=== Pass 1 ===")
        print(summarize_rule_evaluations(evidence.pass1_rule_evaluations).to_string(index=False))
        print("=== Pass 2 ===")
        print(summarize_rule_evaluations(evidence.pass2_rule_evaluations).to_string(index=False))
        print("=== Target evidence ===")
        print(summarize_target_evidence(evidence).to_string(index=False))
    else:
        print("Workbook not found for smoke test. Adjust __main__ path before testing.")
