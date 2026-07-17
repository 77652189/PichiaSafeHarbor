from __future__ import annotations

# Chinese display labels for the small, stable engine vocabularies (orientation,
# structural tier, rule ids, status fields, ...). Every function falls back to
# the raw value when it isn't in the table below -- these lists are not
# formally closed by the engine's contract, and a future engine value must
# never disappear from the UI just because this module hasn't been updated yet.
#
# Format is always "中文（raw_value）": the raw value is ADR-0010's required
# traceability anchor (never hide/rewrite real status), the Chinese text is
# the actual UX fix. Free-text engine-generated prose (the `reason` field on
# excluded regions/rule flags, or `known_limitations` sourced from NCBI
# metadata) is deliberately NOT covered here -- translating dynamically
# generated or externally-sourced sentences risks silently drifting from the
# authoritative source text, and is out of scope for an app-layer-only pass.

ORIENTATION_LABELS = {
    "tandem": "串联",
    "convergent": "收敛",
    "divergent": "背离",
    "unknown": "未知",
}

STRUCTURAL_TIER_LABELS = {
    "convergent_clean": "收敛无标记",
    "flagged": "已标记待复核",
}

EVIDENCE_LEVEL_LABELS = {
    "predicted": "预测",
}

BOUNDARY_CONFIDENCE_LABELS = {
    "high": "高",
    "uncertain": "存疑",
    "unknown": "未知",
}

RULE_ID_LABELS = {
    "neighbor_orientation_tandem": "邻近基因方向为串联",
    "neighbor_orientation_divergent": "邻近基因方向为背离",
    "unusually_long_interval": "基因间区异常过长",
    "boundary_confidence_insufficient": "边界置信度不足",
    "below_minimum_window_length": "缓冲后长度不足最小窗口",
}

COLLINEARITY_STATUS_LABELS = {
    "unavailable": "不可用",
}

EXECUTION_STATUS_LABELS = {
    "complete": "已完成",
}

VERIFICATION_STATUS_LABELS = {
    "not_run": "未验证",
    "passed": "已通过",
}

SCIENTIFIC_ACCEPTANCE_STATUS_LABELS = {
    "blocked": "未通过",
    "accepted": "已通过验收",
}

RUN_PURPOSE_LABELS = {
    "engine_readiness_test_not_scientific_output": "工程自测运行，非科学结论",
}

MISSING_RISK_TRACK_LABELS = {
    "centromere track": "着丝粒轨道",
    "telomere track": "端粒轨道",
    "repeat/mobile-element track": "重复序列/移动元件轨道",
    "replication-origin track": "复制起点轨道",
    "essential-gene track": "必需基因轨道",
    "Strain-T-specific variation track": "Strain-T 特异变异轨道",
}

# Short fixed status codes (not free-generated prose) -- safe to translate,
# unlike `reason`/`known_limitations` which are longer sentences that may
# quote an external source verbatim.
SCIENTIFIC_ACCEPTANCE_BLOCKER_LABELS = {
    "systematic_annotation_boundary_uncertainty_requires_authoritative_review": (
        "系统性注释边界不确定性，需要权威复核"
    ),
    "acceptance_manifest_not_yet_generated": "尚未生成验收清单",
    "rule_parameters_are_illustrative_placeholders_pending_threshold_adr": (
        "规则参数为示意占位值，正式阈值待相关 ADR 冻结"
    ),
    "input_annotation_is_strain-b_genbank_2016_release_not_an_adr_qualified_boundary_track": (
        "输入注释为 Strain-B GenBank 2016 版本，不是 ADR 认定的边界轨道"
    ),
    "no_target_strain_native_sequencing_data_used": "未使用 Strain-T 原生测序数据",
    "collinearity_with_strain-c_not_yet_computed": "尚未计算与 Strain-C 的共线性",
}


def _label(value: str | None, table: dict[str, str]) -> str:
    if value is None:
        return "n/a"
    return f"{table[value]}（{value}）" if value in table else value


def orientation_label(value: str | None) -> str:
    return _label(value, ORIENTATION_LABELS)


def structural_tier_label(value: str | None) -> str:
    return _label(value, STRUCTURAL_TIER_LABELS)


def evidence_level_label(value: str | None) -> str:
    return _label(value, EVIDENCE_LEVEL_LABELS)


def boundary_confidence_label(value: str | None) -> str:
    return _label(value, BOUNDARY_CONFIDENCE_LABELS)


def rule_id_label(value: str | None) -> str:
    return _label(value, RULE_ID_LABELS)


def collinearity_status_label(value: str | None) -> str:
    return _label(value, COLLINEARITY_STATUS_LABELS)


def execution_status_label(value: str | None) -> str:
    return _label(value, EXECUTION_STATUS_LABELS)


def verification_status_label(value: str | None) -> str:
    return _label(value, VERIFICATION_STATUS_LABELS)


def scientific_acceptance_status_label(value: str | None) -> str:
    return _label(value, SCIENTIFIC_ACCEPTANCE_STATUS_LABELS)


def run_purpose_label(value: str | None) -> str:
    return _label(value, RUN_PURPOSE_LABELS)


def missing_risk_track_label(value: str | None) -> str:
    return _label(value, MISSING_RISK_TRACK_LABELS)


def scientific_acceptance_blocker_label(value: str | None) -> str:
    return _label(value, SCIENTIFIC_ACCEPTANCE_BLOCKER_LABELS)
