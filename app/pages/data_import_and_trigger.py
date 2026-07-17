from __future__ import annotations

from datetime import datetime
from pathlib import Path

import streamlit as st

import state as app_state
from common import catalog_root, clear_catalog_cache, load_catalog
from services.catalog import RUN_TYPE_BASELINE, displayable_entries
from services.trigger import trigger_baseline, trigger_candidate_windows

app_state.ensure_defaults()

st.title("数据处理与分析运行")
st.caption("按顺序完成下面两步，就可以从参考基因组数据生成候选安全港位点窗口。")
with st.expander("技术说明（面向开发者）"):
    st.write(
        "这里的运行按钮调用的是与命令行完全相同的核心函数"
        "（pichia_safe_harbor.pipeline.run_baseline / pichia_safe_harbor.slice1.run_slice1），"
        "页面本身不会重新实现参考校验、注释解析、区间计算或候选规则判定（ADR-0010）。"
    )

repo_root = Path(__file__).resolve().parents[2]


def _suggest_output_dir(prefix: str) -> str:
    return str(catalog_root() / f"{prefix}_{datetime.now():%Y%m%d_%H%M%S}")


BASELINE_OUTPUT_KEY = "psh_trigger.baseline_output_dir"
CANDIDATE_OUTPUT_KEY = "psh_trigger.candidate_output_dir"
st.session_state.setdefault(BASELINE_OUTPUT_KEY, _suggest_output_dir("baseline"))
st.session_state.setdefault(CANDIDATE_OUTPUT_KEY, _suggest_output_dir("candidate_windows"))

OUTPUT_DIR_HELP = (
    "每次运行都会新建一个目录保存结果，用来完整保留历史记录、不会覆盖旧结果。"
    "已经默认加上了当前时间，一般不需要修改；如果提示目录已存在，改一个新名字即可。"
)

st.header("步骤一：处理参考基因组，生成基线数据")
st.write("这一步读取参考基因组序列和注释文件，计算基因间区等基础信息，是下一步生成候选窗口的必要输入。")
with st.form("baseline_form"):
    manifest_path = st.text_input(
        "参考基因组清单文件路径",
        value=str(repo_root / "reference" / "manifest.v1.json"),
        help="记录了有哪些可用的参考基因组（比如 Strain-B）以及对应文件的校验信息。一般不需要修改，使用默认路径即可。",
    )
    data_dir = st.text_input(
        "参考基因组文件所在目录",
        value=str(repo_root / "reference" / "data"),
        help="该目录下需要包含清单文件里登记的基因组序列（FASTA）和注释（GFF3）文件。一般不需要修改，使用默认路径即可。",
    )
    baseline_output_dir = st.text_input(
        "运行结果保存目录（新建目录，不能是已经存在的文件夹）",
        key=BASELINE_OUTPUT_KEY,
        help=OUTPUT_DIR_HELP,
    )
    submitted = st.form_submit_button("开始处理，生成基线数据")

if submitted:
    output_path = Path(baseline_output_dir)
    if output_path.exists():
        st.error(
            "运行结果保存目录已经存在。为了完整保留每次运行的历史记录，系统不会写入到已存在的目录。"
            "请把上面的「运行结果保存目录」改成一个新名字后重新点击运行，例如："
        )
        st.code(_suggest_output_dir("baseline"), language=None)
    else:
        with st.spinner("正在处理参考基因组，请稍候..."):
            result = trigger_baseline(Path(manifest_path), "strain-b", Path(data_dir), output_path)
        st.session_state[app_state.KEY_LAST_TRIGGER_RESULT] = result
        if result["status"] == "success":
            manifest = result["run_manifest"]
            st.success(f"基线数据已生成（run_id = {manifest['run_id']}）")
            st.caption(
                f"执行状态 execution_status = {manifest['execution_status']}， "
                f"科学验收状态 scientific_acceptance_status = {manifest['scientific_acceptance_status']}"
            )
            clear_catalog_cache()
        else:
            st.error("处理失败，请检查上面填写的路径是否正确。")
            with st.expander("详细错误信息（供排查问题参考）"):
                st.code(result["error"])

st.divider()

st.header("步骤二：从基线数据生成候选安全港窗口")
st.write("这一步会在步骤一生成的基线数据基础上，按照下面的规则参数，从基因间区里划分出候选窗口。")
entries = load_catalog()
baseline_entries = displayable_entries(entries, RUN_TYPE_BASELINE)
if not baseline_entries:
    st.info("目前还没有可用的基线数据——请先完成上面的步骤一，或者到 **概览 Overview** 页面选择一个已有的运行。")
else:
    # Keyed by directory name, not run_id: two directories can legitimately share
    # the same content-addressed run_id (e.g. baseline has no varying parameters,
    # so re-running it against the same inputs reproduces the same run_id in a
    # different directory) -- run_dir.name is what's actually unique per catalog entry.
    baseline_options = {entry.run_dir.name: entry for entry in baseline_entries}
    with st.form("candidate_windows_form"):
        chosen_baseline_name = st.selectbox(
            "选择要使用的基线数据",
            options=list(baseline_options),
            help="只列出已经成功生成、可以使用的基线数据；如果找不到想要的结果，请到「概览」页面确认。",
        )
        candidate_data_dir = st.text_input(
            "参考基因组文件所在目录",
            value=str(repo_root / "reference" / "data"),
            help=(
                "需要提取候选窗口的实际序列，所以这一步也要用到参考基因组文件；"
                "必须和上面选定的基线数据来自同一份参考基因组，否则会被拒绝运行。一般不需要修改。"
            ),
        )
        buffer_text = st.text_input(
            "缓冲距离（bp）",
            value="",
            placeholder="请输入整数（不含单位）",
            help=(
                "候选窗口与两侧相邻基因之间至少要保留的距离，单位是碱基对（bp）。"
                "目前还没有正式冻结的科学阈值（真实阈值需要单独立项评审），请先输入一个数值用于工程测试，"
                "这个数值本身不代表科学结论。（对应运行结果里的 rule_params.buffer_distance_bp）"
            ),
        )
        min_window_text = st.text_input(
            "最小候选窗口长度（bp）",
            value="",
            placeholder="请输入整数（不含单位）",
            help=(
                "候选窗口至少要有多长才会被保留，单位是碱基对（bp）。同样目前没有正式冻结的科学阈值，"
                "请先输入一个数值用于工程测试。（对应运行结果里的 rule_params.min_candidate_window_bp）"
            ),
        )
        long_interval_percentile = st.number_input(
            "异常过长基因间区判定分位数",
            min_value=0.5,
            max_value=0.999,
            value=0.95,
            step=0.01,
            help=(
                "基因间区长度超过这个分位数时，会被标记为「异常过长」（仅标记、不会被排除）。"
                "默认 0.95（即最长的 5%），一般不需要修改。"
            ),
        )
        candidate_output_dir = st.text_input(
            "运行结果保存目录（新建目录，不能是已经存在的文件夹）",
            key=CANDIDATE_OUTPUT_KEY,
            help=OUTPUT_DIR_HELP,
        )
        submitted_cw = st.form_submit_button("开始生成候选窗口")

    if submitted_cw:
        errors = []
        try:
            buffer_value = int(buffer_text)
            if buffer_value < 0:
                errors.append("「缓冲距离」不能是负数")
                buffer_value = None
        except ValueError:
            buffer_value = None
            errors.append("「缓冲距离」必须填写一个整数")
        try:
            min_window_value = int(min_window_text)
            if min_window_value < 0:
                errors.append("「最小候选窗口长度」不能是负数")
                min_window_value = None
        except ValueError:
            min_window_value = None
            errors.append("「最小候选窗口长度」必须填写一个整数")

        output_path = Path(candidate_output_dir)
        output_exists = output_path.exists()
        if output_exists:
            errors.append(f"运行结果保存目录已经存在：{output_path}")

        if errors:
            for message in errors:
                st.error(message)
            if output_exists:
                st.write("请把上面的「运行结果保存目录」改成一个新名字后重新点击运行，例如：")
                st.code(_suggest_output_dir("candidate_windows"), language=None)
        else:
            baseline_dir = baseline_options[chosen_baseline_name].run_dir
            with st.spinner("正在生成候选安全港窗口，请稍候..."):
                result = trigger_candidate_windows(
                    baseline_dir,
                    Path(str(repo_root / "reference" / "manifest.v1.json")),
                    "strain-b",
                    Path(candidate_data_dir),
                    buffer_value,
                    min_window_value,
                    output_path,
                    long_interval_percentile,
                )
            st.session_state[app_state.KEY_LAST_TRIGGER_RESULT] = result
            if result["status"] == "success":
                manifest = result["run_manifest"]
                st.success(f"候选窗口已生成（run_id = {manifest['run_id']}）")
                st.caption(
                    f"执行状态 execution_status = {manifest['execution_status']}， "
                    f"科学验收状态 scientific_acceptance_status = {manifest['scientific_acceptance_status']}， "
                    f"运行用途 run_purpose = {manifest['run_purpose']}"
                )
                clear_catalog_cache()
            else:
                st.error("生成失败，请检查上面填写的参数是否正确。")
                with st.expander("详细错误信息（供排查问题参考）"):
                    st.code(result["error"])
