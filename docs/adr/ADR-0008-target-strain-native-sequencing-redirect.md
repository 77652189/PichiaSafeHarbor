# ADR-0008：Strain-T 原生测序获批后的证据路线重定向

- 状态：Accepted
- 日期：2026-07-16
- 关系：取代 [ADR-0006](ADR-0006-capability-scoped-transcript-evidence.md) 与 [ADR-0007](ADR-0007-srr31989028-controlled-probe-outcome.md) 中"依赖公开 RNA-seq 数据建立转录证据"的下一步安排；呼应 [ADR-0001](ADR-0001-reference-strategy.md)"未来获得 Strain-T 测序数据后可用新 ADR 替代本决策"的预留条款；[ADR-0005](ADR-0005-evidence-integrated-boundary-track.md) 关于"证据整合边界轨道需先通过数据资格门禁"的原则继续有效，仅证据来源改变。

## 背景

ADR-0006 建立了按能力（coordinate coverage / splice support / strand support / boundary support / strain-specific support）资格化转录证据的路线，并授权对 `PRJNA1210090` 的 `SRR31989028` 单重复执行受控探针。ADR-0007 记录了该探针的正式结果：`execution_status=complete`、`verification_status=passed`，但 `probe_acceptance_status=failed`——核染色体高深度异常 bin 比例 0.0713（阈值 0.05），且 `alignment.bam`/`alignment.bam.bai` 与完整 SAM record stream 在两次运行间不一致，项目保持 `scientific_acceptance_status=blocked`，下一步被限定为只读的 D0 alignment 确定性诊断。

2026-07-16 复核 `PRJNA1210090` 的 NCBI BioProject 元数据（`local_runs/independent_transcript_sources/prjna1210090/source_files/bioproject_1210090.xml`）确认：该项目标题为 *"Ferulic acid tolerance mechanisms in a vanillin-producing Komagataella phaffii"*，中国农业科学院 2025-01-14 提交，研究对象是一株**工程改造的产香兰素细胞工厂菌株在阿魏酸胁迫下的转录组**，与安全港位点研究无关；探针使用的三个"野生型、无胁迫"run 的未改造 Strain-B 身份本身也从未获得确认（见 `local_runs/independent_transcript_sources/prjna1210090/qualification_v4_run1/metadata_qualification.md`）。这意味着即使 D0 诊断顺利解决了可重复性问题，该数据源能提供的证据价值也天然有限——它既不是 Strain-T，也不确定是干净的 Strain-B。

同日，项目负责人在晨会中授权安排对 **Strain-T 本身的真实测序**，作为长期项目推进。这是比继续投入诊断/修复一个数据源存疑的公开数据集更直接、更贴近目标菌株的证据路径。

## 决策

1. 暂停 ADR-0006/ADR-0007 路线定义的下一步工作：不执行 D0 alignment 确定性诊断，不下载 `SRR31989016`/`SRR31989027`，不为修复 `SRR31989028` 探针的 BAM/BAI 可重复性问题投入新预算。这不是因为该路线被证明是错误决策——探针的执行、验证和失败判定本身是严谨且成立的——而是因为它服务的问题（"没有 Strain-T 数据时如何从公开数据中争取证据"）已经被更直接的方案取代。
2. 已产生的探针证据、独立核验与验收记录（`local_runs/controlled_probe/prjna1210090_srr31989028/` 及相关 `local_runs/independent_transcript_sources/`）保留作为历史/审计存档，不删除、不追溯修改。
3. 项目下一次科学证据里程碑改为绑定"Strain-T 原生测序数据到位并通过独立评估"，不再是"修复 BAM 可重复性问题"或"三重复验证"。测序类型（转录组、基因组重测序或从头组装）与交付时间表尚未确定；数据到位后必须由新 ADR 评估其对主坐标、注释边界置信度与候选阈值推断的实际影响，不得直接假定其解除现有 `blocked` 状态。
4. 在等待期间，Strain-B `GCA_001746955.1` 继续按 ADR-0001 作为唯一主坐标空间；允许在此坐标上继续**工程实现与自测**——具体是 Slice 1 候选窗口引擎（`src/pichia_safe_harbor/candidates.py`、`slice1.py`）——但仅限于验证引擎本身可运行、规则语义正确，产物必须显式标记 `run_purpose=engine_readiness_test_not_scientific_output` 且 `scientific_acceptance_status=blocked`。这类工程自测不等同于 ADR-0005 第6条所指"Slice 1 不得开始"中的正式候选生成；正式候选生成、阈值冻结、共线性核查与 Streamlit 仍需等待合格边界证据。
5. 不再对公开 RNA-seq 数据集展开新的来源调查或下载，除非未来出现明确的新理由（例如 Strain-T 真实数据被证明无法覆盖某项特定能力，需要公开数据补充）。

## 放弃的方案

### 继续投入 D0 诊断并按原计划推进三重复验证

即使确定性问题被修复，数据源本身（工程改造菌株胁迫实验、Strain-B 身份未确认）对安全港边界判断的证据价值仍然有限；在真实 Strain-T 测序已授权的情况下，这笔投入的边际价值低于等待真实数据。

### 因为探针失败就整体否定 ADR-0006 的能力分级方法论

按能力而非单一二元门禁资格化证据的方法论本身没有问题，未来 Strain-T 真实数据到位后，仍应沿用"覆盖 / 剪接 / 链方向 / 边界 / 菌株特异性"分别资格化的框架，只是换一个更合适的数据源。

### 立即假定 Strain-T 测序到位即解除科学阻塞

测序尚未交付，类型和时间表未定；提前解除 `blocked` 或跳过独立评估会重复本项目一贯警惕的"把尚不具备的证据当作已具备"的错误。

## 后果

- `HANDOFF.md` 与 `docs/EXECUTION_PLAN.md` 的"下一步授权"从 D0 诊断改为"等待 Strain-T 测序安排；同期可维护/扩展 Slice 1 引擎（占位数据，非科学产物）"。
- `docs/adr/README.md` 中 ADR-0006、ADR-0007 的状态标注为"Superseded by ADR-0008"，决策描述保留不变，作为历史记录。
- 项目下一次规划门禁绑定"Strain-T 测序数据到位后的评估"，不是"三重复是否值得"。
- `scientific_acceptance_status` 继续保持 `blocked`，不因本 ADR 而改变。
