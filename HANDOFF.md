# PichiaSafeHarbor Handoff

## 当前状态

- 工作目录已是 Git 仓库（`C:\Users\63097\Documents\PichiaSafeHarbor`，2026-07-16 与旧的非 Git 工作区 `CursorProject\PichiaSafeHarbor` 合并为单一权威位置）；尚未创建任何 commit。
- SRR31989028 单重复受控探针已正式关闭：Run ID `srr31989028-probe-f439d29b82a7f8d6`，`execution_status=complete`、`verification_status=passed`、`probe_acceptance_status=failed`、`scientific_acceptance_status=blocked`。能力结论：coordinate coverage、splice、unique/multi/unmapped 分类通过；strand 与 strain-specific support 不可用。失败原因：核染色体高深度异常比例超过门槛；BAM、BAI 和完整 SAM record stream 双运行不一致。
- 2026-07-16 复核确认：该探针使用的 `PRJNA1210090`/`SRR31989028` 数据集实际是一项无关的工程改造产香兰素 *K. phaffii* 菌株胁迫转录组研究（NCBI，中国农业科学院 2025-01-14 提交），并非 X-33，其"野生型 GS115"身份本身也未获确认。
- 同日，项目负责人已授权安排 **X-33 本身的真实测序**，作为长期项目推进；具体测序类型与时间表未定。见 [ADR-0008](docs/adr/ADR-0008-x33-native-sequencing-redirect.md)：该 ADR 暂停了 ADR-0006/0007 定义的公开数据证据路线（不再执行 D0 诊断，不再下载 `SRR31989016`/`SRR31989027`），下一次科学里程碑改为绑定真实测序数据到位。
- 同日，[ADR-0009](docs/adr/ADR-0009-deprioritize-proxy-data-compensation.md) 把范围从"暂停 RNA-seq 路线"扩大到"暂停整条代理数据补偿路线"：GS115 主坐标重新定位为开发/测试占位，ADR-0004 的注释来源调查延续工作与 ADR-0005/Slice 0B 的证据整合边界轨道均予以暂停；与数据来源无关的框架代码继续投入。
- 同日，[ADR-0010](docs/adr/ADR-0010-streamlit-import-and-trigger-workflow.md) 修订了 Streamlit 的架构设计：不再要求纯只读，允许导入新参考数据并触发核心引擎生成新运行，但页面仍不得重新实现科学规则或编辑已有结果。这只是架构设计更新，Slice 4 的实际开发仍未获授权。
- 同日新增 Slice 1 候选窗口引擎（`src/pichia_safe_harbor/candidates.py`、`slice1.py`，CLI 子命令 `candidate-windows`），并以现有 GS115 注释作为占位输入验证端到端可运行性（104/104 测试通过）。这属于 ADR-0008/0009 允许的"引擎工程自测"：产物固定标记 `run_purpose=engine_readiness_test_not_scientific_output`、`scientific_acceptance_status=blocked`，buffer/min-window 参数是调用方显式传入的示意占位值，不是冻结阈值，不构成正式候选。
- 同日，[ADR-0011](docs/adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 把框架建设授权范围细化到各切片：**Slice 1 继续**、**Slice 3 的验收/一致性工具部分获准现在建设**、**Slice 4（Streamlit）正式获准开始实现**；同时**明确排除 Slice 2**（GS115-CBS7435 共线性核查）——其原始理由本身也是代理数据补偿，真实数据到位后必要性会大幅下降，因此与 Slice 0A/0B 一起暂停。
- 未提交、未推送、未创建 PR。

## 已授权事项

- 等待 X-33 真实测序安排；数据到位后需要新 ADR 评估其对主坐标、注释边界和候选阈值的影响，不得提前假定解除 `blocked`。
- 继续维护/扩展 Slice 1 引擎（规则逻辑、区间切分逻辑、acceptance/独立核验流程、测试、CLI），仅限于工程自测——不得把当前 GS115 占位产物包装、展示或引用为正式候选或科学结论。
- Streamlit（Slice 4）已实现（`app/`）并通过本地浏览器验收：触发过一次真实 baseline + candidate-windows 运行，六个页面都正确展示三态状态和 `collinearity_status = unavailable`。
- 维护 `docs/` 架构与规划文档，使其反映最新 ADR 决策。
- 不再对公开 RNA-seq 数据源进行新的调查或下载；不再继续 Slice 0A/0B 式的注释来源或证据整合调查；不实现 Slice 2。

## 未授权事项

- D0 alignment 确定性诊断、修复 SRR31989028 探针可重复性问题、下载另外两个重复——均由 ADR-0008 暂停。
- 继续寻找/资格化新的公开注释或转录组数据源——由 ADR-0009 暂停。
- 实现或维护 Slice 2（CBS7435 共线性核查）代码——由 ADR-0011 排除。
- 正式 Slice 1 候选生成、阈值冻结、Slice 3 中依赖可信边界数据的部分（冻结默认阈值、全基因组正式运行、人工抽查、文献回顾）。
- 提交、推送、创建 PR、付费云资源。

## 必读材料

1. [需求](docs/REQUIREMENTS.md)
2. [架构](docs/ARCHITECTURE.md)
3. [执行计划](docs/EXECUTION_PLAN.md)
4. [ADR 索引](docs/adr/README.md)
5. [ADR-0001 双参考策略](docs/adr/ADR-0001-reference-strategy.md)（现重定位为开发/测试占位，见 ADR-0009；Slice 2 部分见 ADR-0011）
6. [ADR-0003 状态模型](docs/adr/ADR-0003-reference-completeness-and-acceptance.md)
7. [ADR-0008 X-33 测序路线重定向](docs/adr/ADR-0008-x33-native-sequencing-redirect.md)（取代 ADR-0006/0007 的下一步安排，历史决策保留存档）
8. [ADR-0009 停用代理数据补偿路线](docs/adr/ADR-0009-deprioritize-proxy-data-compensation.md)（扩大 ADR-0008 范围，暂停 ADR-0004/0005 延续工作）
9. [ADR-0010 Streamlit 导入与触发架构](docs/adr/ADR-0010-streamlit-import-and-trigger-workflow.md)（修订 ADR-0002 的只读限制）
10. [ADR-0011 逐切片授权范围](docs/adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md)（Slice 1/3工具/4 授权，Slice 2 排除）

## 权威证据

- [SRR31989028 失败型验收（历史存档）](local_runs/controlled_probe/prjna1210090_srr31989028/probe_v1_run1/acceptance_manifest.json)
- [SRR31989028 重复性证据（历史存档）](local_runs/controlled_probe/prjna1210090_srr31989028/probe_v1_repeatability.json)
- [PRJNA1210090 元数据资格评估（记录了"未改造身份未确认"的限制）](local_runs/independent_transcript_sources/prjna1210090/qualification_v4_run1/metadata_qualification.md)
- [PRJNA1210090 原始 BioProject 元数据快照（确认实际研究主题）](local_runs/independent_transcript_sources/prjna1210090/source_files/bioproject_1210090.xml)
- Slice 1 引擎：`src/pichia_safe_harbor/candidates.py`、`slice1.py`、`tests/test_slice1.py`（无正式 `local_runs/` 产物；如需查看示例输出，用 CLI 对着已有的 Slice 0 baseline run 目录跑一次 `candidate-windows`）

## 验证与停止线

- 决策出口、资源预算和停止条件以 [执行计划](docs/EXECUTION_PLAN.md) 为准。
- 必须递归保护既有权威 evidence chain（`local_runs/controlled_probe/`、`local_runs/independent_transcript_sources/` 等），运行相关测试和独立 verifier。
- 不得重跑 SRR31989028 探针、解释其高深度异常、下载其他重复或修改其 repeatability contract。
- Slice 1 引擎产物必须保持 `run_purpose=engine_readiness_test_not_scientific_output` 标记；不得移除该标记或以此产物冒充正式候选。
- 未经明确要求，不提交、不推送、不创建 PR。
