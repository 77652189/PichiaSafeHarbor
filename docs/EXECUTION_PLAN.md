# PichiaSafeHarbor 项目执行计划

## 1. 决策摘要

- 项目仍未进入正式候选生成阶段，科学状态仍为 `blocked`。
- 依赖公开 RNA-seq 数据建立转录证据的路线（ADR-0006/ADR-0007）已于 2026-07-16 由 [ADR-0008](adr/ADR-0008-target-strain-native-sequencing-redirect.md) 正式暂停：`SRR31989028` 探针已执行完毕并正式失败（核染色体高深度异常比例超阈值、BAM/BAI 与完整 SAM record stream 双运行不可重复）；复核确认该数据集（`PRJNA1210090`）实际来自一项无关的工程改造产香兰素菌株胁迫转录组研究，并非 Strain-T，其"野生型 Strain-B"身份也从未获确认。继续投入诊断/修复这一数据源的边际价值低。
- 同日，项目负责人已授权安排 **Strain-T 本身的真实测序**，作为长期项目推进——这是解决"没有目标菌株原生数据"这一根本缺口的直接路径，而不是继续在公开数据上打补丁。具体测序类型（转录组、基因组重测序或从头组装）和交付时间表尚未确定。
- 同日，[ADR-0009](adr/ADR-0009-deprioritize-proxy-data-compensation.md) 把范围从"暂停 RNA-seq 路线"扩大到"暂停整条代理数据补偿路线"：Strain-B 主坐标重新定位为开发/测试占位，ADR-0004 的注释来源调查延续工作与 ADR-0005/Slice 0B 的证据整合边界轨道均予以暂停；与数据来源无关的框架代码继续投入。
- 同日，[ADR-0010](adr/ADR-0010-streamlit-import-and-trigger-workflow.md) 修订了 Streamlit 的架构设计：不再要求纯只读，允许导入新参考数据并触发核心引擎生成新运行，但仍禁止页面重新实现科学规则或编辑已有结果。
- 同日，[ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 把"框架代码可继续投入"的原则逐切片细化：**Slice 1**（候选窗口引擎）继续；**Slice 3** 中与数据来源无关的验收/一致性工具部分获准现在建设；**Slice 4**（Streamlit）正式获准开始实现。同时**明确排除 Slice 2**（Strain-B-Strain-C 共线性核查）——其原始理由（借 Strain-C 为代理菌株 Strain-B 佐证）本身也是一种代理数据补偿，真实 Strain-T 数据到位后这个交叉印证的必要性会大幅下降甚至消失，因此与 Slice 0A/0B 一样暂停，不作为占位框架代码现在实现。
- 不授权：继续 D0 alignment 确定性诊断、下载 `SRR31989016`/`SRR31989027`、继续 Slice 0A/0B 式的注释来源或证据整合调查、实现 Slice 2（Strain-C 共线性核查）、Slice 3 中依赖可信边界数据的部分（冻结阈值、全基因组正式运行、人工抽查、文献回顾）。

## 2. 证据边界与当前仓库状态

`C:\Users\63097\Documents\PichiaSafeHarbor` 现在是唯一权威 Git 仓库：2026-07-16 已将此前散落在非 Git 工作区 `C:\Users\63097\Documents\CursorProject\PichiaSafeHarbor` 的全部源码、测试、文档和 `local_runs` 证据合并进来（robocopy 校验通过，0 差异）。旧工作区原样保留作为未删除的备份，尚未提交任何 Git commit。

本计划的事实依据：

- [docs/REQUIREMENTS.md](REQUIREMENTS.md)、[docs/ARCHITECTURE.md](ARCHITECTURE.md)
- [docs/adr/README.md](adr/README.md) 及 ADR-0001 至 ADR-0011，尤其 [ADR-0008](adr/ADR-0008-target-strain-native-sequencing-redirect.md)、[ADR-0009](adr/ADR-0009-deprioritize-proxy-data-compensation.md)（路线重定向的权威记录）、[ADR-0010](adr/ADR-0010-streamlit-import-and-trigger-workflow.md)（Streamlit 架构修订）与 [ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md)（逐切片授权范围，排除 Slice 2）
- `local_runs/controlled_probe/prjna1210090_srr31989028/`（SRR31989028 探针的历史验收与重复性证据）
- `local_runs/independent_transcript_sources/prjna1210090/`（元数据资格评估与原始 BioProject 快照，记录了数据源与安全港研究无关的发现）
- `local_runs/slice0b/`（证据整合边界轨道数据资格门禁的历史评估，5 轮结论均为证据不足）
- `src/pichia_safe_harbor/candidates.py`、`slice1.py`、`tests/test_slice1.py`（Slice 1 引擎及其测试）

## 3. 距离原始 MVP 的主要缺口

### 科学缺口

1. Strain-B 当前公开注释存在系统性边界不确定性，不能作为正式阈值和候选窗口的唯一边界来源（高置信双边界区间仅约 0.66%）。
2. 尚无通过验收的全基因组证据整合边界轨道；`scientific_acceptance_status` 仍为 `blocked`。
3. 唯一执行过的单重复 RNA-seq 探针已正式失败（可重复性与高深度异常门槛均未通过），且其数据源本身与目标菌株无关；链方向、生物重复一致性、菌株特异性支持均不可用。
4. 目标菌株为 Strain-T，当前主坐标仍是 Strain-B `GCA_001746955.1`；`exact_target_strain_coordinates=false`。Strain-T 原生测序已获授权但尚未交付，交付后仍需独立 ADR 评估其实际影响，不能提前假定解除限制。
5. 尚未形成默认基因缓冲距离和最小候选窗口长度的科学依据及 ADR；Slice 1 引擎目前只能接受调用方显式传入的占位参数。
6. Slice 2（Strain-B 到 Strain-C 的候选级共线性核查）已暂停实现，其原始理由本身依赖"缺失 Strain-T 数据"这一前提（ADR-0011）。
7. 尚未完成文献整合位点的回顾性检查，更没有湿实验安全港验证。

### 产品缺口

1. 尚未进入正式 Slice 1 候选生成阶段；已实现的 Slice 1 引擎仅通过占位 Strain-B 数据验证了工程可运行性，未生成任何具备科学效力的 `candidate_window`、候选清单或排除清单。
2. 尚未完成多格式候选产物、规则版本、证据等级和独立验收的一致性闭环（Slice 1 自身的 acceptance/独立核验流程正在补齐，见 §7）。
3. 尚未完成完整 MVP 验收运行（冻结阈值、全基因组正式运行、人工抽查、文献回顾部分依赖真实数据）。
4. Streamlit（Slice 4）已实现并通过本地浏览器验收（`app/`）：数据导入/运行触发、Overview、Candidates、Excluded regions、Genome statistics、Data and methods 六个页面均可用；候选浏览页面在 Slice 2 恢复前如实展示 `collinearity_status = unavailable`。

因此，已完成的是参考与证据资格基础设施 + Slice 1 引擎工程实现，仍不是原始 MVP 本身、也不是任何具备科学效力的候选结果。

## 4. 当前路线的预期决策价值

**等待 Strain-T 真实测序数据**：决策价值高——直接解决项目最大的结构性缺口（没有目标菌株原生数据），且不再依赖存在身份疑问的公开数据。代价是交付时间和数据类型未定，无法现在排期依赖真实数据的部分（Slice 2、Slice 3 剩余部分）。

**建设 Slice 1/3 工具/4（占位数据）**：决策价值中高——不产生科学结论，但让工程管线（候选引擎、验收工具、展示与触发层）提前就绪；一旦真实数据到位，可以直接替换输入并很快跑出初步结果。成本可控，纯工程投入。

**暂停 Slice 2**：决策价值低——其核心价值主张（用 Strain-C 给代理菌株 Strain-B 佐证）会在真实数据到位后大幅贬值，现在投入工程实现的沉没成本风险高于收益。

**继续投入 SRR31989028/D0 诊断、Slice 0A/0B 延续调查**：决策价值低，已分别通过 ADR-0008/ADR-0009 暂停。

## 5. 已投入与后续预算

### 已投入的可观察证据

- Slice 0、Slice 0A、Slice 0B 已完成工程与验收，科学结论保持 `blocked`。
- `PRJNA604658` 与 `PRJNA1210090` 元数据/来源资格调查已关闭；`PRJNA1210090` 进一步核实为无关的工程菌株胁迫实验数据。
- `SRR31989028` 单重复受控探针已完整执行并正式失败（`probe_acceptance_status=failed`），验收与重复性证据齐备，存档于 `local_runs/controlled_probe/prjna1210090_srr31989028/`。
- Slice 1 候选窗口引擎已实现并测试（2026-07-16）：`src/pichia_safe_harbor/candidates.py`、`slice1.py`、CLI `candidate-windows`、11 项新测试，全部 104/104 通过；已用真实 Strain-B baseline 数据跑通一次端到端演示（4220 输入区间 → 4 个占位候选、4216 个排除）。

### 本轮授权：Slice 1/3 工具/4 建设（等待期工作）

- 允许继续完善 Slice 1 候选规则引擎（新规则、边界情况、区间切分逻辑、测试覆盖），前提是产物始终标记为 `engine_readiness_test_not_scientific_output`。
- 允许建设 Slice 1 自身的 acceptance/独立核验流程（比照 Slice 0/0A 已有模式），作为 Slice 3"验证多格式一致性和可重复性"目标中与数据来源无关的部分。
- 允许开始 Streamlit（Slice 4）实现：数据导入/运行触发页面 + 只读展示页面，遵循 ADR-0010 的架构约束。
- 不允许下载任何新的公开或私有测序数据。
- 不允许实现或维护 Slice 2（Strain-B-Strain-C 共线性核查）代码。
- 不允许把占位输出写入 `local_runs/` 并暗示其为正式运行，除非明确标注占位性质。
- 人力/机器预算：以正常开发节奏为准，无需专门预算表——这是纯代码工作，不涉及大规模计算或数据下载。

### 未授权预算

- 三重复 RNA-seq 验证：已通过 ADR-0008 撤销授权前提（原前提是"单重复探针通过"，但探针已失败），不再评估此预算。
- 继续寻找/资格化新的公开注释或转录组数据源（ADR-0004/0005/Slice 0A/0B 的延续工作）：已通过 ADR-0009 撤销，不再评估此预算。
- Slice 2 共线性核查实现：已通过 ADR-0011 排除，需要新 ADR 以不同理由重新论证才会重新评估。
- 正式证据整合边界轨道、Slice 1 正式候选生成、Slice 3 剩余部分（冻结阈值、全基因组正式运行、人工抽查、文献回顾）：均未分配预算，需在 Strain-T 测序数据到位并经 ADR 评估后重新估算。
- 付费云资源、外部采购：0，未授权。

## 6. 路线比较

| 路线 | 决策价值 | 成本与风险 | 结论 |
|---|---|---|---|
| 等待 Strain-T 真实测序数据 | 高；直接解决目标菌株数据缺口 | 时间不确定，当前无法排期依赖真实数据的部分 | 推荐（已授权安排测序） |
| 建设 Slice 1/3 工具/4，仅用占位 Strain-B 数据 | 中高；为真实数据到位后的快速产出做准备 | 成本可控，需持续标注非科学性质 | 推荐（ADR-0011） |
| 实现 Slice 2 共线性核查 | 低；核心价值主张随真实数据到来而贬值 | 沉没成本风险高于收益 | 不授权（ADR-0011） |
| 继续 D0 诊断 / 三重复验证 | 低；数据源本身证据价值有限 | 已产生的探针证据无法解决数据源的根本适用性问题 | 不授权（ADR-0008） |
| 继续寻找新的公开注释/转录组来源（Slice 0A/0B 延续） | 低；已系统性尝试 4 个来源均不合格 | 边际价值低，且会拖延对真实数据的等待 | 不授权（ADR-0009） |
| 立即用 Strain-B 占位数据产出"正式"候选 | 看似有产出，实则用不可靠边界证据冒充科学结论 | 违反 ADR-0004/0005 关于阈值/候选前置条件的规定 | 不授权 |
| 当前停止一切工作 | 避免新增投入 | 浪费等待测序数据的窗口期，框架代码本可以提前就绪 | 不采纳 |

## 7. 推荐路线与阶段授权

### 已授权：Gate P2，框架建设（Slice 1/3 工具/4）+ 测序等待期

授权内容：

- 跟踪 Strain-T 真实测序的安排进度（不由本项目直接执行，由项目负责人/上级安排）。
- 维护、扩展 Slice 1 候选窗口引擎的规则逻辑、测试与 CLI，包括补齐 acceptance/独立核验流程和区间切分逻辑；输入仍为现有 Strain-B Slice 0 baseline，产物必须保持 `engine_readiness_test_not_scientific_output` 标记。
- 开始 Streamlit（Slice 4）实现：应用入口、导航注册、Run Catalog、数据导入/运行触发页面、展示页面（概览、候选、排除记录、统计、数据与方法），遵循 ADR-0010 的架构约束。
- 保持既有权威证据链（Slice 0/0A/0B、SRR31989028 探针、PRJNA1210090 调查）完整、只读、不重跑不篡改。
- 维护 `docs/` 下的架构与规划文档，使其反映最新 ADR 决策（本次更新已完成）。

### 未授权事项

- D0 alignment 确定性诊断、修复 SRR31989028 可重复性问题；
- 下载 `SRR31989016`、`SRR31989027` 或任何新的公开/私有测序数据（在 Strain-T 真实数据到位前）；
- 继续寻找/资格化新的公开注释或转录组数据源（ADR-0009 已暂停的 Slice 0A/0B 延续工作）；
- 实现或维护 Slice 2（Strain-C 共线性核查）代码（ADR-0011 排除）；
- 正式 Slice 1 候选生成、阈值 ADR、候选排序；
- Slice 3 中依赖可信边界数据的部分（冻结默认阈值、运行完整基因组、人工浏览器抽查、文献位点回顾）；
- 付费云资源、外部采购、提交、推送或远端可见性变更（除非用户另行明确要求）。

## 8. 成功、失败与预算停止条件

### Gate P2 成功标志

- Slice 1 引擎（含 acceptance 流程、区间切分逻辑）与 Slice 4 Streamlit 代码可运行、有测试覆盖，且所有产物明确标记非科学性质。
- Streamlit 页面不重新实现科学规则、不提供编辑已有结果的功能；候选页面在 Slice 2 恢复前如实显示 `collinearity_status = unavailable`。
- 既有权威证据链未被意外修改或重跑；未出现 Slice 2 相关代码。
- 没有把占位候选清单误用为正式结论、对外沟通或写入 `local_runs/` 而未标注。

### Gate P2 失败/需要停下重新评估的情况

- 有人（包括 AI 会话自身）试图把 Slice 1/4 占位输出包装为正式候选或科学结论；
- 未经授权重新下载或分析公开 RNA-seq 数据，或未经授权实现 Slice 2；
- Streamlit 页面绕过核心引擎、自行实现规则，或提供编辑已有运行结果的功能。

### 预算停止线

- 出现上述任一失败情况时立即停止并回到本计划重新确认授权范围；
- 除工程维护外，本阶段不应产生显著人力/机器/存储/网络成本；一旦出现明显超出正常开发节奏的资源消耗，应停下检查原因。

## 9. 下一次规划门禁

下一次实质性规划决策绑定 **Strain-T 真实测序数据到位**：

1. 数据到位后，先由新 ADR 评估其类型、质量与对主坐标、注释边界、候选阈值推断的实际影响，不得跳过评估直接假定解除 `blocked`。
2. 评估通过后，再决定是否/如何将 Slice 1 引擎的输入从 Strain-B 占位数据切换为真实数据，是否恢复 Slice 2（以新理由论证），并决定候选长度、缓冲阈值等参数的正式取值（走独立阈值 ADR）。
3. 若测序长期未能交付或数据不足以解决核心缺口，回到"降级产品"或"停止路线"的项目级选择（不属于本计划授权范围）。

在真实数据到位并完成上述评估前，本计划不授权 Slice 2 恢复或 Slice 3 剩余部分的工作。
