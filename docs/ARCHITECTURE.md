# PichiaSafeHarbor 当前架构

## 1. 架构目标

系统以可追溯、可重复的离线分析管线处理完整 Strain-B 核基因组，生成面向 Strain-T 的推定安全港候选，并使用 Strain-C 做跨参考一致性核查。

第一阶段采用"离线分析核心 + Streamlit 展示与触发层"的架构（[ADR-0010](adr/ADR-0010-streamlit-import-and-trigger-workflow.md) 修订自最初的纯只读设计）。分析核心是唯一的科学规则实现；Streamlit 可以导入新参考数据并触发核心引擎生成新的版本化运行，也可以读取并展示已经完成的运行，但不得在页面代码中另行实现或绕过核心规则。架构必须允许后续加入 Strain-T 原生数据，而不改变现有坐标和证据语义；当前使用的 Strain-B 数据是开发/测试占位（[ADR-0009](adr/ADR-0009-deprioritize-proxy-data-compensation.md)），不是长期科学坐标策略。

## 2. 核心不变量

1. Strain-B `GCA_001746955.1` 是当前唯一主坐标空间，现阶段定位为开发/测试占位，等待真实 Strain-T 数据到位后由新 ADR 重新评估（ADR-0009）。
2. Strain-C `GCA_000223565.1` 只提供辅助证据，不覆盖主坐标。
3. Strain-T 是目标菌株，但第一阶段结果不是 Strain-T 精确坐标。
4. 原始基因间区与候选窗口是不同实体。
5. `exclude / flag / annotate` 是规则动作，不是证据等级。
6. 安全风险、结构稳定性、表达环境和证据完整度分别输出。
7. 缺失数据必须表示为 `unavailable`，不能自动解释为安全。
8. 任一结果都能追溯到输入文件、校验值、规则版本和参数快照。
9. 未定位 scaffold、线粒体和坐标不可靠区域默认不产生保守推荐候选。
10. 文献整合成功不等同于安全港验证。
11. Streamlit 页面不得拥有或复制科学规则；页面可以调用核心引擎的顶层编排函数触发新的版本化运行，但只能消费/展示版本化结果契约，不得在页面代码中另行实现规则，也不得提供直接编辑结果的功能。
12. 页面状态和缓存必须绑定运行 ID，切换运行时不得复用旧运行数据。
13. 计算完成、验证通过和科学验收是三个独立状态。
14. 参考序列完整性和注释边界不确定性必须逐记录传播，不能压缩成一个全局备注。

## 3. 逻辑数据流

```text
Reference bundle
  -> input validation
  -> annotation normalization
  -> protected-feature union
  -> complete-genome intergenic inventory
  -> rule-based interval subtraction
  -> candidate windows
  -> Strain-B/Strain-C concordance evidence
  -> risk and evidence records
  -> TSV / JSON / report / run manifest
  -> validated run index
  -> Streamlit presentation (read existing runs, or import data and trigger a new run)
```

## 4. 逻辑模块

### 4.1 Reference Bundle

负责声明和锁定：

- FASTA、GFF3/GTF 和可选风险轨道；
- 主参考与辅助参考身份；
- 下载来源、版本、校验值和许可证；
- 序列名称映射；
- 已知组装缺口、未定位 scaffold 和端粒缺失。
- 每条序列的5′/3′端点完整性和证据来源；
- 注释边界限定属性及其解释策略。

该模块只描述数据，不包含候选筛选逻辑。

### 4.2 Input Validation

在任何候选计算前验证：

- 文件可读取且校验值匹配；
- FASTA 序列与注释序列名称可一一映射；
- 坐标不越界；
- 注释版本与 assembly 契约一致；
- 核染色体、线粒体和未定位 scaffold 分类明确。
- 实际文件哈希与 reference manifest 匹配；
- 序列完整性分类覆盖所有输入序列。

验证失败时停止分析，不生成部分候选结果。

### 4.3 Annotation Normalization

将不同来源的基因、转录本、ncRNA、伪基因、重复序列和结构轨道转换为统一的内部区间模型。内部坐标约定必须唯一；导入和导出时显式转换。

同一功能实体的父子注释需要归并，避免把外显子间隔错误识别成基因间区。

规范化层必须保留 `partial`、`start_range`、`end_range` 和等价限定属性，并产生统一的边界置信模型。不得因为属性在大多数记录中出现而静默忽略。

### 4.3.1 Annotation Evidence Qualification

科学边界来源与坐标 assembly 是两个独立身份。主坐标保持 Strain-B `GCA_001746955.1`，但进入正式候选分析前，必须对可用注释来源进行资格评估。

评估层负责：

- 锁定候选注释的数据来源、版本、许可证、校验值和原始坐标空间；
- 区分独立注释与同一 INSDC 注释的镜像或再发布；
- 统计 gene、转录本、ncRNA、边界限定属性和缺失类型；
- 评估与主坐标的直接兼容性或显式映射需求；
- 比较 gene 身份、方向、边界和邻接关系的一致、冲突及不可映射情况；
- 记录 RNA-seq、蛋白组、人工整理或其他实验支持；
- 输出 `qualified / conditionally-qualified / unqualified / unavailable`，但不在评估代码中自行选择最终来源。

只有后续 ADR 接受的注释或共识边界轨道才能进入阈值推断和 `Candidate Window Engine`。

### 4.4 Protected Feature Model

根据版本化配置，将注释映射为：

- 硬排除区；
- 风险缓冲区；
- 仅注释区。

保护区由区间集合运算产生，不允许在报告层重新推断。

### 4.5 Intergenic Inventory

基于受保护功能实体的外边界，按染色体排序并生成完整的 `intergenic_region` 清单，记录邻近基因、方向、距离和染色体端点状态。

每个区间还必须记录左右边界置信度。描述性基线可以使用提交坐标，但高置信阈值统计和正式候选资格必须能排除或单独处理不确定边界。

该阶段只描述基因组，不应用候选阈值。其统计结果是后续冻结默认阈值的输入。

### 4.6 Candidate Window Engine

从每个 `intergenic_region` 中减去硬排除区和保护缓冲区。每个剩余连续片段生成独立 `candidate_window`。

引擎必须保留：

- 父原始区间；
- 每次扣除的规则和来源；
- 切分前后坐标；
- 被排除片段及理由。

### 4.7 Concordance Evidence

将 Strain-B 候选周边锚点与 Strain-C 对应区域进行核查，输出：

- `confirmed`：一对一共线，邻近基因及方向一致；
- `conflicting`：存在结构冲突、邻近关系变化或无法可靠映射；
- `unavailable`：缺少映射或证据不足。

`conflicting` 默认进入风险标记；是否升级为硬排除由后续 ADR 决定。

### 4.8 Evidence and Reporting

报告层只组合已有事实，不重新计算科学规则。输出包括候选、排除记录、统计摘要、数据缺口和运行清单。

TSV、JSON 和人类可读报告必须由同一规范化结果对象生成，避免格式间漂移。

### 4.9 Run Catalog

运行目录只有在以下条件全部满足后才能进入可展示目录：

- `RunManifest` 存在且 schema 有效；
- 所有必需输出存在且校验值匹配；
- TSV、JSON 和报告的关键数量与坐标一致；
- `execution_status = complete`；
- `verification_status = passed`；
- 当前页面所需的科学阶段已达到相应验收状态；
- 未发现未知或不支持的 schema 版本。

Run Catalog 只建立只读索引，不修改运行产物；触发新运行只会在目录中新增一条记录，不会改写已有记录。`blocked` 运行可以进入审计或基线诊断页面，但不得进入正式候选页面，也不得被误标为可供后续科学决策使用——不管这次运行是命令行触发还是网页触发。

### 4.10 Streamlit Presentation

Streamlit 通过稳定的结果读取服务访问 Run Catalog 和规范化结果来展示已完成的运行；它也可以调用核心引擎的顶层编排函数（`reference.py` 的校验函数、`pipeline.run_baseline`、`slice1.run_slice1`）导入新参考数据并触发新的运行（[ADR-0010](adr/ADR-0010-streamlit-import-and-trigger-workflow.md)）。无论哪种路径，页面都不直接读取内部分析临时文件，也不得绕开顶层编排函数、直接调用或重新实现候选枚举、区间相减、规则判定等内部逻辑。

首版页面结构：

```text
app/streamlit_app.py              应用入口与全局运行选择
pages / navigation registry
  Overview                       运行摘要与限制声明
  Data import & run trigger      导入参考数据、触发 baseline/candidate-windows 运行、查看执行状态
  Candidates                     候选筛选与详情
  Excluded regions               排除记录与规则理由
  Genome statistics              完整基因组基线图表
  Data and methods               数据、规则、坐标和校验信息
```

项目可以采用 Streamlit 原生多页或显式导航注册，但必须只有一个权威页面注册表。页面入口、导航名称、页面模块和测试必须同步维护。

展示类页面只调用展示服务；导入/触发页面调用核心编排函数：

```text
Streamlit page (展示)
  -> presentation service
  -> validated result models / Run Catalog
  -> immutable run artifacts

Streamlit page (导入/触发)
  -> core orchestration functions (pipeline.run_baseline / slice1.run_slice1 / ...)
  -> new versioned run directory + RunManifest
  -> Run Catalog (只读索引更新)
```

禁止页面绕开顶层编排函数，直接导入候选枚举器、区间相减引擎或规则执行器来临时拼装结果；禁止页面提供编辑已有运行候选、排除理由或证据等级的功能。

### 4.11 Streamlit 状态与缓存

所有 `st.session_state` 键必须集中声明、提供安全默认值，并至少包括当前运行 ID、候选筛选条件和当前候选选择。状态键需要命名空间，避免与未来页面冲突。

切换运行时必须清理或重置：

- 当前候选选择；
- 依赖旧运行字段的筛选条件；
- 分页或详情状态；
- 与旧运行相关的派生视图。

缓存只用于读取不可变运行文件和构建只读视图。缓存键必须包含运行 ID、manifest 校验值或等价版本标识。不得缓存可变筛选状态，不得用无版本的全局缓存跨运行复用结果。

### 4.12 Streamlit 错误边界

- 无可用运行时展示操作说明，不伪造示例候选。
- 运行不完整、schema 不支持或校验失败时阻止进入正常页面。
- 单条候选缺少非关键证据时显示 `unavailable`，不得让整个页面崩溃。
- 页面不得静默回退到其他运行或其他 assembly。
- 下载内容必须来自已验证结果或当前筛选视图，并标记运行 ID。
- 触发的运行执行失败、验证失败或科学状态被阻塞时，页面必须如实展示 `execution_status`/`verification_status`/`scientific_acceptance_status`，不得展示为成功或默认跳转到候选页面。
- 长时间运行的触发任务必须有明确的"运行中/完成/失败"状态反馈，不得让页面表现为假死或假成功。

## 5. 核心数据实体

### `ReferenceIdentity`

描述目标菌株、主参考、辅助参考、assembly、annotation 和适用性。

### `GenomicInterval`

统一表示染色体、半开区间坐标、链、来源和坐标空间。

### `IntergenicRegion`

描述原始基因间区、左右邻近功能实体、左右边界置信度和序列端点完整性。

### `CandidateWindow`

描述可供后续实验设计使用的连续候选区间，并引用其父区间。

### `RuleDecision`

包含规则 ID、版本、动作、参数、理由、来源和命中的区间。

### `EvidenceRecord`

描述证据类型、来源、适用菌株、等级及 `available / unavailable` 状态。

### `RunManifest`

锁定一次计算运行的实际输入文件、校验值、序列映射、序列分类、软件版本、规则快照、参数、输出清单和 `execution_status`。它不单独证明科学验收。

### `AcceptanceManifest`

绑定运行 ID、实际输入和产物哈希、自动化测试、独立核验、完整性检查、科学验收状态、阻塞原因和验收版本。该文件是切片验收的权威机器可读记录。

### `RunCatalogEntry`

描述一个可供网页读取的完成运行，包括运行 ID、完成状态、manifest 校验值、schema 版本和产物位置。

### `CandidateView`

面向展示层的只读候选模型。它可以格式化已有字段，但不得产生新的科学判断。

## 6. 坐标契约

- 内部统一使用 0-based、半开区间 `[start, end)`。
- GFF3/GTF 导入时从 1-based 闭区间转换。
- 面向生物学用户的报告可以显示 1-based 闭区间，但必须同时声明约定。
- TSV 和 JSON 字段名必须明确包含坐标约定或由 schema 固定声明。
- 不同 assembly 的坐标不得直接比较，必须经过显式映射并保留映射证据。

## 7. 规则与参数治理

- 每条规则有稳定 ID、版本、动作和理由模板。
- 默认阈值必须来自完整基因组基线统计和独立 ADR。
- 阈值统计必须区分全部提交坐标与高置信边界子集。
- 用户覆盖默认阈值时，覆盖值进入 `RunManifest`。
- 规则升级不得静默改变旧结果语义；结果必须记录所用规则版本。
- 第一阶段不产生不可解释的综合安全总分。

## 8. 运行产物边界

建议边界：

```text
docs/          权威需求、架构、ADR 和 handoff
src/           后续实现代码
tests/         自动化测试
testdata/      小型、固定、可分发的测试数据
reference/     数据清单与下载说明，不默认提交大型原始数据
local_runs/    完整基因组运行产物、缓存和人工核查证据
app/           Streamlit 入口、页面、展示服务和数据导入/运行触发服务
```

大型参考数据和运行产物不应混入源代码或固定测试数据。

## 9. 实施切片

### Slice 0：参考数据与基线统计

- 锁定 Strain-B/Strain-C 数据文件和校验值；
- 验证坐标与序列契约；
- 输出完整核基因组的基因间区和方向分布；
- 传播逐序列端点完整性和注释边界置信度；
- 生成三层状态和正式 acceptance manifest；
- 不产生最终推荐阈值。

Slice 0 已完成工程与验证目标，其科学结果为：当前 GenBank 注释不适合作为唯一阈值边界来源，运行保持 `scientific_acceptance_status = blocked`。

### Slice 0A：注释来源资格评估

- 调查并锁定可获取的独立 Strain-B/Strain-C 注释资源；
- 优先调查2016年 PichiaGenome RNA-seq/蛋白组整理注释及其持久归档；
- 评估旧 Strain-B RefSeq 注释，但不切换主 assembly；
- 证明 Ensembl 或其他资源是否独立于当前 INSDC 注释；
- 对候选来源统计 partial 边界、功能实体和坐标空间；
- 评估映射到 `GCA_001746955.1` 的可行性和冲突；
- 输出注释资格报告和下一份来源选择 ADR；
- 不冻结阈值，不生成正式候选。

Slice 0A 已完成工程与验证目标；结论是当前可获取的注释来源（旧 Strain-B RefSeq、Ensembl、KEGG 等）均不能作为独立边界证据，促成 ADR-0005 转向"实验数据支持的边界轨道"路线。

### Slice 0B：实验数据资格与重注释可行性门禁

- 依据 [ADR-0005](adr/ADR-0005-evidence-integrated-boundary-track.md) 的要求，在建立全基因组证据整合边界轨道前，先评估可获取的实验（RNA-seq/蛋白组）证据是否合格；
- 先后评估 `PRJNA311606`（与当前 assembly/注释生成过程重叠，非独立证据）、`PRJNA604658`（eGFP 工程背景、单次运行、无链特异性声明）、`PRJNA1210090`（后经 [ADR-0008](adr/ADR-0008-target-strain-native-sequencing-redirect.md) 核实为无关的工程菌株胁迫转录组、身份未确认）、`SRR3201093`/`PRJNA304976`（ATCC 亚株编号与文献不一致、无逐文件许可、链方向约定不可用）共 4 个候选来源；
- 跑了 5 轮资格评估（`local_runs/slice0b/qualification_v1_run1` 至 `qualification_v5_run2`）；
- 不修改 Slice 0/0A 既有证据，不生成候选，不冻结阈值，不实现 Streamlit。

Slice 0B 已完成工程与验证目标，正式结论（`local_runs/slice0b/qualification_v5_run1/slice0b_recommendation.md`）为："evidence is insufficient to authorize a full-genome boundary track... Do not start full-genome reannotation or Slice 1."——已系统性评估过的公开数据源均不足以支撑全基因组边界轨道。[ADR-0009](adr/ADR-0009-deprioritize-proxy-data-compensation.md) 据此暂停继续寻找/资格化新公开数据源的延续工作；已完成的 5 轮评估结果作为历史记录保留，结论继续有效，等待真实 Strain-T 数据到位后再重新评估证据整合边界轨道是否值得建设。

### Slice 1：原始区间与候选窗口

- 建立统一注释模型；
- 枚举原始基因间区；
- 实现 `exclude / flag / annotate`；
- 输出候选和排除记录。

引擎已实现（2026-07-16，`src/pichia_safe_harbor/candidates.py` + `slice1.py`，CLI 子命令 `candidate-windows`），并以现有 Strain-B 注释作为占位输入验证端到端可运行性：边界置信度非 `high` 或缓冲后窗口不足最小长度的区间硬排除；方向非 `convergent` 或区间异常长仅标记、不排除。`buffer_distance_bp`/`min_candidate_window_bp` 目前仍是调用方显式传入的示意占位值，不是 ADR 冻结阈值；产物统一标记 `run_purpose=engine_readiness_test_not_scientific_output` 且 `scientific_acceptance_status=blocked`。这属于 [ADR-0008](adr/ADR-0008-target-strain-native-sequencing-redirect.md)/[ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 允许的"引擎工程自测"，不是 ADR-0005 所指的正式 Slice 1 候选生成——正式候选、阈值冻结仍需等待合格边界证据（见 ADR-0008）。引擎已支持把一个基因间区切分为多个独立候选窗口（`extra_exclusion_zones` 参数，`CandidateWindow.split_index`/`split_count` 字段），已用合成测试数据验证；目前没有真实的重复序列/着丝粒等风险轨道数据可以接入，接入方式留待相应数据到位后再定。补齐了与 Slice 0/0A 同构的 acceptance/独立核验流程（`create_slice1_acceptance`，CLI `accept-candidate-windows`）。

**候选窗口序列（2026-07-16 补充）：** REQUIREMENTS.md 4.7 要求每个候选包含"候选窗口长度**和序列**"，早期实现只有长度。`run_slice1` 现在多接收一个 `fasta_path` 参数（`parsers.py` 新增 `read_fasta_sequences`，独立于 `read_fasta_index`，避免改动已测试的 baseline 路径），按 `CandidateWindow.seqid/start/end` 切出实际序列，写入 `CandidateWindow.sequence` 字段（`candidate_windows.json`/`.tsv` 新增 `sequence` 列）以及新增的 `candidate_windows.fasta` 产物（与其余产物一样纳入 `run_manifest.json` 的哈希校验）。为防止候选窗口坐标和序列来自不一致的参考基因组版本，`run_slice1` 会先比对传入 fasta 的 sha256 与父 baseline 运行记录的 `inputs.fasta.sha256`，不一致则拒绝执行（`ContractError`）。CLI `candidate-windows` 子命令和 Streamlit 触发页相应新增 `--data-dir`/参考基因组目录输入。这填补的是"候选=坐标+序列"这一层，必需基因、重复序列、着丝粒/端粒、复制起点等风险轨道仍然缺失（见本节前段与 `missing_risk_tracks`），候选窗口依旧不是已验证的安全港。

### Slice 2：双参考核查（暂停实现）

- 建立 Strain-B 到 Strain-C 的共线性证据；
- 输出 `confirmed / conflicting / unavailable`；
- 复核候选排序和风险解释。

[ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 指出：本切片的原始理由（用 Strain-C 为代理菌株 Strain-B 佐证）本身是另一种"缺失 Strain-T 数据的补偿工作"，一旦真实 Strain-T 数据到位，这个交叉印证的必要性会大幅下降甚至消失。因此明确暂停实现，不作为当前可用占位数据构建的框架代码；如果未来要以"与具体菌株无关的结构稳定性核查"这一新理由重新论证其价值，需要另立 ADR。

### Slice 3：完整 MVP 验收（部分授权）

依赖可信边界数据、须等待真实 Strain-T 数据到位并经独立 ADR 评估的部分：

- 冻结默认阈值 ADR；
- 运行完整基因组；
- 完成人工浏览器抽查和文献位点回顾。

与数据来源无关、已获 [ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 授权现在建设的部分：

- 验证多格式一致性和可重复性（TSV/JSON/报告之间的一致性核验、重复运行确定性核验），可比照 Slice 0/0A 已有的 acceptance/独立核验模式在 Slice 1 上先行落地。

### Slice 4：Streamlit 数据导入、触发与展示（已授权）

[ADR-0011](adr/ADR-0011-broaden-framework-authorization-exclude-slice2.md) 正式授权开始实现，遵循 [ADR-0010](adr/ADR-0010-streamlit-import-and-trigger-workflow.md) 的架构约束：

- 建立唯一应用入口和权威导航注册；
- 建立 Run Catalog 和展示服务；
- 建立数据导入/运行触发页面，调用核心编排函数生成新的版本化运行（ADR-0010）；
- 实现概览、候选、排除记录、统计与数据方法页面；
- 集中管理 `st.session_state` 和版本化缓存；
- 验证页面展示与 TSV/JSON 一致；
- 通过本地启动、健康检查和浏览器验收确认真实页面行为，包括触发运行的执行状态反馈。

页面展示的候选/共线性状态在 Slice 2 恢复实现前应显示为 `unavailable`，不得伪造或跳过该字段。

实现已完成（2026-07-16）：`app/streamlit_app.py`（`st.navigation` 单一权威页面注册表）、`app/services/{catalog,presentation,trigger}.py`、`app/state.py`、`app/common.py`、`app/pages/*.py`。已通过本地浏览器验收：触发一次真实 Slice 0 baseline 和一次 Slice 1 candidate-windows 运行，六个页面（Overview、Data import & run trigger、Candidates、Excluded regions、Genome statistics、Data and methods）均正确展示三态状态、`collinearity_status=unavailable`、`run_purpose` 等字段；触发流程只调用 `pipeline.run_baseline`/`slice1.run_slice1`，未复制或绕过规则；未提供编辑已有结果的功能。Run Catalog 默认扫描 `local_runs/streamlit_runs/`（可在 Overview 页面覆盖）。`app/` 依赖 `streamlit`、`pandas`（见 `pyproject.toml` 的 `[project.optional-dependencies].app`），`src/pichia_safe_harbor` 本身继续保持零依赖。

## 10. 暂缓决策

- 默认基因缓冲距离和最小候选窗口长度；
- 系统性 `partial` 注释最终采用替代注释、跨参考协调还是保守排除；
- 共线性冲突是否一律硬排除；
- Strain-T 原始重测序数据的纳入方式；
- 表达环境证据的具体数据源；
- 是否提供远程部署或 API。

（"是否允许网页触发分析任务或编辑规则"已由 ADR-0010 决定：允许触发核心引擎生成新运行，不允许编辑已有结果。）

这些事项不得在实现中以隐含默认值先行决定。
