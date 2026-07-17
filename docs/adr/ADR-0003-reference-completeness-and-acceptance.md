# ADR-0003：参考完整性、注释边界与验收状态

- 状态：Accepted
- 日期：2026-07-14

## 背景

Slice 0 真实数据运行发现：

- Strain-B 染色体1、3、4在 FASTA/GenBank 中标记为 `partial sequence`；染色体2标记为完整。
- assembly 元数据声明存在缺失的 3′ 端粒区域，但现有运行没有逐染色体、逐端点表达完整性。
- 当前 GFF3 的5288个 gene中有5032个带 `partial=true`，其中5031个同时带 `start_range` 和 `end_range`。
- 原始基线计算和独立区间核验已通过，但 `run_manifest.status=complete` 只能证明产物生成完成，不能证明数据适合阈值推断或 Slice 0 已科学验收。

NCBI 的 GFF3 说明明确指出：`partial=true` 表示 feature 在内部或一个/两个端点不完整；`start_range` 和 `end_range` 分别表示第4列和第5列边界为 partial。依据：<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/>。

## 决策

### 1. 组装端点完整性

- 每条序列分别记录5′和3′端点的 `complete / partial / unknown` 状态、来源和原始声明。
- 不从全局 assembly 备注推导未经证实的逐染色体端点事实。
- 与 `partial` 或 `unknown` 端点相邻的 terminal region 不得解释为真实端粒距离，也不得进入保守候选集合。

### 2. 注释边界不确定性

- 解析并保留 `partial`、`start_range`、`end_range` 及等价属性。
- `IntergenicRegion` 继承左右功能实体的边界置信状态。
- 提交坐标可以用于可重复的描述性基线统计，但必须与高置信边界子集分开报告。
- 任一边界不确定的区间在获得独立注释支持或后续 ADR 明确处理前，不得进入保守候选集合，也不得直接用于冻结默认缓冲阈值。
- 由于当前不确定性具有系统性，不在本 ADR 中武断选择替代注释；Slice 0 必须先形成注释适用性报告，再决定更换、补充、跨参考协调或保守排除。

### 3. 状态模型

运行状态拆分为：

- `execution_status`：`pending / running / complete / failed`；
- `verification_status`：`not_run / passed / failed`；
- `scientific_acceptance_status`：`pending / accepted / blocked / rejected`。

不得继续用单一 `status=complete` 表示科学完成。当前两次 Strain-B 运行应解释为：计算完成、独立验证通过、科学验收阻塞。

### 4. 验收清单

Slice 验收必须生成独立 `AcceptanceManifest`，绑定：

- 实际 FASTA/GFF 哈希；
- 实际序列映射和序列分类；
- 实现哈希、运行 ID 和全部产物哈希；
- 测试和独立核验结果；
- 注释与组装完整性摘要；
- 科学验收状态、阻塞原因和验收版本。

Run ID 必须由实际输入哈希、参考身份、序列映射、序列分类、实现版本和影响结果的配置共同决定。

## 后果

- 现有4220个原始区间仍是有效的“按提交坐标计算的描述性基线”，但不是阈值 ADR 或正式候选资格的充分依据。
- Slice 0 保持未验收状态，先补完元数据、边界传播、状态模型和 acceptance manifest。
- Streamlit 可展示被阻塞运行的诊断信息，但不能把它列为可用于正式候选浏览的已验收运行。
- 如果高置信边界覆盖不足，项目必须评估替代或补充注释来源，而不是通过默认规则掩盖问题。
- 当前高置信覆盖不足后的路线由 [ADR-0004](ADR-0004-annotation-qualification-route.md) 决定。

## 验收

- 四条核染色体均有逐端点完整性记录和来源。
- 5288个 gene 的边界限定属性被保留并汇总。
- 每个原始区间具有左右边界置信状态。
- 基线报告分别给出全部区间和高置信边界子集统计。
- 运行清单使用三层状态，旧 `status` 不再承担科学验收语义。
- 核心入口强制校验实际输入哈希。
- Run ID 包含实际输入、映射和分类。
- 独立核验进入正式 acceptance manifest。
- 参考下载以整个 bundle 原子替换，失败不会留下混合版本目录。
- 重跑测试、双运行一致性和独立核验后，才能重新判断 Slice 0 是否可验收。
