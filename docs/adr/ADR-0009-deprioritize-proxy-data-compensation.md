# ADR-0009：停用代理数据补偿路线，围绕真实 Strain-T 数据重新规划

- 状态：Accepted
- 日期：2026-07-16
- 关系：在 [ADR-0008](ADR-0008-target-strain-native-sequencing-redirect.md)（暂停 ADR-0006/0007 的公开 RNA-seq 路线）基础上扩大范围；将 [ADR-0001](ADR-0001-reference-strategy.md) 的 Strain-B 主坐标策略重新定位为开发/测试占位；暂停 [ADR-0004](ADR-0004-annotation-qualification-route.md) 的注释来源资格评估延续工作与 [ADR-0005](ADR-0005-evidence-integrated-boundary-track.md)/Slice 0B 的证据整合边界轨道路线。三者的历史决策与既有证据继续作为审计存档保留，不删除、不改写。

## 背景

ADR-0008 已经暂停了 ADR-0006/ADR-0007 定义的、依赖 `PRJNA1210090` 单一公开数据源的证据资格化路线。但复盘发现，ADR-0001、ADR-0004、ADR-0005 三者本质上是同一类工作——"没有 Strain-T 原生数据，只能在 Strain-B 这个近缘代理上想办法凑证据"——不只是 ADR-0006/0007 这一条 RNA-seq 子路径。

回顾 Slice 0B 的实际产出：先后评估了 `PRJNA311606`（与当前 assembly/注释生成过程重叠，非独立证据）、`PRJNA604658`（单次运行、eGFP 工程背景、无链特异性声明）、`PRJNA1210090`（经 ADR-0008 核实，实际是无关的工程菌株胁迫转录组，身份未确认）、`SRR3201093`/`PRJNA304976`（ATCC 亚株编号与文献不一致、无逐文件许可、链方向约定不可用）共 4 个候选来源，跑了 5 轮资格评估，最终正式结论是：

> "evidence is insufficient to authorize a full-genome boundary track... Do not start full-genome reannotation or Slice 1."

这说明"通过公开数据在 Strain-B 上拼凑证据"这条路线已经被相当充分地尝试过，且没有成功，不是偶然一次没做好。

ADR-0004 的注释来源调查（旧 Strain-B RefSeq、Ensembl、KEGG、2016 年 PichiaGenome 整理注释）同样是为了在 Strain-B 上找到更可信的边界来源，目前仍停留在"PichiaGenome 原网站不可访问、持久来源未验证"的状态，继续推进的边际价值同样有限。

2026-07-16，Strain-T 原生测序已获项目负责人授权（长期项目）。真实数据到位后，"如何让 Strain-B 代理更可信"这整套问题会被"直接用真实数据重新做注释/边界判断"取代，不需要再追问"这份代理数据到底够不够"。

## 决策

1. 将 ADR-0001 的 Strain-B `GCA_001746955.1` 主坐标策略重新定位为**开发与测试占位**，不再是长期科学坐标策略——它继续用于支撑 Slice 0/Slice 1 等框架代码的开发、集成测试和端到端演示，但项目不再投入精力"提升"这份代理数据本身的可信度。
2. 暂停 ADR-0004 定义的注释来源资格评估延续工作（继续寻找 PichiaGenome 原始文件、深入比较 Ensembl/KEGG 等）。已完成的 Slice 0A 调查结果作为历史记录保留，结论不变。
3. 暂停 ADR-0005 定义的证据整合边界轨道目标与 Slice 0B 后续工作（继续评估新的公开转录组/蛋白组来源）。已完成的 Slice 0B 五轮资格评估作为历史记录保留，其"证据不足"结论继续有效，不需要用新的公开数据源去推翻或补充。
4. 明确以下工作**不因本 ADR 而停止**，因为它们与具体数据来源无关，是可直接复用的框架能力：参考数据校验管线（`reference.py`）、注释归一化与基因间区清单（`analysis.py`/`pipeline.py`）、三态验收模型与 AcceptanceManifest 规范（ADR-0003）、Slice 1 候选窗口规则引擎（`candidates.py`/`slice1.py`）。这些代码换成真实 Strain-T 数据后应当可以直接复用或只需小幅调整。
5. 真实 Strain-T 数据到位后，须由新 ADR 评估：是否切换主坐标、是否需要新的注释归一化规则、是否需要重新设计边界置信度模型；不得默认沿用 Strain-B 时期的具体数值结论（如高置信区间数量、占比等）。

## 放弃的方案

### 继续追加投入寻找第 5、第 6 个公开转录组/蛋白组来源

Slice 0B 已经系统性尝试了 4 个来源并给出"证据不足"的正式结论，没有新信息表明下一个来源会不一样。

### 现在就切换主坐标到旧 Strain-B RefSeq 或其他 assembly

会引入新的坐标迁移成本，而且预期不久会被真实 Strain-T 数据取代，不值得中途切换。

### 保留 ADR-0001/0004/0005 原状不动，假装它们仍是当前最优路线

会让文档继续引导团队追加投入一个已经被证明收益很低的方向，与项目一贯的证据纪律相悖。

## 后果

- `docs/adr/README.md` 中 ADR-0004、ADR-0005 状态标注为"Superseded by ADR-0009"；ADR-0001 保留 Accepted（其"未来用新 ADR 替代"的预留条款正是本 ADR 触发的情形）。
- `docs/REQUIREMENTS.md`、`docs/ARCHITECTURE.md` 中涉及 Strain-B 长期策略、注释来源改善、证据整合边界轨道的表述需要更新为"开发占位、等待真实数据"的框架。
- `docs/EXECUTION_PLAN.md` 的路线比较和预算部分不再为"继续寻找/资格化新公开数据源"分配预算。
- `scientific_acceptance_status` 继续保持 `blocked`，不因本 ADR 而改变。
