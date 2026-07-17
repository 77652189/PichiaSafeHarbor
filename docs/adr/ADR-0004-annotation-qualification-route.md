# ADR-0004：注释来源资格评估路线

- 状态：Accepted
- 日期：2026-07-14
- 后续修订：[ADR-0009](ADR-0009-deprioritize-proxy-data-compensation.md) 因 Strain-T 原生测序获得授权而暂停本 ADR 定义的注释来源资格评估延续工作（继续寻找 PichiaGenome 原始文件、深入比较 Ensembl/KEGG 等）；已完成的 Slice 0A 调查结果作为历史记录保留。

## 背景

Slice 0 补完已经通过工程验证：4220个原始区间可重复，三层状态、边界传播、AcceptanceManifest、实际输入哈希和独立核验均已建立。

科学结果显示：

- 当前 Strain-B GenBank GFF3 中5032/5288个 gene 为 partial；
- 高置信双边界区间仅28/4220（0.006635）；
- 28个区间中只有4个 convergent、7个 divergent、17个 tandem，不足以代表完整基因组的边界与方向分布；
- 当前本地 Strain-C GenBank GFF3 同样有5525/5527个 gene 为 partial；
- Ensembl Fungi 声明其 Strain-B 注释来自 `GCA_001746955.1` 的 INSDC 提交，因此不是当前问题的独立修复来源。

2016年的 PichiaGenome 整理研究使用 RNA-seq 和蛋白组改进 Strain-C 注释，报告5325个 ORF、492个新增 ORF、341项错误 ORF 修正和175项删除，具有优先调查价值。论文：<https://pubmed.ncbi.nlm.nih.gov/27388471/>。其原网站当前不可访问，因此数据可获得性和持久来源尚需验证。

## 决策

1. 不接受当前28个高置信区间作为缓冲距离、最小窗口长度或正式候选规则的依据。
2. 保留 Strain-B `GCA_001746955.1` 作为主坐标空间，不切换到2009年的旧 Strain-B assembly。
3. 在 Slice 1 前新增 Slice 0A“注释来源资格评估”。
4. 优先寻找 PichiaGenome 整理注释的论文附件、机构仓储、公共归档或其他可持久获取副本。
5. 将旧 Strain-B RefSeq `GCF_000027005.1` 作为独立候选注释和坐标比较来源，而不是新的主坐标。
6. 当前 Strain-B/Strain-C GenBank GFF 和 Ensembl 镜像仅作为基线或同源来源，不计为相互独立证据。
7. 如果没有单一注释能够直接用于主 assembly，则评估建立版本化“共识边界轨道”：将合格注释显式映射到主坐标，并保留映射、冲突和不可用证据。
8. Slice 0A 只产出资格报告和来源选择 ADR，不实现候选窗口或阈值。

## 放弃的方案

### 接受28个高置信区间

样本量和方向组成严重偏斜，会使阈值只反映少量非蛋白编码或特殊边界，不能代表完整基因组。

### 切换到旧 Strain-B RefSeq assembly

会放弃已锁定的长读长主坐标，引入坐标迁移和旧组装完整性风险。旧 RefSeq 应用于比较，不应直接替代主 assembly。

### 直接使用 Strain-C GenBank GFF

其 partial 问题与当前 Strain-B 注释同样系统性，且坐标属于另一 assembly，不能直接作为主边界。

### 立即进行全基因组从头注释

范围和依赖显著扩大，且在现有整理注释尚未完成资格调查前没有必要。仅在公共整理数据不可获得或不合格时重新决策。

## Slice 0A 产出

- `annotation_sources.json`：来源、版本、许可证、哈希、坐标空间和独立性；
- `annotation_qualification.tsv/json`：实体数量、partial 比例、边界类型和证据等级；
- `annotation_mapping_probe.tsv/json`：到主坐标的映射覆盖、方向、边界偏移、冲突和不可映射记录；
- 注释资格报告：比较候选路线并提出推荐；
- 新 ADR：选择正式边界来源或共识边界方案。

## 验收

- 至少评估当前 INSDC 注释、旧 Strain-B RefSeq 和一个 RNA-seq/人工整理候选；若第三项不可获得，需提供可复现的获取失败证据。
- 明确识别同源镜像，不能把同一注释重复计算为独立支持。
- 每个来源记录许可证、版本、哈希和原始坐标空间。
- 映射探针覆盖全部核染色体，并输出一致、冲突、不可映射和边界不确定类别。
- 不修改主 assembly，不冻结阈值，不生成正式候选。
- 若没有来源通过资格评估，科学状态继续 `blocked`，并单独规划公共 RNA-seq 支持的重新注释路线。
