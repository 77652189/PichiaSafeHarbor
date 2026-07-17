# ADR-0007：SRR31989028 normalized SRA acquisition 与受控探针结论

- 状态：Accepted
- 日期：2026-07-15
- 关系：修订 [ADR-0006](ADR-0006-capability-scoped-transcript-evidence.md) 中仅以 ENA paired FASTQ 作为 acquisition 门禁的实现方式；ADR-0006 的单重复范围、能力边界和停止线继续有效。
- 后续修订：[ADR-0008](ADR-0008-target-strain-native-sequencing-redirect.md) 因 Strain-T 原生测序获得授权而暂停本 ADR 原定的 D0 alignment 确定性诊断下一步，改为等待真实测序数据。

## 背景

ADR-0006 授权只处理 `SRR31989028`。ENA/Globus paired FASTQ 路径在当前网络下无法在合理时间内可靠完成，且多次中断会增加残留与身份混淆风险。NCBI 同时提供完整 normalized SRA，声明大小 `3931837205` bytes、MD5 `8c8a58880254746890e10c68b402c875`，不是质量值简化的 SRA Lite。

## 决策

1. acquisition 改为下载 NCBI AWS 的完整 normalized SRA，绑定 NCBI 元数据快照、URL、大小、MD5 和本地 SHA-256。
2. 必须先通过 SRA Toolkit `3.4.1` 的 `vdb-validate`，再以本地绝对路径和固定参数执行 `fasterq-dump --split-files`；禁止用裸 accession 静默联网，也禁止 SRA Lite。
3. 派生 FASTQ 固定 read defline、线程、pair/read/base 数和本地 SHA-256；压缩使用固定 Python 3.12 `gzip.GzipFile` / zlib、`compresslevel=6`、`mtime=0`、空文件名头。
4. acceptance 和独立 verifier 必须重新核对 acquisition evidence、SRA、元数据、工具/参数、FASTQ mate、pair 数与总碱基数。
5. 不下载 `SRR31989016` 或 `SRR31989027`；递归发现这两个 accession 的 raw 文件即失败。

## 正式探针结果

- Run ID：`srr31989028-probe-f439d29b82a7f8d6`。
- 权威验收：`local_runs/controlled_probe/prjna1210090_srr31989028/probe_v1_run1/acceptance_manifest.json`。
- 重复性证据：`local_runs/controlled_probe/prjna1210090_srr31989028/probe_v1_repeatability.json`。
- 状态：`execution_status=complete`、`verification_status=passed`、`scientific_acceptance_status=blocked`、`probe_acceptance_status=failed`。
- 21,957,417 对 reads 中，unique 21,382,912、multi 79,611、unmapped 79,057；四条核染色体非零覆盖和剪接支持通过，链方向保持 `unavailable`。
- UniVec 全局信号未触发，但核染色体高深度异常 bin 比例为 `0.0713141881071848`，超过既定 `0.05` gate，668个局部区域不得用于菌株特异边界推断。
- 16项主分析产物中14项逐文件 SHA-256 一致；`alignment.bam` 与 `alignment.bam.bai` 不一致，完整 SAM 流 SHA-256 也不同，因此逐文件重复性失败。

## 后果

- 本次单重复受控探针正式失败，项目继续 `blocked`。
- 不得自动下载另外两个重复，不得建立正式边界轨道、生成候选、冻结阈值或进入 Slice 1。
- 高深度异常的生物学/工程来源，以及 BAM 同坐标记录排序的确定性修复，交由规划侧决定是否形成后续切片；本 ADR 不隐式授权修复后重跑或扩大数据范围。
