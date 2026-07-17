# 架构决策索引

仅记录长期、高影响且存在真实权衡的决策。当前有效决策如下：

| ADR | 状态 | 决策 |
| --- | --- | --- |
| [ADR-0001](ADR-0001-reference-strategy.md) | Accepted（重定位见 ADR-0009） | 以 Strain-B 为主坐标参考，以 Strain-C 为辅助核查参考，结果面向 Strain-T 但不视为 Strain-T 精确坐标；现重新定位为开发/测试占位，长期策略等真实 Strain-T 数据到位后另行评估 |
| [ADR-0002](ADR-0002-streamlit-read-only-presentation.md) | Superseded by ADR-0010 | Streamlit 作为只读展示层消费已验证运行产物，不承载科学计算或结果编辑 |
| [ADR-0003](ADR-0003-reference-completeness-and-acceptance.md) | Accepted | 逐序列记录组装端点完整性，传播注释边界不确定性，并拆分计算、验证和科学验收状态 |
| [ADR-0004](ADR-0004-annotation-qualification-route.md) | Superseded by ADR-0009 | 不接受28个高置信区间作为阈值依据；保留长读长 Strain-B 主坐标并先评估独立补充注释 |
| [ADR-0005](ADR-0005-evidence-integrated-boundary-track.md) | Superseded by ADR-0009 | 不直接采用现有不合格注释；在锁定主坐标上建立实验数据支持的版本化边界轨道，并先执行 Slice 0B 数据资格门禁 |
| [ADR-0006](ADR-0006-capability-scoped-transcript-evidence.md) | Superseded by ADR-0008 | 转录来源按覆盖、剪接、方向、边界和菌株特异能力分别资格化；先对 PRJNA1210090 的单个重复执行受控探针 |
| [ADR-0007](ADR-0007-srr31989028-controlled-probe-outcome.md) | Superseded by ADR-0008 | 以完整 normalized SRA 执行 SRR31989028 受控探针；探针因广泛高深度覆盖和 BAM/BAI 字节重复性失败而保持 blocked |
| [ADR-0008](ADR-0008-target-strain-native-sequencing-redirect.md) | Accepted | Strain-T 原生测序已获授权（长期项目）；暂停 ADR-0006/0007 的公开数据证据路线（该数据源经核实是无关的工程菌株胁迫实验），下一次科学里程碑改为绑定真实测序数据到位；期间允许在 Strain-B 坐标上继续 Slice 1 引擎工程自测，但不产生正式候选 |
| [ADR-0009](ADR-0009-deprioritize-proxy-data-compensation.md) | Accepted | 扩大 ADR-0008 的范围：Strain-B 主坐标重定位为开发/测试占位，暂停 ADR-0004 注释来源调查延续工作与 ADR-0005/Slice 0B 证据整合边界轨道路线；与具体数据来源无关的框架代码（参考校验、注释归一化、三态验收模型、Slice 1 规则引擎）继续投入 |
| [ADR-0010](ADR-0010-streamlit-import-and-trigger-workflow.md) | Accepted | 修订 ADR-0002：Streamlit 允许导入新参考数据并触发 Slice 0/Slice 1 等核心引擎运行，但页面代码不得重新实现科学规则、不得提供直接编辑结果的功能；Run Catalog 准入规则和三态状态模型不变 |
| [ADR-0011](ADR-0011-broaden-framework-authorization-exclude-slice2.md) | Accepted | 把 ADR-0009 的"数据来源无关框架代码继续投入"原则细化到各切片：Slice 1（继续）、Slice 3 的验收/一致性工具部分、Slice 4（Streamlit）均获授权用占位数据实现；Slice 2（Strain-B-Strain-C 共线性核查）本质是另一种代理数据补偿，明确排除、不现在实现 |

受控转录证据、边界轨道、共线性核查与阈值推断均需等待真实 Strain-T 数据到位并经独立 ADR 评估后才能重启；不得提前冻结阈值。ADR-0004/ADR-0005/ADR-0006/ADR-0007 的历史决策与既有证据继续作为审计存档保留，未被删除或改写。
