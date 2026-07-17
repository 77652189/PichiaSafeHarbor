# ADR-0001：Strain-T 候选筛选的双参考策略

- 状态：Accepted
- 日期：2026-07-14
- 后续修订：[ADR-0008](ADR-0008-target-strain-native-sequencing-redirect.md) 记录本 ADR 预留条款（"未来获得 Strain-T 测序数据后可用新 ADR 替代本决策"）的触发情况——Strain-T 原生测序已获授权但尚未交付，本 ADR 的 Strain-B 主坐标策略在数据到位并经评估前继续有效。[ADR-0009](ADR-0009-deprioritize-proxy-data-compensation.md) 进一步把 Strain-B 主坐标重新定位为开发/测试占位，项目不再投入精力提升其可信度，等待真实数据到位。[ADR-0011](ADR-0011-broaden-framework-authorization-exclude-slice2.md) 指出本 ADR 提出的 Strain-C 交叉印证理由（Slice 2）本身也是代理数据补偿，随本 ADR 一并暂停实现。

## 背景

项目目标菌株是 Strain-T，但当前没有确认可用的、同等成熟的 Strain-T 原生染色体级参考组装。公开比较基因组研究表明 Strain-T 与 Strain-B 谱系接近，但并不完全相同。

可用公共资源中，Strain-B `GCA_001746955.1` 是长读长精修的染色体级组装；Strain-C `GCA_000223565.1` 具有成熟的染色体级参考和系统注释。两套组装均存在不能忽略的版本和完整性限制。

## 决策

第一阶段：

1. 使用 Strain-B `GCA_001746955.1` 作为唯一主坐标空间和候选枚举骨架。
2. 使用 Strain-C `GCA_000223565.1` 提供共线性、邻近基因关系和结构一致性证据。
3. 目标菌株仍记录为 Strain-T，但所有结果必须声明 `exact_target_strain_coordinates = false`。
4. 未定位 scaffold、缺失端粒和坐标不可靠区域不进入保守推荐集合。
5. 不允许把辅助参考坐标直接覆盖或替换主坐标。

## 理由

- Strain-B 比 Strain-C 更接近 Strain-T 的已知谱系背景。
- Strain-B 精修组装包含 PacBio 和 Illumina 数据，适合作为主要区间计算骨架。
- Strain-C 可提供独立参考，帮助识别结构冲突和邻近关系不稳定的候选。
- 双参考证据比任一单参考更适合表达当前不确定性。

## 放弃的方案

### 仅使用 Strain-C

参考成熟，但与 Strain-T 的谱系距离大于 Strain-B，作为唯一代理不够贴近目标菌株。

### 仅使用 Strain-B

缺少独立结构核查，且组装含未定位 scaffold 和缺失端粒，容易把参考特有问题误认为低风险区域。

### 构建 Strain-T 原生参考后再启动项目

科学上最直接，但当前没有相应数据，会阻塞可验证的第一阶段。未来获得 Strain-T 测序数据后可用新 ADR 替代本决策。

## 后果

- 输出坐标属于 Strain-B 坐标空间。
- 候选需附带 Strain-C 核查状态。
- 产品文案只能称为“面向 Strain-T 的近缘参考代理预测”。
- 后续若引入 Strain-T 原生组装，需要提供坐标迁移、结果重算和旧版本追溯方案。
- Strain-B 组装与注释完整性的具体处理由 [ADR-0003](ADR-0003-reference-completeness-and-acceptance.md) 约束。

## 验证

- 运行清单准确记录两套 assembly 和注释版本。
- 所有候选均能追溯到 Strain-B 坐标。
- 报告显著展示 Strain-T 非精确坐标声明。
- 未定位 scaffold 和缺失端粒区域无保守推荐候选。
- Strain-C 核查结果使用 `confirmed / conflicting / unavailable`，不伪造缺失证据。
