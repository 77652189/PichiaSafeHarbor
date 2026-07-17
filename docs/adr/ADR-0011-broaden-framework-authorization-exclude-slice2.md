# ADR-0011：扩大框架建设授权范围至 Slice 1/3/4，明确排除 Slice 2

- 状态：Accepted
- 日期：2026-07-16
- 关系：应用并细化 [ADR-0009](ADR-0009-deprioritize-proxy-data-compensation.md)"数据来源无关框架代码可继续投入"的原则到 Slice 1 之外的其余切片；不改变 [ADR-0001](ADR-0001-reference-strategy.md) 关于 Strain-C 辅助参考的定位，但暂停 Slice 2 的当前实现工作；不改变 [ADR-0010](ADR-0010-streamlit-import-and-trigger-workflow.md) 关于 Streamlit 架构的决定，本 ADR 只是正式授权其开工。

## 背景

ADR-0009 建立的原则是：与具体数据来源无关的框架代码（参考校验、注释归一化、三态验收模型、Slice 1 规则引擎）应继续投入，而"补偿缺失 Strain-T 数据"的工作（Strain-B 主坐标改善、注释来源调查、证据整合边界轨道）应当暂停。此前的 `EXECUTION_PLAN.md` 只把这一原则应用到了 Slice 1，把 Slice 2/3/4 整体列为"未授权"，没有逐一区分它们各自哪些部分属于哪一类。

逐一复核发现：

- **Slice 2**（Strain-B 与 Strain-C 共线性核查）按 ADR-0001 的原始理由——"Strain-C 可提供独立参考，帮助识别结构冲突和邻近关系不稳定的候选"——本质上是"Strain-B 只是代理，需要借助另一近缘菌株 Strain-C 交叉印证"，这与 ADR-0004/ADR-0005 是同一类"补偿缺失 Strain-T 数据"的工作，不是数据来源无关的框架代码。一旦真实 Strain-T 数据到位，这个交叉印证的必要性会大幅下降甚至消失——直接使用 Strain-T 数据判断候选安全性即可，不需要再靠 Strain-C 为 Strain-B 佐证。
- **Slice 3**（完整 MVP 验收）实际上是混合的：其中"冻结默认阈值 ADR""运行完整基因组""人工浏览器抽查和文献位点回顾"依赖可信边界数据，属于必须等待真实数据的部分；但"验证多格式一致性和可重复性"是通用验收/核验工具（类似已有的 `create_acceptance_manifest`/`create_slice0a_acceptance` 模式），与具体数据来源无关，属于可以现在建设的框架代码。
- **Slice 4**（Streamlit）按 ADR-0010 的设计，是消费版本化结果契约、触发核心引擎的展示/触发层，同样与具体数据来源无关。

## 决策

1. 授权将"引擎工程自测"的范围从"仅 Slice 1"扩大为：Slice 1（候选窗口规则引擎，继续补强）+ Slice 3 中与数据来源无关的验收/一致性工具部分（例如 Slice 1 自身的 acceptance/独立核验流程）+ Slice 4（Streamlit 全部）。三者均可现在使用 Strain-B/Strain-C 占位数据进行工程实现与自测，产物统一标记 `run_purpose=engine_readiness_test_not_scientific_output` 且 `scientific_acceptance_status=blocked`。
2. **明确排除 Slice 2**：不在当前阶段实现或维护 Strain-B-Strain-C 共线性核查代码。如果未来希望以"独立于具体菌株的结构稳定性核查"这一新理由重新论证 Slice 2 的价值（而不是"Strain-B 需要佐证"这个旧理由），需要另立 ADR 重新评估，不得默认沿用 ADR-0001 的旧理由直接继续实现。
3. Slice 3 中依赖可信边界数据的部分（冻结阈值、全基因组正式运行、人工抽查、文献位点回顾）继续等待真实 Strain-T 数据到位并经独立 ADR 评估，不在本 ADR 授权范围内。
4. Streamlit（Slice 4）正式获准开始实现，遵循 ADR-0010 的架构约束（可触发核心引擎、不得复制或绕过规则、不得提供编辑结果的功能）。

## 放弃的方案

### 把 Slice 2 也当作框架代码现在实现

Slice 2 的核心价值主张（"用 Strain-C 给 Strain-B 撑腰"）会在真实数据到位后大幅贬值甚至失效，现在投入工程实现的沉没成本风险高于收益。

### Slice 3 整体推迟到真实数据到位

会让通用的验收/一致性工具建设也跟着被无谓推迟，浪费等待期。

### Streamlit 继续整体搁置

ADR-0010 已经确认其架构与具体菌株数据无关；继续搁置只是延续旧的"没有可信数据就什么都不做"的思维，而不是基于实际依赖关系判断。

## 后果

- `docs/EXECUTION_PLAN.md` 的门禁更新为：Slice 1（继续）、Slice 3 验收工具、Slice 4（Streamlit）均已授权；Slice 2 与 Slice 3 剩余部分继续等待真实数据。
- `docs/ARCHITECTURE.md` §9 的 Slice 2/3/4 描述相应更新。
- `docs/adr/README.md` 新增本 ADR 条目。
- `scientific_acceptance_status` 继续保持 `blocked`，不因本 ADR 而改变。
