# Alpha 整合与稳定版发布门槛

依据 2026-10-08 的合并审查方案，工程整合和正式 v0.5 发布分别验收。
当前版本仍为 `0.5.0-alpha.2`、schema `1.2.0`。

## 堆叠 PR 的合并顺序

原始依赖：`main` → PR #1 的 `agent/tevox-p0-p1`（`3ba38b5`）→
PR #2 的 `agent/v05-integration-20261005`（原审查基线 `9db4724`）。
两层分别包含 6、9 个原始提交；后续修复追加在 PR #2，保留其祖先关系。

1. 审阅 PR #2 的最新提交，确认其 CI 和原生 Mac 打包通过。
2. PR #2 使用 **merge commit** 合入现有目标 `agent/tevox-p0-p1`。
3. 对随之更新的 PR #1 HEAD 重新审阅，等待对应 HEAD 的 CI。
4. PR #1 使用 **merge commit** 合入 `main`。
5. 检查两个原始 HEAD 与修复 HEAD 都是 `main` 的祖先，最终文件树与
   审阅通过的整合树一致；若有其他并行改动，重新审查差异。

本轮方案要求合并前完成独立 review。仓库没有强制保护规则，并不代表
独立 review 已完成；作者自己的检查也不能记录为另一位审阅者的批准。
准备好供审阅的 PR 不等于发布版本。这里不增加自动合并、管理员绕过、
自动批准或打 tag 的操作。

使用当前 HEAD 的完整 SHA 执行 GitHub CLI 操作，防止合并未审阅的新提交：

```bash
# 前提：PR #2 最新 HEAD 的独立 review 和全部 CI 已完成。
gh pr view 2 --repo Wzhennan-icmm/TEvoX \
  --json baseRefName,headRefOid,reviews,statusCheckRollup
gh pr merge 2 --repo Wzhennan-icmm/TEvoX --merge \
  --match-head-commit '<reviewed-full-PR2-SHA>'

# 等待更新后的 PR #1 完成独立 review 和全部 CI。
gh pr view 1 --repo Wzhennan-icmm/TEvoX \
  --json baseRefName,headRefOid,reviews,statusCheckRollup
gh pr ready 1 --repo Wzhennan-icmm/TEvoX
gh pr merge 1 --repo Wzhennan-icmm/TEvoX --merge \
  --match-head-commit '<reviewed-full-PR1-SHA>'

git fetch origin refs/heads/main:refs/remotes/origin/main
git merge-base --is-ancestor 3ba38b598ef3f1dacf7de2462276c683293735da origin/main
git merge-base --is-ancestor 9db4724c11043aef39f390aea927fb7776696375 origin/main
git merge-base --is-ancestor '<reviewed-full-PR2-SHA>' origin/main
git diff --exit-code '<reviewed-full-PR2-SHA>^{tree}' 'origin/main^{tree}'
```

尖括号参数必须替换为实际审阅过的 SHA。`--match-head-commit` 只检查
HEAD 没有变化，不替代 review、CI 或 base 分支核查。不要 squash/rebase
基础 PR 后原样合并堆叠 PR，也不要删去仍被依赖的开发分支。

## 本轮输出修复的验收

同一前缀普通 TSV 与 gzip TSV 并存会让 Python 读取器报歧义。
输出预检现要求：只要 17 张表中任意一张存在相反模式的路径，就在
分析与写出前拒绝运行；第二次写出前仍复查。既有结果与 `.run.json`
保留，不自动删除。部分输出、目录及悬空符号链接也属于冲突。
同模式重跑继续遵循原有的普通文件覆盖规则。

`tests/check_output_modes.py` 覆盖两个切换方向、全部 17 张部分表、
目录／悬空符号链接和同模式重跑，并在拒绝时比较全部既有文件的
内容、路径、inode 和修改时间。测试已纳入常规 `make check`。
这不提供多个进程同时使用同一前缀的锁；并行任务仍需使用不同前缀。

## 功能声明的边界

- Hungarian 精确匹配只覆盖两侧节点数都不超过阈值的匹配块；超限使用
  明确标注的确定性回退。
- 分量 `OPTIMAL` 只对截断及匹配后的保留图成立。后续约束阻断匹配边时，
  不会重新考虑被 Hungarian 淘汰的替代边；不能宣称所有候选、匹配和
  分区联合求解的全局最优。
- 分数仍为 `BUILTIN_UNCALIBRATED_V1`，不是经过独立验证的后验概率。
- 五物种结果、人类 GIAB 子集、小鼠 PCR 子集、合成拷贝上下文和规模测试
  分别见[整合报告](INTEGRATION_VALIDATION.zh-CN.md)。它们不能相互替代。
- 10 月 6 日冻结的结果及二进制哈希继续描述当时的运行。之后的输出预检
  修复不意味着五物种重新运行过，也不改写既有实验结果。

## 稳定版 v0.5 仍未满足的科学门槛

1. 独立候选关系真值：具有可追溯抽样设计的同位点／不同位点正负标签；
   不把未收录、未判定或 TEvoX 自己的输出当作负例或真值。
2. 正式概率校准：开发、校准、外部测试集隔离；祖先位点、同源组和验证
   批次无泄漏；符合候选生成约定；训练与模型加载实现完成；报告外部
   测试的可靠性曲线、Brier、对数损失及不确定性和适用域。
3. 真实多倍体独立验证：材料、组装、repeat/gene 注释和实验标签对应，
   含有竞争同源拷贝及亚基因组误配评估。合成数据或少量结构变异 PCR
   案例不等于全基因组 TE 拷贝归属验收。
4. 规模结论与适用范围对应：现有稀疏图测试和五个真实 pair 不能证明
   密集图或群体规模均已验证；相应产品声明须有相应资源与准确性证据。

这些门槛通过前，可审阅和整合 alpha 工程代码；不得仅因为 CI 通过就
创建稳定版 v0.5 发布，或将未校准分数改名为概率。
