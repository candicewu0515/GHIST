# GHIST 2026 · Bottleneck

[提交状态和后续结果](STATUS.md) · [完整英文方法稿 v2](submission/GHIST_Bottleneck_Writeup_v2.md) · [English reproduction guide](REPRODUCIBILITY.md)

本仓库已设为 **public**，代码与英文方法稿可公开访问。2026 年 11 月 19 日再次检查时，Synapse writeup 页面仍显示 **“Writeups Currently Not Accepted”**；方法稿尚未正式提交，consortium-level authorship 尚未确认。当前仅整理文档。

## Final 参赛文件

**`submission/bottleneck.final.yaml`** 是当前选定的 Final 预测。标准 VCF 组的两个参数估计为 **1,961.3 代**和 **0.765468**。种群在这次收缩后保留了约 76.55% 的祖先有效群体大小，相当于减少了 23.45%。已按个人身份正式提交。官方提交号为 **9785337**，RRMSE 为 **0.7227318443058701**。

Final 预测已经评分并锁定，不再修改或重新提交。申请 consortium-level authorship 仍需在与 Final 相同的个人身份下，按官方要求另行提交公开的 Project Wiki；GitHub 文档本身不等于正式 writeup 提交。完整英文稿见 `submission/GHIST_Bottleneck_Writeup_v2.md`，原短稿 `submission/GHIST_Bottleneck_Writeup.md` 保留。赛事主页列出的提交截止时间为 2026 年 11 月 20 日 23:59 UTC；writeup 指南未单独列出截止时间或开放日期。[赛程](https://www.synapse.org/Synapse:syn74803798/wiki/) · [writeup 指南](https://www.synapse.org/Synapse:syn74803798/wiki/640753)

## 本次独立复核

40 条连续的 100 Mb 模拟染色体覆盖四种新的瓶颈情景，每种模拟 10 次。每次拟合均使用同一组输入，同时比较祖先等位基因方向、频数合并和统计权重。下表的百分比是已知模拟真值下的相对 RMSE。

| 估计方法 | 两项参数合计 | 发生时间 | 保留比例 |
|---|---:|---:|---:|
| Folded Poisson | 8.30% | 11.64% | 1.53% |
| Polarized Poisson，113 个频数组 | 5.91% | 8.24% | 1.40% |
| Polarized Poisson，14 个频数组 | 5.92% | 8.26% | 1.40% |
| 固定收缩 GLS，正式模型 | 6.32% | 8.83% | 1.45% |
| 自适应收缩 GLS，候选模型 | 6.29% | 8.78% | 1.44% |

在这四种情景中，利用等位基因方向降低了总体误差。加入分段协方差权重比单纯分组合并的误差稍高。按事先固定的规则，只有自适应 GLS 参与正式模型的替换。它在四种情景中改善两种，40 条染色体中改善 18 条。按情景分层配对 bootstrap 10,000 次，综合误差相对降低 0.53%，95% 区间为 −1.71% 至 +3.03%。两项更换条件均未满足。scrm 采用 10 kb 历史剪枝近似；这些误差只适用于这批模拟。

复核还发现梯度优化器有时在较差的位置报告成功，尤其是在较轻的瓶颈下。现在所有方法都从参数网格和梯度搜索的候选点分别运行独立的无梯度优化，重算所有比较。发布流程从当前频谱重算正式预测并核对已有结果。单独的 ALT 计数用于核对 bcftools，使验证也适用于祖先等位基因为 ALT 的数据。

## 官方数据拟合

| 数据 | 方法 | 发生时间（代） | 收缩后保留比例 | 综合相对误差 |
|---|---|---:|---:|---:|
| 官方 Testing VCF | 旧方法，folded SFS | 4,585.6 | 52.2420% | 2.57% |
| 官方 Testing VCF | 新方法，分段加权 polarized SFS | 4,598.6 | 52.1876% | 2.37% |
| 官方 Testing 真值 | 主办方公开值 | 4,757 | 52.0000% | — |
| 官方 Final VCF | 旧方法，folded SFS | 1,829.7 | 75.7763% | 未知 |
| 官方 Final VCF | 新方法，分段加权 polarized SFS | **1,961.3** | **76.5468%** | **72.27%（官方 RRMSE）** |

误差等于两个参数各自相对误差的均方根。例如，Testing 新方法的时间估计偏低 3.33%，保留比例偏高约 0.19 个百分点。Testing 误差在本地计算。Final 的 RRMSE 来自主办方评分邮件。Python 在原始 VCF 上复现了 R 基线，独立的 R 数值优化也复现了新方法，Final 时间估计相差不到 0.01 代。

## 目录

| 目录 | 内容 |
|---|---|
| `submission/` | 已锁定预测、原英文短稿和完整英文稿 v2 |
| `data/` | GitHub 仅含官方输入清单；VCF 与频谱在本地复现时获取或生成 |
| `scripts/` | 数据处理、拟合、模拟与独立复核脚本 |
| `results/optimization/` | 已归档汇总表；频谱矩阵、逐次模拟结果与验证图由本地分析生成 |
| `results/review/` | 五种方法的汇总结果与替换决定；逐次模拟输入需本地生成 |
| `results/baseline/` | 本地运行基线拟合后生成，不包含在 GitHub 归档中 |

## 为什么选择新方法

两种方法使用相同的双时期人口史模型、突变率和长度。旧方法按 minor allele count 合并后，以 Poisson composite likelihood 拟合。新方法利用经核查的祖先等位基因方向，分开保留低频与高频 derived allele 信息。它把 100 Mb 基因组分成 100 个 1 Mb 区段，用不同频数类别在各区段中共同波动的程度来决定拟合权重。

具体计算将 113 个多态频数类别合并为 14 组，上界依次为 1、2、3、4、6、9、14、21、32、46、60、79、99、113。把奇数区段和偶数区段分别取均值，用另一组区段的协方差逆矩阵加权。协方差矩阵固定向其对角线收缩 50%。两组共用一个模型。固定大小变化比例和时间时，突变尺度可以直接估计，再由 μ 和长度换算为有效种群大小及代数。

选择规则在比较结果前固定。新方法须同时改善官方原始 Testing VCF 和模拟总体误差，在四种模拟情景中至少改善三种，且在官方 Testing 数据上得到搜索范围内部的解。本次满足全部条件。公开 Testing 真值已知，真正的盲评对象是 Final。

| 模拟情景 | 祖先 Ne | 真值（代） | 真保留比例 | 旧方法相对 RMSE | 新方法相对 RMSE |
|---|---:|---:|---:|---:|---:|
| 1 | 20,000 | 1,600 | 0.75 | 14.86% | 9.85% |
| 2 | 20,000 | 5,000 | 0.50 | 5.42% | 3.80% |
| 3 | 25,000 | 2,500 | 0.30 | 4.45% | 3.17% |
| 4 | 15,000 | 8,000 | 0.70 | 17.08% | 7.85% |
| 全部 80 次 | | | | **11.85%** | **6.77%** |

每种情景用 Hudson’s ms 独立模拟 20 次，样本大小为 114 条单倍体染色体，每次使用独立随机种子。每次模拟包含 1,000 个独立的 25 kb 片段，片段内按题目给定的重组率 1.87 × 10⁻⁸ /bp/generation 模拟重组，每十个片段合为一个统计区段。总长为 25 Mb。官方数据覆盖一条连续的 100 Mb 染色体，所以模拟用于比较参数恢复能力，不能据此推算 Final 的实际误差。两种方法在每次模拟中使用相同的输入。新方法在 80 次中有 52 次单次误差更低。各情景汇总时先汇总平方相对误差，再开平方。

## 原始数据核查和结果稳定性

- **Testing** 使用 syn74933321.2。原始 VCF 有 1,194,350 条通过过滤的 SNP 记录，主办方示例应用附带的频谱有 1,186,635 条。按坐标仅保留最后一条后，113 个频数类别的计数全部与示例频谱一致。上一版报告的 4,550.4 代来自该示例频谱。上表两个方法均使用完整 VCF。
- **Final** 使用 syn74933318.3。57 个二倍体样本含 1,588,748 条可用 SNP 记录。按坐标仅保留最后一条后剩余 1,575,522 条，新方法估计 1,931.1 代、76.58%。提交值采用完整记录。两份 VCF 的 MD5 均与各自的 Synapse manifest 一致。
- 原始 VCF 的 AA 元数据行结束符不完整。bcftools 解析时发出警告。独立 Python 解析直接从 GT 计数并读取 AA。两个数据集的逐频数计数都与 bcftools 一致，全部保留的 Final 记录均为 AA = REF。
- 顺次删除连续 10 Mb 并重新估计协方差及全部参数，10 次时间估计在 1,853–2,073 代之间，保留比例在 75.92%–77.07% 之间。对应 block jackknife 标准误约为 204 代和 1.12 个百分点。删除 2 Mb 的 50 次检查见结果表。这里的范围描述删段敏感性。
- 祖先有效大小估计约 23,350，当前约 17,874。两者是模型的有效群体大小。题目给定 μ = 3.42 × 10⁻⁸ /bp/generation 和测序长度 100,000,000 bp，用于绝对定标。

[模拟验证结果表](results/review/summary.tsv) · [Final 删段敏感性](results/optimization/optimized_jackknife_summary.tsv)。完整项目的验证图保存在 `results/optimization/validation.png`，运行最终拟合脚本时也会生成。

## 复现

详细英文步骤见 [REPRODUCIBILITY.md](REPRODUCIBILITY.md)。在单独的工作副本中运行，以免覆盖归档结果表。历史 `scripts/finalize.py` 会替换已锁定的预测文件，因此已从下列复现步骤移除；不要对本次参赛记录运行该脚本。

本次使用 R 4.5.2、phyclust 0.1-34、bcftools 1.23.1、scrm 1.7.5 和带 NumPy、SciPy、Matplotlib、scikit-learn 的 Python 3.13。优化脚本需要 Python 3.11 或更新版本。以下命令从项目根目录运行。本地完整包包含 VCF、频谱与模拟矩阵。GitHub 版保存代码、最终预测、方法稿、输入清单及汇总表。复算时，先下载下方两个官方 VCF 和主办方示例频谱到 `data/`，使用脚本中对应的文件名，再执行下列命令。

```sh
# 重算两份原始 VCF 的频谱，基线模型和独立 GT 核查。
python scripts/vcf_to_fs.py data/bottleneck.testing.vcf.gz data/bottleneck.testing.fs
python scripts/vcf_to_fs.py data/bottleneck.final.vcf.gz data/bottleneck.final.fs
Rscript --vanilla scripts/fit.R data/bottleneck.testing.fs results/baseline/testing/raw_folded folded
Rscript --vanilla scripts/fit.R data/bottleneck.final.fs results/baseline/final folded
python scripts/validate_final.py
python scripts/make_blocks.py

# 项目含每次模拟的频谱；从 GitHub 复算时运行以下两条 R 模拟命令；本地完整包已包含这些输入。
Rscript --vanilla scripts/simulate.R results/optimization/simulated_blocks
Rscript --vanilla scripts/simulate_review.R results/review/chromosomes
python scripts/compare.py results/optimization/simulated_blocks
python scripts/review.py
Rscript --vanilla scripts/check_gls.R results/optimization/final_blocks.tsv results/review/official_fits.tsv final block_gls
```

基线参数以 TSV 保存。上述步骤重新计算拟合与模型选择表，不写入已锁定的 YAML，也不向赛事平台提交。坐标敏感性、jackknife 与验证图属于历史 `finalize.py` 的输出，本流程保留其已有汇总表，不重新生成。记录的历史 Python 依赖没有完整版本锁定，跨环境复现需记录实际软件版本。

## 来源

- [GHIST 2026 官方赛题与 Testing 真值](https://www.synapse.org/Synapse:syn74803798/wiki/640740) · [参赛规则](https://www.synapse.org/Synapse:syn74803798/wiki/640738)
- [官方 Testing VCF](https://www.synapse.org/Synapse:syn74933321.2) · [官方 Final VCF](https://www.synapse.org/Synapse:syn74933318.3) · [官方 YAML 模板](https://www.synapse.org/Synapse:syn74933357)
- [主办方示例频谱](https://huggingface.co/spaces/RyanGutenkunst/dadi_two_epoch/blob/425751dad7abd64203c8db30309ec88568787178/GHIST_2026_bottleneck_testing.fs) · [示例应用输出格式](https://huggingface.co/spaces/RyanGutenkunst/dadi_two_epoch/blob/425751dad7abd64203c8db30309ec88568787178/app_natparams.py)
- [Jouganous et al.，频谱人口史推断](https://doi.org/10.1534/genetics.117.200493) · [Hudson，ms 模拟器](https://doi.org/10.1093/bioinformatics/18.2.337) · [phyclust 的 ms 接口](https://stat.ethz.ch/CRAN/web/packages/phyclust/phyclust.pdf)
