# GHIST 2026 · Bottleneck

[Submission status and outcome](STATUS.md) · [Full method write-up (v2)](submission/GHIST_Bottleneck_Writeup_v2.md) · [Reproduction guide](REPRODUCIBILITY.md)

This repository is public. The code, the MIT license and the English method write-up are openly accessible, and the tag `v1.0-ghist-final` marks the state of the entry that produced the scored Final submission. The method write-up has also been published to the main Wiki of the Synapse project [syn78752926](https://www.synapse.org/Synapse:syn78752926), which is publicly viewable. As of 20 November 2026 the Synapse writeup page still shows **"Writeups Currently Not Accepted"**, so the write-up has not yet been formally submitted and consortium-level authorship is not yet confirmed.

## Final entry file

**`submission/bottleneck.final.yaml`** holds the selected Final prediction. For the standard VCF group the two parameters are estimated as **1,961.3 generations** and **0.765468**. The population retained about 76.55% of its ancestral effective size after the decline, a reduction of 23.45%. The entry was submitted as an individual participant. The official submission ID is **9785337** and the RRMSE is **0.7227318443058701**.

The Final prediction is scored and locked; it will not be modified or resubmitted. Applying for consortium-level authorship additionally requires submitting a publicly accessible Project Wiki under the same individual identity as the Final, as the organizers specify; this GitHub documentation does not itself constitute a formal write-up submission. The full English write-up is `submission/GHIST_Bottleneck_Writeup_v2.md`; the original short draft `submission/GHIST_Bottleneck_Writeup.md` is retained. The challenge home page lists 20 November 2026, 23:59 UTC as the prediction deadline; the writeup instructions page gives no separate deadline or opening date. [Schedule](https://www.synapse.org/Synapse:syn74803798/wiki/) · [Writeup instructions](https://www.synapse.org/Synapse:syn74803798/wiki/640753)

## Independent review

Forty continuous 100 Mb simulated chromosomes cover four new bottleneck scenarios, ten replicates each. Every fit used the same set of inputs, while comparing ancestral-allele polarization, frequency binning and statistical weighting. The percentages below are relative RMSE against known simulated truth.

| Estimator | Both parameters | Time of decline | Retained fraction |
|---|---:|---:|---:|
| Folded Poisson | 8.30% | 11.64% | 1.53% |
| Polarized Poisson, 113 frequency classes | 5.91% | 8.24% | 1.40% |
| Polarized Poisson, 14 frequency bins | 5.92% | 8.26% | 1.40% |
| Fixed-shrinkage GLS, incumbent model | 6.32% | 8.83% | 1.45% |
| Adaptive-shrinkage GLS, candidate model | 6.29% | 8.78% | 1.44% |

Across these four scenarios, using allele polarization lowered the overall error. Adding block covariance weighting gave a slightly higher error than binning alone. Under a rule fixed in advance, only adaptive GLS was eligible to replace the incumbent model. It improved two of the four scenarios and 18 of the 40 chromosomes. A paired bootstrap stratified by scenario over 10,000 resamples put the relative reduction in combined error at 0.53%, with a 95% interval of −1.71% to +3.03%. Neither replacement condition was met. scrm uses a 10 kb history-pruning approximation; these errors apply only to this set of simulations.

The review also found that the gradient optimizer sometimes reported success at a poor optimum, most often under milder bottlenecks. All methods now run an independent gradient-free optimization from both grid-search and gradient-search candidate points, and every comparison was recomputed. The release procedure recomputes the official prediction from the current spectrum and checks it against the existing results. A separate ALT count is used to cross-check bcftools, so the validation also applies to data whose ancestral allele is ALT.

## Fits to the official data

| Data | Method | Generations | Retained fraction | Combined relative error |
|---|---|---:|---:|---:|
| Official Testing VCF | Previous method, folded SFS | 4,585.6 | 52.2420% | 2.57% |
| Official Testing VCF | New method, block-weighted polarized SFS | 4,598.6 | 52.1876% | 2.37% |
| Official Testing truth | Published by the organizers | 4,757 | 52.0000% | — |
| Official Final VCF | Previous method, folded SFS | 1,829.7 | 75.7763% | unknown |
| Official Final VCF | New method, block-weighted polarized SFS | **1,961.3** | **76.5468%** | **72.27% (official RRMSE)** |

The error is the root mean square of each parameter's relative error. For the new method on Testing, for example, the time estimate is 3.33% low and the retained fraction about 0.19 percentage points high. The Testing errors were computed locally. The Final RRMSE comes from the organizers' scoring email. Python reproduced the R baseline on the raw VCFs, and an independent R numerical optimization reproduced the new method, with the Final time estimate agreeing to within 0.01 generations.

## Directory layout

| Directory | Contents |
|---|---|
| `submission/` | The locked prediction, the original short English draft, and the full write-up v2 |
| `data/` | GitHub holds only the official input manifests; VCFs and spectra are obtained or generated during local reproduction |
| `scripts/` | Data processing, fitting, simulation and independent review scripts |
| `results/optimization/` | Archived summary tables; spectrum matrices, per-replicate results and validation figures are produced by local analysis |
| `results/review/` | Summary results for five methods and the replacement decision; per-replicate simulation inputs must be generated locally |
| `results/baseline/` | Produced by running the baseline fits locally; not included in the GitHub archive |

## Why the new method was selected

Both methods use the same two-epoch demographic model, mutation rate and sequence length. The previous method binned by minor allele count and fitted a Poisson composite likelihood. The new method uses the verified ancestral-allele polarization, keeping low-frequency and high-frequency derived-allele information separate. It divides the 100 Mb genome into 100 blocks of 1 Mb and derives the fitting weights from how strongly the different frequency classes co-vary across blocks.

Concretely, the 113 polymorphic frequency classes are pooled into 14 bins with upper bounds 1, 2, 3, 4, 6, 9, 14, 21, 32, 46, 60, 79, 99 and 113. The odd-numbered and even-numbered blocks are averaged separately, and each set is weighted by the inverse covariance of the other. The covariance matrix is shrunk toward its diagonal by a fixed 50%. Both sets share one model. Holding the size-change ratio and the time fixed makes the mutation scale directly estimable, which μ and the sequence length then convert into effective population sizes and generations.

The selection rule was fixed before the comparison. The new method had to improve both the official raw Testing VCF and the pooled simulation error, improve at least three of the four simulated scenarios, and return an interior solution on the official Testing data. All conditions were met. The public Testing truth is known, so the genuinely blind target is the Final.

| Scenario | Ancestral Ne | True generations | True retained fraction | Previous method rRMSE | New method rRMSE |
|---|---:|---:|---:|---:|---:|
| 1 | 20,000 | 1,600 | 0.75 | 14.86% | 9.85% |
| 2 | 20,000 | 5,000 | 0.50 | 5.42% | 3.80% |
| 3 | 25,000 | 2,500 | 0.30 | 4.45% | 3.17% |
| 4 | 15,000 | 8,000 | 0.70 | 17.08% | 7.85% |
| All 80 replicates | | | | **11.85%** | **6.77%** |

Each scenario was simulated 20 times independently with Hudson's ms, at a sample size of 114 haploid chromosomes and an independent random seed per replicate. Each replicate contains 1,000 independent 25 kb segments, with within-segment recombination at the stated rate of 1.87 × 10⁻⁸ per base per generation, and every ten segments form one statistical block. The total length is 25 Mb. The official data cover a single continuous 100 Mb chromosome, so the simulations serve to compare parameter recovery and cannot be used to infer the actual Final error. Both methods used identical inputs in every replicate. The new method had the lower single-replicate error in 52 of 80. Scenario summaries pool squared relative errors before taking the square root.

## Raw-data checks and stability of the results

- **Testing** uses syn74933321.2. The raw VCF has 1,194,350 filter-passing SNP records, while the spectrum distributed with the organizers' example application has 1,186,635. Keeping only the last record per coordinate reproduces all 113 frequency-class counts of the example spectrum exactly. The 4,550.4 generations reported in an earlier version came from that example spectrum. Both methods in the table above use the complete VCF.
- **Final** uses syn74933318.3. The 57 diploid samples contain 1,588,748 usable SNP records. Keeping only the last record per coordinate leaves 1,575,522, for which the new method estimates 1,931.1 generations and 76.58%. The submitted value uses the complete set of records. The MD5 of both VCFs matches their respective Synapse manifests.
- The AA metadata lines in the raw VCFs have incomplete line terminators, which makes bcftools emit a parsing warning. An independent Python parser counts genotypes directly and reads AA. The per-frequency counts agree with bcftools for both datasets, and all retained Final records have AA = REF.
- Deleting successive 10 Mb segments and re-estimating the covariance and all parameters gives ten time estimates between 1,853 and 2,073 generations and retained fractions between 75.92% and 77.07%. The corresponding block jackknife standard errors are about 204 generations and 1.12 percentage points. Fifty checks deleting 2 Mb are in the results tables. These ranges describe deletion sensitivity.
- The ancestral effective size is estimated at about 23,350 and the current size at about 17,874. Both are model effective population sizes. The stated μ = 3.42 × 10⁻⁸ per base per generation and a sequenced length of 100,000,000 bp provide the absolute scaling.

[Simulation validation table](results/review/summary.tsv) · [Final deletion sensitivity](results/optimization/optimized_jackknife_summary.tsv). The validation figure for the complete project is saved at `results/optimization/validation.png` and is regenerated when the final fitting script is run.

## Reproduction

Detailed steps are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). Run them in a separate working copy so that archived result tables are not overwritten. The historical `scripts/finalize.py` replaces the locked prediction file and has therefore been removed from the steps below; do not run it against this entry.

This work used R 4.5.2, phyclust 0.1-34, bcftools 1.23.1, scrm 1.7.5, and Python 3.13 with NumPy, SciPy, Matplotlib and scikit-learn. The optimization scripts require Python 3.11 or newer. The commands below run from the repository root. The complete local package contains the VCFs, spectra and simulation matrices. The GitHub version keeps the code, the final prediction, the write-up, the input manifests and the summary tables. To recompute, first download the two official VCFs and the organizers' example spectrum into `data/` under the filenames used by the scripts, then run:

```sh
# Recompute the spectra of both raw VCFs, the baseline models and an independent GT check.
python scripts/vcf_to_fs.py data/bottleneck.testing.vcf.gz data/bottleneck.testing.fs
python scripts/vcf_to_fs.py data/bottleneck.final.vcf.gz data/bottleneck.final.fs
Rscript --vanilla scripts/fit.R data/bottleneck.testing.fs results/baseline/testing/raw_folded folded
Rscript --vanilla scripts/fit.R data/bottleneck.final.fs results/baseline/final folded
python scripts/validate_final.py
python scripts/make_blocks.py

# The project contains the per-replicate spectra; when recomputing from GitHub, run the two R simulation commands below.
Rscript --vanilla scripts/simulate.R results/optimization/simulated_blocks
Rscript --vanilla scripts/simulate_review.R results/review/chromosomes
python scripts/compare.py results/optimization/simulated_blocks
python scripts/review.py
Rscript --vanilla scripts/check_gls.R results/optimization/final_blocks.tsv results/review/official_fits.tsv final block_gls
```

Baseline parameters are saved as TSV. These steps recompute the fitting and model-selection tables; they do not write the locked YAML and do not submit anything to the challenge platform. Coordinate sensitivity, the jackknife and the validation figure are outputs of the historical `finalize.py`; this procedure keeps their existing summary tables rather than regenerating them. The recorded historical Python dependencies are not fully version-locked, so record the actual software versions for any cross-environment reproduction.

## Sources

- [GHIST 2026 challenge description and Testing truth](https://www.synapse.org/Synapse:syn74803798/wiki/640740) · [Challenge rules](https://www.synapse.org/Synapse:syn74803798/wiki/640738)
- [Official Testing VCF](https://www.synapse.org/Synapse:syn74933321.2) · [Official Final VCF](https://www.synapse.org/Synapse:syn74933318.3) · [Official YAML template](https://www.synapse.org/Synapse:syn74933357)
- [Organizers' example spectrum](https://huggingface.co/spaces/RyanGutenkunst/dadi_two_epoch/blob/425751dad7abd64203c8db30309ec88568787178/GHIST_2026_bottleneck_testing.fs) · [Example application output format](https://huggingface.co/spaces/RyanGutenkunst/dadi_two_epoch/blob/425751dad7abd64203c8db30309ec88568787178/app_natparams.py)
- [Jouganous et al., allele frequency spectrum inference](https://doi.org/10.1534/genetics.117.200493) · [Hudson, the ms coalescent simulator](https://doi.org/10.1093/bioinformatics/18.2.337) · [The ms interface in phyclust](https://stat.ethz.ch/CRAN/web/packages/phyclust/phyclust.pdf)
