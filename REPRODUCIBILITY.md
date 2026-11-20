# Reproducing the GHIST 2026 Bottleneck analysis

These instructions describe local analysis of the already scored individual entry. They do not submit to Synapse and do not regenerate the locked prediction file. Use a separate working copy for reruns: analysis scripts overwrite derived tables under data/ and results/.

## Inputs

Obtain the following inputs through their original sources and place them under data/ with these exact filenames:

- [Testing VCF, syn74933321 version 2](https://www.synapse.org/Synapse:syn74933321.2): `data/bottleneck.testing.vcf.gz`; expected MD5 `8f12cba8d8f65076d3abe242d56b8a86`.
- [Final VCF, syn74933318 version 3](https://www.synapse.org/Synapse:syn74933318.3): `data/bottleneck.final.vcf.gz`; expected MD5 `3d38a487825c2641ca7def0d6b311ffb`.
- [Organizer example SFS at the recorded commit](https://huggingface.co/spaces/RyanGutenkunst/dadi_two_epoch/blob/425751dad7abd64203c8db30309ec88568787178/GHIST_2026_bottleneck_testing.fs): `data/GHIST_2026_bottleneck_testing.fs`.

The two Synapse manifests are tracked in data/. VCFs, per-replicate simulation matrices, and other raw inputs are not distributed in this repository. Input acquisition remains subject to the organizer's access conditions.

## Software

The existing analysis record lists Python 3.13, R 4.5.2, bcftools 1.23.1, phyclust 0.1-34, and scrm 1.7.5. Python dependencies are NumPy, SciPy, Matplotlib, and scikit-learn. The Python code requires Python 3.11 or newer because it uses hashlib.file_digest. Historical Python package versions were not recorded, so this is not a fully locked environment. Record exact versions for any new rerun.

The initial simulation script uses parallel::mclapply with three cores; run that script on macOS or Linux.

## Analysis sequence

From the repository root, with the inputs and software already available:

```sh
# Generate independent ALT-count spectra from the official VCFs.
python scripts/vcf_to_fs.py data/bottleneck.testing.vcf.gz data/bottleneck.testing.fs
python scripts/vcf_to_fs.py data/bottleneck.final.vcf.gz data/bottleneck.final.fs

# Fit the folded baselines required by the comparison checks.
Rscript --vanilla scripts/fit.R data/bottleneck.testing.fs results/baseline/testing/raw_folded folded
Rscript --vanilla scripts/fit.R data/bottleneck.final.fs results/baseline/final folded

# Verify genotype counts and create the polarized 1 Mb block spectra.
python scripts/validate_final.py
python scripts/make_blocks.py

# Generate the two simulation panels using the seeds in the scripts.
Rscript --vanilla scripts/simulate.R results/optimization/simulated_blocks
Rscript --vanilla scripts/simulate_review.R results/review/chromosomes

# Fit, compare, and record the selection decisions.
python scripts/compare.py results/optimization/simulated_blocks
python scripts/review.py

# Independently check the selected GLS estimates using the R implementation.
Rscript --vanilla scripts/check_gls.R results/optimization/testing_blocks.tsv results/review/official_fits.tsv testing block_gls
Rscript --vanilla scripts/check_gls.R results/optimization/final_blocks.tsv results/review/official_fits.tsv final block_gls
```

The historical `scripts/finalize.py` is intentionally excluded: it replaces `submission/bottleneck.final.yaml`. Do not run it for this locked entry. It is retained as provenance for the archived coordinate-sensitivity, jackknife, and figure-generation procedure. The commands above regenerate fitting and selection tables but do not regenerate those finalize-only outputs.

## What to compare

The selected method should remain `block_gls`. The archived full-precision estimates in `results/review/decision.tsv` are 1961.326267785825 generations and 0.7654675823454838 for the size ratio. The submitted rounding is 1,961.3 and 0.765468. Do not replace the historical prediction if a different environment produces a discrepancy; investigate it in the separate working copy.

The tracked evidence includes:

- `results/optimization/input_qc.tsv`: input checksums, retained records, and ancestral-state checks.
- `results/optimization/simulation_summary.tsv`: initial 80-replicate comparison.
- `results/review/summary.tsv` and `results/review/decision.tsv`: 40-chromosome comparison and retention of the incumbent.
- `results/optimization/coordinate_sensitivity.tsv`: all-record versus last-record-per-position fits.
- `results/optimization/optimized_jackknife_summary.tsv`: archived deletion sensitivity and jackknife standard errors.

The official Final RRMSE, 0.7227318443058701, is an organizer-provided score. It cannot be recomputed from the archived predicted parameters without the hidden Final truth. Local simulation errors and Testing error are separate results.

## Scope of the November 19 documentation review

The methods draft and instructions were checked against the existing source code and tracked summary tables. No VCFs were downloaded or uploaded, and the full fitting pipeline was not rerun during this documentation review. The historical prediction, original method draft, and analytical code were preserved.
