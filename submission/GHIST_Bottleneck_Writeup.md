# GHIST 2026 · Bottleneck

## Authors and affiliations

Candice Wu (Synapse username: candicewu), University of Iowa.

## Bottleneck Challenge

The estimates for the standard VCF challenge are **1,961.3 generations since the decline** and **0.765468 for the ratio of current to ancestral effective population size**. The analysis used the organizer's Final VCF (syn74933318, version 3). It contains 57 diploid individuals and 1,588,748 polymorphic records on a 100 Mb chromosome. We counted derived alleles using the ancestral allele in each VCF record, retained records sharing a position, and pooled the 113 derived-frequency classes into 14 bins. We fitted a single instantaneous decline in effective size to the expected site frequency spectrum. We partitioned the chromosome into 100 consecutive 1 Mb blocks. For each alternating set of blocks we fitted its mean spectrum using the other set's covariance, shrinking the off-diagonal elements by 50%. We estimated the mutation scale and converted the fitted time to generations using the supplied mutation rate of 3.42 × 10⁻⁸ per base per generation. Python performed the primary fit and an independent R implementation verified the released estimates.

## Validation and discussion

The public Testing data gave 4,598.6 generations and a retained fraction of 0.521876, compared with the published values of 4,757 and 0.52. We compared the estimator with a folded-spectrum baseline on 80 simulated genomes and evaluated adaptive covariance shrinkage on 40 independently simulated 100 Mb chromosomes. In the latter comparison, the fixed and adaptive fits had pooled relative RMSEs of 6.32% and 6.29%. The paired 95% interval for the relative error reduction crossed zero, so we retained fixed shrinkage. Using only one record per Final coordinate shifted the inferred date to 1,931.1 generations. A statement on the challenge page about how to handle repeated coordinates would make it easier for participants to reproduce the example spectrum.

## Data, code and references

The challenge VCFs were the only data inputs for the official predictions. The accompanying source and reproduction commands are in the scripts directory and README. The local validation used Hudson's ms through phyclust and scrm 1.7.5 with a 10 kb pruning window.

- [GHIST 2026 demographic history challenge and public Testing answer](https://www.synapse.org/Synapse:syn74803798/wiki/640740)
- [Testing VCF](https://www.synapse.org/Synapse:syn74933321.2) and [Final VCF](https://www.synapse.org/Synapse:syn74933318.3)
- [Jouganous et al. (2017), allele frequency spectrum inference](https://doi.org/10.1534/genetics.117.200493)
- [Hudson (2002), the ms coalescent simulator](https://doi.org/10.1093/bioinformatics/18.2.337)
- [scrm package and documentation](https://cran.r-project.org/package=scrm)

## Author contributions

This was an individual entry under the Synapse account candicewu.
