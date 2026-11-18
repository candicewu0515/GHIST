#!/usr/bin/env python3
"""Count ALT copies in complete biallelic diploid SNPs using bcftools."""
import argparse
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
import subprocess


def read_spectrum(vcf):
    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(vcf)], text=True
    ).splitlines()
    if len(samples) < 3 or len(set(samples)) != len(samples):
        raise ValueError("Expected at least three distinct diploid samples.")
    n = 2 * len(samples)
    counts = [0] * (n + 1)
    qc = Counter()
    # Recompute AC and AN from GT, since existing INFO tags can be stale.
    # Use FILTER unchanged for QC and discard non-PASS sites locally.
    commands = [
        ["bcftools", "view", "-m2", "-M2", "-v", "snps", "-Ou", str(vcf)],
        ["bcftools", "+fill-tags", "-Ou", "--", "-t", "AC,AN"],
        ["bcftools", "query", "-f", "%FILTER\t%AN\t%AC\n"],
    ]
    with ExitStack() as stack:
        previous = None
        jobs = []
        for command in commands:
            proc = stack.enter_context(subprocess.Popen(
                command, stdin=previous, stdout=subprocess.PIPE
            ))
            if previous is not None:
                previous.close()
            previous = proc.stdout
            jobs.append(proc)
        for row in jobs[-1].stdout:
            filt, called, alt = row.rstrip(b"\n").split(b"\t")
            qc["biallelic_snp_records"] += 1
            if filt not in (b".", b"PASS"):
                qc["excluded_filter"] += 1
                continue
            if called == b"." or int(called) != n:
                qc["excluded_incomplete_diploid"] += 1
                continue
            ac = int(alt)
            if not 0 <= ac <= n:
                raise ValueError("Inconsistent allele count")
            if ac in (0, n):
                qc["excluded_monomorphic"] += 1
                continue
            counts[ac] += 1
        for proc in jobs:
            if proc.wait() != 0:
                raise RuntimeError("bcftools failed at " + " ".join(proc.args))
    if not sum(counts):
        raise ValueError("Zero complete, polymorphic biallelic SNPs.")
    qc["retained"] = sum(counts)
    qc["haploid_n"] = n
    return counts, qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", type=Path)
    parser.add_argument("fs", type=Path)
    args = parser.parse_args()
    counts, qc = read_spectrum(args.vcf)
    mask = [1] + [0] * (len(counts) - 2) + [1]
    args.fs.parent.mkdir(parents=True, exist_ok=True)
    args.fs.write_text('{} unfolded "ALT_count"\n{}\n{}\n'.format(
        len(counts), " ".join(map(str, counts)), " ".join(map(str, mask))
    ))
    for name, value in sorted(qc.items()):
        print("{}\t{}".format(name, value))


if __name__ == "__main__":
    main()
