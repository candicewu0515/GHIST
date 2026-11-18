#!/usr/bin/env python3
"""Verify the official VCF and estimate variability across twenty 5 Mb blocks."""
import csv
import gzip
import hashlib
import math
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results/baseline/final"
VCF = DATA / "bottleneck.final.vcf.gz"
N_BLOCKS = 20
LENGTH = 100000000
BLOCK_BP = LENGTH // N_BLOCKS
EXPECTED_MD5 = "3d38a487825c2641ca7def0d6b311ffb"


def write_fs(path, counts):
    mask = [1] + [0] * (len(counts) - 2) + [1]
    path.write_text('{} unfolded "sample"\n{}\n{}\n'.format(
        len(counts), " ".join(map(str, counts)), " ".join(map(str, mask))
    ))


def fit(path, destination, mode, length):
    subprocess.run(["Rscript", "--vanilla", str(ROOT / "scripts" / "fit.R"), str(path),
                    str(destination), mode, str(length)], check=True,
                   stdout=subprocess.DEVNULL)
    with (destination / "fit.tsv").open() as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def main():
    with VCF.open("rb") as stream:
        digest = hashlib.md5()
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != EXPECTED_MD5:
        raise ValueError("Input differs from Synapse syn74933318.3")
    aa_ref, aa_alt, missing_aa = 0, 0, 0
    with gzip.open(VCF, "rb") as stream:
        for line in stream:
            if line.startswith(b"#CHROM\t"):
                break
        n_diploid = len(line.rstrip().split(b"\t")) - 9
        n = 2 * n_diploid
        if n_diploid != 57:
            raise ValueError("Unexpected sample count")
        gt_pattern = re.compile(rb"[01][|/][01](?:\t[01][|/][01]){%d}" % (n_diploid - 1))
        blocks = [[0] * (n + 1) for _ in range(N_BLOCKS)]
        polarized = [0] * (n + 1)
        for line in stream:
            chrom, pos, _, ref, alt, _, filt, info, fmt, gt = line.rstrip().split(b"\t", 9)
            if chrom != b"1" or not 1 <= int(pos) <= LENGTH:
                raise ValueError("Record outside the challenge region")
            if filt not in (b".", b"PASS") or ref not in (b"A",b"C",b"G",b"T") or alt not in (b"A",b"C",b"G",b"T"):
                continue
            if fmt != b"GT":
                raise ValueError("Expected GT-only FORMAT")
            if not gt_pattern.fullmatch(gt):
                if b"." in gt:
                    continue
                raise ValueError("Unexpected diploid genotype")
            count = gt.count(b"1")
            if count in (0, n):
                continue
            blocks[(int(pos)-1) // BLOCK_BP][count] += 1
            tags = set(info.split(b";"))
            if b"AA=" + ref in tags:
                polarized[count] += 1
                aa_ref += 1
            elif b"AA=" + alt in tags:
                polarized[n-count] += 1
                aa_alt += 1
            else:
                missing_aa += 1

    counts = list(map(sum, zip(*blocks)))
    saved = list(map(int, (DATA / "bottleneck.final.fs").read_text().splitlines()[1].split()))
    if counts != saved:
        raise ValueError("Independent GT counts disagree with bcftools")
    if missing_aa:
        raise ValueError("Review incomplete polarization before applying full length scaling")
    write_fs(DATA / "bottleneck.final.polarized.fs", polarized)
    fit(DATA / "bottleneck.final.polarized.fs", RESULTS / "polarized", "polarized", LENGTH)
    with (RESULTS / "allele_qc.tsv").open("w") as handle:
        handle.write("observed_SNPs\tAA_equals_REF\tAA_equals_ALT\tother_AA\tvcf_md5\n")
        handle.write(f"{sum(counts)}\t{aa_ref}\t{aa_alt}\t{missing_aa}\t{EXPECTED_MD5}\n")
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        tmp = Path(temp)
        for b, block in enumerate(blocks):
            rest = [x-y for x,y in zip(counts, block)]
            write_fs(tmp/"leave_out.fs", rest)
            result = fit(tmp/"leave_out.fs", tmp/"fit", "folded", LENGTH - BLOCK_BP)
            rows.append([b*BLOCK_BP + 1,(b+1)*BLOCK_BP,sum(block),float(result["generations"]),
                         float(result["post_decline_fraction"])])
    with (RESULTS/"leave_one_block_out.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["excluded_start_1based","excluded_end_1based","SNPs_in_excluded_block",
                         "generations","post_decline_fraction"])
        writer.writerows(rows)
    for j, name in [(3,"generations"),(4,"post_decline_fraction")]:
        values = [row[j] for row in rows]
        print(name, "range",min(values),max(values),"jackknife_SE",statistics.stdev(values)*(N_BLOCKS-1)/math.sqrt(N_BLOCKS))


if __name__ == "__main__":
    main()
