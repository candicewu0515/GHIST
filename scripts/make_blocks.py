#!/usr/bin/env python3
"""Create 100 consecutive 1-Mb derived allele count spectra for the two releases."""
import csv
import gzip
import hashlib
from pathlib import Path
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results/optimization"
BASES = (b"A", b"C", b"G", b"T")


def make_blocks(phase):
    vcf = DATA / f"bottleneck.{phase}.vcf.gz"
    with (DATA / f"{phase}_synapse_manifest.csv").open() as handle:
        meta = next(csv.DictReader(handle))
    with vcf.open("rb") as handle:
        md5 = hashlib.file_digest(handle, "md5").hexdigest()
    if md5 != meta["dataFileMD5Hex"]:
        raise ValueError(f"{phase}: MD5 mismatch")
    counts = np.zeros((100, 113), dtype=np.int64)
    alt_spectrum = np.zeros(113, dtype=np.int64)
    aa_is_alt = 0
    with gzip.open(vcf, "rb") as handle:
        for line in handle:
            if line.startswith(b"#CHROM\t"):
                break
        samples = len(line.rstrip().split(b"\t")) - 9
        if samples != 57:
            raise ValueError(f"{phase}: unexpected ploidy/sample dimensions")
        gt_pattern = re.compile(rb"[01][|/][01](?:\t[01][|/][01]){56}")
        total = excluded = 0
        previous_pos = previous_ac = None
        overwritten = np.zeros_like(counts)
        for line in handle:
            chrom, pos, _, ref, alt, _, filt, info, fmt, gt = line.rstrip().split(b"\t", 9)
            total += 1
            pos = int(pos)
            if chrom != b"1" or not 1 <= pos <= 100_000_000:
                raise ValueError("Unexpected sequence/coordinate")
            if alt not in BASES or ref not in BASES or filt not in (b".", b"PASS"):
                excluded += 1
                continue
            if fmt != b"GT" or not gt_pattern.fullmatch(gt):
                raise ValueError("A complete, diploid GT is required")
            alt_copies = gt.count(b"1")
            ac = alt_copies
            tags = info.split(b";")
            if b"AA=" + alt in tags:
                ac = 2 * samples - ac
                if 0 < ac < 2 * samples:
                    aa_is_alt += 1
            elif b"AA=" + ref not in tags:
                raise ValueError("Missing ancestral state")
            if 0 < ac < 2 * samples:
                if previous_pos is not None and pos < previous_pos:
                    raise ValueError("Unsorted VCF")
                if pos == previous_pos:
                    overwritten[(pos-1)//1_000_000, previous_ac-1] += 1
                previous_pos, previous_ac = pos, ac
                counts[(pos-1)//1_000_000, ac-1] += 1
                alt_spectrum[alt_copies-1] += 1
            else:
                excluded += 1
    fs = DATA / f"bottleneck.{phase}.fs"
    reference = np.fromstring(fs.read_text().splitlines()[1], dtype=int, sep=" ")[1:-1]
    # bcftools reports ALT counts, including when ALT is ancestral.
    np.testing.assert_array_equal(alt_spectrum, reference)
    last_blocks = counts - overwritten
    last_counts = last_blocks.sum(axis=0)
    app_matches = ""
    if phase == "testing":
        app = DATA / "GHIST_2026_bottleneck_testing.fs"
        app_counts = np.fromstring(app.read_text().splitlines()[1], dtype=int, sep=" ")[1:-1]
        app_matches = np.array_equal(app_counts, last_counts)
        np.savetxt(OUT / "testing_source_comparison.tsv", np.column_stack(
            (np.arange(1,114),counts.sum(axis=0),last_counts,app_counts)), fmt="%d", delimiter="\t",
            header="derived_count\tVCF_records\tlast_record_per_position\torganizer_app",comments="")
    np.savetxt(OUT / f"{phase}_last_per_position_blocks.tsv", last_blocks, fmt="%d", delimiter="\t")
    np.savetxt(OUT / f"{phase}_blocks.tsv", counts, fmt="%d", delimiter="\t")
    return [phase, meta["ID"], meta["versionNumber"], md5, total, excluded, int(counts.sum()), int(overwritten.sum()), app_matches, aa_is_alt]


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    result = [make_blocks(p) for p in ("testing", "final")]
    with (OUT / "input_qc.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["round","synapse_id","version","md5","total_records","excluded","retained","extra_records_at_shared_positions","last_record_counts_equal_app","retained_AA_equals_ALT"])
        writer.writerows(result)
