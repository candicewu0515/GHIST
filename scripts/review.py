#!/usr/bin/env python3
"""Score fresh linked chromosomes, run matched ablations, and verify the release fit."""
import csv
import hashlib
from pathlib import Path
import numpy as np
from compare import fit, OUT, ROOT, prediction_error, table

REVIEW = ROOT / "results/review"
HOLDOUT = REVIEW / "chromosomes"
METHODS = ("folded", "polarized", "grouped", "block_gls", "adaptive_gls")
N_REPS = 10
LENGTH = 100_000_000


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_release():
    # The published Testing answers are used only to report the review metrics.
    # The incumbent was selected earlier. New promotion is based on holdout errors.
    decision = read_rows(OUT/"decision.tsv")
    if len(decision) != 1 or decision[0]["selected"] not in ("folded", "block_gls"):
        raise ValueError("Invalid first-stage decision")
    incumbent = decision[0]["selected"]
    cases = read_rows(HOLDOUT/"scenarios.tsv")
    if [int(s["scenario"]) for s in cases] != [1, 2, 3, 4]:
        raise ValueError("Expected four independent holdout scenarios")

    files_seen = set()
    rows = []
    scores = {m:np.empty((len(cases),N_REPS)) for m in METHODS}
    for k,case in enumerate(cases):
        sid = int(case["scenario"])
        truth_t,truth_f = float(case["T"]),float(case["fraction"])
        for rep in range(1,N_REPS+1):
            path = HOLDOUT/f"sim_{sid}_{rep:02d}.tsv"
            digest = hashlib.sha256(path.read_bytes()).digest()
            if digest in files_seen:
                raise ValueError("Duplicate holdout genome")
            files_seen.add(digest)
            blocks = np.loadtxt(path, dtype=int)
            if blocks.shape != (100,113):
                raise ValueError("Expected 100 one-megabase blocks")
            for method in METHODS:
                est = fit(blocks, LENGTH, method)
                loss = prediction_error(est, truth_t, truth_f)
                scores[method][k,rep-1] = loss
                rows.append(dict(scenario=sid,replicate=rep,seed=465201+100*sid+rep,
                    method=method,truth_generations=truth_t,truth_fraction=truth_f,
                    generations=est["generations"],fraction=est["fraction"],relative_rmse=loss,
                    boundary=est["boundary"]))
    table(REVIEW/"fits.tsv",rows)
    summary = []
    for method,errors in scores.items():
        if not np.isfinite(errors).all():
            raise ValueError("A holdout fit returned a non-finite error")
        for j,label in enumerate(["1","2","3","4","all"]):
            e = errors.ravel() if label=="all" else errors[j,:]
            rr = [r for r in rows if r["method"]==method and (label=="all" or str(r["scenario"])==label)]
            te = [r["generations"]/r["truth_generations"]-1 for r in rr]
            fe = [r["fraction"]/r["truth_fraction"]-1 for r in rr]
            summary.append(dict(method=method,scenario=label,n=len(e),
                relative_rmse=np.sqrt(np.mean(e**2)),time_relative_rmse=np.sqrt(np.mean(np.square(te))),
                fraction_relative_rmse=np.sqrt(np.mean(np.square(fe))),boundary_fits=sum(r["boundary"] for r in rr)))
    table(REVIEW/"summary.tsv",summary)
    challenger = "adaptive_gls"
    old, new = scores[incumbent], scores[challenger]
    improvement = 1-np.sqrt(np.mean(new**2)/np.mean(old**2))
    # Paired resampling by chromosome, stratified by its simulated history.
    rng = np.random.default_rng(106741)
    picks = rng.integers(0,N_REPS,(10000,len(cases),N_REPS))
    strata = np.arange(len(cases))[None,:,None]
    delta = 1-np.sqrt((new[strata,picks]**2).mean(axis=(1,2))/(old[strata,picks]**2).mean(axis=(1,2)))
    lo,hi = np.quantile(delta,[0.025,0.975])
    wins = int((np.mean(new**2,axis=1)<np.mean(old**2,axis=1)).sum())
    chosen = challenger if lo>0 and wins>=3 else incumbent

    # Compute both official predictions from the actual block files. Detect stale
    # first-stage tables before they can be used to issue a submission.
    old_official = read_rows(OUT/"official_fits.tsv")
    official = []
    for phase in ("testing","final"):
        blocks = np.loadtxt(OUT/f"{phase}_blocks.tsv",dtype=int)
        if blocks.shape!=(100,113):
            raise ValueError(f"Invalid block layout for {phase}")
        for method in METHODS:
            est = fit(blocks,LENGTH,method)
            if method in ("folded","block_gls"):
                existing = [r for r in old_official if r["method"]==method and r["phase"]==phase]
                if len(existing)!=1:
                    raise ValueError("Missing or duplicated first-stage prediction")
                np.testing.assert_allclose([est["generations"],est["fraction"]],
                    [float(existing[0]["generations"]),float(existing[0]["fraction"])],rtol=2e-6)
            official.append(dict(phase=phase,method=method,**est,
                testing_relative_rmse=prediction_error(est,4757,.52) if phase=="testing" else ""))
    table(REVIEW/"official_fits.tsv",official)
    released = next(r for r in official if r["phase"]=="final" and r["method"]==chosen)
    if released["boundary"]:
        raise ValueError("Final prediction is at a model boundary; time is unresolved")
    decision_row = dict(incumbent=incumbent,challenger=challenger,selected=chosen,replicates=40,
        scenarios_improved=wins,paired_wins=int((new<old).sum()),relative_gain=improvement,
        gain_interval_low=lo,gain_interval_high=hi,bootstrap_draws=10000,
        generations=released["generations"],fraction=released["fraction"])
    table(REVIEW/"decision.tsv",[decision_row])
    return chosen, official, summary


if __name__=="__main__":
    print(score_release()[0])
