#!/usr/bin/env python3
"""Compare a folded-SFS fit with a cross-fitted covariance-weighted polarized SFS."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/optimization"
MU = 3.42e-8
N = 114
I = np.arange(1, N)
# Retain rare allele resolution and pool progressively noisier high-count bins.
ENDS = np.array([1, 2, 3, 4, 6, 9, 14, 21, 32, 46, 60, 79, 99, 113])
FOLDED = np.eye(N // 2)[np.minimum(I, N-I)-1].T
GROUPED = np.eye(len(ENDS))[np.searchsorted(ENDS, I)].T
log_bounds = np.log([[0.0001, 1.0], [1e-5, 20.0]])
drift = I * (N-I)/2
eigenvalues, vectors = eigh_tridiagonal(-2*drift, np.sqrt(drift[:-1]*drift[1:]))
old_basis = (vectors / np.sqrt(drift)[:, None]) * (vectors.T @ (np.sqrt(drift)/I))
METHODS = ("folded", "block_gls")
POISSON_TRANSFORMS = {"folded": FOLDED, "polarized": np.eye(N-1), "grouped": GROUPED}
GLS_METHODS = ("block_gls", "adaptive_gls")


def precision(blocks, method):
    cov = np.cov(blocks, rowvar=False)
    var = np.diag(cov)
    if np.any(var <= 0):
        raise ValueError("A grouped frequency class has no block variance")
    if method == "adaptive_gls":
        from sklearn.covariance import ledoit_wolf
        # Standardize before shrinkage so the target preserves each bin's units.
        sd = np.sqrt(var)
        z = (blocks - blocks.mean(axis=0)) / sd
        shrunk, _ = ledoit_wolf(z, assume_centered=True)
        cov = shrunk * np.outer(sd, sd)
    else:
        cov = 0.5*cov + 0.5*np.diag(var)
    return np.linalg.solve(cov, np.eye(len(var)))


def fit(blocks, length, method):
    blocks = np.asarray(blocks)
    if (blocks.ndim != 2 or blocks.shape[1] != N-1 or len(blocks) == 0 or
            not np.isfinite(blocks).all() or (blocks < 0).any() or blocks.sum() <= 0):
        raise ValueError("Expected nonnegative block spectra with 113 bins and positive total")
    if not np.isfinite(length) or length <= 0:
        raise ValueError("The sequence length must be positive and finite")
    if method not in POISSON_TRANSFORMS and method not in GLS_METHODS:
        raise ValueError(f"Unknown fitting method: {method}")
    poisson = method in POISSON_TRANSFORMS
    transform = POISSON_TRANSFORMS.get(method, GROUPED)
    eq = transform @ (1/I)
    basis = transform @ old_basis
    x = blocks @ transform.T
    if poisson:
        y = x.sum(axis=0)
        L = length
    else:
        a, b = x[::2], x[1::2]
        if len(a) != len(b) or len(a) < 2:
            raise ValueError("Use at least four equal-length genomic blocks in even numbers")
        means = (a.mean(axis=0), b.mean(axis=0))
        weights = (precision(b, method), precision(a, method))
        W = weights[0] + weights[1]
        z = weights[0] @ means[0] + weights[1] @ means[1]
        L = length/len(x)

    def evaluate(p, detailed=False):
        frac, tau = np.exp(p)
        e = np.exp(eigenvalues*tau)
        old = basis @ e
        s = frac*eq + (1-frac)*old
        jac = np.stack([frac*(eq-old), (1-frac)*(basis @ (e*eigenvalues*tau))])
        if poisson:
            theta = y.sum()/s.sum()
            expected = theta*s
            nz = y > 0
            loss = 2*(y[nz] @ np.log(y[nz]/expected[nz]))
            grad = 2*jac @ (theta-y/s)
        else:
            theta = (s @ z)/(s @ W @ s)
            if theta <= 0:
                return (1e10, np.zeros(2)) if not detailed else None
            residuals = [u-theta*s for u in means]
            loss = sum(r @ w @ r for r,w in zip(residuals,weights))
            grad = 2*theta*(jac @ (theta*(W @ s)-z))
        if detailed:
            na = theta/(4*MU*L)
            return dict(generations=2*na*frac*tau, fraction=frac, ancestral_Ne=na,
                        current_Ne=na*frac, objective=loss)
        return loss, grad

    starts = [(0.5, 0.03), (0.6, 0.2), (0.45, 1.2)]
    fits = [minimize(evaluate, np.log(start), jac=True, bounds=log_bounds,
                     method="L-BFGS-B", options={"ftol":1e-13,"gtol":1e-9,"maxiter":350})
            for start in starts]
    # A line search can report success on the flat no-decline boundary or
    # reject a much better interior result because of roundoff. Independently
    # refine the best interior grid point and the best gradient-search result.
    grid = [np.log([f,t]) for f in (0.08,0.35,0.65,0.9)
            for t in (0.003,0.03,0.3,3.0)]
    initial = min(grid,key=lambda p:evaluate(p)[0])
    finite = [r for r in fits if np.isfinite(r.fun) and np.isfinite(r.x).all()
              and evaluate(r.x,True) is not None]
    candidates = [initial]
    if finite:
        candidates.append(min(finite,key=lambda r:r.fun).x)
    verified = []
    for point in candidates:
        opt = minimize(lambda p:evaluate(p)[0],point,method="Nelder-Mead",bounds=log_bounds,
                       options={"xatol":1e-8,"fatol":1e-8,"maxiter":1200})
        if opt.success and np.isfinite(opt.fun) and evaluate(opt.x,True) is not None:
            verified.append(opt)
    if not verified:
        raise RuntimeError("Unable to independently verify a finite minimum")
    best = min(verified,key=lambda r:r.fun)
    if finite and best.fun > min(r.fun for r in finite) + 1e-7*max(1,abs(best.fun)):
        raise RuntimeError("Independent optimizer could not reproduce the lower objective")
    result = evaluate(best.x, True)
    result["boundary"] = bool(np.any(np.minimum(best.x-log_bounds[:,0], log_bounds[:,1]-best.x) < 1e-5))
    return result


def prediction_error(pred, t, f):
    return np.sqrt(((pred["generations"]/t-1)**2 + (pred["fraction"]/f-1)**2)/2)


def table(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulations", type=Path)
    args = parser.parse_args()
    with (args.simulations/"scenarios.tsv").open() as f:
        scenarios = {int(r["scenario"]):r for r in csv.DictReader(f, delimiter="\t")}
    sim_rows = []
    scores = {}
    for sid, scenario in scenarios.items():
        truth_t, truth_f = float(scenario["T"]), float(scenario["fraction"])
        group = []
        for rep in range(1, 21):
            blocks = np.loadtxt(args.simulations/f"sim_{sid}_{rep:02d}.tsv", dtype=int)
            if blocks.shape != (100,113):
                raise ValueError("Incomplete simulation")
            for method in METHODS:
                pred = fit(blocks, 25_000_000, method)
                error = prediction_error(pred, truth_t, truth_f)
                group.append([rep,method,error])
                sim_rows.append(dict(scenario=sid, replicate=rep, method=method,
                                     true_generations=truth_t, true_fraction=truth_f,
                                     **pred, relative_rmse=error))
        scores[sid] = {m: np.array([row[2] for row in group if row[1]==m]) for m in METHODS}
    table(OUT/"simulation_fits.tsv", sim_rows)
    summary = []
    for sid in [*scenarios,"pooled"]:
        for m in METHODS:
            errors = np.concatenate([scores[i][m] for i in scenarios]) if sid=="pooled" else scores[sid][m]
            summary.append(dict(scenario=sid, method=m, replicates=len(errors),
                                pooled_relative_rmse=np.sqrt(np.mean(errors**2)),
                                median_relative_rmse=np.median(errors)))
    table(OUT/"simulation_summary.tsv",summary)
    official = []
    for phase in ["testing","final"]:
        blocks = np.loadtxt(OUT/f"{phase}_blocks.tsv",dtype=int)
        for method in METHODS:
            pred = fit(blocks,100_000_000,method)
            err = prediction_error(pred,4757,0.52) if phase=="testing" else ""
            official.append(dict(phase=phase,method=method,**pred,testing_relative_rmse=err))
    table(OUT/"official_fits.tsv",official)
    for idx,phase in enumerate(("testing","final")):
        original_path = ROOT/("results/baseline/testing/raw_folded/fit.tsv" if phase=="testing" else "results/baseline/final/fit.tsv")
        with original_path.open() as f:
            original = next(csv.DictReader(f, delimiter="\t"))
        matched = official[2*idx]
        np.testing.assert_allclose([matched["generations"],matched["fraction"]],
            [float(original["generations"]),float(original["post_decline_fraction"])],rtol=2e-5)
    wins = sum(np.mean(scores[i]["block_gls"]**2) < np.mean(scores[i]["folded"]**2) for i in scenarios)
    sim_improved = summary[-1]["pooled_relative_rmse"] < summary[-2]["pooled_relative_rmse"]
    test_improved = official[1]["testing_relative_rmse"] < official[0]["testing_relative_rmse"]
    interior = not official[1]["boundary"]  # The Final input is reserved for prediction.
    selected = "block_gls" if (wins >=3 and sim_improved and test_improved and interior) else "folded"
    primary = next(row for row in official if row["phase"]=="final" and row["method"]==selected)
    table(OUT/"decision.tsv",[dict(selected=selected, improved_scenarios=wins, scenarios=len(scenarios),
        simulation_improved=sim_improved, testing_improved=test_improved, testing_interior=interior,
        generations=primary["generations"], fraction=primary["fraction"])])
    print("selected", selected, "scenarios improved",wins)
    for r in official: print(r)
    print("pooled simulations",summary[-2:])
    print("simulation boundary count",sum(r["boundary"] for r in sim_rows))


if __name__=="__main__":
    main()
