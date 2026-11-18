#!/usr/bin/env python3
"""Write the selected submission, spatial sensitivity results, and validation figure."""
from pathlib import Path
from tempfile import NamedTemporaryFile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compare import fit, ROOT, OUT, table
from review import score_release


def main():
    method, official, review_summary = score_release()
    final = next(r for r in official if r["phase"]=="final" and r["method"]==method)
    blocks = np.loadtxt(OUT/"final_blocks.tsv",dtype=int)
    sensitivity = []
    for phase in ("testing", "final"):
        for scheme, stem in (("all_records","blocks"), ("last_record_per_position","last_per_position_blocks")):
            b = np.loadtxt(OUT/f"{phase}_{stem}.tsv",dtype=int)
            p = fit(b,100_000_000,method)
            sensitivity.append(dict(phase=phase, records=scheme, variants=int(b.sum()),
                                    generations=p["generations"], fraction=p["fraction"]))
    table(OUT/"coordinate_sensitivity.tsv",sensitivity)
    jack = []
    for width in (2,10):
        values = []
        for start in range(0,100,width):
            b = np.concatenate((blocks[:start],blocks[start+width:]))
            p = fit(b,1_000_000*len(b),method)
            values.append([p["generations"],p["fraction"]])
            jack.append(dict(excluded_Mb=width, first_base=start*1_000_000+1,
                             last_base=(start+width)*1_000_000,generations=p["generations"],
                             fraction=p["fraction"]))
        values = np.array(values)
        if not np.isfinite(values).all():
            raise ValueError("Non-finite jackknife replicate")
    table(OUT/"optimized_jackknife.tsv",jack)
    summary = []
    for w in (2,10):
        rows = [r for r in jack if r["excluded_Mb"]==w]
        for field in ("generations","fraction"):
            a = np.array([r[field] for r in rows])
            summary.append(dict(block_Mb=w,parameter=field,n=len(a),leave_out_min=a.min(),
                                leave_out_max=a.max(),jackknife_SE=a.std(ddof=1)*(len(a)-1)/np.sqrt(len(a))))
    table(OUT/"optimized_jackknife_summary.tsv",summary)
    # Write a single canonical submission only after the validation has succeeded.
    dest = ROOT/"submission"
    dest.mkdir(exist_ok=True)
    yaml = ("description: Single-population bottleneck demography.\n"
            "time_units: generations\n"
            "post_decline_fraction: ratio_of_modern_individuals_to_ancient_individuals\n"
            "parameters:\n"
            f"  generations: {final['generations']:.1f}\n"
            f"  post_decline_fraction: {final['fraction']:.6f}\n")
    fig, (left,right) = plt.subplots(1, 2, figsize=(11,4.7))
    labels = {"folded":"Folded Poisson", "polarized":"Polarized Poisson",
              "grouped":"Grouped Poisson", "block_gls":"Fixed GLS", "adaptive_gls":"Adaptive GLS"}
    selected_color = "#007A78"
    for j, row in enumerate(r for r in review_summary if r["scenario"]=="all"):
        key = row["method"]
        val = row["relative_rmse"]*100
        left.barh(j,val,color=selected_color if key==method else "#7891A6",height=0.65)
        left.text(val+0.12,j,f"{val:.2f}",va="center",fontsize=9)
    left.set_yticks(range(len(labels)),labels.values())
    left.invert_yaxis()
    left.margins(x=0.18)
    left.set_xlabel("Pooled relative RMSE (%)")
    left.set_title("A   40 independent 100 Mb chromosomes",fontsize=11,loc="left")
    two = [r for r in jack if r["excluded_Mb"]==2]
    right.scatter([float(r["generations"]) for r in two],[float(r["fraction"])*100 for r in two],
                  s=22, color="#77B5AD", alpha=.8, label="Leave out 2 Mb")
    for key,marker in (("folded","s"),("block_gls","D"),("adaptive_gls","^")):
        rec = next(r for r in official if r["phase"]=="final" and r["method"]==key)
        label = labels[key] + (" (selected)" if key==method else "")
        right.scatter(rec["generations"],100*rec["fraction"],s=65,marker=marker,
                      color=selected_color if key==method else "#534D65",label=label,zorder=4)
    right.set_xlabel("Generations since decline")
    right.set_ylabel("Population retained (%)")
    right.set_title("B   Final predictions and block sensitivity",fontsize=12,loc="left")
    right.legend(frameon=False,fontsize=9,loc="upper left")
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(length=3)
    fig.tight_layout(pad=1.6)
    fig.savefig(OUT/"validation.png",dpi=170)
    plt.close(fig)
    with NamedTemporaryFile(mode="w",dir=dest,delete=False) as stage:
        stage.write(yaml)
    Path(stage.name).replace(dest/"bottleneck.final.yaml")


if __name__=="__main__":
    main()
