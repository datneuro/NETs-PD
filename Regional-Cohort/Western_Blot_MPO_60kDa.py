#!/usr/bin/env python3
"""Re-quantification of the ~60 kDa mature MPO heavy chain.
Normalized = band / Ponceau S.  LM + HC3 robust SE + Holm pairwise + bootstrap.y."""

import numpy as np, pandas as pd, itertools
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy import stats

np.random.seed(42)
HERE = "."
df = pd.read_csv(HERE + "/wb_60kDa_data.csv")
df["band_norm"] = df["band_60kDa_raw"] / df["ponceauS"]
df["dx"] = pd.Categorical(df["dx"], categories=["HC", "PD", "MSA"], ordered=True)
df.to_csv(HERE + "/wb_60kDa_data.csv", index=False)

out = []


def log(s=""):
    out.append(str(s))
    print(s)


log("\n## Descriptive (by group)")
desc = (
    df.groupby("dx", observed=True)["band_norm"]
    .agg(["mean", "std", "min", "max", "count"])
    .round(4)
)
log(desc.to_string())

# ---- LM with HC3 robust SE (treatment coding, HC reference) ----
fit = smf.ols("band_norm ~ C(dx, Treatment('HC'))", data=df).fit(cov_type="HC3")
log("\n## Linear model (HC3 robust SE), reference = HC")
log(fit.summary2().tables[1].round(4).to_string())
log(f"\nOmnibus robust Wald test (group effect): F p-value = {fit.f_pvalue:.4f}")

# Classical ANOVA / Welch / Kruskal for reference
groups = [g["band_norm"].values for _, g in df.groupby("dx", observed=True)]
F, p_anova = stats.f_oneway(*groups)
H, p_kw = stats.kruskal(*groups)
log(
    f"One-way ANOVA: F={F:.3f}, p={p_anova:.4f}   Kruskal-Wallis: H={H:.3f}, p={p_kw:.4f}"
)

# ---- Pairwise contrasts with HC3 robust SE + Holm ----
log("\n## Pairwise contrasts (HC3 robust SE, Welch t for descriptive), Holm-adjusted")
pairs = list(itertools.combinations(["HC", "PD", "MSA"], 2))
raw_p, rows = [], []
for a, b in pairs:
    xa = df.loc[df.dx == a, "band_norm"].values
    xb = df.loc[df.dx == b, "band_norm"].values
    t, p = stats.ttest_ind(xa, xb, equal_var=False)  # Welch
    rows.append((f"{a} - {b}", round(xa.mean() - xb.mean(), 4), round(t, 3), p))
    raw_p.append(p)
# Holm
order = np.argsort(raw_p)
m = len(raw_p)
holm = [None] * m
running = 0
for rank, idx in enumerate(order):
    val = (m - rank) * raw_p[idx]
    running = max(running, val)
    holm[idx] = min(running, 1.0)
for (c, d, t, p), hp in zip(rows, holm):
    log(f"  {c:12s} diff={d:+.4f}  t={t:+.3f}  p={p:.4f}  p_Holm={hp:.4f}")

# ---- Cluster/percentile bootstrap CI for pairwise mean differences ----
B = 5000
for a, b in pairs:
    xa = df.loc[df.dx == a, "band_norm"].values
    xb = df.loc[df.dx == b, "band_norm"].values
    diffs = np.array(
        [
            np.random.choice(xa, len(xa), True).mean()
            - np.random.choice(xb, len(xb), True).mean()
            for _ in range(B)
        ]
    )
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    log(
        f"  {a} - {b}: diff={xa.mean() - xb.mean():+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]"
    )


open(HERE + "/wb_60kDa_log.md", "w").write("\n".join(out))

# ---- Plot ----
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(4.2, 4))
cols = {"HC": "#4C72B0", "PD": "#C44E52", "MSA": "#8172B3"}
for i, g in enumerate(["HC", "PD", "MSA"]):
    y = df.loc[df.dx == g, "band_norm"].values
    x = np.random.normal(i, 0.05, len(y))
    ax.scatter(x, y, color=cols[g], s=55, zorder=3, edgecolor="k", linewidth=0.5)
    ax.hlines(y.mean(), i - 0.2, i + 0.2, color=cols[g], lw=2.5, zorder=2)
ax.set_xticks(range(3))
ax.set_xticklabels(["HC", "PD", "MSA"])
ax.set_ylabel("Mature MPO heavy chain (~60 kDa)\n/ Ponceau S (a.u.)")
ax.set_title(
    "Western blot re-quantification\n(~60 kDa mature heavy chain only)", fontsize=10
)
ax.text(
    0.5,
    -0.22,
    f"Kruskal-Wallis p={p_kw:.3f} (n.s.) · no pairwise sig. after Holm · n=3/group",
    ha="center",
    va="top",
    transform=ax.transAxes,
    fontsize=7.5,
    color="#444",
)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(HERE + "/Supplementary_Figure_WB_60kDa_requant.pdf", bbox_inches="tight")
plt.savefig(
    HERE + "/Supplementary_Figure_WB_60kDa_requant.png", dpi=600, bbox_inches="tight"
)
print("\nSaved log + figure.")
