#!/usr/bin/env python3
"""Brain NET-complex ELISA, 3 NET-specific markers + composite (R2 Minor#2, R1#1).
Markers: MPO-DNA, NE-DNA, CitH3-DNA (new). n=3/group (Control, PD, MSA).
Robust frequentist (Kruskal-Wallis + Welch pairwise + Holm) + percentile bootstrap CI
+ Bayesian-bootstrap posterior P(MSA>PD>HC). Small n reported transparently."""
import numpy as np, pandas as pd, itertools, openpyxl
from scipy import stats
np.random.seed(42)

SRC = "data/net.brain.xlsx"
HERE = "."

wb = openpyxl.load_workbook(SRC, data_only=True)
ws = wb["Sheet1"]
rows = [r for r in ws.iter_rows(values_only=True)][1:10]
df = pd.DataFrame(rows, columns=["id","dx","mpo_dna","ne_dna","cith3_dna","net","_a","_b","_c","_d"]).iloc[:, :6]
df["dx"] = df["dx"].astype(str).str.strip().replace({"Control":"HC"})
df["dx"] = pd.Categorical(df["dx"], categories=["HC","PD","MSA"], ordered=True)
df = df.rename(columns={"net":"composite"})
df.to_csv(HERE + "/brain_elisa_3markers_data.csv", index=False)

markers = {"mpo_dna":"MPO-DNA","ne_dna":"NE-DNA","cith3_dna":"CitH3-DNA","composite":"Composite NET index"}
out=[]
def log(s=""): out.append(str(s)); print(s)
log("# Brain NET-complex ELISA — 3 markers + composite (n=3/group)\n")
log(df[["id","dx","mpo_dna","ne_dna","cith3_dna","composite"]].round(4).to_string(index=False))

def holm(ps):
    order=np.argsort(ps); m=len(ps); adj=[None]*m; run=0
    for rank,idx in enumerate(order):
        run=max(run,(m-rank)*ps[idx]); adj[idx]=min(run,1.0)
    return adj

def bayes_bootstrap_order(g):
    # Dirichlet(1,...) weights -> posterior of weighted group means; P(MSA>PD>HC)
    B=4000; hc,pd,ms=g["HC"],g["PD"],g["MSA"]
    def draw(x):
        w=np.random.dirichlet(np.ones(len(x)),B); return (w*x).sum(1)
    mh,mp,mm=draw(hc),draw(pd),draw(ms)
    return float(np.mean((mm>mp)&(mp>mh))), float(np.mean(mp>mh)), float(np.mean(mm>mh))

summary=[]
for col,name in markers.items():
    log("\n" + "="*60); log(f"## {name}")
    g={lv:df.loc[df.dx==lv,col].astype(float).values for lv in ["HC","PD","MSA"]}
    desc=df.groupby("dx",observed=True)[col].agg(["mean","std","min","max"]).round(4)
    log(desc.to_string())
    H,pkw=stats.kruskal(*g.values()); F,pan=stats.f_oneway(*g.values())
    log(f"Kruskal-Wallis p={pkw:.4f} | one-way ANOVA p={pan:.4f}")
    pairs=list(itertools.combinations(["HC","PD","MSA"],2)); raw=[]; rws=[]
    for a,b in pairs:
        t,p=stats.ttest_ind(g[a],g[b],equal_var=False); raw.append(p)
        rws.append((f"{a}-{b}", g[a].mean()-g[b].mean(), t, p))
    hp=holm(raw)
    for (c,d,t,p),h in zip(rws,hp):
        # percentile bootstrap CI on diff
        diffs=np.array([np.random.choice(g[c.split('-')[0]],3,True).mean()-np.random.choice(g[c.split('-')[1]],3,True).mean() for _ in range(5000)])
        lo,hi=np.percentile(diffs,[2.5,97.5])
        log(f"  {c:7s} diff={d:+.4f} Welch t={t:+.2f} p={p:.4f} p_Holm={h:.4f} | boot95%CI[{lo:+.3f},{hi:+.3f}]")
    p_ord,p_pd,p_msa=bayes_bootstrap_order(g)
    log(f"  Bayesian bootstrap: P(MSA>PD>HC)={p_ord:.3f} | P(PD>HC)={p_pd:.3f} | P(MSA>HC)={p_msa:.3f}")
    summary.append({"marker":name,"HC":round(g['HC'].mean(),3),"PD":round(g['PD'].mean(),3),
                    "MSA":round(g['MSA'].mean(),3),"KW_p":round(pkw,4),
                    "P_ordered":round(p_ord,3),"P_PDgtHC":round(p_pd,3),"P_MSAgtHC":round(p_msa,3)})

pd.DataFrame(summary).to_csv(HERE + "/brain_elisa_summary.csv", index=False)
log("\n## SUMMARY"); log(pd.DataFrame(summary).to_string(index=False))
log("\nNote: n=3/group; frequentist pairwise tests underpowered. Bayesian-bootstrap")
log("ordered-gradient probabilities provide an honest, distribution-free summary.")
open(HERE + "/brain_elisa_log.md","w").write("\n".join(out))

# ---- 4-panel figure ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
cols={"HC":"#4C72B0","PD":"#C44E52","MSA":"#8172B3"}
fig,axes=plt.subplots(1,4,figsize=(13,3.6))
for ax,(col,name) in zip(axes,markers.items()):
    g={lv:df.loc[df.dx==lv,col].astype(float).values for lv in ["HC","PD","MSA"]}
    for i,lv in enumerate(["HC","PD","MSA"]):
        y=g[lv]; ax.scatter(np.random.normal(i,0.05,len(y)),y,color=cols[lv],s=45,zorder=3,edgecolor="k",lw=0.4)
        ax.hlines(y.mean(),i-0.22,i+0.22,color=cols[lv],lw=2.3,zorder=2)
    H,pkw=stats.kruskal(*g.values())
    ax.set_xticks(range(3)); ax.set_xticklabels(["HC","PD","MSA"])
    ax.set_title(name,fontsize=10); ax.spines[["top","right"]].set_visible(False)
    ax.text(0.5,0.99,f"KW p={pkw:.3f}",ha="center",va="top",transform=ax.transAxes,fontsize=7.5,color="#555")
axes[0].set_ylabel("NET complex (relative units)")
fig.suptitle("Brain NET-complex ELISA (n=3/group) — individual markers + composite",fontsize=11,y=1.04)
plt.tight_layout()
plt.savefig(HERE + "/Figure4a_Brain_ELISA_3markers.pdf",bbox_inches="tight")
plt.savefig(HERE + "/Figure4a_Brain_ELISA_3markers.png",dpi=300,bbox_inches="tight")
print("\nSaved log, summary, figure.")
