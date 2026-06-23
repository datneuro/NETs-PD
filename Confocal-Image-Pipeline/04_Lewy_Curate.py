#!/usr/bin/env python3
"""
04_Lewy_Curate.py — Human-in-the-loop CURATION of Lewy bodies (no manual tracing).
The machine proposes candidates -> exports a numbered montage + CSV -> the human just ticks correct/incorrect.

Output:
  results_stats/lewy_curation/<image_id>_montage.png   <- look at this
  results_stats/lewy_candidates.csv                    <- fill the 'confirm' column (1=Lewy, 0/blank=not)

auto_suggest column: 1 = machine is confident it is Lewy (pre-ticked); you only confirm.
Does NOT touch 02 / NET. Runs independently.
"""
import os, glob, argparse
import numpy as np, pandas as pd, tifffile
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, white_tophat
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, label, ball
from skimage.measure import regionprops

BASE = "."
TIFF_DIR = os.path.join(BASE, "tiff_3d")
MANIFEST = os.path.join(BASE, "manifest_final.csv")
OUTDIR = os.path.join(BASE, "results_stats", "lewy_curation")
CSV = os.path.join(BASE, "results_stats", "lewy_candidates.csv")
os.makedirs(OUTDIR, exist_ok=True)

# ── candidate Lewy criteria (permissive, to NOT miss weak Lewy) ──
CAND_MIN_UM3, CAND_MAX_UM3 = 5.0, 8000.0  
CAND_MIN_SOLIDITY = 0.20                    
CAND_MIN_ZSPAN = 3                         
TOPN_PER_FOV = 24                          
AUTO_CONTRAST, AUTO_SOLIDITY = 4.0, 0.30

def norm01(img, bit):
    img = img.astype(np.float32); return img / {8:255.,12:4095.,16:65535.}.get(bit, max(1.0,img.max()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contains", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--thresh-factor", type=float, default=1.0, help="Otsu x factor for the candidate mask")
    args = ap.parse_args()

    m = pd.read_csv(MANIFEST)
    m = m[(m.scene_type=="3D_stack") & (m.combo=="Combo3") & (m.condition=="PD")]
    if args.contains: m = m[m.image_id.apply(lambda s: any(c in s for c in args.contains))]
    if args.limit: m = m.head(args.limit)
    print(f"Processing {len(m)} PD FOVs\n", flush=True)

    all_rows = []
    for _, r in m.iterrows():
        tp = os.path.join(TIFF_DIR, os.path.basename(str(r.tiff_path)))
        ai = r.get("ch_asyn_idx", -1)
        if not os.path.exists(tp) or pd.isna(ai) or ai < 0:
            print(f"  skip {r.image_id}"); continue
        vox = (float(r.voxel_z_um), float(r.voxel_y_um), float(r.voxel_x_um)); vvol=float(np.prod(vox))
        data = tifffile.imread(tp);  data = data[0] if data.ndim==5 else data
        asyn = norm01(data[int(ai)], int(r.bit_depth))
        dapi = norm01(data[int(r.ch_dapi_idx)], int(r.bit_depth))
        tissue = dapi > dapi.mean()
        # white top-hat (enhances compact ~Lewy blobs, removes diffuse background) — radius ~ Lewy size
        th_rad = max(1, int(8.0 / vox[1]))   # 8µm structuring element (2D per z for speed)
        asyn_th = np.stack([white_tophat(asyn[z], size=th_rad) for z in range(asyn.shape[0])])
        try: thr = threshold_otsu(asyn_th[tissue]) * args.thresh_factor
        except: thr = asyn_th.mean() + asyn_th.std()
        cand = remove_small_objects(asyn_th > thr, min_size=max(1,int(CAND_MIN_UM3/vvol)))
        lab = label(cand)
        bg = float(np.median(asyn[tissue])) or 1e-6
        props = regionprops(lab, intensity_image=asyn)
        cands = []
        for p in props:
            area = p.area*vvol
            if not (CAND_MIN_UM3 <= area <= CAND_MAX_UM3): continue
            zspan = p.bbox[3]-p.bbox[0]
            if zspan < CAND_MIN_ZSPAN: continue
            sol = float(p.solidity)
            if sol < CAND_MIN_SOLIDITY: continue
            contrast = float(p.intensity_mean)/bg
            score = sol * (zspan/asyn.shape[0]) * np.log1p(area) * contrast
            cands.append(dict(image_id=r.image_id, case_id=r.case_id, area_um3=round(area,1),
                z_span=int(zspan), solidity=round(sol,2), contrast=round(contrast,2),
                cz=float(p.centroid[0]), cy=float(p.centroid[1]), cx=float(p.centroid[2]),
                score=round(float(score),2),
                auto_suggest=int(contrast>=AUTO_CONTRAST and sol>=AUTO_SOLIDITY)))
        cands = sorted(cands, key=lambda d:-d["score"])[:TOPN_PER_FOV]
        for i,c in enumerate(cands): c["cand_id"]=i+1
        print(f"  {r.image_id:22s} candidates={len(cands)}  auto_suggest={sum(c['auto_suggest'] for c in cands)}", flush=True)

        # montage
        if cands:
            n=len(cands); ncol=6; nrow=int(np.ceil(n/ncol))
            half=int(25.0/vox[1])  # crop ±25µm
            amax=asyn.max(axis=0)
            fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.2,nrow*2.4)); axes=np.array(axes).reshape(-1)
            for k,c in enumerate(cands):
                ax=axes[k]; cy,cx=int(c["cy"]),int(c["cx"])
                y0,y1=max(0,cy-half),cy+half; x0,x1=max(0,cx-half),cx+half
                ax.imshow(amax[y0:y1,x0:x1],cmap="magma"); ax.axis("off")
                col="lime" if c["auto_suggest"] else "white"
                ax.set_title(f"#{c['cand_id']} a{c['area_um3']:.0f} z{c['z_span']} s{c['solidity']:.2f} c{c['contrast']:.1f}",
                             fontsize=7,color=col)
                for s in ax.spines.values(): s.set_visible(False)
            for k in range(n,len(axes)): axes[k].axis("off")
            fig.suptitle(f"{r.image_id}  (blue outline=machine-suggested Lewy)  -> write the real cand_id into the confirm column",fontsize=10)
            plt.tight_layout(); plt.savefig(os.path.join(OUTDIR,f"{r.image_id}_montage.png"),dpi=110); plt.close(fig)
        all_rows += cands

    df=pd.DataFrame(all_rows)
    if df.empty: print("No candidates."); return
    df["confirm"]=""   # <- column for you to fill 1 (Lewy) / blank (not)
    cols=["image_id","case_id","cand_id","area_um3","z_span","solidity","contrast","score","auto_suggest","confirm","cz","cy","cx"]
    df[cols].to_csv(CSV,index=False)
    print(f"\n-> Montage: {OUTDIR}\n-> CSV to fill in: {CSV}")
    print(f"  Total candidates: {len(df)} | machine-suggested (auto_suggest=1): {int(df.auto_suggest.sum())}")
    print("  Usage: open the montage; for each #id that is a true Lewy -> write 1 in the 'confirm' column of the matching CSV row.")

if __name__=="__main__": main()
