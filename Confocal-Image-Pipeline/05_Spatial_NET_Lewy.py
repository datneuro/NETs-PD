#!/usr/bin/env python3
"""
05_Spatial_NET_Lewy.py — Spatial NET <-> Lewy (manually curated Lewy) + Monte Carlo null.
Input:
  results_stats/master_NET_events_xyz.csv  (NET centroids, from 03_Extract_NET_Centroids.py)
  results_stats/lewy_confirmed.csv         (manually curated Lewy from napari)
Output:
  results_stats/master_spatial_events.csv  (event-level: True_NET + Random_Point, dist_to_nearest_lewy_um)
  -> consumed by 06_NET_Statistics.R (Module 2)
Distances are in µm (scaled by voxel size). Only Combo3-PD FOVs contain Lewy.
"""
import os, argparse
import numpy as np, pandas as pd, tifffile
from skimage.filters import threshold_triangle

BASE = "."
TIFF_DIR = os.path.join(BASE, "tiff_3d")
MANIFEST = os.path.join(BASE, "manifest_final.csv")
NET_XYZ = os.path.join(BASE, "results_stats", "master_NET_events_xyz.csv")
LEWY = os.path.join(BASE, "results_stats", "lewy_confirmed.csv")
OUT = os.path.join(BASE, "results_stats", "master_spatial_events.csv")

N_RANDOM = 1000
SEED = 42

def norm01(img, bit):
    img = img.astype(np.float32); return img / {8:255.,12:4095.,16:65535.}.get(bit, max(1.0, img.max()))

def nearest_um(pts, lewy, vox):
    """pts (N,3), lewy (M,3), voxel -> nearest distance (µm) for each pt."""
    vox = np.asarray(vox)
    out = np.empty(len(pts))
    for i, p in enumerate(pts):
        d = np.sqrt((((lewy - p) * vox)**2).sum(axis=1))
        out[i] = d.min()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cith3-min", type=float, default=0.0, help="use only NETs with mean_cith3 >= threshold")
    ap.add_argument("--vol-min", type=float, default=0.0, help="use only NETs with volume_um3 >= threshold")
    args = ap.parse_args()

    net = pd.read_csv(NET_XYZ)
    n0 = len(net)
    net = net[(net.mean_cith3 >= args.cith3_min) & (net.volume_um3 >= args.vol_min)]
    print(f"Filter NETs: {n0} -> {len(net)} (cith3>={args.cith3_min}, vol>={args.vol_min})", flush=True)
    lewy = pd.read_csv(LEWY)
    man = pd.read_csv(MANIFEST)
    man["tp"] = man.tiff_path.apply(lambda p: os.path.join(TIFF_DIR, os.path.basename(str(p))))
    minfo = man.set_index("image_id").to_dict("index")

    rng = np.random.default_rng(SEED)
    rows = []
    fov_rows = []
    images = sorted(lewy.image_id.unique())   # only FOVs with curated Lewy (all Combo3-PD)
    print(f"{len(images)} FOVs with manually curated Lewy\n", flush=True)

    for iid in images:
        ldf = lewy[lewy.image_id == iid]
        lcoords = ldf[["cz","cy","cx"]].to_numpy()
        info = minfo.get(iid)
        if info is None: print(f"  skip {iid} (not found in manifest)"); continue
        vox = (float(info["voxel_z_um"]), float(info["voxel_y_um"]), float(info["voxel_x_um"]))
        case_id = info.get("case_id"); cond = info.get("condition")

        # NET centroids for this FOV
        ndf = net[net.image_id == iid]
        ncoords = ndf[["cz","cy","cx"]].to_numpy()

        # True NET → nearest Lewy
        if len(ncoords):
            d_true = nearest_um(ncoords, lcoords, vox)
            for d in d_true:
                rows.append(dict(image_id=iid, case_id=case_id, condition=cond,
                                 event_type="True_NET", dist_to_nearest_lewy_um=float(d)))
        else:
            d_true = np.array([])

        # Monte Carlo: random points in tissue (reload DAPI -> tissue mask)
        d_rand = np.array([])
        if os.path.exists(info["tp"]):
            data = tifffile.imread(info["tp"]); data = data[0] if data.ndim==5 else data
            dapi = norm01(data[int(info["ch_dapi_idx"])], int(info["bit_depth"]))
            try: tthr = threshold_triangle(dapi)
            except: tthr = dapi.mean()
            tidx = np.argwhere(dapi > tthr)
            if len(tidx):
                sel = rng.choice(len(tidx), size=min(N_RANDOM, len(tidx)), replace=False)
                rpts = tidx[sel]
                d_rand = nearest_um(rpts, lcoords, vox)
                for d in d_rand:
                    rows.append(dict(image_id=iid, case_id=case_id, condition=cond,
                                     event_type="Random_Point", dist_to_nearest_lewy_um=float(d)))

        fov_rows.append(dict(image_id=iid, case_id=case_id, n_lewy=len(lcoords), n_net=len(ncoords),
            mean_true_um=float(np.mean(d_true)) if len(d_true) else np.nan,
            mean_random_um=float(np.mean(d_rand)) if len(d_rand) else np.nan))
        print(f"  {iid:25s} Lewy={len(lcoords):2d} NET={len(ncoords):3d} "
              f"mean_true={np.mean(d_true):.1f}µm  mean_rand={np.mean(d_rand):.1f}µm" if len(d_true) and len(d_rand)
              else f"  {iid:25s} Lewy={len(lcoords)} NET={len(ncoords)} (missing NET or tissue)", flush=True)

    pd.DataFrame(rows).to_csv(OUT, index=False)
    fov = pd.DataFrame(fov_rows)
    print(f"\n→ {len(rows)} event → {OUT}\n")

    # FOV-level summary (primary, avoids pseudo-replication)
    ok = fov.dropna(subset=["mean_true_um","mean_random_um"])
    print(f"=== FOV-LEVEL (n={len(ok)} FOVs with both NET & random) ===")
    if len(ok) >= 2:
        from scipy.stats import wilcoxon
        diff = ok.mean_true_um - ok.mean_random_um
        print(f"  mean True (NET→Lewy):   {ok.mean_true_um.mean():.2f} µm")
        print(f"  mean Random (null):     {ok.mean_random_um.mean():.2f} µm")
        print(f"  Δ (true-random):        {diff.mean():.2f} µm  ({'NET CLOSER to Lewy' if diff.mean()<0 else 'farther'})")
        try:
            w = wilcoxon(ok.mean_true_um, ok.mean_random_um)
            print(f"  Wilcoxon paired (FOV): p = {w.pvalue:.4g}")
        except Exception as e:
            print(f"  (wilcoxon error: {e})")
    else:
        print("  not enough FOVs to test.")
    print("\n(formal stats + Bayesian are run in 06_NET_Statistics.R)")

if __name__=="__main__": main()
