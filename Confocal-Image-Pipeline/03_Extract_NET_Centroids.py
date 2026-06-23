#!/usr/bin/env python3
"""
03_Extract_NET_Centroids.py — Extract NET-event CENTROIDS (cz,cy,cx) WITHOUT re-running full 02.
Mirrors the EXACT NET segmentation logic of 02 (same parameters), but:
  - SKIPS Costes (the slowest part), SKIPS Cellpose (reuses the nuclei cache in labels/),
  - only exports coordinates + per-event morphology.
Output: results_stats/master_NET_events_xyz.csv
"""
import os, argparse
import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, zoom
from skimage import morphology
from skimage.filters import threshold_otsu, threshold_triangle, threshold_li
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops

BASE = "."
TIFF_DIR = os.path.join(BASE, "tiff_3d")
MANIFEST = os.path.join(BASE, "manifest_final.csv")
LABELS_DIR = os.path.join(BASE, "results_stats", "labels")
OUT = os.path.join(BASE, "results_stats", "master_NET_events_xyz.csv")
OUT_FOV = os.path.join(BASE, "results_stats", "master_fov_xyz.csv")

# ── MUST MATCH 02 ──
NET_BG_RADIUS_UM = 20.0
NET_THRESH_FACTOR = 0.7
NET_USE_LI = False
COLOC_TOLERANCE_VOX = 1
NET_CLOSING_RADIUS_UM = 1.0
NET_MIN_VOL_UM3 = 3.0

def norm01(img, bit):
    img = img.astype(np.float32); return img / {8:255.,12:4095.,16:65535.}.get(bit, max(1.0, img.max()))

def highpass(img, radius_um, vox):
    """High-pass = img - gaussian background. Background computed on a DOWNSAMPLED image then upscaled
    -> ~10x faster for large sigma, nearly identical (background is low-frequency)."""
    img = img.astype(np.float32)
    sig = np.array([radius_um/v if v>0 else 0 for v in vox], float)
    fac = np.maximum(1, np.floor(sig/4)).astype(int)        # downsample sao cho sigma_ds ~4px
    fac = np.minimum(fac, np.maximum(1, np.array(img.shape)//4))
    small = img[::fac[0], ::fac[1], ::fac[2]]
    bg_s = gaussian_filter(small, sigma=sig/fac, mode="reflect")
    bg = zoom(bg_s, np.array(img.shape)/np.array(bg_s.shape), order=1)
    if bg.shape != img.shape:                                # correct rounding offsets
        bg = bg[:img.shape[0], :img.shape[1], :img.shape[2]]
        pad = [(0, img.shape[i]-bg.shape[i]) for i in range(3)]
        bg = np.pad(bg, pad, mode="edge")
    out = img - bg; out[out<0]=0
    return out.astype(np.float32)

def seg_net(dapi, cith3, mpo, vox):
    vvol = float(np.prod(vox))
    try: tthr = threshold_triangle(dapi)
    except: tthr = dapi.mean()
    tissue = dapi > tthr
    try: dna_thr = threshold_otsu(dapi[tissue])
    except: dna_thr = dapi.mean()+dapi.std()
    mask_dna = dapi > dna_thr
    cith3_bg = highpass(cith3, NET_BG_RADIUS_UM, vox)
    mpo_bg = highpass(mpo, NET_BG_RADIUS_UM, vox)
    def thr_of(x):
        try:
            b = threshold_li(x[tissue]) if NET_USE_LI else threshold_otsu(x[tissue])
            return b*NET_THRESH_FACTOR
        except: return (x.mean()+2*x.std())*NET_THRESH_FACTOR
    mc = cith3_bg > thr_of(cith3_bg); mm = mpo_bg > thr_of(mpo_bg)
    if COLOC_TOLERANCE_VOX>0:
        mc = ndi.binary_dilation(mc, iterations=COLOC_TOLERANCE_VOX)
        mm = ndi.binary_dilation(mm, iterations=COLOC_TOLERANCE_VOX)
    net = mc & mm & mask_dna
    net = remove_small_objects(net, min_size=max(1,int(NET_MIN_VOL_UM3/vvol)))
    if NET_CLOSING_RADIUS_UM>0:
        cv = max(1,int(NET_CLOSING_RADIUS_UM/np.mean(vox)))
        net = morphology.binary_closing(net, morphology.ball(cv))
    return morphology.label(net), cith3_bg, mpo_bg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filter-combo", default=None)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    m = pd.read_csv(MANIFEST)
    m["tp"] = m.tiff_path.apply(lambda p: os.path.join(TIFF_DIR, os.path.basename(str(p))))
    m = m[(m.scene_type=="3D_stack") & (m.tp.apply(os.path.exists))]
    if args.filter_combo: m = m[m.combo==args.filter_combo]
    if args.limit: m = m.head(args.limit)
    print(f"{len(m)} FOV\n", flush=True)

    rows = []
    for _, r in m.iterrows():
        if pd.isna(r.ch_cith3_idx) or pd.isna(r.ch_mpo_idx) or r.ch_cith3_idx<0 or r.ch_mpo_idx<0:
            continue
        vox = (float(r.voxel_z_um), float(r.voxel_y_um), float(r.voxel_x_um)); vvol=float(np.prod(vox))
        data = tifffile.imread(r.tp); data = data[0] if data.ndim==5 else data
        dapi = norm01(data[int(r.ch_dapi_idx)], int(r.bit_depth))
        cith3 = norm01(data[int(r.ch_cith3_idx)], int(r.bit_depth))
        mpo = norm01(data[int(r.ch_mpo_idx)], int(r.bit_depth))
        net_lab, cith3_bg, mpo_bg = seg_net(dapi, cith3, mpo, vox)

        lp = os.path.join(LABELS_DIR, f"{r.image_id}_nuclei.tif")
        nuc_mask = (tifffile.imread(lp) > 0) if os.path.exists(lp) else None

        nz_total = dapi.shape[0]
        props = regionprops(net_lab, intensity_image=cith3_bg)
        n=0
        for p in props:
            c = p.coords  # (N,3) voxels — O(object), does NOT scan the volume
            zc = c[:,0]
            z_span = int(zc.max()-zc.min()+1)           # thickness along z (slices)
            n_zplanes = int(np.unique(zc).size)          # number of z-planes the object actually occupies
            # brightness uniformity along z (CitH3): low cv = even across z
            cith3_vals = cith3_bg[c[:,0],c[:,1],c[:,2]]
            zmean = pd.Series(cith3_vals).groupby(zc).mean()
            z_cv = float(zmean.std()/zmean.mean()) if len(zmean)>1 and zmean.mean()>0 else 0.0
            nuc_frac = float(nuc_mask[c[:,0],c[:,1],c[:,2]].mean()) if nuc_mask is not None else np.nan
            loc = "intracellular/perinuclear" if (not np.isnan(nuc_frac) and nuc_frac>0.5) else "extracellular"
            rows.append(dict(image_id=r.image_id, combo=r.combo, condition=r.condition, case_id=r.case_id,
                net_id=int(p.label), cz=float(p.centroid[0]), cy=float(p.centroid[1]), cx=float(p.centroid[2]),
                volume_um3=float(p.area*vvol), solidity=float(p.solidity),
                z_span=z_span, n_zplanes=n_zplanes, nz_total=nz_total, z_cv=round(z_cv,3),
                mean_cith3=float(p.intensity_mean),
                mean_mpo=float(mpo_bg[c[:,0],c[:,1],c[:,2]].mean()),
                localization=loc))
            n+=1
        print(f"  {r.image_id:25s} NET={n}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\n-> {len(df)} NET events (with centroid) -> {OUT}")
    if len(df):
        print(df.groupby(['combo','condition']).agg(nFOV=('image_id','nunique'),nNET=('net_id','count')).to_string())

if __name__=="__main__": main()
