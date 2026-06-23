#!/usr/bin/env python3
"""
02_Segment_Quantify_NETs.py
───────────────────────────────────────────────────────────────────────────────
Confocal analysis pipeline for NETs.
Updated per spec: NETs = CitH3 + MPO + DAPI (regardless of intra- or extracellular).

Input: manifest_final.csv
Output:
  - per_image_csv/ (event-level data and null distribution)
  - master_NET_events.csv
  - master_fov.csv
  - qc_overlays/
"""

import os
import sys
import time
import glob
import warnings
import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr as _pearsonr, spearmanr as _spearmanr

from skimage import morphology
from skimage.filters import threshold_otsu, threshold_triangle, threshold_li, gaussian
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops

from cellpose import models as cp_models
from csbdeep.utils import normalize as cd_normalize

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "."
MANIFEST_FILENAME = "manifest_final.csv"
TIFF_SUBFOLDER = "tiff_3d"

MANIFEST_PATH = os.path.join(BASE_DIR, MANIFEST_FILENAME)
TIFF_DIR = os.path.join(BASE_DIR, TIFF_SUBFOLDER)

OUT_DIR = os.path.join(BASE_DIR, "results_stats")
TMP_DIR = os.path.join(BASE_DIR, "temp_processing")

for sub in ["labels", "qc_overlays", "per_image_csv"]:
    os.makedirs(os.path.join(OUT_DIR, sub), exist_ok=True)
    os.makedirs(os.path.join(TMP_DIR, sub), exist_ok=True)

OUT_DIRS = {
    "labels": os.path.join(OUT_DIR, "labels"),
    "qc": os.path.join(OUT_DIR, "qc_overlays"),
    "per_image_csv": os.path.join(OUT_DIR, "per_image_csv"),
}

# Cellpose
CELLPOSE_DIAMETER_UM = 15.0
CELLPOSE_STITCH_THRESHOLD = 0.25

# Lewy Body Settings
LEWY_MIN_VOL_UM3 = 10.0 

# NET Detection Settings
NET_BG_RADIUS_UM      = 20.0   
NET_THRESH_FACTOR     = 0.7    
NET_USE_LI            = False     
COLOC_TOLERANCE_VOX   = 1      
NET_CLOSING_RADIUS_UM = 1.0    
NET_MIN_VOL_UM3       = 3.0   
# Monte Carlo Settings
N_RANDOM_POINTS_PER_FOV = 1000

# ═══════════════════════════════════════════════════════════════════════════════
# UTILS
# ═══════════════════════════════════════════════════════════════════════════════

def normalize_to_float01(img, bit_depth):
    img = img.astype(np.float32)
    if bit_depth == 8:
        return img / 255.0
    if bit_depth == 12:
        return img / 4095.0
    if bit_depth == 16:
        return img / 65535.0
    a_min, a_max = float(img.min()), float(img.max())
    if a_max - a_min < 1e-9:
        return np.zeros_like(img, dtype=np.float32)
    return ((img - a_min) / (a_max - a_min)).astype(np.float32)

def load_ome_tiff(tiff_path, ch_dapi, ch_cith3, ch_mpo, ch_asyn, bit_depth):
    with tifffile.TiffFile(str(tiff_path)) as tf:
        data = tf.asarray()
    if data.ndim == 5:
        data = data[0]
    
    result = {}
    ch_mapping = {
        "dapi": ch_dapi,
        "cith3": ch_cith3,
        "mpo": ch_mpo,
        "asyn": ch_asyn
    }
    for ch_name, ch_idx in ch_mapping.items():
        if pd.isna(ch_idx) or ch_idx < 0:
            result[ch_name] = None
        else:
            result[ch_name] = normalize_to_float01(data[int(ch_idx)], bit_depth)
    return result

def gaussian_high_pass_subtract(img, radius_um, voxel_size_um):
    sigma_z = radius_um / voxel_size_um[0] if voxel_size_um[0] > 0 else 0
    sigma_y = radius_um / voxel_size_um[1] if voxel_size_um[1] > 0 else 0
    sigma_x = radius_um / voxel_size_um[2] if voxel_size_um[2] > 0 else 0
    sigmas = (sigma_z, sigma_y, sigma_x) if img.ndim == 3 else (sigma_y, sigma_x)
    bg = gaussian_filter(img.astype(np.float32), sigma=sigmas, mode="reflect")
    out = img - bg
    out[out < 0] = 0
    return out.astype(np.float32)

_cellpose_model = None
def get_cellpose_model():
    global _cellpose_model
    if _cellpose_model is None:
        import torch
        device = "mps" if torch.backends.mps.is_available() else "cpu"
        _cellpose_model = cp_models.CellposeModel(gpu=(device != "cpu"))
    return _cellpose_model

def segment_nuclei_stardist(dapi_data, voxel_size_um):
    model = get_cellpose_model()
    diam_px = CELLPOSE_DIAMETER_UM / voxel_size_um[1]
    dapi_norm = cd_normalize(dapi_data, 1, 99.8)
    masks, _, _ = model.eval(
        dapi_norm,
        diameter=diam_px,
        do_3D=False,
        stitch_threshold=CELLPOSE_STITCH_THRESHOLD,
        z_axis=0,
    )
    if masks is None:
        return np.zeros_like(dapi_data, dtype=np.uint16)
    return masks.astype(np.uint16)

# ═══════════════════════════════════════════════════════════════════════════════
# COLOCALIZATION (Manders & Costes per FOV)
# ═══════════════════════════════════════════════════════════════════════════════

def _block_shuffle_2d(arr_2d, block_size, rng):
    h, w = arr_2d.shape
    nh, nw = h // block_size, w // block_size
    if nh == 0 or nw == 0:
        flat = arr_2d.flatten()
        rng.shuffle(flat)
        return flat.reshape(arr_2d.shape)
    blocks = []
    for i in range(nh):
        for j in range(nw):
            blocks.append(arr_2d[i * block_size : (i + 1) * block_size, j * block_size : (j + 1) * block_size].copy())
    idx = list(range(len(blocks)))
    rng.shuffle(idx)
    blocks = [blocks[i] for i in idx]
    out = arr_2d.copy()
    k = 0
    for i in range(nh):
        for j in range(nw):
            out[i * block_size : (i + 1) * block_size, j * block_size : (j + 1) * block_size] = blocks[k]
            k += 1
    return out

def coloc_costes(ch_a, ch_b, mask, observed_r, n_permutations=200, block_size=16):
    if mask.sum() < 100 or np.isnan(observed_r): return np.nan, np.nan
    coords = np.where(mask)
    if len(coords[0]) == 0: return np.nan, np.nan
    bbox = tuple(slice(c.min(), c.max() + 1) for c in coords)
    a_crop = ch_a[bbox].copy()
    b_crop = ch_b[bbox].copy()
    mask_crop = mask[bbox]
    rng = np.random.default_rng(seed=42)
    random_rs = []
    for _ in range(n_permutations):
        b_shuf = b_crop.copy()
        if b_shuf.ndim == 3:
            for z in range(b_shuf.shape[0]):
                b_shuf[z] = _block_shuffle_2d(b_shuf[z], block_size, rng)
        else:
            b_shuf = _block_shuffle_2d(b_shuf, block_size, rng)
        a_in, b_in = a_crop[mask_crop], b_shuf[mask_crop]
        if np.std(a_in) > 1e-9 and np.std(b_in) > 1e-9:
            r, _ = _pearsonr(a_in, b_in)
            random_rs.append(r)
    if len(random_rs) < 10: return np.nan, np.nan
    random_rs = np.array(random_rs)
    p_value = float((random_rs >= observed_r).sum() / len(random_rs))
    return p_value, float(random_rs.mean())

def fov_coloc_metrics(ch_a, ch_b, tissue_mask):
    if tissue_mask.sum() < 100:
        return np.nan, np.nan, np.nan, np.nan
    a_vals = ch_a[tissue_mask].astype(np.float64)
    b_vals = ch_b[tissue_mask].astype(np.float64)
    if np.std(a_vals) < 1e-9 or np.std(b_vals) < 1e-9:
        return np.nan, np.nan, np.nan, np.nan
    
    # Pearson
    r, _ = _pearsonr(a_vals, b_vals)
    
    # Manders
    try:
        t_a = threshold_otsu(a_vals)
        t_b = threshold_otsu(b_vals)
    except:
        t_a, t_b = a_vals.mean(), b_vals.mean()
    a_above, b_above = a_vals > t_a, b_vals > t_b
    sum_a, sum_b = a_vals.sum(), b_vals.sum()
    m1 = float(a_vals[b_above].sum() / sum_a) if sum_a > 0 else np.nan
    m2 = float(b_vals[a_above].sum() / sum_b) if sum_b > 0 else np.nan
    
    # Costes
    costes_p, _ = coloc_costes(ch_a, ch_b, tissue_mask, r)
    return r, m1, m2, costes_p

# ═══════════════════════════════════════════════════════════════════════════════
# PROCESS ONE
# ═══════════════════════════════════════════════════════════════════════════════

def process_one_image(row, out_dirs, overwrite):
    image_id = row["image_id"]
    csv_out_path = os.path.join(out_dirs["per_image_csv"], f"{image_id}_nets.csv")
    json_out_path = os.path.join(out_dirs["per_image_csv"], f"{image_id}_summary.json")
    
    if not overwrite and os.path.exists(csv_out_path) and os.path.exists(json_out_path):
        return None, None, "skipped_existing"

    voxel_size = (float(row["voxel_z_um"]), float(row["voxel_y_um"]), float(row["voxel_x_um"]))
    voxel_vol_um3 = np.prod(voxel_size)
    
    # Load TIFF
    channels = load_ome_tiff(
        row["tiff_path"], 
        row.get("ch_dapi_idx", -1), 
        row.get("ch_cith3_idx", -1), 
        row.get("ch_mpo_idx", -1), 
        row.get("ch_asyn_idx", -1),
        int(row["bit_depth"])
    )
    
    dapi = channels["dapi"]
    cith3 = channels["cith3"]
    mpo = channels["mpo"]
    asyn = channels["asyn"]
    condition = row.get("condition", "")
    
    if dapi is None or cith3 is None or mpo is None:
        return None, None, "Missing DAPI, CitH3, or MPO channel"
        
    # Tissue mask (very low thresh to exclude empty space)
    try:
        tissue_thresh = threshold_triangle(dapi)
    except:
        tissue_thresh = dapi.mean()
    tissue_mask = dapi > tissue_thresh
    tissue_vol_mm3 = tissue_mask.sum() * voxel_vol_um3 / 1e9
    
    # DAPI mask for actual DNA signal (higher threshold than tissue mask)
    try:
        # Use Otsu on the tissue area to find actual DNA signal
        dna_thresh = threshold_otsu(dapi[tissue_mask])
    except:
        dna_thresh = dapi.mean() + dapi.std()
    mask_dna = dapi > dna_thresh

    # BG sub
    cith3_bg = gaussian_high_pass_subtract(cith3, NET_BG_RADIUS_UM, voxel_size)
    mpo_bg = gaussian_high_pass_subtract(mpo, NET_BG_RADIUS_UM, voxel_size)
    
    # Thresholds for NET markers
    try:
        base_thresh = threshold_li(cith3_bg[tissue_mask]) if NET_USE_LI else threshold_otsu(cith3_bg[tissue_mask])
        cith3_thresh = base_thresh * NET_THRESH_FACTOR
    except:
        cith3_thresh = (cith3_bg.mean() + 2*cith3_bg.std()) * NET_THRESH_FACTOR
        
    try:
        base_thresh = threshold_li(mpo_bg[tissue_mask]) if NET_USE_LI else threshold_otsu(mpo_bg[tissue_mask])
        mpo_thresh = base_thresh * NET_THRESH_FACTOR
    except:
        mpo_thresh = (mpo_bg.mean() + 2*mpo_bg.std()) * NET_THRESH_FACTOR
        
    mask_cith3 = cith3_bg > cith3_thresh
    mask_mpo = mpo_bg > mpo_thresh
    
    # Dung sai co-loc: dilate mask
    if COLOC_TOLERANCE_VOX > 0:
        mask_cith3_d = ndi.binary_dilation(mask_cith3, iterations=COLOC_TOLERANCE_VOX)
        mask_mpo_d = ndi.binary_dilation(mask_mpo, iterations=COLOC_TOLERANCE_VOX)
    else:
        mask_cith3_d = mask_cith3
        mask_mpo_d = mask_mpo
        
    # NET Mask = CitH3 ∩ MPO ∩ DNA
    net_raw = mask_cith3_d & mask_mpo_d & mask_dna
    
    # Minimum physical size
    min_net_voxels = max(1, int(NET_MIN_VOL_UM3 / voxel_vol_um3))
    net_clean = remove_small_objects(net_raw, min_size=min_net_voxels)
    
    # Closing to connect pieces
    if NET_CLOSING_RADIUS_UM > 0:
        closing_vox = max(1, int(NET_CLOSING_RADIUS_UM / np.mean(voxel_size)))
        net_clean = morphology.binary_closing(net_clean, morphology.ball(closing_vox))
        
    net_labels = morphology.label(net_clean)
    
    # Cellpose Nuclei for localization flag
    label_path = os.path.join(out_dirs["labels"], f"{image_id}_nuclei.tif")
    if os.path.exists(label_path):
        nuclei_labels = tifffile.imread(label_path)
    else:
        dapi_smooth = gaussian(dapi, sigma=1.5)
        nuclei_labels = segment_nuclei_stardist(dapi_smooth, voxel_size)
        tifffile.imwrite(label_path, nuclei_labels, compression="zlib")
        
    # Lewy body segmentation (Only for Combo3 / PD cases)
    lewy_labels = None
    lewy_coords = []
    if asyn is not None and row.get("combo","") == "Combo3" and condition == "PD":
        asyn_bg = gaussian_high_pass_subtract(asyn, 3.0, voxel_size)
        try: asyn_thresh = threshold_otsu(asyn_bg[tissue_mask])
        except: asyn_thresh = asyn_bg.mean() + 3*asyn_bg.std()
        
        # Physical size filter for Lewy bodies
        min_lewy_voxels = int(LEWY_MIN_VOL_UM3 / voxel_vol_um3)
        lewy_clean = remove_small_objects(asyn_bg > asyn_thresh, min_size=min_lewy_voxels)
        lewy_labels = morphology.label(lewy_clean)
        
        if lewy_labels.max() > 0:
            l_props = regionprops(lewy_labels)
            lewy_coords = np.array([p.centroid for p in l_props])
        
    # Coloc stats FOV
    r_val, m1, m2, c_p = fov_coloc_metrics(cith3_bg, mpo_bg, tissue_mask)
    
    # Extract NET events
    events = []
    props = regionprops(net_labels, intensity_image=cith3_bg)
    
    # True NET events
    for p in props:
        mask_obj = net_labels == p.label
        vol_voxels = mask_obj.sum()
        vol_um3 = vol_voxels * voxel_vol_um3
        
        # Intra vs Extra
        nuc_overlap = (mask_obj & (nuclei_labels > 0)).sum()
        loc = "intracellular/perinuclear" if (nuc_overlap / vol_voxels) > 0.5 else "extracellular"
        
        # Morphology (solidity as inverse fragmentation)
        solidity = float(p.solidity) if p.solidity is not None else np.nan
        
        # Nearest Lewy body
        dist_to_lewy = np.nan
        if len(lewy_coords) > 0:
            c = np.array(p.centroid)
            dists = np.sqrt(np.sum(((lewy_coords - c) * voxel_size)**2, axis=1))
            dist_to_lewy = float(np.min(dists))
            
        events.append({
            "image_id": image_id,
            "event_type": "True_NET",
            "net_id": p.label,
            "volume_um3": vol_um3,
            "solidity": solidity,
            "localization": loc,
            "mean_cith3": float(cith3_bg[mask_obj].mean()),
            "mean_mpo": float(mpo_bg[mask_obj].mean()),
            "dist_to_nearest_lewy_um": dist_to_lewy
        })
        
    # Monte Carlo Random Events (Null Distribution for Spatial)
    if len(lewy_coords) > 0:
        tissue_indices = np.argwhere(tissue_mask)
        if len(tissue_indices) > 0:
            np.random.seed(42) # For reproducibility
            # Sample N random points
            rand_idx = np.random.choice(len(tissue_indices), size=min(N_RANDOM_POINTS_PER_FOV, len(tissue_indices)), replace=False)
            rand_points = tissue_indices[rand_idx]
            
            for i, rand_c in enumerate(rand_points):
                rand_dists = np.sqrt(np.sum(((lewy_coords - rand_c) * voxel_size)**2, axis=1))
                dist_to_lewy = float(np.min(rand_dists))
                
                events.append({
                    "image_id": image_id,
                    "event_type": "Random_Point",
                    "net_id": f"rand_{i}",
                    "volume_um3": np.nan,
                    "solidity": np.nan,
                    "localization": "random",
                    "mean_cith3": np.nan,
                    "mean_mpo": np.nan,
                    "dist_to_nearest_lewy_um": dist_to_lewy
                })

    events_df = pd.DataFrame(events)
    if len(events_df) > 0:
        events_df.to_csv(csv_out_path, index=False)
        
    # FOV Summary
    n_true_nets = len([e for e in events if e["event_type"] == "True_NET"])
    fov_summary = {
        "image_id": image_id,
        "combo": row.get("combo", ""),
        "condition": condition,
        "case_id": row.get("case_id", row.get("replicate", "")),
        "tissue_vol_mm3": tissue_vol_mm3,
        "n_net_events": n_true_nets,
        "net_density_per_mm3": n_true_nets / tissue_vol_mm3 if tissue_vol_mm3 > 0 else 0,
        "net_vol_fraction_pct": net_clean.sum() / tissue_mask.sum() * 100 if tissue_mask.sum() > 0 else 0,
        "cith3_mpo_pearson_r": r_val,
        "cith3_mpo_manders_M1": m1,
        "cith3_mpo_manders_M2": m2,
        "cith3_mpo_costes_p": c_p,
    }
    
    with open(json_out_path, "w") as f:
        json.dump(fov_summary, f)
    
    # QC Overlay
    try:
        has_lewy = asyn is not None and row.get("combo","") == "Combo3"
        n_panels = 4 if has_lewy else 3
        fig, axes = plt.subplots(1, n_panels, figsize=(5*n_panels, 5))
        if n_panels == 3:
            axes_arr = axes
        else:
            axes_arr = axes
            
        axes_arr[0].imshow(dapi.max(axis=0), cmap="gray")
        axes_arr[0].set_title("DAPI")
        axes_arr[0].axis("off")
        
        rgb = np.zeros((*dapi.shape[1:], 3))
        rgb[..., 0] = mpo_bg.max(axis=0) / (mpo_bg.max() + 1e-9)
        rgb[..., 1] = cith3_bg.max(axis=0) / (cith3_bg.max() + 1e-9)
        rgb[..., 2] = dapi.max(axis=0) / (dapi.max() + 1e-9)
        rgb = np.clip(rgb * 2.5, 0, 1)
        axes_arr[1].imshow(rgb)
        axes_arr[1].set_title("Merged (R=MPO, G=CitH3, B=DAPI)")
        axes_arr[1].axis("off")
        
        axes_arr[2].imshow(net_labels.max(axis=0) > 0, cmap="magma")
        axes_arr[2].set_title(f"NET Objects (n={n_true_nets})")
        axes_arr[2].axis("off")
        
        if has_lewy:
            asyn_max = asyn.max(axis=0)
            axes_arr[3].imshow(asyn_max, cmap="viridis")
            axes_arr[3].set_title("αSyn (max-proj)")
            axes_arr[3].axis("off")
            if lewy_labels is not None and lewy_labels.max() > 0:
                axes_arr[3].contour(lewy_labels.max(axis=0) > 0, colors='red', linewidths=0.5)
                
        plt.tight_layout()
        plt.savefig(os.path.join(out_dirs["qc"], f"{image_id}_qc.png"), dpi=120)
        plt.close(fig)
    except Exception as e:
        print(f"QC overlay failed: {e}")

    return events_df, fov_summary, None

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, help="Process only first N images for testing")
    parser.add_argument("--filter-combo", type=str, help="Filter by combo (e.g., 'Combo3')")
    parser.add_argument("--filter-image", type=str, help="Filter by image name")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing results")
    args = parser.parse_args()

    if not os.path.exists(MANIFEST_PATH):
        print(f"Error: Manifest not found at {MANIFEST_PATH}")
        return
        
    df = pd.read_csv(MANIFEST_PATH)
    df["tiff_path"] = df["tiff_path"].apply(lambda p: os.path.join(TIFF_DIR, os.path.basename(str(p))) if pd.notna(p) and p else "")
    df = df[df["scene_type"] == "3D_stack"]
    if "export_status" in df.columns:
        df = df[df["export_status"] != "excluded"]
    df = df[df["tiff_path"].apply(lambda p: os.path.exists(p))]
    
    if getattr(args, "filter_combo", None):
        df = df[df["combo"] == args.filter_combo]
    if getattr(args, "filter_image", None):
        df = df[df["source_lif"].str.contains(args.filter_image)]
    if getattr(args, "limit", None):
        df = df.head(args.limit)
        
    print(f"Found {len(df)} 3D scenes to process.")
    
    all_events = []
    all_fovs = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        events_df, summary, err = process_one_image(row, OUT_DIRS, args.overwrite)
        if err == "skipped_existing":
            # Load existing
            image_id = row["image_id"]
            csv_out_path = os.path.join(OUT_DIRS["per_image_csv"], f"{image_id}_nets.csv")
            json_out_path = os.path.join(OUT_DIRS["per_image_csv"], f"{image_id}_summary.json")
            try:
                events_df = pd.read_csv(csv_out_path)
                for c in ["combo", "condition"]:
                    events_df[c] = row.get(c, "")
                events_df["case_id"] = row.get("case_id", row.get("replicate", ""))
                all_events.append(events_df)
                
                if os.path.exists(json_out_path):
                    with open(json_out_path, "r") as f:
                        summary = json.load(f)
                    summary["combo"] = row.get("combo", "")
                    summary["condition"] = row.get("condition", "")
                    summary["case_id"] = row.get("case_id", row.get("replicate", ""))
                    all_fovs.append(summary)
            except: pass
            continue
        elif err:
            print(f"Failed {row['image_id']}: {err}")
            continue
            
        if events_df is not None and len(events_df) > 0:
            for c in ["combo", "condition"]:
                events_df[c] = row.get(c, "")
            events_df["case_id"] = row.get("case_id", row.get("replicate", ""))
            all_events.append(events_df)
        all_fovs.append(summary)
        
    if all_events:
        pd.concat(all_events, ignore_index=True).to_csv(os.path.join(OUT_DIR, "master_NET_events.csv"), index=False)
    if all_fovs:
        pd.DataFrame(all_fovs).to_csv(os.path.join(OUT_DIR, "master_fov.csv"), index=False)
        
    print(f"✅ Pipeline complete. Outputs in: {OUT_DIR}")

if __name__ == "__main__":
    main()
