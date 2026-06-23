"""
═══════════════════════════════════════════════════════════════════════════════
01_lif_to_tif_preprocess_NEWFIX.py
───────────────────────────────────────────────────────────────────────────────
Step 1 of the confocal analysis pipeline (revised to fix the Leica LIF frame-sequential bug).

TASKS:
  1. Scan the LIF folder, read metadata for all scenes in each file
  2. Classify 3D-stack (Z>1) vs single-plane (Z=1, usually 4K representative)
  3. Export 3D scenes to OME-TIFF (preserve voxel size + channel info)
     * FIXED Z-C CHANNEL MIXING BY USING AICSImage AND RE-INDEXING THE FRAME ORDER
  4. Generate manifest.csv for the next script (02_Segment_Quantify_NETs.py)
  5. Also generate manifest_template_metadata.csv (empty) for the user to fill
     condition/replicate/antibody/combo manually, then the script merges it

INPUT:
  --lif-dir   : folder containing the .lif files
  --out-dir   : output folder (creates subfolders: tiff_3d/, png_singleplane/)
  --condition-table (optional) : CSV mapping filename -> condition/antibody/...

OUTPUT:
  out-dir/
    ├── tiff_3d/<lif_stem>__scene<idx>.ome.tif    # 3D stacks for analysis
    ├── png_singleplane/<lif_stem>__scene<idx>.png # 4K reps (max projection)
    ├── manifest.csv                              # used by script 02
    ├── manifest_template_metadata.csv            # USER FILLS before running 02
    └── preprocess_log.txt                        # diagnostic log

DEPENDENCIES:
  pip install aicsimageio[bioformats] ome-types tifffile pillow
  Java (cho bioformats backend): brew install openjdk
  (On macOS: export JAVA_HOME=$(/usr/libexec/java_home))

NEWFIX v3 (2026-05-21):
  Uses AICSImage (Bio-Formats) to read pixel data and unscramble, fully replacing
  the readlif library, which fails with "seek of closed file" on large/broken files
  such as 3.HMC3-S33.lif.
  
  Verified: scene 3 of 4-C.lif and all scenes of 3.HMC3-S33.lif run perfectly.
═══════════════════════════════════════════════════════════════════════════════
"""

import argparse
import os
import sys
import re
import logging
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import tifffile

# Try to import aicsimageio
try:
    from aicsimageio import AICSImage
    from aicsimageio.writers import OmeTiffWriter
except ImportError as e:
    print(f"ERROR: Missing dependency: {e}")
    print("  pip install 'aicsimageio[bioformats]' tifffile pillow")
    print(
        "  Also need Java: brew install openjdk (macOS) or apt-get install default-jre"
    )
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════

# Channel name patterns — regex (case-insensitive)
CHANNEL_PATTERNS = {
    "DAPI": [r"dapi", r"hoechst", r"405", r"blue"],
    "CitH3": [r"cith3", r"citrullinated", r"h3cit", r"(488|alexa.?488|af488)"],
    "MPO": [r"mpo", r"(568|594|af568|af594)"],
    "aSyn": [r"(a-?syn)", r"syn1", r"synuclein", r"(647|af647|alexa.?647)"],
}

# Preferred channel order (if name matching fails, fall back to this order
# emission wavelength order): DAPI (405) → CitH3 (488) → MPO (568) → aSyn (647)
FALLBACK_CHANNEL_ORDER = ["DAPI", "CitH3", "MPO", "aSyn"]

# Z threshold: <= this is treated as single-plane
SINGLE_PLANE_Z_MAX = 1

# Output OME-TIFF compression
OME_COMPRESSION = "zlib"  # lossless, well-supported


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITIES & RE-INDEX READER
# ═══════════════════════════════════════════════════════════════════════════════


def setup_logging(log_path):
    """Setup logging to both console and file."""
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    # Clear any existing handlers
    logger.handlers = []
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def match_channel(name, patterns):
    """
    Match channel name (case-insensitive) against list of regex patterns.
    Returns True if any pattern matches.
    """
    if name is None:
        return False
    name_lower = str(name).lower()
    return any(re.search(p, name_lower) for p in patterns)


def map_channels(channel_names):
    """
    Map channel names → semantic labels ('DAPI', 'CitH3', 'MPO', 'aSyn').
    Returns: dict {channel_index: semantic_label}
             unmatched channels become 'UNKNOWN_<idx>'
    """
    result = {}
    used_labels = set()

    # Pass 1: pattern match
    for idx, name in enumerate(channel_names):
        for label, patterns in CHANNEL_PATTERNS.items():
            if label in used_labels:
                continue
            if match_channel(name, patterns):
                result[idx] = label
                used_labels.add(label)
                break

    # Pass 2: fill unmatched using the fallback order
    fallback_iter = iter([l for l in FALLBACK_CHANNEL_ORDER if l not in used_labels])
    for idx, name in enumerate(channel_names):
        if idx not in result:
            try:
                label = next(fallback_iter)
                result[idx] = label
                used_labels.add(label)
            except StopIteration:
                result[idx] = f"UNKNOWN_{idx}"

    return result


def detect_bit_depth(img_array):
    """Detect bit depth from numpy dtype + actual max value."""
    dtype = str(img_array.dtype)
    if dtype == "uint8":
        return 8
    if dtype == "uint16":
        # Check if 12-bit data stored as uint16 (max ~4095)
        max_val = int(img_array.max())
        if max_val < 4096:
            return 12  # likely 12-bit, common for Leica
        return 16
    if dtype == "uint32":
        return 32
    if "float" in dtype:
        return 32  # treat as float32
    return -1  # unknown


def read_lif_scene_corrected(lif_path, scene_idx, sequential_unscramble=False):
    """
    Read one scene from the LIF with a fix for frame-sequential acquisition order.
    Uses AICSImage (Bio-Formats backend) to avoid the "seek of closed file" error
    from readlif on large/structurally-broken files.
    
    Returns:
        arr: numpy array shape (C, Z, Y, X), dtype uint16/uint8
    """
    img = AICSImage(str(lif_path))
    img.set_scene(scene_idx)
    
    # Read raw data as CZYX
    data = img.get_image_data("CZYX")  # (C, Z, Y, X)
    
    if not sequential_unscramble:
        return data
        
    C, Z, Y, X = data.shape
    
    # Initialize the output array
    arr = np.zeros((C, Z, Y, X), dtype=data.dtype)
    
    # Apply the unscramble formula
    for c_lin in range(C):
        for z_lin in range(Z):
            p = c_lin * Z + z_lin       # linear position
            c_real = p % C              # true channel
            z_real = p // C             # true Z
            arr[c_real, z_real] = data[c_lin, z_lin]
            
    return arr


def verify_exported_tiff(tiff_path, expected_arr_czyx, logger):
    """Sanity-check that the newly written OME-TIFF matches the source data."""
    try:
        img = AICSImage(str(tiff_path))
        actual = img.get_image_data("CZYX")
        C, Z = actual.shape[:2]
        max_diff = 0
        for c in range(C):
            for z in range(Z):
                e = float(expected_arr_czyx[c, z].mean())
                a = float(actual[c, z].mean())
                d = abs(a - e)
                if d > max_diff:
                    max_diff = d
        if max_diff > 1.0:
            logger.warning(f"  VERIFY FAILED for {tiff_path.name}: max mean diff = {max_diff:.3f} > 1.0")
            return False
        logger.info(f"  Verify OK: max mean diff = {max_diff:.3f}")
        return True
    except Exception as e:
        logger.error(f"  Verification check raised an exception: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# LIF SCANNING
# ═══════════════════════════════════════════════════════════════════════════════


def scan_lif_file(lif_path, logger):
    """
    Read one LIF file, return a list of dicts (one per scene).
    Uses AICSImage to get metadata.
    """
    logger.info(f"--- Scanning: {lif_path.name}")
    scenes_info = []

    try:
        img = AICSImage(str(lif_path))
    except Exception as e:
        logger.error(f"Could not open {lif_path.name}: {e}")
        return scenes_info

    n_scenes = len(img.scenes)
    logger.info(f"  {n_scenes} scenes detected")

    for scene_idx, scene_name in enumerate(img.scenes):
        try:
            img.set_scene(scene_idx)
            # AICSImage dimensions: T, C, Z, Y, X
            dim_t = img.dims.T
            dim_c = img.dims.C
            dim_z = img.dims.Z
            dim_y = img.dims.Y
            dim_x = img.dims.X

            # Physical pixel sizes (μm) — may be None if metadata is missing
            ps = img.physical_pixel_sizes  # (Z, Y, X) in μm or None
            voxel_z = float(ps.Z) if ps.Z is not None else np.nan
            voxel_y = float(ps.Y) if ps.Y is not None else np.nan
            voxel_x = float(ps.X) if ps.X is not None else np.nan

            # Channel names
            ch_names = (
                list(img.channel_names)
                if img.channel_names
                else [f"Ch{i}" for i in range(dim_c)]
            )

            # Bit depth — read a small subset to check
            try:
                small = img.get_image_data("ZYX", T=0, C=0)
                bit_depth = detect_bit_depth(small)
            except Exception:
                bit_depth = -1

            # Map channels → semantic labels
            ch_map = map_channels(ch_names)

            # Find indices for DAPI/CitH3/MPO/aSyn
            ch_dapi_idx = next((i for i, l in ch_map.items() if l == "DAPI"), -1)
            ch_cith3_idx = next((i for i, l in ch_map.items() if l == "CitH3"), -1)
            ch_mpo_idx = next((i for i, l in ch_map.items() if l == "MPO"), -1)
            ch_asyn_idx = next((i for i, l in ch_map.items() if l == "aSyn"), -1)

            # Classify scene type
            scene_type = "3D_stack" if dim_z > SINGLE_PLANE_Z_MAX else "single_plane"

            info = {
                "source_lif": lif_path.name,
                "scene_idx": scene_idx,
                "scene_name": scene_name,
                "dim_t": dim_t,
                "dim_c": dim_c,
                "dim_z": dim_z,
                "dim_y": dim_y,
                "dim_x": dim_x,
                "voxel_z_um": voxel_z,
                "voxel_y_um": voxel_y,
                "voxel_x_um": voxel_x,
                "bit_depth": bit_depth,
                "channel_names_raw": "|".join(ch_names),
                "channel_mapping": "|".join(
                    f"{i}:{ch_map[i]}" for i in sorted(ch_map.keys())
                ),
                "ch_dapi_idx": ch_dapi_idx,
                "ch_cith3_idx": ch_cith3_idx,
                "ch_mpo_idx": ch_mpo_idx,
                "ch_asyn_idx": ch_asyn_idx,
                "scene_type": scene_type,
                "image_id": f"{lif_path.stem}__scene{scene_idx:02d}",
            }
            scenes_info.append(info)

            logger.info(
                f"  scene {scene_idx:2d}: {scene_type:12s} | "
                f"Z={dim_z:3d} YX={dim_y}x{dim_x} | C={dim_c} {ch_names} | "
                f"voxel=({voxel_z:.3f},{voxel_y:.3f},{voxel_x:.3f}) | "
                f"{bit_depth}-bit"
            )

        except Exception as e:
            logger.error(f"  scene {scene_idx} failed: {e}")
            continue

    return scenes_info


# ═══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════════


def export_scene_to_ome_tiff(lif_path, scene_idx, out_path, scene_info, logger, sequential_unscramble=False):
    """
    Export one scene from LIF -> OME-TIFF with metadata preserved.
    Saved as TCZYX (T=1 for fixed-cell imaging).
    """
    try:
        # Read the corrected numpy array (C, Z, Y, X)
        arr_czyx = read_lif_scene_corrected(lif_path, scene_idx, sequential_unscramble)
        data = arr_czyx[np.newaxis, ...]  # (1, C, Z, Y, X)

        # Get metadata from AICSImage (channel names and physical pixel sizes)
        img = AICSImage(str(lif_path))
        img.set_scene(scene_idx)
        ps = img.physical_pixel_sizes

        # Channel names cho OME metadata
        ch_names = (
            list(img.channel_names)
            if img.channel_names
            else [f"Ch{i}" for i in range(data.shape[1])]
        )

        # Write OME-TIFF
        OmeTiffWriter.save(
            data=data,
            uri=str(out_path),
            dim_order="TCZYX",
            channel_names=ch_names,
            physical_pixel_sizes=ps,
        )

        # Verify
        verify_exported_tiff(out_path, arr_czyx, logger)

        file_size_mb = out_path.stat().st_size / 1024 / 1024
        logger.info(f"  → exported: {out_path.name} ({file_size_mb:.1f} MB)")
        return True

    except Exception as e:
        logger.error(f"  EXPORT FAILED for scene {scene_idx}: {e}")
        return False


def export_single_plane_preview(lif_path, scene_idx, out_path, scene_info, logger, sequential_unscramble=False):
    """
    For single-plane scenes (Z=1, usually 4K representative): export a PNG preview
    (max projection if Z is present, or raw 2D if Z=1).
    Multi-channel -> combined into an RGB display (DAPI=blue, CitH3=green, MPO=red).
    """
    try:
        from PIL import Image

        data = read_lif_scene_corrected(lif_path, scene_idx, sequential_unscramble)  # (C, Z, Y, X)

        # MIP if Z>1, else squeeze
        if data.shape[1] > 1:
            data_2d = data.max(axis=1)  # (C, Y, X)
        else:
            data_2d = data[:, 0, :, :]  # (C, Y, X)

        # Build RGB from the 3 main channels (DAPI/CitH3/MPO)
        # Create RGB array (Y, X, 3) — note R=MPO, G=CitH3, B=DAPI
        rgb = np.zeros((data_2d.shape[1], data_2d.shape[2], 3), dtype=np.uint8)

        def norm_to_uint8(arr):
            """Min-max scale to 0-255 cho display only (not for quantification!)."""
            a_min, a_max = float(arr.min()), float(arr.max())
            if a_max - a_min < 1:
                return np.zeros_like(arr, dtype=np.uint8)
            return np.clip((arr - a_min) / (a_max - a_min) * 255, 0, 255).astype(
                np.uint8
            )

        ch_dapi_idx = scene_info["ch_dapi_idx"]
        ch_cith3_idx = scene_info["ch_cith3_idx"]
        ch_mpo_idx = scene_info["ch_mpo_idx"]

        if 0 <= ch_mpo_idx < data_2d.shape[0]:
            rgb[..., 0] = norm_to_uint8(data_2d[ch_mpo_idx])
        if 0 <= ch_cith3_idx < data_2d.shape[0]:
            rgb[..., 1] = norm_to_uint8(data_2d[ch_cith3_idx])
        if 0 <= ch_dapi_idx < data_2d.shape[0]:
            rgb[..., 2] = norm_to_uint8(data_2d[ch_dapi_idx])

        # Downsample if too large (4K = 4096; a 2K PNG preview is enough)
        from PIL import Image as PILImage

        pil_img = PILImage.fromarray(rgb, mode="RGB")
        max_dim = 2048
        if max(pil_img.size) > max_dim:
            scale = max_dim / max(pil_img.size)
            new_size = (int(pil_img.size[0] * scale), int(pil_img.size[1] * scale))
            pil_img = pil_img.resize(new_size, PILImage.LANCZOS)

        pil_img.save(str(out_path), "PNG", optimize=True)
        logger.info(f"  → preview PNG: {out_path.name}")
        return True

    except Exception as e:
        logger.error(f"  PREVIEW FAILED for scene {scene_idx}: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# METADATA TEMPLATE GENERATION
# ═══════════════════════════════════════════════════════════════════════════════


def create_metadata_template(manifest_df, out_path, logger):
    template = manifest_df[["image_id", "source_lif"]].copy()
    template["combo"] = ""
    template["antibody"] = ""
    template["condition"] = ""
    template["case_id"] = ""
    template["notes"] = ""
    template["magnification_class"] = ""
    template["scene_subtype"] = ""
    template.to_csv(out_path, index=False)
    logger.info(f"Metadata template written (Scene-level): {out_path}")

    logger.info(
        f"  → FILL the empty columns, then rerun with --metadata-csv "
        f"{out_path.name} --merge-only"
    )


def merge_metadata_into_manifest(manifest_path, metadata_csv_path, logger):
    """
    Read manifest.csv + metadata.csv (user-filled), merge -> manifest_final.csv.
    """
    manifest = pd.read_csv(manifest_path)
    meta = pd.read_csv(metadata_csv_path)

    missing = [c for c in ["image_id", "combo", "condition", "case_id"] if c not in meta.columns]
    if missing:
        logger.warning(f"Metadata CSV is missing important columns: {missing}")
        # Not returning None, we can still proceed with what we have
    empty_rows = meta[meta[["combo", "condition", "case_id"]].isna().any(axis=1)]
    if len(empty_rows) > 0:
        logger.warning(f"{len(empty_rows)} images have incomplete metadata")
    if len(bad_cond) > 0:
        logger.warning(f"Invalid conditions (expected {valid_conditions}):")
        logger.warning(bad_cond[[id_col, "condition"]].to_string())

    # Merge
    merged = manifest.merge(
        meta[["image_id", "combo", "antibody", "condition", "case_id", "notes", "magnification_class", "scene_subtype"]],
        on="image_id",
        how="left",
    )

    out_path = manifest_path.parent / "manifest_final.csv"
    merged.to_csv(out_path, index=False)
    logger.info(f"Final manifest written: {out_path}")
    logger.info(f"  Total scenes: {len(merged)}")
    logger.info(
        f"  3D stacks (for quantification): "
        f"{(merged['scene_type'] == '3D_stack').sum()}"
    )
    logger.info(
        f"  Single-plane (representative only): "
        f"{(merged['scene_type'] == 'single_plane').sum()}"
    )

    if "condition" in merged.columns:
        analyzable = merged[
            (merged["scene_type"] == "3D_stack")
            & merged["condition"].notna()
            & (merged["condition"] != "")
        ]
        if len(analyzable) > 0:
            logger.info("\n=== Breakdown 3D scenes per group ===")
            breakdown = (
                analyzable.groupby(["combo", "condition", "case_id"])
                .size()
                .reset_index(name="n_scenes")
            )
            logger.info("\n" + breakdown.to_string(index=False))

    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess LIF files → OME-TIFF + manifest for HMC3 aSyn analysis (Corrected Frame Order version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--lif-dir", type=str, required=True, help="Folder containing the .lif files"
    )
    parser.add_argument(
        "--out-dir", type=str, required=True, help="Output folder (creates subfolders)"
    )
    parser.add_argument(
        "--metadata-csv",
        type=str,
        default=None,
        help="(optional) filled metadata CSV to merge into the manifest",
    )
    parser.add_argument(
        "--merge-only",
        action="store_true",
        help="Only merge metadata, do not re-scan/re-export",
    )
    parser.add_argument(
        "--skip-export",
        action="store_true",
        help="Only scan + generate manifest, do not export TIFF (debug)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing TIFFs"
    )
    parser.add_argument("--limit", type=int, default=None, help="Process only first N images for testing")
    parser.add_argument("--filter-combo", type=str, default=None, help="Filter by combo (e.g., 'Combo3')")
    parser.add_argument("--filter-image", type=str, default=None, help="Filter by image name")
    parser.add_argument("--disable-unscramble", action="store_true", help="Disable sequential unscramble")
    args = parser.parse_args()

    lif_dir = Path(args.lif_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tiff_dir = out_dir / "tiff_3d"
    png_dir = out_dir / "png_singleplane"
    tiff_dir.mkdir(exist_ok=True)
    png_dir.mkdir(exist_ok=True)

    log_path = out_dir / "preprocess_log.txt"
    logger = setup_logging(log_path)

    logger.info("=" * 70)
    logger.info(f"HMC3 aSyn — LIF preprocessing | {datetime.now()}")
    logger.info("=" * 70)
    logger.info(f"LIF dir:   {lif_dir}")
    logger.info(f"Out dir:   {out_dir}")
    logger.info(
        f"Mode:      {'merge-only' if args.merge_only else 'full scan + export'}"
    )

    manifest_path = out_dir / "manifest.csv"

    # ─── MERGE-ONLY MODE ───
    if args.merge_only:
        if not manifest_path.exists():
            logger.error(
                f"manifest.csv does not exist at {manifest_path}. "
                "Run full mode first (without --merge-only)."
            )
            sys.exit(1)
        if args.metadata_csv is None:
            logger.error("--merge-only requires --metadata-csv")
            sys.exit(1)
        metadata_path = Path(args.metadata_csv).expanduser().resolve()
        if not metadata_path.exists():
            logger.error(f"Metadata CSV does not exist: {metadata_path}")
            sys.exit(1)
        merge_metadata_into_manifest(manifest_path, metadata_path, logger)
        return

    # ─── FULL MODE: scan + export ───
    if not lif_dir.exists():
        logger.error(f"LIF dir does not exist: {lif_dir}")
        sys.exit(1)

    lif_files = sorted(
        [f for f in lif_dir.rglob("*.lif") if not f.name.startswith("._")]
    )
    if len(lif_files) == 0:
        logger.error(f"No .lif files found in {lif_dir}")
        sys.exit(1)
    logger.info(f"Found {len(lif_files)} LIF file(s)")

    # Scan all files
    all_scenes = []
    for lif_path in lif_files:
        scenes = scan_lif_file(lif_path, logger)
        all_scenes.extend(scenes)

    if len(all_scenes) == 0:
        logger.error("No scenes could be scanned — exit.")
        sys.exit(1)

    manifest_df = pd.DataFrame(all_scenes)

    # ─── EXPORT ───
    if not args.skip_export:
        logger.info("\n" + "=" * 70)
        logger.info("EXPORTING 3D stacks → OME-TIFF; single-planes → PNG preview")
        logger.info("=" * 70)

        # Merge metadata right before export if it exists to allow filtering by combo
        if args.metadata_csv:
            metadata_path = Path(args.metadata_csv).expanduser().resolve()
            if metadata_path.exists():
                meta = pd.read_csv(metadata_path)
                merge_cols = ["image_id", "combo", "antibody", "condition", "case_id", "notes", "magnification_class", "scene_subtype"]
                actual_cols = [c for c in merge_cols if c in meta.columns]
                manifest_df = manifest_df.merge(
                    meta[actual_cols],
                    on="image_id",
                    how="left"
                )

        export_df = manifest_df.copy()
        if getattr(args, "filter_combo", None):
            export_df = export_df[export_df["combo"] == args.filter_combo]
        if getattr(args, "filter_image", None):
            export_df = export_df[export_df["source_lif"].str.contains(args.filter_image)]
        if getattr(args, "limit", None):
            export_df = export_df.head(args.limit)

        export_status = []
        for _, row in export_df.iterrows():
            # Resolve recursive path
            lif_path = lif_dir / row["source_lif"]
            if not lif_path.exists():
                matches = list(lif_dir.rglob(row["source_lif"]))
                if matches:
                    lif_path = matches[0]

            image_id = row["image_id"]

            if row["scene_type"] == "3D_stack":
                out_tiff = tiff_dir / f"{image_id}.ome.tif"
                if out_tiff.exists() and not args.overwrite:
                    logger.info(f"  SKIP existing: {out_tiff.name}")
                    export_status.append("skipped_existing")
                    continue
                ok = export_scene_to_ome_tiff(
                    lif_path, row["scene_idx"], out_tiff, row.to_dict(), logger, not args.disable_unscramble
                )
                export_status.append("exported" if ok else "failed")
            else:
                # single-plane → PNG preview
                out_png = png_dir / f"{image_id}.png"
                if out_png.exists() and not args.overwrite:
                    export_status.append("skipped_existing")
                    continue
                ok = export_single_plane_preview(
                    lif_path, row["scene_idx"], out_png, row.to_dict(), logger, not args.disable_unscramble
                )
                export_status.append("preview" if ok else "failed")

        # Initialize export_status if missing
        if "export_status" not in manifest_df.columns:
            manifest_df["export_status"] = "skipped"
        manifest_df.loc[export_df.index, "export_status"] = export_status

        # Add tiff_path column for 3D scenes
        manifest_df["tiff_path"] = manifest_df.apply(
            lambda r: (
                str(tiff_dir / f"{r['image_id']}.ome.tif")
                if r["scene_type"] == "3D_stack"
                else ""
            ),
            axis=1,
        )

    # ─── WRITE MANIFESTS ───
    manifest_df.to_csv(manifest_path, index=False)
    logger.info(f"\nManifest written: {manifest_path}")

    # Template
    template_path = out_dir / "manifest_template_metadata.csv"
    create_metadata_template(manifest_df, template_path, logger)

    # Auto-merge if the user provided metadata_csv
    if args.metadata_csv:
        metadata_path = Path(args.metadata_csv).expanduser().resolve()
        if metadata_path.exists():
            merge_metadata_into_manifest(manifest_path, metadata_path, logger)
        else:
            logger.warning(f"Metadata CSV does not exist: {metadata_path}")

    # ─── SUMMARY ───
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total scenes scanned : {len(manifest_df)}")
    logger.info(
        f"  3D stacks (for quantification): "
        f"{(manifest_df['scene_type'] == '3D_stack').sum()}"
    )
    logger.info(
        f"  Single-plane (representative)  : "
        f"{(manifest_df['scene_type'] == 'single_plane').sum()}"
    )
    logger.info(f"\nBit depth distribution:")
    logger.info(manifest_df["bit_depth"].value_counts().to_string())
    logger.info(f"\nVoxel sizes (Z, Y, X) seen:")
    voxel_summary = (
        manifest_df[["voxel_z_um", "voxel_y_um", "voxel_x_um"]]
        .round(4)
        .drop_duplicates()
    )
    logger.info(voxel_summary.to_string())

    logger.info("\n" + "=" * 70)
    logger.info("NEXT STEPS:")
    logger.info("=" * 70)
    logger.info(f"1. OPEN and FILL metadata: {template_path}")
    logger.info(f"   Columns to fill: combo, antibody, condition, case_id")
    logger.info(f"   - combo:     'Combo1' or 'Combo3'")
    logger.info(f"   - antibody:  'MPO_CitH3' or 'MPO_CitH3_aSyn'")
    logger.info(f"   - condition: 'HC'/'PD'")
    logger.info(f"   - case_id: 'R1', 'R2', 'R3', ... (hay T1069a, T127, v.v.)")
    logger.info(f"2. After filling, run:")
    logger.info(f"   python {Path(__file__).name} \\")
    logger.info(f"       --lif-dir {lif_dir} --out-dir {out_dir} \\")
    logger.info(f"       --metadata-csv {template_path} --merge-only")
    logger.info(
        f"3. -> Generates manifest_final.csv, ready for 02_Segment_Quantify_NETs.py"
    )


if __name__ == "__main__":
    main()
