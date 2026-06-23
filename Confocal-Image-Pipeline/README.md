# Confocal NETs image-analysis pipeline (Figure 5)

Quantitative 3D confocal pipeline for neutrophil extracellular traps (NETs) in post-mortem
Parkinson's disease cortex and substantia nigra, with α-synuclein (Lewy) spatial analysis.
A **NET event** is defined as voxel-level co-localisation of background-subtracted **CitH3**
and **MPO** signal overlapping **DAPI⁺** DNA. Image processing is in Python; statistics in R.

## Run order

| Step | Script | Role | Lang |
|------|--------|------|------|
| 01 | `01_Preprocess_LIF_to_TIFF.py` | Leica LIF → OME-TIFF (fixes frame-sequential Z/C order); builds `manifest.csv` | Py |
| 02 | `02_Segment_Quantify_NETs.py` | Cellpose nucleus segmentation; CitH3∩MPO∩DAPI NET segmentation; per-NET volume/solidity; per-FOV density; Manders'/Costes' co-localisation | Py |
| 03 | `03_Extract_NET_Centroids.py` | Fast NET centroid extraction (mirrors step 02, skips Costes/Cellpose) for spatial analysis | Py |
| 04 | `04_Lewy_Curate.py` | Candidate Lewy-body detection on the α-synuclein channel + human-in-the-loop confirmation (montage CSV) | Py |
| 05 | `05_Spatial_NET_Lewy.py` | NET↔Lewy nearest-neighbour distance vs Monte-Carlo null (per FOV) | Py |
| 06 | `06_NET_Statistics.R` | Primary statistics: Poisson GLMM (density) + LMM (log-volume), each with a Bayesian (brms) model; posterior probability of direction at field and donor levels | R |
| 06b | `06b_Sensitivity_Random_Effects.R` | Field-level sensitivity (no case random effect); credible intervals for plot annotation | R |

## Inputs / outputs

- **Input:** Leica `.lif` z-stacks (DAPI, CitH3, MPO; + α-synuclein/Syn1 for substantia nigra),
  placed under a local working directory; metadata in `manifest_final.csv`.
- **Output:** `results_stats/` (per-image CSVs, master NET-event tables, statistics, QC overlays).

Source images and confirmed-Lewy coordinates are not distributed here (available from the
corresponding author on reasonable request). Acquisition: Leica STELLARIS 5, sequential scan,
40× objective; substantia nigra 1× zoom (0.284 µm/px), cortex 2× zoom (0.142 µm/px).
