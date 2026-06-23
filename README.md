<div align="center">

# **Longitudinal systemic transcriptomic profiling and neuropathological deposition of neutrophil extracellular traps in Parkinson’s Disease**

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-DA291C)](https://www.gnu.org/licenses/gpl-3.0)
![R](https://img.shields.io/badge/R-4.4-DA291C?logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.11-DA291C?logo=python&logoColor=white)

</div>

## **Authors**
**Huu Dat Nguyen**<sup>1,2,3*</sup>, Seungmin Lee<sup>4</sup>,  Hyeo Il Ma<sup>1,2,3</sup>, Yun Joong Kim<sup>5</sup>, Han-Joon Kim<sup>4</sup>, **Young Eun Kim**<sup>1,2,3,†</sup>

<sup>†</sup>Corresponding author  
<sup>*</sup>First author, Lead contact   

## **Affiliations**  
<sup>1</sup> Department of Neurology, Hallym University Sacred Heart Hospital, Hallym University, Anyang, Gyeonggi, Republic of Korea  
<sup>2</sup> Laboratory of Parkinson’s Disease and Neurodegenerative disease, Hallym Institute for Translational Medicine, Anyang, Gyeonggi, Republic of Korea  
<sup>3</sup> Hallym Neurological Institute, Hallym University, Anyang, Gyeonggi, Republic of Korea  
<sup>4</sup> Department of Neurology, Seoul National University Hospital, Seoul, Republic of Korea  
<sup>5</sup> Department of Neurology, Yongin Severance Hospital, Yonsei University College of Medicine, Yongin, Gyeonggi, Republic of Korea  

## **Contacts**  
****Professor Young Eun Kim, MD., PhD.**** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; – [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7182--6569-green)](https://orcid.org/0000-0002-7182-6569) | ✉️ [Email](mailto:yekneurology@hallym.or.kr)  
**Huu Dat Nguyen, Engr., MMSc., PhD.** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; – [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--2491--5566-green)](https://orcid.org/0000-0003-2491-5566) | ✉️ [Email](mailto:datneurosci@gmail.com)

---

## Repository structure

Analyses are organised by data domain. Statistics are in **R**; confocal image
processing and the post-mortem brain-ELISA / Western-blot quantification are in **Python**.

```
NETs-PD/
├── PPMI-Baseline/            # PPMI baseline cohort (Figure 1, Suppl. S1–S3)
│   ├── Fig1bc_S1_S3a.r                       # baseline differential expression
│   ├── PPMI-BL-visit_Volcano.ipynb           # baseline volcano
│   ├── Baseline_Neutrophil_Deconvolution.R   # CIBERSORTx + per-cell PADI4 (Suppl. S2)
│   ├── Baseline_ROC_Biomarkers.R             # diagnostic ROC panel (Suppl. S1)
│   └── Baseline_Clinical_Correlations.R      # DaT / UPDRS correlations (Suppl. S3)
├── PPMI-Longitudinal/        # PPMI longitudinal cohort (Figure 2, Suppl. S5)
│   ├── Fig2_S3.r                             # longitudinal mixed-model trajectories
│   ├── PPMI-all-visit_NETs.ipynb             # all-visit preparation
│   ├── Longitudinal_Trajectories.R           # random-slope LMM + Bayesian (Figure 2)
│   └── Longitudinal_Stability.R              # ICC / variance partition / Bayesian ICC (Suppl. S5)
├── Regional-Cohort/          # serum + post-mortem brain cohorts (Figures 3–4, Suppl. WB)
│   ├── Serum_MPO-DNA_ELISA.R                 # serum MPO-DNA (Figure 3a)
│   ├── Serum_CitH3-DNA_ELISA.R               # serum CitH3-DNA (Figure 3b)
│   ├── Serum_Biomarker_ROC.R                 # serum biomarker ROC (Figure 3)
│   ├── Brain_ELISA_3Markers.py               # brain MPO/NE/CitH3-DNA + composite (Figure 4)
│   └── Western_Blot_MPO_60kDa.py             # mature MPO ~60 kDa re-quantification (Suppl. WB)
└── Confocal-Image-Pipeline/  # post-mortem confocal NETs in cortex + substantia nigra (Figure 5)
    └── (see Confocal-Image-Pipeline/README.md for the ordered pipeline)
```

## Data

Scripts read their inputs from a local `data/` directory (or a path passed as the first
command-line argument) and write to `results/`. **Source data are not distributed here.**
PPMI transcriptomic and clinical data are controlled-access and available from the
[PPMI](https://www.ppmi-info.org/) upon application; human serum and post-mortem brain data
are available from the corresponding author on reasonable request.

## Software

- **R ≥ 4.4** — `tidyverse`, `lme4`/`lmerTest`, `emmeans`, `sandwich`, `boot`, `brms`,
  `performance`, `pROC`, `glmnet`, `limma`/`edgeR`, `WRS2`.
- **Python ≥ 3.11** — `numpy`, `pandas`, `scipy`, `statsmodels`, `scikit-image`,
  `cellpose`, `napari`, `aicsimageio`/Bio-Formats, `tifffile`.

Statistical framework (applied throughout): covariate-adjusted (mixed-effects) models with
HC3 robust standard errors, estimated marginal means with Holm adjustment, stratified/cluster
bootstrap confidence intervals, and Bayesian sensitivity analyses with weakly informative
priors (primary inference for the small post-mortem cohorts).
