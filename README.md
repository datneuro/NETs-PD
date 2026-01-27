<div align="center">

# **Longitudinal systemic transcriptomic profiling and neuropathological infiltration of neutrophil extracellular traps in Parkinson’s Disease**

``[![EBioMedicine](https://img.shields.io/badge/DOI-10.1016/j.ebiom.2025.105990-DA291C?logo=elsevier)](https://doi.org/10.1016/j.ebiom.2025.105990)"
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-DA291C)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Python-3.9-DA291C?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.2-DA291C?logo=r&logoColor=white)

</div>

## **🧬 Overview**
Parkinson’s disease (PD) is a neurodegenerative disorder characterized by **α-synuclein accumulation and immune dysregulation.** This repository contains **Jupyter Notebooks (Python) and R scripts** used for analyzing **single-cell RNA sequencing (scRNA-seq), and PBMCs related data** to understand the role of **PSMB8 immunoproteasome** in synucleinopathies.

Our findings highlight:
- **Increased PSMB8 immunoproteasome expression** in **PD, MSA, and RBD**.
- **α-Synuclein pathology links immunoproteasome activation** in multiple experimental models.
- **PSMB8 inhibition reduces α-synuclein accumulation, alters immune responses, and enhances neuronal survival**.

## **📁 Repository structure** 
```
📁 Immunoproteasome_PD
 ┣ 📜 Fig.1_PBMC.R                 # R script for PBMC analysis (Demographic, RT-qPCR, Western blot data)
 ┣ 📜 Fig.3_sc-RNAseq.ipynb        # Jupyter Notebook for single-cell RNA-seq
 ┣ 📜 Fig.6_ONX0914-Treatment.R    # R script for ONX-0914 immunoproteasome inhibition
 ┣ 📜 README.md                    # Project documentation
```

## **🔍 Key analyses**
### **🧬 Single-cell RNA sequencing (scRNA-seq)**
- **Dataset:** iPSC-derived dopaminergic neurons (wild-type, rotenone-treated, SNCA-A53T mutant)
- **Preprocessing:** Cell filtering, batch correction using **scVI**, clustering via **Leiden algorithm**
- **Findings:** **PSMB8 significantly upregulated in synucleinopathy models (p < 10<sup>-300</sup>)**

### **🧪 PBMCs and NEVs analysis**
- **Peripheral blood mononuclear cells (PBMCs) from HC, RBD, PD & MSA participants**
- **Neuronal Extracellular vesicles (NEVs) from HC, PD, MSA participants**
- **RT-qPCR & Western blot:** Higher **PSMB8 immunoproteasome level in patient samples**
- **PSMB8 activity:** **Enhanced in NEVs** of patients with PD, MSA 

### **⚙️ PSMB8 Immunoproteasome inhibition with ONX-0914**
- **Treatment:** ONX-0914 (PSMB8 inhibitor) in **Human dopaminergic neurons & PBMCs**
- **Effects:**
  - **Reduced α-synuclein accumulation**
  - **Decreased HLA-I expression & cytotoxic T cells**
  - **Restored neuronal survival & decreased oxidative stress**
  - **Non-proteasomal trypsin-like activity enhances α-synuclein clearance**``

## **👥 Authors**
**Huu Dat Nguyen**<sup>1,2,3*</sup>, Seungmin Lee <sup>4</sup>,  Hyeo Il Ma<sup>1,2,3</sup>, , Yun Joong Kim<sup>5</sup>, Han-Joon Kim<sup>4</sup>, **Young Eun Kim**<sup>1,2,3,†</sup>

<sup>†</sup>Corresponding author  
<sup>*</sup>First author, Lead contact   

## **🏢 Affiliations**  
<sup>1</sup> Department of Neurology, Hallym University Sacred Heart Hospital, Hallym University, Anyang, Gyeonggi, Republic of Korea 
<sup>2</sup> Laboratory of Parkinson’s Disease and Neurodegenerative disease, Hallym Institute for Translational Medicine, Anyang, Gyeonggi, Republic of Korea
<sup>3</sup> Hallym Neurological Institute, Hallym University, Anyang, Gyeonggi, Republic of Korea  
<sup>4</sup> Department of Neurology, Seoul National University Hospital, Seoul, Republic of Korea 
<sup>5</sup> Department of Neurology, Yongin Severance Hospital, Yonsei University College of Medicine, Yongin, Gyeonggi, Republic of Korea 

## **✉️ Contacts**  
🔗 ****Professor Young Eun Kim**** – [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7182--6569-green)](https://orcid.org/0000-0002-7182-6569) | ✉️ [Email](mailto:yekneurology@hallym.or.kr)  
🔗 **Huu Dat Nguyen ENG., MMSc., PhD.** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; – [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--2491--5566-green)](https://orcid.org/0000-0003-2491-5566) | ✉️ [Email](mailto:datneurosci@gmail.com)



