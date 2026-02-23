#README – Supplementary Data**

## Integrating evolutionary signals and protein structure reveals localized adaptive divergence in *Cereus* lineages

---

## Overview

This repository contains the supplementary datasets associated with the manuscript:

**Integrating evolutionary signals and protein structure reveals localized adaptive divergence in *Cereus* lineages**

**Authors:** Danilo T. Amaral, João Alfredo Teodoro, Maria Izadora O. Cardoso, Fernando F. Franco, Isabel A. S. Bonatelli

These files support the evolutionary, structural, and biophysical analyses presented in the main text and supplementary materials.

The dataset enables full reproducibility of the sequence variability, structural comparison, pocket detection, and stability analyses conducted in this study.

---

## Contents

### 1. Positive Selection and Reference Data

**Directory:** `01_psg_sites/`

* `BEB_sites_prevpaper.csv`
* `PSG_master_prevpaper.csv`

**Supplementary Tables:**

* Table S1 – Complete list of positively selected sites
* Table S2 – Summary per orthogroup

**Directory:** `00_inputs/selection_prevpaper/`

* `425_2024_4442_MOESM1_ESM.xlsx`

**Supplementary Table:**

* Table S4 – Reference data from previous study (Amaral et al., 2024)

---

### 2. Quality Control

**Directory:** `01_qc_tables/`

* `og_summary.tsv`

**Supplementary Table:**

* Table S3 – Orthogroup quality statistics

---

### 3. Orthogroup Datasets

**Directory:** `02_selected_ogs/`

* `main/` – Primary datasets
* `method_control/` – Control datasets

Contains curated protein sequences for selected orthogroups:

* OG0028003
* OG0028976
* OG0029756
* OG0030099
* OG0031271

---

### 4. Multiple Sequence Alignments

**Directory:** `03_alignments/raw_nr/`

Includes non-redundant multiple sequence alignments used for entropy analyses, for example:

* `OG0028003.fa.nr.fa.aln.fa`
* `OG0028976.fa.nr.fa.aln.fa`
* `OG0029756.fa.nr.fa.aln.fa`
* `OG0030099.fa.nr.fa.aln.fa`

These files were used to compute Shannon entropy and identify variability hotspots.

---

### 5. Structural Models

**Directory:** `04_structures/`

Subdirectories include:

* `alphafold_best/` – Best-confidence AlphaFold models
* Additional folders (when available)

**Supplementary File S2**
Representative AlphaFold models (one per species per orthogroup) are provided in CIF format:

```
04_structures/alphafold_best/OG*/fold_*.cif
```

---

### 6. Structural Comparison Data

Includes results from structural alignments:

* Pairwise RMSD matrices
* Pairwise TM-score matrices
* Heatmaps and visualization files

Used to generate Figure 2 and related supplementary figures.

---

### 7. Pocket Detection Results

**Directory:** `fpocket/`

Contains output files from fpocket analyses, including:

* Predicted cavities
* Pocket scores
* Pocket coordinates

Used for hotspot–pocket enrichment analyses (Table 2).

---

### 8. Stability Analysis Data

**Directory:** `foldx/`

Contains FoldX output files from mutational scanning:

* Individual mutation results
* Replicate runs
* ΔΔG summaries

Used to generate Table 3 and related analyses.

---

### 9. Supplementary Figures

Includes additional visualization files:

* Heatmaps
* Scatterplots
* ΔΔG plots
* Entropy profiles

Correspond to Figures S1–S5 in the Supplementary Material.

---

### 10. Scripts and Pipelines

Includes custom analysis scripts:

* `*.py` – Python scripts
* `*.sh` – Shell pipelines

Used for entropy calculation, data integration, visualization, and statistical analyses.

---

## Methods Summary

### Sequence Analysis

* Orthologs retrieved from OrthoFinder outputs
* Multiple sequence alignments generated with MAFFT
* Shannon entropy calculated per alignment position
* Gap-rich columns (>20%) excluded
* Top 5% entropy positions defined as hotspots

### Structural Modeling

* Structures predicted using AlphaFold
* Best-confidence models selected based on mean pLDDT
* Regions with pLDDT < 70 interpreted conservatively

### Structural Comparison

* Pairwise alignments using TM-align and Foldseek
* RMSD and TM-score extracted
* Results visualized as heatmaps

### Pocket Detection

* Cavities predicted using fpocket
* Hotspot proximity assessed within 5 Å
* Enrichment tested using Fisher’s exact tests

### Stability Analysis

* Mutations introduced using FoldX
* RepairPDB and BuildModel workflows applied
* Five replicates per mutation
* Mean ΔΔG values reported

Detailed descriptions are provided in the Methods section of the manuscript.

---

## File Organization

```
/
├── 00_inputs/
├── 01_psg_sites/
├── 01_qc_tables/
├── 02_selected_ogs/
├── 03_alignments/
├── 04_structures/
├── foldx/
├── fpocket/
├── Figures/
├── Scripts/
└── README.md
```

(Folder names may vary slightly depending on the final Figshare upload.)

---

## Usage Notes

* All data are provided for transparency and reproducibility.
* Structural matrices and tabular data can be imported into R, Python, or MATLAB.
* Alignments are compatible with standard bioinformatics tools.
* AlphaFold models can be visualized using PyMOL, ChimeraX, or UCSF Chimera.

---

## Reproducibility

All analyses were conducted using standardized and reproducible pipelines.

Random seeds, preprocessing parameters, and filtering criteria are described in the main manuscript.

Whenever possible, software versions are reported in the Methods section.

---

## Contact

For questions regarding these data, please contact:

**Danilo Trabuco do Amaral**
Centro de Ciências Naturais e Humanas, UFABC
Email: [danilo.trabuco@ufabc.edu.br](mailto:danilo.trabuco@ufabc.edu.br)
