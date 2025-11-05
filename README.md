# Somatic Variant Calling and Annotation in High-Grade Serous Ovarian Cancer (HGSOC)

## Overview

This project performs **somatic variant calling and functional annotation** on paired tumor–normal **Whole Exome Sequencing (WES)** samples from
**High-Grade Serous Ovarian Cancer (HGSOC)**, one of the most lethal cancers in women.
Data was obtained from **NCBI BioProject [PRJNA1068328](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1068328)**.

The study investigates molecular characteristics of HGSOC, focusing on **BRCA1/2 mutation status**, **homologous recombination repair (HRR)** activity, and **epithelial-mesenchymal transition (EMT)** phenotypes.

---

## Sample Information

| Sample ID       | Type   | Description                                      |
| --------------- | ------ | ------------------------------------------------ |
| **SRR27718988** | Tumor  | Chemotherapy-naive primary ovarian cancer tissue |
| **SRR27718966** | Normal | Matched blood sample from the same patient       |

**Reference genome:** GRCh38.p14
**Study type:** Whole Exome Sequencing (WES)
**Variant type analyzed:** Somatic (Tumor vs Normal)

---

## Pipeline Overview

The full analysis was performed on a **Conda environment** using standard bioinformatics tools.
Below is the summarized pipeline:

1. **Quality Control** → `FastQC`, `MultiQC`
2. **Reference Preparation** → `samtools`, `bwa`
3. **Alignment** → `bwa mem`
4. **Mark Duplicates & Index** → `gatk MarkDuplicates`
5. **Somatic Variant Calling** → `gatk Mutect2`
6. **Filtering Variants** → `gatk FilterMutectCalls`
7. **Annotation** → `Ensembl VEP`
8. **Visualization** → `Python (matplotlib, pandas, seaborn)`

---

## Reproducibility

### Environment setup

To reproduce the exact environment:

```bash
conda env create -f environment.yml
conda activate somatic_variant_calling
```

### Running the entire pipeline

A single bash script automates all steps from QC to annotation:

```bash
bash somatic_variant_calling_clean.sh
```

This script:

* Runs FastQC on all FASTQ files
* Indexes the reference genome
* Aligns reads using BWA and sorts with samtools
* Performs variant calling with GATK Mutect2
* Filters low-confidence calls
* Annotates the variants with Ensembl VEP
* Outputs annotated `.vcf` and `.vep.txt` results for downstream visualization

---

## Tools Used

| Tool                             | Purpose                      |
| -------------------------------- | ---------------------------- |
| **FastQC / MultiQC**             | Quality control of raw reads |
| **BWA-MEM**                      | Alignment of reads to GRCh38 |
| **Samtools**                     | File indexing & BAM handling |
| **GATK4 (Mutect2)**              | Somatic variant calling      |
| **Ensembl VEP**                  | Variant effect annotation    |
| **Python (matplotlib, seaborn)** | Visualization of results     |

---

## Key Observations

* The **protein-coding biotype** harbored the highest number of variants.
* Most variants showed **MODIFIER or MODERATE** impact levels.
* **Chromosomes 1, 2, 19, and 11** showed higher mutation densities.
* **Transition mutations (A→C, T→G, C→A)** were predominant.
* Key high-impact genes detected include **DGLUCY**, **KANK1**, **RANGAP1**, and **TP53I11**, indicating potential driver events.

---

## Repository Structure

```
ngs/
├── results/
│   ├── fastqc_reports/
│   ├── plots/
│   └── variant_calling_results/
├── scripts/
│   ├── plot_variants1.py
│   └── somatic_variantcalling.sh
├── .gitignore
├── environment.yml
└── README.md

```

---

## Biological Context

The original study revealed two distinct molecular subtypes of HGSOC:

* **HRR-activated type:** Enriched with BRCA1/2 mutations and high genomic instability
* **Mesenchymal type:** Exhibits EMT activation, low genomic alterations, and poor prognosis

Our somatic variant calling highlights several genes and variant types consistent with **genomic instability and HRR pathway activation**, supporting the findings of the published study.

---

## Reproducibility Checklist

* `environment.yml` – reproducible Conda setup
* `somatic_variant_calling_clean.sh` – full automation
* `plots/` – all visual outputs
* `results/` – filtered & annotated VCFs

---

## References

* Kim et al., *Integrative multi-omics analysis reveals HRR and EMT-based subtypes of high-grade serous ovarian cancer*, BioProject: [PRJNA1068328](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1068328)
* Broad Institute GATK Best Practices for Somatic Short Variant Discovery
* Ensembl Variant Effect Predictor (VEP)

---

## Author

**Karthik R**
MSc Bioinformatics – Garden City University
Focused on drug and vaccine design, genome annotation, and variant analysis.
