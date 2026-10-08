# DNA methylation profile of plasma cell-free DNA in patients with Alzheimer's disease

**Trabajo Fin de Máster — 2026**
Máster en Bioinformática · Universidad Internacional de Valencia (VIU)

Autora: Johana Álvarez-Jiménez


---

## Overview

This repository contains the bioinformatics pipelines and analysis scripts for an
integrative study of plasma cell-free DNA (cfDNA) from Alzheimer's disease (AD)
patients and healthy controls (CT), profiled with **Enzymatic Methyl-seq (EM-seq)**.
The project characterises disease-associated epigenetic differences and infers the
tissue of origin of circulating cfDNA.

Analytical components:

| # | Component | Directory |
|---|-----------|-----------|
| 1 | Preprocessing: QC, trimming, alignment, methylation calling | [`PrePro/`](PrePro/) |
| 2 | Genome-wide differentially methylated loci and regions (DMLs / DMRs) | [`DMA/`](DMA/) |
| 3 | Genomic distribution, CpG context and functional enrichment of DMLs/DMRs | [`DMA/`](DMA/) |
| 4 | Tissue-of-origin deconvolution (reference-based and reference-free) | [`Deconvolution/`](Deconvolution/) |
| 5 | Statistical analysis of clinical covariates and group comparisons | [`statistics/`](statistics/) |

> **Data availability.** No sequencing data or clinical metadata are stored in this
> repository. Raw EM-seq data contain identifiable patient information and are
> available only under a data access agreement.

---

## Directory structure

```
DMA_ADvsCT/
├── PrePro/               # Preprocessing with nf-core/methylseq (FastQC, Bismark, Trim Galore, MultiQC)
├── DMA/                  # Differential Methylation Analysis (BSseq, DSS)
├── Deconvolution/        # Tissue-of-origin estimation (UXM)
├── pData/                # Clinical and phenotypic metadata + exploratory statistics
├── references/           # Genome reference files (GRCh38, SNPs, annotations)
├── Suplementary_data.xlsx  # Supplementary tables summarizing key results (DMLs, DMRs, cell-type composition)
└── README.md             # Main project documentation
```

- `bismark/`: Contains BAM files aligned with Bismark and methylation extraction outputs.
- `fastqc/`: Contains FastQC quality control reports for raw and trimmed reads.
- `multiqc_data/`: Raw data used to generate the MultiQC report.
- `multiqc_report.html`: Summary report aggregating metrics from all steps.
- `pipeline_info/`: Metadata about the pipeline run, including software versions and logs.
- `trimgalore/`: Output of adapter and quality trimming performed by Trim Galore.
---

## Installation & dependencies

### 1. Conda environment

Install `mamba` once (into a dedicated environment, not `base`):

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -n base -c conda-forge mamba
```

Create the project environment:

```bash
mamba create -n tfm_AD -c bioconda -c conda-forge -c defaults \
samtools picard  r-ggplot2 r-tidyverse python nextflow nf-core singularity
conda activate tfm_AD
```

The environment file pins every R/Bioconductor and Python dependency used in this
project, including `BSseq`, `DSS`, `ChIPseeker`, `annotatr`.
Installing Bioconductor packages through `install.packages()` will fail —
they are resolved from the `bioconda` channel, or inside R with `BiocManager`:
**R Libraries**:
```R :
BiocManager::install(c("bsseq", "DSS", "ChIPseeker", "annotatr", "tidyverse" ))
```
**Python Libraries**:
```bash
# Python packages
pip install pandas numpy seaborn 
```

### 2. Tools installed outside conda

| Tool | Purpose | Installation |
|------|---------|--------------|
| `wgbstools` | BAM → PAT conversion for UXM | https://github.com/nloyfer/wgbs_tools |
| `UXM` | Tissue-of-origin deconvolution | https://github.com/nloyfer/UXM_deconv |

Both are cloned and built manually; neither is available on `bioconda`. See
[`Deconvolution/README.md`](Deconvolution/README.md).

### 3. Reference genome

GRCh38 and the Bismark bisulfite index are downloaded into `references/`.
See [`references/README.md`](references/README.md) for the exact source and commands.
These files are large and are deliberately excluded from version control.

---

## Quick start

### Preprocessing — nf-core/methylseq

Raw sequencing data were preprocessed with the **nf-core/methylseq** pipeline, run in
EM-seq mode (`--em_seq true`), which adapts the workflow to Enzymatic Methyl-seq
chemistry (notably by disabling the bisulfite-specific trimming defaults). Alignment
and methylation extraction used the human reference genome **GRCh38**. The workflow was
executed with **Singularity** for containerised reproducibility, with all parameters
declared in an external YAML file rather than on the command line.

```bash
cd PrePro/
nextflow run nf-core/methylseq -r 2.7.1 \
  -profile singularity \
  -params-file params.yaml \
  -resume
```

Pinning the pipeline revision with `-r` is what makes the run reproducible — without
it, `nextflow run` pulls the latest release and the results may change.

Output subdirectories:

- `bismark/` — alignment and methylation extraction results
- `fastqc/` — per-sample quality control of raw reads
- `trimgalore/` — adapter and low-quality base trimming reports
- `multiqc/` — aggregated QC report across all samples

Technical details: [`PrePro/README.md`](PrePro/README.md).

### Differential methylation
Genome-wide differential methylation analysis was conducted on EM-seq cfDNA data using the DSS package. Two comparisons were explored:

- AD vs Control samples

After rigorous quality control (SNP exclusion, coverage filtering, the following results were obtained:
|Comparison	|Significant DMLs	|Significant DMRs|
|-----------|-----------------|----------------|
|AD vs Control|	3384	|96|

Functional enrichment (GO, KEGG) highlighted neuronal and AD-related pathways. Scripts and outputs are organized in the [`DMA/`](DMA) folder and include annotated DMLs/DMRs, PCA, and enrichment summaries.

## 🧬 Deconv

### Deconvolution

```bash
cd Deconvolution/
```
---

## Reproducibility notes

- No absolute paths appear inside the scripts. All paths and analysis parameters are
  read from `config.yaml` (see `config.example.yaml`), so the pipeline runs unchanged
  on another machine.
- Scripts are numbered by execution order within each module.
- `sessionInfo()` output for every R analysis is written to `results/tables/`.
