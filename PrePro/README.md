# Preprocessing and Methylation Extraction

This step involves the preprocessing and methylation extraction of cfDNA samples using the [nf-core/methylseq](https://github.com/nf-core/methylseq) pipeline. The alignment tool used is **Bismark**, and the data corresponds to **Enzymatic Methyl-seq (EM-seq)** libraries.

## ⚙️ Pipeline Characteristics

The configuration for this step is defined in the [params.yaml](./params.yaml) file. The main characteristics are:

- **Alignment Tool**: Bismark
- **Reference Genome**: GRCh38 (downloaded from the official source)
- **EM-seq Flag**: `em_seq = true`  
  This is recommended by the Bismark User Guide for EM-seq data and activates the following adapter trimming parameters:

  ```
  --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10
  ```

## 🛠️ Environment Setup and `nf-core/methylseq` Execution

To execute the pipeline, a computational environment was set up using the **mamba** package manager, due to its performance advantages over conda. The following steps were followed:

```bash
# Add conda-forge as a priority channel
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install mamba
conda install mamba

# Create a dedicated environment for the analysis
mamba create -c conda-forge -c bioconda -c defaults --name nf25 \\
    python nextflow nf-core singularity python-keycloak

# Activate the environment
conda activate nf25
```

Run the test profile to verify installation:

```bash
nextflow run nf-core/methylseq -profile test_full,singularity --outdir testing -resume
```

Then, execute the actual pipeline with the sample sheet and configuration:

```bash
nextflow run nf-core/methylseq -profile singularity \
  -params-file params.yaml \
  -resume
```

## 🔁 Workflow Overview
Below is the official workflow diagram for nf-core/methylseq:

![nf-core methylseq workflow](https://github.com/nf-core/methylseq/blob/3.0.0/docs/images/3.0.0_metromap.png?raw=true)

## 📋 Sample Sheet Format

The input samples are listed in a CSV file [samplesheet.csv](./samplesheet.csv) with the following format:

```
sample,fastq_1,fastq_2,genome
BACW_44,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-44_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-44_R2.fastq.gz,
BACW_46,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-46_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-46_R2.fastq.gz,
BACW_49,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-49_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-49_R2.fastq.gz,
BACW_51,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-51_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-51_R2.fastq.gz,
BACW_52,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-52_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-52_R2.fastq.gz,
BACW_54,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-54_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-54_R2.fastq.gz,
BACW_55,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-55_R1.fastq.gz,/home/projects/NAVA1/Data/241122_A00902_B_L1-2_BACW-55_R2.fastq.gz,
```
## 📂 Output Structure

After the pipeline execution, the following directory structure was generated within the output directory:

```
PrePro/
├── bismark
│   ├── alignments
│   ├── deduplicated
│   ├── methylation_calls
│   ├── reports
│   └── summary
├── fastqc
│   ├── *_1_fastqc.html
│   ├── *_2_fastqc.html
│   └── zips
├── multiqc
│   └── bismark
└── trimgalore
    ├── fastqc
    └── logs
```

- `bismark/`: Contains BAM files aligned with Bismark and methylation extraction outputs.
- `fastqc/`: Contains FastQC quality control reports for raw and trimmed reads.
- `multiqc_data/`: Raw data used to generate the MultiQC report.
- `multiqc_report.html`: Summary report aggregating metrics from all steps.
- `pipeline_info/`: Metadata about the pipeline run, including software versions and logs.
- `trimgalore/`: Output of adapter and quality trimming performed by Trim Galore.
