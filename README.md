# Asterix Project

## Overview

The Asterix project is designed to facilitate pharmacogenetic (PGx) analysis by converting individual genetic variants into medication-level advice using translation tables.

## Pharmacogenetic Data Analysis

Pharmacogenetic data analysis is performed using a series of scripts and an imputation pipeline. These are initiated using a parameter file. The workflow starts with filtering genetic data for reliably typed variants and samples that meet certain quality criteria. Next, CNV determination is performed for the CYP2D6 gene, the typed data is statistically imputed based on reference datasets, and genotypes at the variant level are translated to haplotype-level genotypes per gene (star alleles) using predefined tables. Finally, genotypes per gene are translated to a predicted phenotype per gene and per medication using predefined tables.

## Data Analysis

The data analysis uses the following GitHub repositories:

- [PGx-pipeline](https://github.com/molgenis/PGx-pipeline) for overall data analysis. 
- [Asterix](https://github.com/molgenis/asterix) for CNV determination in CYP2D6 and translating and logging genotyped variants to predicted phenotypes.
- [GAP](https://github.com/molgenis/GAP-QC) for processing data from raw .idat files to final report files and Oxford .gen and .sample files.


## Repositories

### Asterix

The Asterix repository contains a Java tool (built with Maven) that converts individual genetic variants into medication-level advice. It also includes Python software for CYP2D6 copy number (CN) calling and R scripts for visualizing the results.

#### Structure

- **src/main/java/org/molgenis/asterix**: Contains the Java source code for the Asterix tool.
- **src/main/python/cnvcaller**: Python scripts for running the CYP2D6 CNV calling
- **src/main/R/cnvcaller**: R scripts for visualizing results.

### PGx-pipeline

The PGx-pipeline repository encapsulates a series of scripts and a Nextflow pipeline for comprehensive PGx analysis. This includes variant and sample processing, variant calling, CYP2D6 CN calling (using Asterix), phasing, imputation, annotation, and running the Asterix Java tool to convert variants to medication-level advice.

#### Structure

- **data**: Contains input data files.
- **pgx-imputation-pipeline**: Nextflow pipeline for phasing and imputation.
- **protocols**: Job scripts not yet in the nextflow format
- **scripts**: Various scripts for processing and analysis.
- **templates**: Template files for configuration and reporting.
- **workflow**: Scripts for starting jobs and the nextflow pipeline.
- **parameters_template.sh**: Template parameter file for the pipeline.

## Getting Started

### Prerequisites

- Java 8 or higher
- Python (tested with 3.9)
- R (tested with 4.0.3)
- PLINK (both 1.9 and 2.0)
- Maven
- Nextflow

### Installation

1. Clone the repositories:
   ```sh
   git clone https://github.com/molgenis/asterix.git
   git clone https://github.com/molgenis/PGx-pipeline.git

## Steps

1. Create a combined samplesheet with background/reference samples and the samples of the dataset to be analyzed. Both need to be genotyped on the same genotyping chip:

2. Create a directory and link the folders with each of the SentrixBarcodes with .idat or .gtc files from a background/reference dataset and the dataset to be analyzed.

3. If .idat files are used, convert all .idat files to .gtc files.

4. Use the GAP pipeline to convert the GTC files to final report files and Oxford .gen and .sample files.

5. Fill in the parameter file with the parameters for the dataset to be analyzed: `PGx-pipeline/parameters_template.sh`

6. Start the following steps sequentially, ensuring each step is completed before starting the next. Replace `<parameters_file>` with the path to the parameter file.

   - Convert .GTC files to final report files:
     ```sh
     bash PGx-pipeline/workflow/gtc_to_final_report.submit.sh -p <parameters_file>
     ```

   **Comment:** The GAP pipeline also generates final report files, and automatically combines these into one large final report file. Since the step of processing final report files to Python pickle files is now written to process one final report file per slide, this step is included in this workflow as well.

   - Filter variants and samples:
     ```sh
     bash PGx-pipeline/workflow/quality_control.submit.sh -p <parameters_file>
     ```

   - Merge chromosomes into one binary plink dataset:
     ```sh
     bash PGx-pipeline/workflow/concatenate_chromosomes.submit.sh -p <parameters_file>
     ```

   - Process final report files to Python pickle files for CNV determination in CYP2D6:
     ```sh
     bash PGx-pipeline/workflow/stage_intensities.submit.sh -p <parameters_file>
     ```

   - Perform CNV determination in CYP2D6:
     ```sh
     bash PGx-pipeline/workflow/cnv_calling.submit.sh -p <parameters_file>
     ```

   - Integrate variants from CYP2D6 with the other data:
     ```sh
     bash PGx-pipeline/workflow/integrate_cnv_calls.submit.sh -p <parameters_file>
     ```

   - Start the imputation pipeline:
     ```sh
     bash PGx-pipeline/workflow/impute_pgx_genes.submit.sh -p <parameters_file>
     ```

   - Start the imputation pipeline for CYP2D6:
     ```sh
     bash PGx-pipeline/workflow/impute_cyp2d6.submit.sh -p <parameters_file>
     ```

   - Translate genotypes at the variant level to haplotypes at the gene level, and translate these haplotypes to predicted phenotypes:
     ```sh
     bash PGx-pipeline/workflow/asterix.submit.sh -p <parameters_file>
     ```
