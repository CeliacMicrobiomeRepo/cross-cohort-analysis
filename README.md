# Celiac Gut Microbiome Sequencing Cross-Cohort Analysis
This repository contains scripts supporting the cross-cohort analysis of celiac gut microbiome sequencing data from the Celiac Microbiome Repository.

## About

### Publication
In this work we investigated the stool and duodenal microbiomes across the progession of celiac disease using multiple 16S rRNA and shotgun sequencing datasets from across the globe. Analysis involved differential abundance analysis, comparison of diversity metrics, machine learning prediction of celiac disease, and more. 

### The Celiac Microbiome Repository
The Celiac Microbiome Repository (CMR) is the best effort to comprehensively combine all high throughput sequencing datasets of the gut microbiome related to celiac disease. This cross-cohort analysis used data from version 1.0 of The Celiac Microbiome Repository (CMR), which was up to date as of **15th July 2025*. CMR is being extended beyond this date and is accessible on the [CMR GitHub Repo](https://github.com/CeliacMicrobiomeRepo/celiac-repository/tree/main/)

### The Celiac Microbiome Repository Web App
Our [R Shiny web application](https://celiac.shinyapps.io/celiac-webapp/) draws on the CMR's data, allowing for visualisation and exploration. The code behind this site is open source: [Webapp GitHub](https://github.com/CeliacMicrobiomeRepo/celiac-webapp/tree/main)

## Contents

### / (root)
 - `all_samples.tsv` - A TSV file containing all samples in version 1.0 of the CMR (and used in the cross-cohort analysis)
 - `low_read_samples.tsv` - A TSV file containing all samples with less than 1000 reads in the CMR
 - `samples.py` - Script to simulate and print out details of sample filtering for the cross-cohort analysis

### preprocessing/
The `preprocessing` directory contains Python and R scripts for processing of sequencing data to prepare for analysis and machine learning. These script draw from the data in the CMR. Downloading and trimming of raw reads as well as DADA2 scripts are found in the CMR's GitHub repository.
 - `truncate_asvs.R` - Script to truncate ASVs to align with no overhangs
 - `format_gtdb_database.py` - Script to format the GTDB database, required before running `combine_trunc_asvs_to_phyloseq.R`
 - `combine_trunc_asvs_to_phyloseq.R` - Script to merge and filter truncated ASVs across datasets into phyloseq objects
 - `get_read_counts.R` - Script for summarising read counts of samples in the phyloseq objects after filtering ASVs
 - `plot_asvs.py` - Script to visualise the ASVs for each analysis group after running the `combine_trunc_asvs_to_phyloseq.R` script

### analysis/
The `analysis` directory contains all Python and R scripts used to conduct the analysis.
 - `shotgun_processing/` - Directory for shotgun pre-processing scripts
 - `differential_abundance/` - Directory for differential abundance analyses
    - Scripts for differential abundance analysis (ASV, Genus, and Species seperately)
 - `diversity/` - Directory for diversity analyses
    - Scripts for alpha diversity analysis (16S and Shotgun seperately)
    - Scripts for beta diversity analysis (16S and Shotgun seperately)
    - Scripts for constrained beta diversity analysis (CAP) (16S and Shotgun seperately)
 - `phyla/` - Directory for phyla-related analyses
    - Scripts to investigate log₂(Bacillota/Bacteroidota) abundance ratios
    - Scripts to investigate phylum compositions
 - `agreement/` - Directory for agreement analyses

### machine_learning/
The `machine_learning` directory contains all Python and R scripts used to conduct the machine learning.
 - `01_phylogseq_to_tsv.R` - Script to extract data from a phyloseq object for ML
 - `02_construct_dataset.py` - Script to construct a specific ASV-level ML dataset
 - `03_kfold_train.py` - Script to perform K-fold training on an ML dataset
 - `03_lodo_xset_train.py` - Script to perform LODO and XSET training on an ML dataset
 - `custom_sklearn.py` - Script that includes a custom MLP model with dropout

## Requirements & Licenses

The code within this repository is licensed under the **GNU Affero General Public License v3.0 (AGPL-3.0)** (see the `LICENSE` file).

The scripts rely on several external tools and libraries, each with its own license:

### Python
[Python 3.12.2](https://www.python.org/downloads/release/python-3122/) was used for many scripts in this repository. Packages utilised include:
 - `pandas` 2.2.3 ([BSD 3-Clause](https://github.com/pandas-dev/pandas/blob/main/LICENSE))
 - `composition_stats` 2.0.0 ([BSD 3-Clause](https://github.com/ntessore/composition_stats/blob/main/LICENSE.txt))
 - `numpy` 1.26.4 ([BSD 3-Clause](https://github.com/numpy/numpy/blob/main/LICENSE.txt))
 - `sklearn` (scikit-learn) 1.6.1 ([BSD 3-Clause](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING))
 - `xgboost` 3.0.2 ([Apache License 2.0](https://github.com/dmlc/xgboost/blob/master/LICENSE))
 - `matplotlib` 3.10.0 ([PSF license](https://github.com/matplotlib/matplotlib/blob/main/LICENSE/LICENSE))
 - `seaborn` 0.13.2 ([BSD 3-Clause](https://github.com/mwaskom/seaborn/blob/master/LICENSE.md))

### R
[R version 4.5.1](https://www.r-project.org/) was used for many scripts in this repository. Packages utilised include:
 - `phyloseq` 1.52.0 ([AGPL-3](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html))
 - `ANCOMBC` 2.10.1 ([Artistic-2.0](https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html))
 - `dplyr` 1.1.4 ([MIT](https://dplyr.tidyverse.org/LICENSE.html))
 - `tidyr` 1.3.1 ([MIT](https://tidyr.tidyverse.org/LICENSE.html))
 - `broom` 1.0.8 ([MIT](https://broom.tidymodels.org/LICENSE.html))
 - `metafor` 4.8.0 ([GPL-2 | GPL-3](https://cran.r-project.org/package=metafor))
 - `tibble` 3.3.0 ([MIT](https://tibble.tidyverse.org/LICENSE.html))
 - `readr` 2.1.5 ([MIT](https://readr.tidyverse.org/LICENSE.html))
 - `stringr` 1.5.1 ([MIT](https://stringr.tidyverse.org/LICENSE.html))
 - `purrr` 1.1.0 ([MIT](https://purrr.tidyverse.org/LICENSE.html))
 - `dada2` 1.36.0 ([LGPL-3.0](https://github.com/benjjneb/dada2/blob/master/LICENSE))
 - `ggplot2` 3.5.2 ([MIT](https://github.com/tidyverse/ggplot2/blob/main/LICENSE.md))
 - `Biostrings` 2.76.0 ([Artistic License 2.0](https://bioconductor.org/packages/release/bioc/html/Biostrings.html))
 - `microbiome` 1.30.0 ([BSD-2-Clause + file LICENSE](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html))
 - `meta` 8.2.0 ([GPL-2 | GPL-3](https://cran.r-project.org/package=meta))
 - `vegan` 2.7.1 ([GPL-2](https://cran.r-project.org/package=vegan))
 - `philentropy` 0.9.0 ([GPL-2](https://cran.r-project.org/package=philentropy))
 - `compositions` 2.0.8 ([GPL-2 | GPL-3](https://cran.r-project.org/package=compositions))
 - `DECIPHER` 3.4.0 ([GPL-3](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html))
 - `BiocManager` 1.30.26 ([Artistic-2.0](https://cran.r-project.org/package=BiocManager))


### SRAtoolkit
[SRAtoolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) was used for downloading raw reads from NCBI SRA. (License: [Public Domain - U.S. Government Work](https://github.com/ncbi/sra-tools/blob/master/LICENSE))

### FastQC
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used for checking for the presence of adapter sequences. (License: [GPL v3.0](https://github.com/s-andrews/FastQC/blob/master/LICENSE.txt))

### Trimmomatic
[Trimmomatic](https://github.com/usadellab/Trimmomatic/releases) was used for trimming adapter sequences from raw reads. (License: [GPL v3.0](https://github.com/usadellab/Trimmomatic/blob/main/LICENSE))

### Cutadapt
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) was used for trimming adapter sequences from raw reads. (License: [MIT](https://github.com/marcelm/cutadapt/blob/main/LICENSE))

### DADA2
[DADA2](https://benjjneb.github.io/dada2/tutorial.html) 1.36.0 was used for processing raw 16S reads into ASVs and abundances. (License: [LGPL-3.0](https://github.com/benjjneb/dada2/blob/master/LICENSE))

### Mothur
[Mothur v.1.47.0](https://github.com/mothur/mothur/releases/tag/v.1.47.0) was used for aligning ASVs to a full length 16S gene reference alignment. (License: [GPL v3.0](https://github.com/mothur/mothur/blob/master/LICENSE.md))


## Authors
- **Haig Bishop**:    haig.bishop@pg.canterbury.ac.nz
- **Peter Prendergast**:    peter.prendergast@pg.canterbury.ac.nz

