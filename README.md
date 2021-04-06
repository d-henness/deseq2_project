# deseq2_project
## Prerequisites
* Python 3
* Miniconda
* Snakemake
### R packages
* DESeq2
* tximportData
* tximport
* readr
* ggplot2
* dplyr
* gage
* gageData
* pathview
* AnnotationDbi
* org.Hs.eg.db
* RColorBrewer
* pheatmap
* EnhancedVolcano
* clusterProfiler
* biomaRt
* edgeR
* tidyverse

## Steps
### Aligning RNA-seq data
A prebuilt bowtie2 index will be needed. Edit ref.yaml to include the path to the index. The scripts_dir/make_config_3.py script can be used to make a configfile from the paired end fastq files
you want to run on for snakemake.  The paired end fastq files for different libraries will need to be in separate directories. For example MF1_1.fq and MF1_2.fq should be in a directory called
MF1, MF2_1.fq and MF2_2.fq should be in a directory called MF2, RG3_1.fq and RG3_2.fq should be in a directory called RG3, etc. It takes one or more arguments which are a strings that match
something that all libraries we want to work on contain. In the previous example the command would be
```
python3 /path/to/scripts_dir/make_config_3.py RG MF > config.yaml
```

Once the config.yaml file has been generated, the alignment can be run with
```
snakemake -s /path/to/rsem.snakefile --use-conda --configfile config.yaml
```
