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
### Aligning and Quantifying Paired-End RNA-seq Data
A prebuilt bowtie2 index will be needed. Edit the 
```
"bowtie_index": "/path/to/bowtie/index/"
```
line in ref.yaml to include the path to the index. The scripts_dir/make_config_3.py script can be used to make a config file from the paired end fastq files
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
Once completed, the results can be found in the rsem directory.

### Differential Expression Analysis
The code written to perform the differential expression analysis is in the R directory and should be run using the RStudio IDE. 
The following two(2) dependencies should be available in the same directory containing this "Xiao_MF_DeSeq2.R" script to avoid renaming file paths.   
1. "patient_info_final.cvs" -- a master file that holds information on sample ID (e.g. MF4_2), lesion type (e.g. TMR), disease stage (e.g. IIB)
2. A folder, named "RSEM_37" to hold all 37 .genes.results files. (e.g. "MF4_2.genes.results")

The variable `SAMPLE_GROUP` designates the comparison group and can hold the value `'TMR'`, `'LSP'` or `'ESP'`. Similarly, `REFERENCE_GROUP` can be either `'ESP'` or `'LSP'`.
The default code is set up to compare tumor (`TMR`, sample group) with early stage plaque (`ESP`, reference group), i.e.  `SAMPLE_GROUP <- 'TMR'` and `REFERENCE_GROUP <- 'ESP'`.
To explore `LSP` (sample) vs `ESP` (reference), please relabel the variables as: `SAMPLE_GROUP <- 'LSP'` vs `REFERENCE_GROUP <- 'ESP'`.
To explore `TMR` (sample) vs `LSP` (reference), please relabel the variables as: `SAMPLE_GROUP <- 'TMR'` vs `REFERENCE_GROUP <- 'LSP'`.
