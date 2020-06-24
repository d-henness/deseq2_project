#Life does not matter 

# Load libraries
library('DESeq2')
library('tximportData')
library('tximport')
library('readr')
library('ggplot2')
library('dplyr')
library('gage')
library('gageData')
library('pathview')
library('AnnotationDbi')
library('org.Hs.eg.db')
library('RColorBrewer')
library("pheatmap")
library('EnhancedVolcano')
library('clusterProfiler')
library('biomaRt')
library('edgeR')
library('tidyverse')


# Set up biomaRt
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


#generate tx2gene
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# SETS THE WORKING DIRECTORY TO THE FOLDER CONTAINING THIS FILE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("~/Maggie/melon")

# coldata
# coldata <- read.csv('patient_info.csv')
# coldata trimmed to clean up RNAseq data



#### DEPRECATED -- Import kallisto ####
# CUSTOM TX2GENE FROM GRCh38.p13 #
# tx2gene <- read.table('tx2gene.txt', stringsAsFactors = F)

# file_list <- list.files('~/Maggie/abundances')

coldata <- read.csv('patient_info_final.csv')
file_list <- sapply(coldata$sample, FUN = function(x) {paste('MF',x,'.tsv', sep = '')})
files <- file.path(getwd(), 'abundances', file_list)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

###EdgeR####
#import RSEM
edgeR_input <- readDGE(files, columns=c("gene_id","TPM"), group=coldata$x, labels=coldata$sample)
keep <- filterByExpr(edgeR_input)
edgeR_input <- edgeR_input[keep,,keep.lib.sizes=FALSE]
edgeR_input <-calcNormFactors(edgeR_input)
design <- model.matrix(~coldata$x)
edgeR_input <- estimateDisp(edgeR_input,design)
logcpm <- cpm(edgeR_input, log=TRUE)
plotMDS(logcpm)

#### Import RSEM ####
dir <- normalizePath('RSEM_37')
coldata <- read.csv('patient_info_final.csv', stringsAsFactors = F)
files <- file.path(dir, paste0('MF', coldata$sample, ".genes.results"))
# finames(files) <- paste0("sample", 1:6)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
# head(txi.rsem$counts)
#countsFromAbundance="scaledTPM" 

#added column names 
colnames(txi.rsem$abundance) [1:37] <- as.character(coldata$sample)
colnames(txi.rsem$counts) [1:37] <- as.character(coldata$sample)
colnames(txi.rsem$length) [1:37] <- as.character(coldata$sample)


# Fix gene length = 0 in txi.rsem
txi.rsem$length[txi.rsem$length == 0] <- 1

# DESeq
VARIABLE_OF_INTEREST <- 'x'
SAMPLE_GROUP <- 'LSP'
REFERENCE_GROUP <- 'ESP'

dds_formula = as.formula(paste('~', VARIABLE_OF_INTEREST))
# dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = coldata, design = ~ x)
# dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = coldata, design = dds_formula)
dds <- DESeqDataSetFromTximport(txi.rsem, colData = coldata, design = dds_formula)
dds <- DESeq(dds)

# Function to use AnnotationDbi and biomaRt for gene ids
convert_ensembl <- function(ensembl_ids) {
  to_query <- sapply(ensembl_ids, FUN = function(x) {unlist(strsplit(x, "(\\s+)|(?=[[:punct:]])", perl = TRUE))[1]})
  
  # Use biomaRt for symbols, mapIds for entrez ids and names
  gene_ids <- biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', 'entrezgene_description'),
    filters = c('ensembl_gene_id'),
    values = to_query,
    mart = ensembl,
    uniqueRows = TRUE)
  
  entrez <- unlist(AnnotationDbi::mapIds(org.Hs.eg.db, keys=gene_ids$hgnc_symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first"))
  gene_ids$annotationdbi_entrezid <- NA
  gene_ids[match(names(entrez), gene_ids$hgnc_symbol), 'annotationdbi_entrezid'] <- entrez
  
  name <- unlist(AnnotationDbi::mapIds(org.Hs.eg.db, keys=gene_ids$hgnc_symbol, column="GENENAME", keytype="SYMBOL", multiVals="first"))
  gene_ids$annotationdbi_name <- NA
  
  gene_ids[match(names(name), gene_ids$hgnc_symbol), 'annotationdbi_name'] <- name
  
  # Set empty strings to NA
  gene_ids[gene_ids == ""] <- NA
  
  # Remove duplicate ensembl ids while keeping correct entrez
  duplicated_rows <- which(duplicated(gene_ids$ensembl_gene_id) | duplicated(gene_ids$ensembl_gene_id, fromLast=TRUE))
  correct_rows <- which(gene_ids[duplicated_rows,"entrezgene_id"] == gene_ids[duplicated_rows, "annotationdbi_entrezid"])
  gene_ids <- gene_ids[-setdiff(duplicated_rows, correct_rows),]   
  gene_ids <- merge(as.data.frame(to_query), gene_ids, by.x = 'to_query', by.y = 'ensembl_gene_id', all.x = TRUE)
  
  # Remove duplicate hgnc_symbols
  symbols <- na.omit(gene_ids$hgnc_symbol)
  duplicated_symbols <- symbols[which(duplicated(symbols))]
  duplicated_symbol_rows <- gene_ids[which(gene_ids$hgnc_symbol %in% duplicated_symbols),]
  indices_to_remove <- rownames(duplicated_symbol_rows)[which(is.na(duplicated_symbol_rows$annotationdbi_entrezid))]
  gene_ids[indices_to_remove, "hgnc_symbol"] <- NA
  
  # Cleanup
  stopifnot(nrow(gene_ids) == length(to_query))
  gene_ids <- gene_ids[,c('to_query', 'hgnc_symbol', 'annotationdbi_entrezid', 'annotationdbi_name')]
  colnames(gene_ids) <- c('ensembl', 'symbol', 'entrez', 'name')
  gene_ids[gene_ids==""] <- NA
  return(gene_ids)
}
gene_ids <- convert_ensembl(rownames(counts(dds)))


# VST for downstream steps i.e. heatmap and gage; WILL TAKE LONG TIME
vsd <- vst(dds, blind = TRUE)

### Vignette plots ####
# PCA plot
pcaData <- plotPCA(vsd, intgroup=c("tp", "x"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=x, shape=tp)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle('PCA of MF dataset') +
  theme(plot.title = element_text(hjust=0.5))
coord_fixed()

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$x, vsd$sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

### Output counts and run gage without results(dds) ####
# Try a VST (base log2 output)
vst_counts <- assay(vsd)
rownames(vst_counts) <- gene_ids$entrez

# gage setup
samp <- grep(SAMPLE_GROUP, coldata[[VARIABLE_OF_INTEREST]])
ref <- grep(REFERENCE_GROUP, coldata[[VARIABLE_OF_INTEREST]])
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# KEGG

vsd.kegg.sigmet.2d.p <- gage(vst_counts, gsets = kegg.sigmet.gs, ref = ref, samp = samp, compare = 'unpaired', same.dir=F)
vsd.kegg.dise.2d.p <- gage(vst_counts, gsets = kegg.dise.gs, ref = ref, samp = samp, compare = 'unpaired', same.dir=F)

vsd.kegg.sigmet.p <- gage(vst_counts, gsets = kegg.sigmet.gs, ref = ref, samp = samp, compare = 'unpaired')
vsd.kegg.dise.p <- gage(vst_counts, gsets = kegg.dise.gs, ref = ref, samp = samp, compare = 'unpaired')

# GO - WARNING: WILL TAKE A LONG TIME
vsd.go.bp.p <- gage(vst_counts, gsets = go.bp.gs, ref = ref, samp = samp, compare = 'unpaired')
vsd.go.mf.p <- gage(vst_counts, gsets = go.mf.gs, ref = ref, samp = samp, compare = 'unpaired')
vsd.go.cc.p <- gage(vst_counts, gsets = go.cc.gs, ref = ref, samp = samp, compare = 'unpaired')

# Essential groups analyses
fast.esset.grp <- function(x, gset, up) {
  return(esset.grp(x,
                   vst_counts,
                   gsets=gset,
                   ref = ref,
                   samp = samp,
                   test4up = up,
                   output = T,
                   outname = paste(deparse(substitute(x)),'.esg', sep=''),
                   make.plot = F,
                   compare = 'unpaired'
  ))
}

vsd.kegg.sigmet.esg.up <- fast.esset.grp(vsd.kegg.sigmet.p$greater, kegg.sigmet.gs, up = T)
vsd.kegg.sigmet.esg.down <- fast.esset.grp(vsd.kegg.sigmet.p$less, kegg.sigmet.gs, up = F)
vsd.kegg.dise.esg.up <- fast.esset.grp(vsd.kegg.dise.p$greater, kegg.dise.gs, up = T)
vsd.kegg.dise.esg.down <- fast.esset.grp(vsd.kegg.dise.p$less, kegg.dise.gs, up = F)

vsd.go.bp.esg.up <- fast.esset.grp(vsd.go.bp.p$greater, go.bp.gs, up = T)
vsd.go.bp.esg.down <- fast.esset.grp(vsd.go.bp.p$less, go.bp.gs, up = F)

vsd.go.mf.esg.up <- fast.esset.grp(vsd.go.mf.p$greater, go.mf.gs, up = T)
vsd.go.mf.esg.down <- fast.esset.grp(vsd.go.mf.p$less, go.mf.gs, up = F)

vsd.go.cc.esg.up <- fast.esset.grp(vsd.go.cc.p$greater, go.cc.gs, up = T)
vsd.go.cc.esg.down <- fast.esset.grp(vsd.go.cc.p$less, go.cc.gs, up = F)



####

vsd.kegg.sigmet.esg.up <- fast.esset.grp(vsd.kegg.sigmet.p$greater, kegg.sigmet.gs, up = T)
#View(vsd.kegg.sigmet.p$greater[which(rownames(vsd.kegg.sigmet.p$greater) %in% vsd.kegg.sigmet.esg.up$essentialSets),]) 

#write.csv(vsd.kegg.sigmet.p$greater[which(rownames(vsd.kegg.sigmet.p$greater) %in% vsd.kegg.sigmet.esg.up$essentialSets),],"DESeq2_ESP_TMR_KEGG_SIGMET.csv", row.names = TRUE)


#vsd.kegg.sigmet.p$greater
#head(kegg.sigmet.gs$hsa04660) 

#View(vst_counts[which(rownames(vst_counts) %in% go.bp.gs[[rownames(vsd.go.bp.p$greater[which(rownames(vsd.go.bp.p$greater) %in% vsd.go.bp.esg.up$essentialSets),])[2]]]),])


vsd.kegg.sigmet.esg.down <- fast.esset.grp(vsd.kegg.sigmet.p$less, kegg.sigmet.gs, up = F)
#View(vsd.kegg.sigmet.p$less[which(rownames(vsd.kegg.sigmet.p$less) %in% vsd.kegg.sigmet.esg.down$essentialSets),]) 
vsd.kegg.dise.esg.up <- fast.esset.grp(vsd.kegg.dise.p$greater, kegg.dise.gs, up = T)
#View(vsd.kegg.dise.p$greater[which(rownames(vsd.kegg.dise.p$greater) %in% vsd.kegg.dise.esg.up$essentialSets),]) 
vsd.kegg.dise.esg.down <- fast.esset.grp(vsd.kegg.dise.p$less, kegg.dise.gs, up = F)
#View(vsd.kegg.dise.p$less[which(rownames(vsd.kegg.dise.p$less) %in% vsd.kegg.dise.esg.down$essentialSets),]) 


vsd.go.bp.esg.up <- fast.esset.grp(vsd.go.bp.p$greater, go.bp.gs, up = T)
#View(vsd.go.bp.p$greater[which(rownames(vsd.go.bp.p$greater) %in% vsd.go.bp.esg.up$essentialSets),])

#write.csv(vsd.go.bp.p$greater[which(rownames(vsd.go.bp.p$greater) %in% vsd.go.bp.esg.up$essentialSets),],"DESeq2_ESP_TMR_GO_BP_UP.csv", row.names = TRUE)



#View(vst_counts[which(vst_counts[,0] %in% go.bp.gs[[rownames(vsd.go.bp.p$greater)[1]]]),])

vsd.go.bp.esg.down <- fast.esset.grp(vsd.go.bp.p$less, go.bp.gs, up = F)
#View(vsd.go.bp.p$less[which(rownames(vsd.go.bp.p$less) %in% vsd.go.bp.esg.down$essentialSets),])

#write.csv(vsd.go.bp.p$less[which(rownames(vsd.go.bp.p$less) %in% vsd.go.bp.esg.down$essentialSets),],"DESeq2_ESP_TMR_GO_BP_DOWN.csv", row.names = TRUE)



vsd.go.mf.esg.up <- fast.esset.grp(vsd.go.mf.p$greater, go.mf.gs, up = T)
#View(vsd.go.mf.p$greater[which(rownames(vsd.go.mf.p$greater) %in% vsd.go.mf.esg.up$essentialSets),])
vsd.go.mf.esg.down <- fast.esset.grp(vsd.go.mf.p$less, go.mf.gs, up = F)
#View(vsd.go.mf.p$less[which(rownames(vsd.go.mf.p$less) %in% vsd.go.mf.esg.down$essentialSets),])

vsd.go.cc.esg.up <- fast.esset.grp(vsd.go.cc.p$greater, go.cc.gs, up = T)
#View(vsd.go.cc.p$greater[which(rownames(vsd.go.cc.p$greater) %in% vsd.go.cc.esg.up$essentialSets),])


vsd.go.cc.esg.down <- fast.esset.grp(vsd.go.cc.p$less, go.cc.gs, up = F)
#View(vsd.go.cc.p$less[which(rownames(vsd.go.cc.p$less) %in% vsd.go.cc.esg.down$essentialSets),])



# Heatmaps of GO stats
n.genes <- 20
go.stat.greater <- vsd.go.bp.p$stats[match(rownames(vsd.go.bp.esg.up$essentialSets)[1:n.genes], rownames(vsd.go.bp.p$stats)),]
go.stat.less <- vsd.go.bp.p$stats[match(rownames(vsd.go.bp.esg.down$essentialSets)[1:n.genes], rownames(vsd.go.bp.p$stats)),]
go.stat.heatmap.df <- rbind(go.stat.greater, go.stat.less)[,2:ncol(go.stat.greater)]
go.stat.heatmap <- pheatmap(go.stat.heatmap.df,
                            cluster_rows = FALSE,
                            cluster_cols = TRUE,
                            color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(100)
) 

n.genes <- 25

go.esg.stat.greater <- vsd.go.bp.p$stats[match(vsd.go.bp.esg.up$essentialSets[1:n.genes], rownames(vsd.go.bp.p$stats)),]
go.esg.stat.less<- vsd.go.bp.p$stats[match(vsd.go.bp.esg.down$essentialSets[1:n.genes], rownames(vsd.go.bp.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)

n.genes <- 10
go.esg.stat.greater <- vsd.go.cc.p$stats[match(vsd.go.cc.esg.up$essentialSets[1:n.genes], rownames(vsd.go.cc.p$stats)),]
go.esg.stat.less<- vsd.go.cc.p$stats[match(vsd.go.cc.esg.down$essentialSets[1:n.genes], rownames(vsd.go.cc.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df,
                                cluster_rows = FALSE,
                                cluster_cols = TRUE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)

n.genes <- 20
go.esg.stat.greater <- vsd.go.mf.p$stats[match(vsd.go.mf.esg.up$essentialSets[1:n.genes], rownames(vsd.go.mf.p$stats)),]
go.esg.stat.less<- vsd.go.mf.p$stats[match(vsd.go.mf.esg.down$essentialSets[1:n.genes], rownames(vsd.go.mf.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df,
                                cluster_rows = FALSE,
                                cluster_cols = TRUE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)

n.genes <- 15
go.esg.stat.greater <- vsd.kegg.sigmet.p$stats[match(vsd.kegg.sigmet.esg.up$essentialSets[1:n.genes], rownames(vsd.kegg.sigmet.p$stats)),]
go.esg.stat.less <- vsd.kegg.sigmet.p$stats[match(vsd.kegg.sigmet.esg.down$essentialSets[1:n.genes], rownames(vsd.kegg.sigmet.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df,
                                cluster_rows = FALSE,
                                cluster_cols = TRUE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)


n.genes <- 4
go.esg.stat.greater <- vsd.kegg.dise.p$stats[match(vsd.kegg.dise.esg.up$essentialSets[1:n.genes], rownames(vsd.kegg.dise.p$stats)),]
go.esg.stat.less <- vsd.kegg.dise.p$stats[match(vsd.kegg.dise.esg.down$essentialSets[1:n.genes], rownames(vsd.kegg.dise.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df,
                                cluster_rows = FALSE,
                                cluster_cols = TRUE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)

### Output results for individual DEGs ####
res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, REFERENCE_GROUP, SAMPLE_GROUP)))
# res <- as.data.frame(results(dds, name = 'tp_Tumor_vs_Plaque'))
res <- res[order(res$padj),]
# res <- res[which(res$padj <0.1),]
# res <- res[which(res$padj < 0.1 & abs(res$log2FoldChange) > 1.5),]
# rownames(res) <- make.names(tx2gene[match(rownames(res), tx2gene$V2),3], unique = TRUE)



###testicular genes###


res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'TUM', 'ESP')))
rownames(res) <- gene_ids$ensembl
testicular_genes <-c('SYCP1', 'PLS3', 'BCL2', 'JUNB', 'SYCP3', 'REC8', 'PPA2', 'SMC1A', 'SMC1B', 'SMC3', 'SPO11', 'GTSF1', 'STAG3', 'SGO2', 'DMC1')


pathway_of_interest <- rownames(vst_counts[which(rownames(vst_counts) %in% go.bp.gs[[rownames(vsd.go.bp.p$greater[which(rownames(vsd.go.bp.p$greater) %in% vsd.go.bp.esg.up$essentialSets),])[7]]]),])
pathway_results <- res[match(gene_ids$ensembl[match(pathway_of_interest, gene_ids$entrez)], rownames(res)),]

rownames(pathway_results) <- mapIds(org.Hs.eg.db,
       keys=sapply(rownames(pathway_results), FUN = function(x) {unlist(strsplit(x, "\\_"))[1]}),
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first")
View(pathway_results)


testicular_results <- res[match(gene_ids$ensembl[match(testicular_genes, gene_ids$symbol)], rownames(res)),]
rownames(testicular_results) <- testicular_genes
View(testicular_results)





res <- res[which(res$padj < 0.1 & abs(res$log2FoldChange) > 1.5),]



genes <- mapIds(org.Hs.eg.db,
                        keys=sapply(rownames(res), FUN = function(x) {unlist(strsplit(x, "\\_"))[1]}),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'LSP', 'ESP')))

testicular_genes <-c('SYCP1', 'CT45A3', 'PLS3', 'BCL2', 'JUNB', 'SYCP3', 'REC8', 'PPA2', 'SMC1A', 'SMC1B', 'SMC3', 'SPO11', 'GTSF1', 'STAG3', 'SGO2', 'DMC1')
testicular_results <- res[match(gene_ids$ensembl[match(testicular_genes, gene_ids$symbol)], rownames(res)),]
rownames(testicular_results) <- testicular_genes
View(testicular_results)

res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'TUM', 'ESP')))
rownames(res) <- gene_ids$ensembl
testicular_genes <-c('FOXP3', 'CTLA4', 'BCL2', 'STAT4')
testicular_results <- res[match(gene_ids$ensembl[match(testicular_genes, gene_ids$symbol)], rownames(res)),]
rownames(testicular_results) <- testicular_genes
View(testicular_results)

res <- res[which(res$padj < 0.005 & abs(res$log2FoldChange) > 1),]
rownames(res) <-sub("^[^_]*_", "", rownames(res))
View(res)


res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'LSP', 'ESP')))

res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]
rownames(res) <-sub("^[^_]*_", "", rownames(res))

write.csv(res,"ESP_LSP_short_005.csv", row.names = TRUE)

#view(res[row.names(res) %in% testicular_results, ])

###Clusterprofiler### 

res_cluster <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'TUM', 'ESP')))

###GO###
# we want the log2 fold change
original_gene_list <- res_cluster$log2FoldChange

# name the vector
names(original_gene_list) <- gene_ids$ensembl

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

##BP
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

##MF
gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")


##CC
gse <- gseGO(geneList=gene_list, 
             ont ="CC", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")


##ALL
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
emapplot(gse, showCategory = 10)

###KEGG### 

# Convert gene IDs for gseKEGG function

# Create a vector of the gene unuiverse
kegg_gene_list <- res_cluster$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- gene_ids$entrez

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

#### EnhancedVolcano ####
labels_to_plot <- rev(levels(as.factor(coldata[[VARIABLE_OF_INTEREST]])))
deseq_results_list = list()
volcano_plot_list = list()
for (i in 1:length(labels_to_plot)) {
  sample_level <- labels_to_plot[i]
  for (j in (i+1):length(labels_to_plot)) {
    reference_level <-labels_to_plot[j]
    if (is.na(reference_level)) {
      break
    }
    print(paste('Generating result and plot object for', sample_level,'versus',reference_level))
    res_temp <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, sample_level, reference_level)))
   
    res_temp <- res_temp[res_temp$log2FoldChange <10,]
    
    #res_temp  %>% filter(res_temp$log2FoldChange > 7)
    deseq_results_list[[paste(sample_level,'_',reference_level,'_res', sep = '')]] <- res_temp
     gene_symbols_for_plot <- mapIds(org.Hs.eg.db,
                                    keys=sapply(rownames(res_temp), FUN = function(x) {unlist(strsplit(x, "\\_"))[1]}),
                                     column="SYMBOL",
                                     keytype="ENSEMBL",
                                     multiVals="first")
     please <-sub("^[^_]*_", "", rownames(res_temp))
     
    plot_temp <- EnhancedVolcano(res_temp,
                                 lab =  please,
                                 selectLab = c("PRSS21", "KIR3DL2"),  
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pCutoff = 0.05, 
                                 FCcutoff = 1,
                                 title = paste(sample_level,'versus',reference_level))
    volcano_plot_list[[paste(sample_level,'_',reference_level,'_volcano_plot', sep = '')]] <- plot_temp
  }
}
#volcano_plot_list$TUM_ESP_volcano_plot
#volcano_plot_list$LSP_ESP_volcano_plot  #885, 775 
#volcano_plot_list$TUM_LSP_volcano_plot


### DEPRECATED - GO and KEGG with a different vignette https://genviz.org/module-04-expression/0004/03/01/pathwayAnalysis/ ####
# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez

fc.kegg.sigmet.p <- gage(foldchanges, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(foldchanges, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs)
fc.go.mf.p <- gage(foldchanges, gsets = go.mf.gs)
fc.go.cc.p <- gage(foldchanges, gsets = go.cc.gs)

# covert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

### DEPRECATED - GO and KEGG ####
indices <- which(res$padj < 0.1 & abs(res$log2FoldChange) > 1.5)
foldchanges <- res$log2FoldChange[indices]
names(foldchanges) <- mapIds(org.Hs.eg.db,
                             keys=rownames(res)[indices],
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
keggdfs <- lapply(keggres[names(keggres)[!names(keggres) %in% 'stats']], FUN = function(x) {as.data.frame(x)[order(as.data.frame(x)['q.val']),]})

data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres <- gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
gobpdfs <- lapply(gobpres[names(gobpres)[!names(gobpres) %in% 'stats']], FUN = function(x) {as.data.frame(x)[order(as.data.frame(x)['q.val'], decreasing = FALSE),]})

