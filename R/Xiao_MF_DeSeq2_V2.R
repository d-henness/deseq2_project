#####Xiao et al. Comparative Transcriptome Analysis of Early- and Advanced-Stage Mycosis Fungoides######
#####Division of Dermatology, University of Alberta, Edmonton, AB, Canada

######################### Load libraries #########################
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

# SETS THE WORKING DIRECTORY TO THE FOLDER CONTAINING THIS FILE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############# Import RSEM for DESeq #############
##### The "RSEM_37" folder holds all 37 .genes.results files.
dir <- normalizePath('RSEM_37') 
##### The "patient_info_final.csv" is a spreadsheet with information on sample ID (MF4_2), lesion type (TMR), stage (IIB)                                
coldata <- read.csv('patient_info_final.csv', stringsAsFactors = F)
coldata$x[coldata$x == "TUM"] <- "TMR"
files <- file.path(dir, paste0('MF', coldata$sample, ".genes.results"))
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

#Adding column names 
coldata$label <- paste(coldata$x, coldata$sample, sep= '')
colnames(txi.rsem$abundance) [1:37] <- as.character(coldata$label)
colnames(txi.rsem$counts) [1:37] <- as.character(coldata$label)
colnames(txi.rsem$length) [1:37] <- as.character(coldata$label)

# Fix gene length = 0 in txi.rsem
txi.rsem$length[txi.rsem$length == 0] <- 1

########PCoA plot using EdgeR ##############
edgeR_input <- readDGE(files, columns=c("gene_id","TPM"), group=coldata$x, labels=coldata$sample)
keep <- filterByExpr(edgeR_input)
edgeR_input <- edgeR_input[keep,,keep.lib.sizes=FALSE]
edgeR_input <-calcNormFactors(edgeR_input)
design <- model.matrix(~coldata$x)
edgeR_input <- estimateDisp(edgeR_input,design)
logcpm <- cpm(edgeR_input, log=TRUE)
logcpm_plot <- plotMDS(logcpm)
mds_data <- data.frame("x" = logcpm_plot$x, "y" = logcpm_plot$y, "Disease_Stage" = coldata$x, "Lesion" = coldata$tp)

mds_plot <- ggplot(mds_data, aes(x, y, color = Disease_Stage, shape = Lesion)) +                                                        
  geom_point(size=3) +                                                                                                    
  xlab("Leading logFC dim 1") +                                                                                            
  ylab("Leading logFC dim 2") +                                                                                            
  theme(axis.text.x = element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size=14), plot.title = element_text(hjust=0.5), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))

mds_plot

####################### DESeq ##########################
VARIABLE_OF_INTEREST <- 'x'

dds_formula = as.formula(paste('~', VARIABLE_OF_INTEREST))
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

### Output counts and run gage without results(dds) ####
# VST (base log2 output)
vst_counts <- assay(vsd)
rownames(vst_counts) <- gene_ids$entrez

###SAMPLE_GROUP designates the comparison group and can hold the value 'TMR', 'LSP' or 'ESP'. 
###Similarly, REFERENCE_GROUP can be either 'ESP' or 'LSP'
###In this analysis, three pairwise comparisons are possible:  
###(1) SAMPLE_GROUP (TMR) vs REFERENCE_GROUP (ESP) 
###(2) SAMPLE_GROUP (LSP) vs REFERENCE_GROUP (ESP)
###(3) SAMPLE_GROUP (TMR) vs REFERENCE_GROUP (LSP)
SAMPLE_GROUP <- 'TMR'
REFERENCE_GROUP <- 'ESP'

###SET up GAGE ###
samp <- grep(SAMPLE_GROUP, coldata[[VARIABLE_OF_INTEREST]])
ref <- grep(REFERENCE_GROUP, coldata[[VARIABLE_OF_INTEREST]])

###SET up KEGG database###
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]

###SET up Gene Ontology database###
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]

# KEGG  - WARNING: WILL TAKE CONSIDERABLE TIME TO RUN ##
vsd.kegg.sigmet.p <- gage(vst_counts, gsets = kegg.sigmet.gs, ref = ref, samp = samp, compare = 'unpaired')

# Gene Ontology - WARNING: WILL TAKE CONSIDERABLE TIME TO RUN ##
vsd.go.bp.p <- gage(vst_counts, gsets = go.bp.gs, ref = ref, samp = samp, compare = 'unpaired')

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


vsd.go.bp.esg.up <- fast.esset.grp(vsd.go.bp.p$greater, go.bp.gs, up = T)
vsd.go.bp.esg.down <- fast.esset.grp(vsd.go.bp.p$less, go.bp.gs, up = F)

####KEGG SIGMET PATHWAYS UPREGULATED IN COMPARISON VS. REFERENCE GROUP#### 
vsd.kegg.sigmet.esg.up <- fast.esset.grp(vsd.kegg.sigmet.p$greater, kegg.sigmet.gs, up = T)
View(vsd.kegg.sigmet.p$greater[which(rownames(vsd.kegg.sigmet.p$greater) %in% vsd.kegg.sigmet.esg.up$essentialSets),]) 

####KEGG SIGMET PATHWAYS DOWNREGULATED IN COMPARISON VS. REFERENCE 
vsd.kegg.sigmet.esg.down <- fast.esset.grp(vsd.kegg.sigmet.p$less, kegg.sigmet.gs, up = F)
View(vsd.kegg.sigmet.p$less[which(rownames(vsd.kegg.sigmet.p$less) %in% vsd.kegg.sigmet.esg.down$essentialSets),]) 

####GENE ONTOLOGY BIOLOGICAL PROCESS PATHWAYS UPREGULATED IN COMPARISON VS. REFERENCE GROUP#### 
vsd.go.bp.esg.up <- fast.esset.grp(vsd.go.bp.p$greater, go.bp.gs, up = T)
View(vsd.go.bp.p$greater[which(rownames(vsd.go.bp.p$greater) %in% vsd.go.bp.esg.up$essentialSets),])

####GENE ONTOLOGY BIOLOGICAL PROCESS PATHWAYS DOWNREGULATED IN COMPARISON VS. REFERENCE GROUP#### 
vsd.go.bp.esg.down <- fast.esset.grp(vsd.go.bp.p$less, go.bp.gs, up = F)
View(vsd.go.bp.p$less[which(rownames(vsd.go.bp.p$less) %in% vsd.go.bp.esg.down$essentialSets),])


#############Heatmaps of GO or KEGG stats #####################

####Example 1 : Generate a heatmap of GO upregulated annotations 3,4,5,12,13,25 

go.esg.stat.greater <- vsd.go.bp.p$stats[match(vsd.go.bp.esg.up$essentialSets[c(3,4,5,12,13,25)], rownames(vsd.go.bp.p$stats)),]
go.esg.stat.heatmap.df <- go.esg.stat.greater[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df ,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                cellwidth = 15, cellheight = 12,
                                fontsize = 12,
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)

####Example 2 : Generate a heatmap of all GO upregulated and GO downregulated annotations 
n.genes <-58
n.genes1 <-7
go.esg.stat.greater <- vsd.go.bp.p$stats[match(vsd.go.bp.esg.up$essentialSets[1:n.genes], rownames(vsd.go.bp.p$stats)),]
go.esg.stat.less<- vsd.go.bp.p$stats[match(vsd.go.bp.esg.down$essentialSets[1:n.genes1], rownames(vsd.go.bp.p$stats)),]
go.esg.stat.heatmap.df <- rbind(go.esg.stat.greater, go.esg.stat.less)[,2:ncol(go.esg.stat.greater)]
go.esg.stat.heatmap <- pheatmap(go.esg.stat.heatmap.df ,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                               gaps_row = 58,
                               cellwidth = 15, cellheight = 12,
                               fontsize = 12,
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
)



####Example 3 : Generate a heatmap of KEGG upregulated annotations 3,6,8 
vsd.kegg.stat.greater <- vsd.kegg.sigmet.p$stats[match(vsd.kegg.sigmet.esg.up$essentialSets[c(3,6,8)], rownames(vsd.kegg.sigmet.p$stats)),]

vsd.kegg.stat.heatmap.df <- vsd.kegg.stat.greater[,2:ncol(vsd.kegg.stat.greater)]
vsd.kegg.stat.heatmap <- pheatmap(vsd.kegg.stat.heatmap.df,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                breaks = seq(from=-5, to=5, length.out = 1000),
                                cellwidth = 15, cellheight = 12,
                                fontsize = 12,
                                color = colorRampPalette(c("red", "black", "green"),space = 'rgb')(1000)
) 



################Recapitulating MeiCT genes#####################
res <- as.data.frame(results(dds, contrast = c(VARIABLE_OF_INTEREST, 'TMR', 'ESP')))
rownames(res) <- gene_ids$ensembl
testicular_genes <-c('SYCP1', 'TP53','VEGFA', 'PLS3', 'BCL2', 'JUNB', 'SYCP3', 'REC8', 'PPA2', 'SMC1A', 'SMC1B', 'SMC3', 'SPO11', 'GTSF1', 'STAG3', 'SGO2', 'DMC1')
testicular_results <- res[match(gene_ids$ensembl[match(testicular_genes, gene_ids$symbol)], rownames(res)),]
rownames(testicular_results) <- testicular_genes
View(testicular_results)


################## EnhancedVolcano Plots ############################
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
    deseq_results_list[[paste(sample_level,'_',reference_level,'_res', sep = '')]] <- res_temp
    gene_symbols_for_plot <- mapIds(org.Hs.eg.db,
                                    keys=sapply(rownames(res_temp), FUN = function(x) {unlist(strsplit(x, "\\_"))[1]}),
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first")
    labels <-sub("^[^_]*_", "", rownames(res_temp))
    
    plot_temp <- EnhancedVolcano(res_temp,
                                 lab =  labels,
                                 #selectLab = c("GNLY", "IL1B","NCR1"), ##### If LSP vs. ESP 
                                 selectLab = c("PRSS21", "ZAP70", "CARD11", "TRAF2", "IL2RB", "IL2RG", "PLCG1", "NFKB2", "IKBKE", "IKBKG", "TNFRSF17", "TNFRSF18", "LTA", "TRAF2", "MMP3", "MMP9", "MMP19", "GNLY", "ELF4", "STAG3", "GTSF1", "CT45A1", "CT45A3", "PAGE5", "PNMA5", "BRDT", "PRAME", "SPAG4", "KIR3DL2", "KIR3DL4"), ##### If TMR vs. ESP  
                                 #selectLab = c("IL1B", "CX3CR1", "LRP2"),   ##### If TMR vs. LSP 
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pointSize = 3.0,
                                 labSize = 5.0,
                                pCutoff = 0.05, 
                                 FCcutoff = 1,
                                #pCutoff = 0, 
                                #FCcutoff = 0,
                                 title = paste(sample_level,'versus',reference_level))
    volcano_plot_list[[paste(sample_level,'_',reference_level,'_volcano_plot', sep = '')]] <- plot_temp
  }
}
volcano_plot_list$TMR_ESP_volcano_plot
volcano_plot_list$LSP_ESP_volcano_plot
volcano_plot_list$TMR_LSP_volcano_plot



