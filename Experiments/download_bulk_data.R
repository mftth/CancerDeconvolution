## Mastherthesis, Melanie Fattohi
## download bulk RNA-seq datasets from different sources
## bulk RNA-seq of common cancers and their meta data (specifically tumor grading)
## pancreatic ductal adenocarcinoma (PDAC)
## invasive ductal carcinoma (IDC) (breast cancer)
## diffuse large B-cell lymphoma (DLBCL)

library(tidyverse)
library(GEOquery)
library(HelpersMG)
library(readxl)

get_max_var_genes <- function(expr_data, gene, genes){ 
  # expr_data = expression in data.frame, gene = gene of question, genes = vector of all genes in expr_data
  gene_idx <- which(genes == gene)
  var_genes <- apply(expr_data[gene_idx,], MARGIN = 1, var)
  max_var_gene <- names(which.max(var_genes)) #as.numeric(names(which.max(var_genes)))
  return(max_var_gene)
}

#################### PDAC ####################

## PACA-CA
# setwd("~/Masterthesis/Data/Bulk/PACA_CA")
# HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PACA-CA/exp_seq.PACA-CA.tsv.gz")
# gunzip("~/Masterthesis/Data/Bulk/PACA_CA/exp_seq.PACA-CA.tsv.gz")
# PACA_CA <- read.table("~/Masterthesis/Data/Bulk/PACA_CA/exp_seq.PACA-CA.tsv",
#                         header = TRUE, sep = "\t")
# 
# HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PACA-CA/simple_somatic_mutation.open.PACA-CA.tsv.gz")
# gunzip("~/Masterthesis/Data/Bulk/PACA_CA/simple_somatic_mutation.open.PACA-CA.tsv.gz")
# PACA_CA <- read.table("~/Masterthesis/Data/Bulk/PACA_CA/simple_somatic_mutation.open.PACA-CA.tsv",
#                         header = TRUE, sep = "\t")
# 
# 
# ## PACA-AU
# setwd("~/Masterthesis/Data/Bulk/PACA_AU")
# HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PACA-AU/exp_array.PACA-AU.tsv.gz")
# gunzip("~/Masterthesis/Data/Bulk/PACA_AU/exp_array.PACA-AU.tsv.gz")
# PACA_AU <- read.table("~/Masterthesis/Data/Bulk/PACA_AU/exp_array.PACA-AU.tsv",
#                         header = TRUE, sep = "\t")


## Yang, 2016, GSE62452, microarray, 69 samples
setwd("~/Masterthesis/Data/Bulk/Yang")
HelpersMG::wget("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/matrix/GSE62452_series_matrix.txt.gz")
gunzip("~/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt.gz")
#Yang_metadata <- read.table("~/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt",
#                            header = FALSE, sep = "\t")
Yang <- getGEO(filename="~/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt")
Yang_metadata <- pData(Yang)
Yang_metadata_tumor <- Yang_metadata[which(Yang_metadata$'tissue:ch1'== 'Pancreatic tumor'),]
## cols of interest: Stage:ch1, grading:ch1, survival months:ch1, survival status:ch1, tissue:ch1
Yang_metadata_tumor <- Yang_metadata_tumor[, c("Stage:ch1", "grading:ch1",
                                               "survival months:ch1", "survival status:ch1", 
                                               "tissue:ch1")]
Yang_bulk <- exprs(Yang)
Yang_bulk_tumor <- Yang_bulk[, match(rownames(Yang_metadata_tumor),colnames(Yang_bulk))]

## annotation of genes
#HelpersMG::wget("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/soft/GSE62452_family.soft.gz")
#gunzip("~/Masterthesis/Data/Bulk/Yang/GSE62452_family.soft.gz")
#Yang_anno <- getGEO(filename="~/Masterthesis/Data/Bulk/Yang/GSE62452_family.soft")

#require(hugene10sttranscriptcluster.db)
#Yang_annotLookup <- select(hugene10sttranscriptcluster.db, keys = rownames(Yang_bulk_tumor),
#  columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
Yang_annotLookUp <- getBM(filters = "affy_hugene_1_0_st_v1", 
                          attributes = c("affy_hugene_1_0_st_v1", "ensembl_gene_id", "hgnc_symbol"), 
                          values = rownames(Yang_bulk_tumor),
                          mart= mart, useCache = FALSE)
yang_intersection <- intersect(rownames(Yang_bulk_tumor), Yang_annotLookUp$affy_hugene_1_0_st_v1)
Yang_sub_annotLookUp <- Yang_annotLookUp[which(yang_intersection %in% Yang_annotLookUp$affy_hugene_1_0_st_v1),]
Yang_sub_annotLookUp$affy_hugene_1_0_st_v1 <- as.character(Yang_sub_annotLookUp$affy_hugene_1_0_st_v1)
Yang_sub_bulk_tumor <- as.data.frame(Yang_bulk_tumor[which(yang_intersection %in% rownames(Yang_bulk_tumor)),])
## rwonames von bulk verschieben zu spalte. no rownames
Yang_sub_bulk_tumor$gene <- rownames(Yang_sub_bulk_tumor)
rownames(Yang_sub_bulk_tumor) <- NULL
## dann diese spalten matchen mit subannot
Yang_sub_bulk_tumor <- Yang_sub_bulk_tumor[match(Yang_sub_annotLookUp$affy_hugene_1_0_st_v1, Yang_sub_bulk_tumor$gene),]
all(Yang_sub_bulk_tumor$gene == Yang_sub_annotLookUp$affy_hugene_1_0_st_v1)
## dann diese spalte gleichsetzen mit hgnc genen
Yang_sub_bulk_tumor$gene <- Yang_sub_annotLookUp$hgnc_symbol
## dann nach duplikaten suchen und maxvargene function anwenden
## remove non-unique genes by taking the one with highest variance
yang_max_var_genes <- sapply(unique(Yang_sub_bulk_tumor$gene), 
                             function(x) get_max_var_genes(Yang_sub_bulk_tumor[,-ncol(Yang_sub_bulk_tumor)], x, 
                                                           Yang_sub_bulk_tumor$gene))
Yang_sub_bulk_tumor <- Yang_sub_bulk_tumor[yang_max_var_genes,]
rownames(Yang_sub_bulk_tumor) <- names(yang_max_var_genes)
Yang_sub_bulk_tumor <- Yang_sub_bulk_tumor[-1,-ncol(Yang_sub_bulk_tumor)]

write.table(Yang_sub_bulk_tumor, file = "~/Masterthesis/Data/Bulk/Yang/Yang_bulk.tsv",
           sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(Yang_metadata_tumor, file = "~/Masterthesis/Data/Bulk/Yang/Yang_metadata.tsv",
           sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(Yang_sub_bulk_tumor, file = "~/Masterthesis/Data/Bulk/Yang/Yang_bulk.RDS")
saveRDS(Yang_metadata_tumor, file = "~/Masterthesis/Data/Bulk/Yang/Yang_metadata.RDS")


## PAAD
#paad <- read.table("~/Masterthesis/Data/Bulk/PAAD/illuminahiseq_rnaseqv2-RSEM_genes/PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
#                   header = TRUE, sep = "\t")
paad_raw <- read.table("~/Masterthesis/Data/Bulk/PAAD/mRNAseq_Preprocess/PAAD.uncv2.mRNAseq_raw_counts.txt",
                   header = TRUE, sep = "\t", row.names = 1)
paad_zscore <- read.table("~/Masterthesis/Data/Bulk/PAAD/mRNAseq_Preprocess/PAAD.uncv2.mRNAseq_RSEM_Z_Score.txt",
                       header = TRUE, sep = "\t", row.names = 1)
paad_meta <- read.delim("~/Masterthesis/Data/Bulk/PAAD/Clinical_Pick_Tier1/All_CDEs.txt",
                        header = TRUE, sep = "\t", row.names = 1)
paad_meta <- as.data.frame(t(paad_meta))
rownames(paad_meta) <- toupper(rownames(paad_meta))
colnames(paad_raw) <- unname(sapply(colnames(paad_raw), function(x) substr(x, 1, nchar(x)-3)))
paad_meta <- paad_meta[match(colnames(paad_raw), rownames(paad_meta)),]
paad_genes <- unname(sapply(rownames(paad_raw), function(x) strsplit(x, "\\|")[[1]][1]))
paad_raw <- paad_raw[-which(paad_genes == "?"),]
paad_genes <- paad_genes[-which(paad_genes == "?")]
paad_max_var_genes <- sapply(unique(paad_genes), function(x) get_max_var_genes(paad_raw, x, paad_genes))
paad_raw <- paad_raw[paad_max_var_genes,]
rownames(paad_raw) <- names(paad_max_var_genes)
colnames(paad_raw)[which(duplicated(colnames(paad_raw)))] <- rownames(paad_meta)[which(duplicated(colnames(paad_raw)))]

write.table(paad_raw, file = "~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(paad_meta, file = "~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(paad_raw, file = "~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.RDS")
saveRDS(paad_meta, file = "~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.RDS")
  

## Guo
guo <- read.table("~/Masterthesis/Data/Bulk/Guo/GSE172356_PDA_gene_expression_matrix.txt",
                  header = TRUE, sep = "\t", row.names = 1)
guo_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/Guo/GSE172356_series_matrix.txt") 
guo_meta <- guo_obj@phenoData@data
guo_meta_supl <- lapply(1:19, function(x) read_excel("~/Masterthesis/Data/Bulk/Guo/42003_2021_2557_MOESM7_ESM.xlsx",
                            sheet = x)) ## survival is #11
guo_survival <- guo_meta_supl[[11]]
guo_meta <- cbind(guo_meta, guo_survival[match(guo_meta$description.1, guo_survival$Sample_ID),])
guo_meta <- data.frame(guo_meta)
colnames(guo) <- rownames(guo_meta)

write.table(guo, file = "~/Masterthesis/Data/Bulk/Guo/Guo_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(guo, file = "~/Masterthesis/Data/Bulk/Guo/Guo_bulk.RDS")
write.table(guo_meta, file = "~/Masterthesis/Data/Bulk/Guo/Guo_metadata.tsv",
           sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(guo_meta, file = "~/Masterthesis/Data/Bulk/Guo/Guo_metadata.RDS")


## Kirby
kirby <- read.table("~/Masterthesis/Data/Bulk/Kirby/GSE79668_51_tumors_sharedgenecounts.txt",
                  header = TRUE, sep = "\t", row.names = 1)
kirby_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/Kirby/GSE79668_series_matrix.txt") 
kirby_meta <- kirby_obj@phenoData@data
kirby_genes <- unname(sapply(rownames(kirby), function(x) strsplit(x, split = "_")[[1]][1]))
kirby_max_var_genes <- sapply(unique(kirby_genes), function(x) get_max_var_genes(kirby, x, kirby_genes))
kirby <- kirby[kirby_max_var_genes,]
rownames(kirby) <- names(kirby_max_var_genes)
colnames(kirby) <- rownames(kirby_meta)

write.table(kirby, file = "~/Masterthesis/Data/Bulk/Kirby/Kirby_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(kirby, file = "~/Masterthesis/Data/Bulk/Kirby/Kirby_bulk.RDS")
write.table(kirby_meta, file = "~/Masterthesis/Data/Bulk/Kirby/Kirby_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(kirby_meta, file = "~/Masterthesis/Data/Bulk/Kirby/Kirby_metadata.RDS")



## Moffitt
## RNA-seq
moffit_seq <- read.table("~/Masterthesis/Data/Bulk/Moffitt/Moffitt_rna_seq.txt", 
                         sep = "\t", header = TRUE, row.names = NULL)
colnames(moffit_seq) <- colnames(moffit_seq)[-1]
moffit_seq <- moffit_seq[,-66]
moffit_seq_genes_meta <- moffit_seq[,1:4]
moffit_seq <- moffit_seq[,-c(2,3,4)]
moffit_seq_meta <- read.delim("~/Masterthesis/Data/Bulk/Moffitt/Moffitt_meta.txt", 
                              sep = "\t", header = FALSE, row.names = 1)
moffit_seq_meta <- t(moffit_seq_meta[,-62])
rownames(moffit_seq_meta) <- sub("-", "_",moffit_seq_meta[,1])
moffit_seq$Gene <- toupper(moffit_seq$Gene)
moffit_seq_max_var_genes <- sapply(unique(moffit_seq$Gene), function(x) get_max_var_genes(moffit_seq[,-1], x, moffit_seq$Gene))
moffit_seq <- moffit_seq[moffit_seq_max_var_genes,]
rownames(moffit_seq) <- names(moffit_seq_max_var_genes)
moffit_seq <- moffit_seq[,-1]
moffit_seq_meta <- as.data.frame(moffit_seq_meta)

write.table(moffit_seq, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_seq_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(moffit_seq, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_seq_bulk.RDS")
write.table(moffit_seq_meta, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_seq_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(moffit_seq_meta, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_seq_metadata.RDS")

## microarray
moffit_array_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/Moffitt/GSE71729_series_matrix.txt")  
moffit_array <- moffit_array_obj@assayData$exprs #getGEO(filename="~/Masterthesis/Data/Bulk/Moffitt/GSE71729_RAW.tar")#
moffit_array_meta <- moffit_array_obj@phenoData@data
moffit_array <- as.data.frame(moffit_array)

write.table(moffit_array, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_array_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(moffit_array, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_array_bulk.RDS")
write.table(moffit_array_meta, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_array_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(moffit_array_meta, file = "~/Masterthesis/Data/Bulk/Moffitt/Moffitt_array_metadata.RDS")


## Janky
janky_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/Janky/GSE62165_series_matrix.txt")  
janky <- janky_obj@assayData$exprs # gen namen sind in feature data (aber auch duplicates) # transcripts need to be filtered
janky_meta <-janky_obj@phenoData@data
janky_genes <- janky_obj@featureData@data
janky_genes <- janky_genes[c("ID", "Gene Symbol")]
janky_genes <- janky_genes[-grep("NA", rownames(janky_genes)),]
janky_genes$`Gene Symbol` <- unname(sapply(janky_genes$`Gene Symbol`, function(x) strsplit(x, " ///")[[1]][1]))
janky <- janky[match(janky_genes$ID, rownames(janky)),]
janky <- as.data.frame(janky)
janky_max_var_genes <- sapply(unique(janky_genes$`Gene Symbol`), function(x) get_max_var_genes(janky, x, janky_genes$`Gene Symbol`))
janky <- janky[unlist(janky_max_var_genes),]
rownames(janky) <- names(unlist(janky_max_var_genes))

write.table(janky, file = "~/Masterthesis/Data/Bulk/Janky/Janky_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(janky, file = "~/Masterthesis/Data/Bulk/Janky/Janky_bulk.RDS")
write.table(janky_meta, file = "~/Masterthesis/Data/Bulk/Janky/Janky_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(janky_meta, file = "~/Masterthesis/Data/Bulk/Janky/Janky_metadata.RDS")
