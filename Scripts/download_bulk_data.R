## Mastherthesis, Melanie Fattohi
## download bulk RNA-seq datasets from different sources
## bulk RNA-seq of common cancers and their meta data (specifically tumor grading)
## pancreatic ductal adenocarcinoma (PDAC)
## invasive ductal carcinoma (IDC) (breast cancer)
## adenocarcinoma of the lung (e.g. NSCLC)
## colorectal adenocarcinoma (CRA)
## head and neck squamous cell carcinoma (HNSCC)

library(tidyverse)
library(GEOquery)
library(HelpersMG)

#################### PDAC ####################

## PACA-CA
setwd("/home/fattohim/Masterthesis/Data/Bulk/PACA_CA")
HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PACA-CA/exp_seq.PACA-CA.tsv.gz")
gunzip("/home/fattohim/Masterthesis/Data/Bulk/PACA_CA/exp_seq.PACA-CA.tsv.gz")
PACA_CA <- read.table("/home/fattohim/Masterthesis/Data/Bulk/PACA_CA/exp_seq.PACA-CA.tsv",
                        header = TRUE, sep = "\t")

HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/release_27/Projects/PACA-CA/simple_somatic_mutation.open.PACA-CA.tsv.gz")
gunzip("/home/fattohim/Masterthesis/Data/Bulk/PACA_CA/simple_somatic_mutation.open.PACA-CA.tsv.gz")
PACA_CA <- read.table("/home/fattohim/Masterthesis/Data/Bulk/PACA_CA/simple_somatic_mutation.open.PACA-CA.tsv",
                        header = TRUE, sep = "\t")


## PACA-AU
setwd("/home/fattohim/Masterthesis/Data/Bulk/PACA_AU")
HelpersMG::wget("https://dcc.icgc.org/api/v1/download?fn=/current/Projects/PACA-AU/exp_array.PACA-AU.tsv.gz")
gunzip("/home/fattohim/Masterthesis/Data/Bulk/PACA_AU/exp_array.PACA-AU.tsv.gz")
PACA_AU <- read.table("/home/fattohim/Masterthesis/Data/Bulk/PACA_AU/exp_array.PACA-AU.tsv",
                        header = TRUE, sep = "\t")


## Yang, 2016, GSE62452, microarray, 69 samples
setwd("/home/fattohim/Masterthesis/Data/Bulk/Yang")
HelpersMG::wget("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/matrix/GSE62452_series_matrix.txt.gz")
gunzip("/home/fattohim/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt.gz")
#Yang_metadata <- read.table("/home/fattohim/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt",
#                            header = FALSE, sep = "\t")
Yang <- getGEO(filename="/home/fattohim/Masterthesis/Data/Bulk/Yang/GSE62452_series_matrix.txt")
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
#gunzip("/home/fattohim/Masterthesis/Data/Bulk/Yang/GSE62452_family.soft.gz")
#Yang_anno <- getGEO(filename="/home/fattohim/Masterthesis/Data/Bulk/Yang/GSE62452_family.soft")

#require(hugene10sttranscriptcluster.db)
#Yang_annotLookup <- select(hugene10sttranscriptcluster.db, keys = rownames(Yang_bulk_tumor),
#  columns = c('PROBEID', 'ENSEMBL', 'SYMBOL'))

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
Yang_annotLookUp <- getBM(filters = "affy_hugene_1_0_st_v1", 
                          attributes = c("affy_hugene_1_0_st_v1", "ensembl_gene_id", "hgnc_symbol"), 
                          values = rownames(Yang_bulk_tumor),
                          mart= mart, useCache = FALSE)
Yang_sub_annotLookUp <- Yang_annotLookUp[match(rownames(Yang_bulk_tumor), Yang_annotLookUp$affy_hugene_1_0_st_v1),]
Yang_sub_bulk_tumor <- Yang_bulk_tumor[which(!is.na(Yang_sub_annotLookUp$affy_hugene_1_0_st_v1)),]
rownames(Yang_sub_bulk_tumor) <- Yang_sub_annotLookUp$hgnc_symbol[which(!is.na(Yang_sub_annotLookUp$hgnc_symbol))]
## remove non-unique genes
Yang_sub_bulk_tumor <- Yang_sub_bulk_tumor[!duplicated(rownames(Yang_sub_bulk_tumor)),]

write.table(Yang_sub_bulk_tumor, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_bulk.tsv", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(Yang_metadata_tumor, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(Yang_sub_bulk_tumor, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_bulk.RDS")
saveRDS(Yang_metadata_tumor, file = "/home/fattohim/Masterthesis/CancerDeconvolution/Data/Bulk/Yang/Yang_metadata.RDS")
