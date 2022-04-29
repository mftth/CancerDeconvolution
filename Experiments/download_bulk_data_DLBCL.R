## Mastherthesis, Melanie Fattohi
## download bulk RNA-seq datasets from different sources
## bulk RNA-seq of common cancers and their meta data (specifically tumor grading)
## diffuse large B-cell lymphoma (DLBCL)

library(tidyverse)
library(GEOquery)
library(HelpersMG)
library(readxl)
library(biomaRt)

get_max_var_genes <- function(expr_data, gene, genes){ 
  # expr_data = expression in data.frame, gene = gene of question, genes = vector of all genes in expr_data
  gene_idx <- which(genes == gene)
  var_genes <- apply(expr_data[gene_idx,], MARGIN = 1, var)
  max_var_gene <- names(which.max(var_genes)) #as.numeric(names(which.max(var_genes)))
  return(max_var_gene)
}

mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))


## Chapuy, GSE98588
## suppl table 1 + 2 + 8
chapuy_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/DLBCL/Chapuy/GSE98588_series_matrix.txt") 
chapuy_m <- chapuy_obj@phenoData@data
chapuy_m$title2 <- gsub("_NULLPAIR", "", chapuy_m$title)
chapuy_bulk <- as.data.frame(chapuy_obj@assayData$exprs)
chapuy_bulk <- chapuy_bulk[,match(chapuy_m$geo_accession, colnames(chapuy_bulk))]
chapuy_bulk_genes <- rownames(chapuy_bulk)
rownames(chapuy_bulk) <- gsub("_at", "", rownames(chapuy_bulk))
#mapping_hgu113plus2_ensemble <- read.table("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/GPL23432_HGU133Plus2_Hs_ENSG_mapping.txt",
#                                           header = TRUE, sep = "\t")
chapuy_annotLookUp <- getBM(filters = "ensembl_gene_id",
                            attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            values = rownames(chapuy_bulk),
                            mart = mart)
#chapuy_gene_idx <- match(chapuy_annotLookUp$ensembl_gene_id, rownames(chapuy_bulk))
#chapuy_bulk <- chapuy_bulk[chapuy_gene_idx,]
#all(rownames(chapuy_bulk) == chapuy_annotLookUp$ensembl_gene_id)
#rownames(chapuy_bulk) <- chapuy_annotLookUp$hgnc_symbol
#all(colnames(chapuy_bulk) == chapuy_m$geo_accession)
chapuy_genes <- intersect(rownames(chapuy_bulk), chapuy_annotLookUp$ensembl_gene_id)
chapuy_gene_idx <- match(chapuy_genes, rownames(chapuy_bulk))
chapuy_bulk <- chapuy_bulk[chapuy_gene_idx,]
chapuy_gene_idx <- match(chapuy_genes, chapuy_annotLookUp$ensembl_gene_id)
chapuy_annotLookUp <- chapuy_annotLookUp[chapuy_gene_idx,]
all(rownames(chapuy_bulk) == chapuy_annotLookUp$ensembl_gene_id)
chapuy_max_var_genes <- mclapply(unique(chapuy_annotLookUp$hgnc_symbol), 
                                   function(x) get_max_var_genes(chapuy_bulk, gene = x, genes = chapuy_annotLookUp$hgnc_symbol),
                                   mc.cores = 5)
chapuy_max_var_genes <- Reduce(c, chapuy_max_var_genes)
chapuy_annotLookUp <- chapuy_annotLookUp[match(chapuy_max_var_genes, chapuy_annotLookUp$ensembl_gene_id),]
chapuy_bulk <- chapuy_bulk[match(chapuy_max_var_genes, rownames(chapuy_bulk)),]
all(rownames(chapuy_bulk) == chapuy_annotLookUp$ensembl_gene_id)
rownames(chapuy_bulk) <- chapuy_annotLookUp$hgnc_symbol
chapuy_suppl1 <- read_excel("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/41591_2018_16_MOESM3_ESM.xlsx", sheet = 1, skip = 1)
chapuy_suppl1$individual_id <- gsub("-", "_", chapuy_suppl1$individual_id)
chapuy_suppl1$individual_id <- toupper(chapuy_suppl1$individual_id)
chapuy_suppl1 <- chapuy_suppl1[match(chapuy_m$title2, chapuy_suppl1$individual_id),]
all(chapuy_suppl1$individual_id == chapuy_m$title2)
chapuy_suppl1$geo_accession <- chapuy_m$geo_accession
chapuy_suppl2 <- read_excel("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/41591_2018_16_MOESM4_ESM.xlsx", sheet = 1, skip = 1)
chapuy_suppl2$individual_id <- gsub("-", "_", chapuy_suppl2$individual_id)
chapuy_suppl2$individual_id <- toupper(chapuy_suppl2$individual_id)
chapuy_suppl2 <- chapuy_suppl2[match(chapuy_m$title2, chapuy_suppl2$individual_id),]
all(chapuy_suppl2$individual_id == chapuy_m$title2)
chapuy_suppl2$geo_accession <- chapuy_m$geo_accession
chapuy_suppl8 <- read_excel("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/41591_2018_16_MOESM10_ESM.xlsx", sheet = 2, skip = 1)
chapuy_suppl8 <- chapuy_suppl8[,c("individual_id", "Cluster...6")]
colnames(chapuy_suppl8) <- c("individual_id", "Cluster")
chapuy_suppl8$individual_id <- gsub("-", "_", chapuy_suppl8$individual_id)
chapuy_suppl8$individual_id <- toupper(chapuy_suppl8$individual_id)
chapuy_suppl8 <- chapuy_suppl8[match(chapuy_m$title2, chapuy_suppl8$individual_id),]
all(chapuy_suppl8$individual_id == chapuy_m$title2)
chapuy_suppl8$geo_accession <- chapuy_m$geo_accession
chapuy_meta <- cbind(chapuy_m[,c("title", "title2", "geo_accession", "tissue type:ch1")],
                     chapuy_suppl1, chapuy_suppl2, chapuy_suppl8[,"Cluster"])
rownames(chapuy_meta) <- chapuy_meta$geo_accession ## COO, OS, OS_stat, Cluster
chapuy_meta <- chapuy_meta[,-duplicated(colnames(chapuy_meta))]
all(colnames(chapuy_bulk) == rownames(chapuy_meta))
write.table(chapuy_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(chapuy_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_bulk.RDS")
write.table(chapuy_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(chapuy_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_metadata.RDS")


## Schleich, GSE134751
schleich_obj <- getGEO(filename="~/Masterthesis/Data/Bulk/DLBCL/Schleich/GSE134751_series_matrix.txt") 
schleich_bulk <- schleich_obj@assayData$exprs  ## gene names need to be converted
schleich_meta <- schleich_obj@phenoData@data
schleich_meta$treatment_response <- rep(NA, nrow(schleich_meta))
schleich_meta$treatment_response[schleich_meta$`response:ch1` == "Never Relapsed"] <- "NR"
schleich_meta$treatment_response[schleich_meta$`response:ch1` == "Relapse Prone"] <- "RP"
schleich_meta$treatment_response[schleich_meta$`response:ch1` == "Resistant"] <- "RES"
schleich_meta$treatment <- rep(NA, nrow(schleich_meta))
schleich_meta$treatment[schleich_meta$`treatment:ch1` == "native"] <- "native"
schleich_meta$treatment[schleich_meta$`treatment:ch1` == "CTX 4h"] <- "CTX"
schleich_bulk <- schleich_bulk[,match(rownames(schleich_meta), colnames(schleich_bulk))]
mart2 <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "genes", mirror = "useast")
schleich_annotLookUp <- getBM(filters = "affy_mogene_1_0_st_v1",
                            attributes = c("affy_mogene_1_0_st_v1", "hgnc_symbol", "hgnc_id", "mgi_symbol", "ensembl_gene_id"),
                            values = rownames(schleich_bulk),
                            mart = mart2)
schleich_genes <- intersect(rownames(schleich_bulk), schleich_annotLookUp$affy_mogene_1_0_st_v1)
schleich_gene_idx <- match(schleich_genes, rownames(schleich_bulk))
schleich_bulk <- schleich_bulk[schleich_gene_idx,]
schleich_gene_idx <- match(schleich_genes, schleich_annotLookUp$affy_mogene_1_0_st_v1)
schleich_annotLookUp <- schleich_annotLookUp[schleich_gene_idx,]
all(rownames(schleich_bulk) == schleich_annotLookUp$affy_mogene_1_0_st_v1)
schleich_bulk <- as.data.frame(schleich_bulk)
schleich_max_var_genes <- mclapply(unique(schleich_annotLookUp$mgi_symbol), 
                                 function(x) get_max_var_genes(schleich_bulk, gene = x, genes = schleich_annotLookUp$mgi_symbol),
                                 mc.cores = 5)
schleich_max_var_genes <- Reduce(c, schleich_max_var_genes)
schleich_annotLookUp <- schleich_annotLookUp[match(schleich_max_var_genes, schleich_annotLookUp$affy_mogene_1_0_st_v1),]
schleich_bulk <- schleich_bulk[match(schleich_max_var_genes, rownames(schleich_bulk)),]
all(rownames(schleich_bulk) == schleich_annotLookUp$affy_mogene_1_0_st_v1)
rownames(schleich_bulk) <- schleich_annotLookUp$mgi_symbol
schleich_bulk <- schleich_bulk[-7,]
all(rownames(schleich_meta) == colnames(schleich_bulk))
write.table(schleich_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(schleich_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_bulk.RDS")
write.table(schleich_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(schleich_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_metadata.RDS")


## Schmitz
schmitz_m <- read.csv("Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/coldata_Schmitz_DLBCL.csv", header = TRUE)
schmitz_suppl <- read_excel("~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/nejmoa1801445_appendix_2.xlsx", sheet = 9)
schmitz_suppl <- schmitz_suppl[match(schmitz_m$ID, schmitz_suppl$`dbGaP subject ID`),]
all(schmitz_m$ID == schmitz_suppl$`dbGaP subject ID`)
schmitz_meta <- cbind(schmitz_m, schmitz_suppl)
rownames(schmitz_meta) <- schmitz_meta$sample
schmitz_meta <- schmitz_meta[,-duplicated(colnames(schmitz_meta))]
schmitz_b <- read.table("Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_562.DLBCL_counts_norm_annotated.csv", sep = ";", header = TRUE)
schmitz_genes <- schmitz_b$gene_symbol
schmitz_bulk <- schmitz_b[,-c(1,2,3)]
schmitz_bulk <- schmitz_bulk[,match(rownames(schmitz_meta), colnames(schmitz_bulk))]
schmitz_b2 <- schmitz_bulk
schmitz_b2 <- sapply(1:ncol(schmitz_b2), function(x) as.numeric(gsub(",", "\\.", schmitz_b2[,x])))
schmitz_b2 <- data.frame(schmitz_b2, row.names = rownames(schmitz_bulk))
colnames(schmitz_b2) <- colnames(schmitz_bulk)
schmitz_bulk <- schmitz_b2
schmitz_max_var_genes <- mclapply(unique(schmitz_genes), 
                                   function(x) get_max_var_genes(schmitz_bulk, gene = x, genes = schmitz_genes),
                                   mc.cores = 10)
schmitz_max_var_genes <- Reduce(c, schmitz_max_var_genes)
schmitz_max_var_genes <- as.numeric(schmitz_max_var_genes)
schmitz_genes <- schmitz_genes[schmitz_max_var_genes]
schmitz_bulk <- schmitz_bulk[schmitz_max_var_genes,]
rownames(schmitz_bulk) <- schmitz_genes
all(rownames(schmitz_meta) == colnames(schmitz_bulk))
write.table(schmitz_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_bulk.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(schmitz_bulk, file = "~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_bulk.RDS")
write.table(schmitz_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_metadata.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
saveRDS(schmitz_meta, file = "~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_metadata.RDS")

save.image("~/Masterthesis/Workspaces/download_bulk_data_DLBCL.RData")
