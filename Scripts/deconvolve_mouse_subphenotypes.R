## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Dorothy (~/SeneSys_scRNA_mouse/senescence_mouse_subpheno.tsv)
## contains 5 phenotypes (ADR_ki67_high, ADR_ki67_low, ADR_OHT_ki67_low, PS_ki67_high, PS_ki67_low)
## Deconvolve Mouse data based on the scRNA-seq data

source("~/Masterthesis/CancerDeconvolution/Scripts/Quality_control.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(readxl)
library(readtext)


## 1) Perform SCDC QC on the scRNA-seq data
senescence <- read.table("~/SeneSys_scRNA_mouse/senescence_mouse_subpheno.tsv", sep = "\t", header = TRUE, row.names = 1)
senescence <- apply(senescence, MARGIN = c(1,2), FUN = as.numeric) # takes much time
rownames(senescence) <- toupper(rownames(senescence))
phenotypes <- sapply(colnames(senescence), function(x) strsplit(x, "_\\d+")[[1]][1])
phenotypes <- as.data.frame(phenotypes)
phenotypes$sample <- rep("donor", nrow(phenotypes))
colnames(phenotypes)[1] <- "cluster"

qc_senescence <- Quality_control(sc_data = senescence, sc_meta = phenotypes, 
                                 sc_path = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_senescence_mouse.RDS",
                                 cell_types = unique(phenotypes$cluster),
                                 multiple_donors = FALSE)


## 2) Import Chapuy DLBCL bulk RNA-seq data (~/SeneSys/Chapuy_Enrichment/chapuy_cohorts_expr.gct)
## Meta data in ~/SeneSys/Chapuy_Enrichment/41591_2018_16_MOESM3_ESM.xlsx

# Chapuy_88_processed <- readRDS("~/SeneSys/Survival_SeneSys/Data/Chapuy_88_processed.RDS")
# Chapuy_series_matrix <- readtext(file ="~/SeneSys/Chapuy_Enrichment/GSE98588_series_matrix.txt")
# Chapuy_series_matrix_split <- strsplit(substr(Chapuy_series_matrix$text, 1, 8000), split = "!")
# chapuy_sample_title <- Chapuy_series_matrix_split[[1]][32]
# chapuy_sample_geo_accession <- Chapuy_series_matrix_split[[1]][33]
# 
# chapuy_sample_title <- strsplit(chapuy_sample_title, split = "\t")[[1]][-1]
# chapuy_sample_title <- unname(sapply(chapuy_sample_title, function(x) strsplit(x, split = "\"")[[1]][2]))
# 
# chapuy_sample_geo_accession <- strsplit(chapuy_sample_geo_accession, split = "\t")[[1]][-1]
# chapuy_sample_geo_accession <- unname(sapply(chapuy_sample_geo_accession, function(x) strsplit(x, split = "\"")[[1]][2]))
# 
# chapuy_meta <- read_excel(path = "~/SeneSys/Chapuy_Enrichment/41591_2018_16_MOESM3_ESM.xlsx", skip = 1, col_names = TRUE)
# chapuy_meta$individual_id <- gsub("-", "_", chapuy_meta$individual_id)
# chapuy_meta$individual_id <- toupper(chapuy_meta$individual_id)
# chapuy_meta$pair_id <- gsub("-", "_", chapuy_meta$pair_id)
# chapuy_meta$pair_id <- toupper(chapuy_meta$pair_id)
# chapuy_meta$individual_id[which(chapuy_meta$pair_id %in% chapuy_sample_title)] <- chapuy_meta$pair_id[which(chapuy_meta$pair_id %in% chapuy_sample_title)]
# 
# chapuy_meta <- chapuy_meta[match(chapuy_sample_title, chapuy_meta$individual_id) ,]
# all(chapuy_meta$individual_id == chapuy_sample_title)
# unclassified <- which(chapuy_meta$COO_byGEP == "Unclassified")
# chapuy_meta <- chapuy_meta[-unclassified,]
# chapuy_sample_geo_accession <- chapuy_sample_geo_accession[-unclassified]
# chapuy_sample_title <- chapuy_sample_title[-unclassified]
# Chapuy_88_processed <- Chapuy_88_processed[, -unclassified]
# all(colnames(Chapuy_88_processed) == chapuy_sample_geo_accession)
# colnames(Chapuy_88_processed) <- chapuy_sample_title
# Chapuy_88_processed_gct <- cbind(rownames(Chapuy_88_processed), rep(NA, nrow(Chapuy_88_processed)), Chapuy_88_processed)
# colnames(Chapuy_88_processed_gct) <- c("Name", "Description", colnames(Chapuy_88_processed))
# 
# #chapuy <- read.table("~/SeneSys/Chapuy_Enrichment/chapuy_cohorts_expr.gct", header = TRUE, row.names = 1,
# #                     sep = "\t")
# chapuy <- Chapuy_88_processed_gct[,-c(1,2)]
# #chapuy_meta <- read_excel(path = "~/SeneSys/Chapuy_Enrichment/41591_2018_16_MOESM3_ESM.xlsx", 
# #                          skip = 1, col_names = TRUE)
# rownames(chapuy_meta) <- chapuy_meta$individual_id
# all(colnames(chapuy) == rownames(chapuy_meta))
# 
## 3) Import Mouse DLBCL bulk RNA-seq data
mouse <- read.table("~/SeneSys/Mouse_data/Data_9461.Counts.HGNC.tsv", header = TRUE, row.names = 1, sep = "\t")
mouse_meta <- read.table("~/SeneSys/Mouse_data/Data_9461.Phenotypes.tsv", header = TRUE, row.names = 1, sep = "\t")
# 
# ## 4) Perform Calculate_pvalue on the data
# chapuy_meta <- as.data.frame(chapuy_meta)
# chapuy <- apply(chapuy, MARGIN = c(1,2), FUN = as.numeric)
# 
# chapuy_decon <- Calculate_pvalue(nrep = 1000, ncores = 10, bulk_data = chapuy,
#                                  bulk_meta = chapuy_meta, sc_data = qc_senescence$sc.eset.qc,
#                                  cell_types = c("PS", "ADR", "ADR_OHT"),
#                                  ensemble = FALSE, multiple_donors = FALSE)


mouse_decon <- Calculate_pvalue(nrep = 1000, ncores = 10, bulk_data = mouse,
                                bulk_meta = mouse_meta, sc_data = qc_senescence$sc.eset.qc,
                                cell_types = unique(phenotypes$cluster),
                                ensemble = FALSE, multiple_donors = FALSE)
pheatmap(mouse_decon$decon_res$prop.est.mvw, 
         annotation_row = as.data.frame(mouse_meta$ResponseCode, 
                                        row.names = rownames(mouse_meta)), 
         show_rownames = FALSE)

## ADR_OHT is completely 1, always
## try normalize scRNA-seq data

#library(edgeR)
#library(limma)
#senescence_norm <- edgeR::DGEList(senescence)
#v = limma::voom(senescence_norm, design = NULL)
#senescence_norm <- v$E

senescence_norm <- (senescence - colMeans(senescence)) / apply(senescence, MARGIN = 2, FUN = sd)

qc_senescence_norm <- Quality_control(sc_data = senescence_norm, sc_meta = phenotypes, 
                                      sc_path = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_senescence_norm.RDS",
                                      cell_types = c("PS", "ADR", "ADR_OHT"),
                                      multiple_donors = FALSE)

chapuy_decon_norm <- Calculate_pvalue(nrep = 1000, ncores = 10, bulk_data = chapuy,
                                      bulk_meta = chapuy_meta, sc_data = qc_senescence_norm$sc.eset.qc,
                                      cell_types = c("PS", "ADR", "ADR_OHT"),
                                      ensemble = FALSE, multiple_donors = FALSE)