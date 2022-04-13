## Mastherthesis, Melanie Fattohi
## ablation study of Tosti scRNA-seq dataset
## bulk RNA-seq dataset: Guo, PAAD

source("~/Masterthesis/CancerDeconvolution/Scripts/ablation_study.R")

## bulk RNA-seq datasets
## 183 samples; survival; tumor grading; RNA-seq
PAAD_bulk <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.RDS")
PAAD_meta <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.RDS")
hayashi_PAAD_meta <- read.table("~/Masterthesis/Data/Bulk/PAAD/hayashi_suppl.txt", 
                                sep = "\t", header = TRUE)
hayashi_PAAD_meta$TCGA_ID <- hayashi_PAAD_meta$Tumor.Sample.ID
hayashi_PAAD_meta$TCGA_ID <- gsub("-", ".", hayashi_PAAD_meta$TCGA_ID)
hayashi_PAAD_meta$TCGA_ID <- sapply(hayashi_PAAD_meta$TCGA_ID, function(x) substr(x, 1, nchar(x)-4))
rownames(hayashi_PAAD_meta) <- hayashi_PAAD_meta$TCGA_ID
hayashi_idx <- match(rownames(PAAD_meta), rownames(hayashi_PAAD_meta))
hayashi_PAAD_meta$tumor_moffitt <- hayashi_PAAD_meta$mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical
hayashi_PAAD_meta$tumor_moffitt[hayashi_PAAD_meta$tumor_moffitt == 1] <- "Basal"
hayashi_PAAD_meta$tumor_moffitt[hayashi_PAAD_meta$tumor_moffitt == 2] <- "Classical"
hayashi_PAAD_meta$tumor_collisson <- hayashi_PAAD_meta$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 1] <- "Classical"
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 2] <- "Exocrine"
hayashi_PAAD_meta$tumor_collisson[hayashi_PAAD_meta$tumor_collisson == 3] <- "QM"
hayashi_PAAD_meta$tumor_bailey <- hayashi_PAAD_meta$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 1] <- "Squamous"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 2] <- "Immunogenic"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 3] <- "Progenitor"
hayashi_PAAD_meta$tumor_bailey[hayashi_PAAD_meta$tumor_bailey == 4] <- "ADEX"
PAAD_meta$tumor_moffitt <- hayashi_PAAD_meta$tumor_moffitt[hayashi_idx]
PAAD_meta$tumor_collisson <- hayashi_PAAD_meta$tumor_collisson[hayashi_idx]
PAAD_meta$tumor_bailey <- hayashi_PAAD_meta$tumor_bailey[hayashi_idx]
PAAD_clinical_char <- data.frame("moffitt" = PAAD_meta$tumor_moffitt,
                                 "collisson" = PAAD_meta$tumor_collisson,
                                 "bailey" = PAAD_meta$tumor_bailey,
                                 "grading" = PAAD_meta$neoplasm_histologic_grade,
                                 row.names = rownames(PAAD_meta))
## 62 samples; survival; tumor subtype (basal, classical, hybrid); RNA-seq
Guo_bulk <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_bulk.RDS")
Guo_meta <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_metadata.RDS")
Guo_clinical_char <- data.frame("tumor_subtype" = Guo_meta$description,
                                row.names = rownames(Guo_meta))
## single-cell RNA-seq datasets
qc_tosti_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_tosti.RDS")

## ablation study
ct_set <- c("alpha", "beta", "gamma", "delta")
sub_ct_set <- c("sacinar", "racinar", "iacinar", "ductal", "mductal")

tosti_guo_ablation <- ablation_study(ct_combi_list = NULL, ct_set = ct_set, sub_ct_set = sub_ct_set, 
                                     res_path = "~/Masterthesis/CancerDeconvolution/Results/tosti_ablation_study_guo", 
                                     clinical_char = Guo_clinical_char, 
                                     bulk_data = Guo_bulk, bulk_meta = Guo_meta, 
                                     sc_data = qc_tosti_sc$sc.eset.qc, ensemble = FALSE, 
                                     multiple_donors = TRUE, nrep = 1000, ncores = 15)

tosti_paad_ablation <- ablation_study(ct_combi_list = NULL, ct_set = ct_set, sub_ct_set = sub_ct_set, 
                                      res_path = "~/Masterthesis/CancerDeconvolution/Results/tosti_ablation_study_paad", 
                                      clinical_char = PAAD_clinical_char, 
                                      bulk_data = PAAD_bulk, bulk_meta = PAAD_meta, 
                                      sc_data = qc_tosti_sc$sc.eset.qc, ensemble = FALSE, 
                                      multiple_donors = TRUE, nrep = 1000, ncores = 15)

save.image("~/Masterthesis/Workspaces/ablation_study_tosti.RData")

