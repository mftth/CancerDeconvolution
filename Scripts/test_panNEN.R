## Mastherthesis, Melanie Fattohi
## test functions pn panNEN bulk RNA-seq and healthy neuroendocrine pancreatic scRNA-seq data
## bulk
### Alvarez/Califano
### Fadista (healthy)
### Missiaglia
### Sadnanandam
### Scarpa
### Riemer
### RepSet
## scRNA
### Baron
### Segerstolpe
### Lawlor (one donor)
### Tosti (ablation study)

source("~/Masterthesis/CancerDeconvolution/Scripts/Quality_control.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")

library(parallel)
library(robustbase)
library(caret)
library(skimr)

set.seed(4)

## load qc'ed scRNA-seq references
qc_baron <- readRDS(file = "~/Praktikum/Data/Baron/qc_baron_exo.RDS")
qc_segerstolpe <- readRDS("~/Praktikum/Deko_SCDC/Training_Data/Segerstolpe_qc_exo.RDS")
qc_lawlor <- readRDS("~/Praktikum/Deko_SCDC/Training_Data/Lawlor_qc_exo.RDS")
pannen_reference <- list("baron" = qc_baron$sc.eset.qc, "segerstolpe" = qc_segerstolpe$sc.eset.qc,
                         "lawlor" = qc_lawlor$sc.eset.qc)
pannen_reference_donors <- list(TRUE, TRUE, FALSE)
qc_tosti <- readRDS("~/Tosti/Tosti_qc_complete.RDS")

## load bulk RNA-seq
## Repset
repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")

## set conditions for deconvolution
cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")
reps <- 1000
nc <- 10

## perform deconvolution + p-value calculation
start_time <- Sys.time()
repset_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                       bulk_meta = repset_meta, sc_data = qc_baron$sc.eset.qc, 
                                       cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 11.93 Minuten
repset_baron_decon$p_value_wy_pearson ## 0.000999001
repset_baron_decon$p_value_wy_spearman ## 0.000999001
repset_baron_decon$p_value_wy_mad ## 0.1428571
repset_baron_decon$p_value_wy_rmsd ## 1
start_time <- Sys.time()
repset_segerstolpe_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                             bulk_meta = repset_meta, 
                                             sc_data = qc_segerstolpe$sc.eset.qc, cell_types = cts, 
                                             ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 11.18751 Minuten
repset_segerstolpe_decon$p_value_wy_pearson ## 0.982018
repset_segerstolpe_decon$p_value_wy_spearman ## 0.000999001
repset_segerstolpe_decon$p_value_wy_mad ## 0.1488511
repset_segerstolpe_decon$p_value_wy_rmsd ## 1
start_time <- Sys.time()
repset_lawlor_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                        bulk_meta = repset_meta, sc_data = qc_lawlor$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = FALSE)
end_time <- Sys.time()
end_time - start_time  ## 13.68682 Minuten
repset_lawlor_decon$p_value_wy_pearson ## 0.000999001
repset_lawlor_decon$p_value_wy_spearman ## 0.000999001
repset_lawlor_decon$p_value_wy_mad ## 0.1128871
repset_lawlor_decon$p_value_wy_rmsd ## 0.6213786
start_time <- Sys.time()
repset_tosti_ac_s_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                            bulk_meta = repset_meta, sc_data = qc_tosti$sc.eset.qc, 
                                            cell_types = c("alpha", "beta", "gamma", "delta", 
                                                           "acinar_s", "ductal"), 
                                            ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 9.839013 Minuten
repset_tosti_ac_s_decon$p_value_wy_pearson ## 0.989011
repset_tosti_ac_s_decon$p_value_wy_spearman ## 0.6643357 
repset_tosti_ac_s_decon$p_value_wy_mad ## 0.5224775
repset_tosti_ac_s_decon$p_value_wy_rmsd ## 0.98002
start_time <- Sys.time()
repset_tosti_ac_i_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                            bulk_meta = repset_meta, sc_data = qc_tosti$sc.eset.qc, 
                                            cell_types = c("alpha", "beta", "gamma", "delta", 
                                                           "acinar_i", "ductal"), 
                                            ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 10.00426 Minuten
repset_tosti_ac_i_decon$p_value_wy_pearson ## 0.9370629 
repset_tosti_ac_i_decon$p_value_wy_spearman ## 0.7012987
repset_tosti_ac_i_decon$p_value_wy_mad ## 0.4215784
repset_tosti_ac_i_decon$p_value_wy_rmsd ## 0.974026
start_time <- Sys.time()
repset_tosti_ac_r_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset, 
                                            bulk_meta = repset_meta, sc_data = qc_tosti$sc.eset.qc, 
                                            cell_types = c("alpha", "beta", "gamma", "delta", 
                                                           "acinar_reg", "ductal"), 
                                            ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 10.73886
repset_tosti_ac_r_decon$p_value_wy_pearson ## 0.989011
repset_tosti_ac_r_decon$p_value_wy_spearman ## 0.6703297
repset_tosti_ac_r_decon$p_value_wy_mad ## 0.5284715
repset_tosti_ac_r_decon$p_value_wy_rmsd ## 0.974026
start_time <- Sys.time()
repset_ensemble_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = repset,
                                          bulk_meta = repset_meta, sc_data = pannen_reference,
                                          cell_types = cts, ensemble = TRUE, 
                                          multiple_donors = pannen_reference_donors)
end_time <- Sys.time()
end_time - start_time ##
repset_ensemble_decon$p_value_wy_pearson ## 
repset_ensemble_decon$p_value_wy_spearman ## 
repset_ensemble_decon$p_value_wy_mad ## 
repset_ensemble_decon$p_value_wy_rmsd ## 
