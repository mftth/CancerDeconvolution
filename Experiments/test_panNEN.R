## Mastherthesis, Melanie Fattohi
## test functions pn panNEN bulk RNA-seq and healthy neuroendocrine pancreatic scRNA-seq data
## bulk
### Alvarez/Califano
### Fadista (healthy)
### Missiaglia
### Diedisheim
### Sadanandam
### Scarpa
### Master
### Groetzinger
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
end_time - start_time ## 5.483678 hours
repset_ensemble_decon$p_value_wy_pearson ## 0.06293706
repset_ensemble_decon$p_value_wy_spearman ## 0.5054945
repset_ensemble_decon$p_value_wy_mad ## 0.8211788
repset_ensemble_decon$p_value_wy_rmsd ## 0.7692308

###################################################
## perform decon on raw bulk panNEN data

## read in Meta data
meta_data <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Meta_information.tsv", header = TRUE, sep = "\t")
## read in bulk panNEN data
alvarez <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Alvarez.S104.HGNC.tsv", header = TRUE, sep = "\t")
alvarez_meta <- meta_data[which(meta_data$Sample %in% colnames(alvarez)),]
rownames(alvarez_meta) <- alvarez_meta$Sample

diedisheim <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Diedisheim.S66.HGNC.tsv", header = TRUE, sep = "\t", row.names = 1)
diedisheim_meta <- meta_data[which(meta_data$Sample %in% colnames(diedisheim)),]
rownames(diedisheim_meta) <- diedisheim_meta$Sample

fadista <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Fadista.S89.tsv", header = TRUE, sep = "\t", row.names = 1)
fadista_meta <- meta_data[which(meta_data$Sample %in% colnames(fadista)),]
rownames(fadista_meta) <- fadista_meta$Sample

groetzinger <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Groetzinger.New.HGNC.S39.tsv", header = TRUE, sep = "\t", check.names = FALSE)
groetzinger_meta <- meta_data[which(meta_data$Sample %in% colnames(groetzinger)),]
rownames(groetzinger_meta) <- groetzinger_meta$Sample

master <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Master.S34.HGNC.tsv", header = TRUE, sep = "\t", check.names = FALSE)
master_meta <- meta_data[which(meta_data$Sample %in% colnames(master)),]
rownames(master_meta) <- master_meta$Sample

missiaglia <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Missaglia.S75.tsv", header = TRUE, sep = "\t", check.names = FALSE)
missiaglia_meta <- meta_data[which(meta_data$Sample %in% colnames(missiaglia)),]
rownames(missiaglia_meta) <- missiaglia_meta$Sample

sadanandam <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Sadanandam.S29.tsv", header = TRUE, sep = "\t", check.names = FALSE)
sadanandam_meta <- meta_data[which(meta_data$Sample %in% colnames(sadanandam)),]
rownames(sadanandam_meta) <- sadanandam_meta$Sample

scarpa <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Scarpa.New.HGNC.tsv", header = TRUE, sep = "\t", check.names = FALSE)
scarpa_meta <- meta_data[which(meta_data$Sample %in% colnames(scarpa)),]
rownames(scarpa_meta) <- scarpa_meta$Sample

## perform decon on raw data
start_time <- Sys.time()
alvarez_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = alvarez, 
                                        bulk_meta = alvarez_meta, sc_data = qc_baron$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 1.197906 hours
alvarez_baron_decon$p_value_wy_pearson ## 0.000999001
alvarez_baron_decon$p_value_wy_spearman ## 0.000999001
alvarez_baron_decon$p_value_wy_mad ## 0.000999001
alvarez_baron_decon$p_value_wy_rmsd ## 1
start_time <- Sys.time()
diedisheim_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = diedisheim, 
                                           bulk_meta = diedisheim_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 5.608439 Minuten
diedisheim_baron_decon$p_value_wy_pearson ## 0.8611389
diedisheim_baron_decon$p_value_wy_spearman ## 0.991009
diedisheim_baron_decon$p_value_wy_mad ## 0.984016
diedisheim_baron_decon$p_value_wy_rmsd ## 0.99001
start_time <- Sys.time()
fadista_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = fadista, 
                                        bulk_meta = fadista_meta, sc_data = qc_baron$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 24.25848 Minuten
fadista_baron_decon$p_value_wy_pearson ## 0.000999001
fadista_baron_decon$p_value_wy_spearman ## 0.000999001
fadista_baron_decon$p_value_wy_mad ## 0.000999001
fadista_baron_decon$p_value_wy_rmsd ## 0.001998002
start_time <- Sys.time()
groetzinger_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = groetzinger, 
                                            bulk_meta = groetzinger_meta, sc_data = qc_baron$sc.eset.qc, 
                                            cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 10.13119 Minuten
groetzinger_baron_decon$p_value_wy_pearson ## 0.000999001
groetzinger_baron_decon$p_value_wy_spearman ## 0.000999001
groetzinger_baron_decon$p_value_wy_mad ## 0.000999001
groetzinger_baron_decon$p_value_wy_rmsd ## 0.999001
start_time <- Sys.time()
master_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = master, 
                                        bulk_meta = master_meta, sc_data = qc_baron$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 6.056476 Minuten
master_baron_decon$p_value_wy_pearson ## 0.03696304
master_baron_decon$p_value_wy_spearman ## 0.4305694
master_baron_decon$p_value_wy_mad ## 0.003996004
master_baron_decon$p_value_wy_rmsd ## 0.982018
start_time <- Sys.time()
missiaglia_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = missiaglia, 
                                           bulk_meta = missiaglia_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 10.82134 Minuten
missiaglia_baron_decon$p_value_wy_pearson ## 0.000999001
missiaglia_baron_decon$p_value_wy_spearman ## 0.7162837
missiaglia_baron_decon$p_value_wy_mad ## 0.000999001
missiaglia_baron_decon$p_value_wy_rmsd ## 0.999001
start_time <- Sys.time()
sadanandam_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = sadanandam, 
                                           bulk_meta = sadanandam_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 6.453358 Minuten
sadanandam_baron_decon$p_value_wy_pearson ## 0.000999001
sadanandam_baron_decon$p_value_wy_spearman ## 0.02197802
sadanandam_baron_decon$p_value_wy_mad ## 0.04095904
sadanandam_baron_decon$p_value_wy_rmsd ## 0.98002
start_time <- Sys.time()
scarpa_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = scarpa, 
                                           bulk_meta = scarpa_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 9.52831 Minuten
scarpa_baron_decon$p_value_wy_pearson ## 0.002997003
scarpa_baron_decon$p_value_wy_spearman ## 0.000999001
scarpa_baron_decon$p_value_wy_mad ## 0.000999001
scarpa_baron_decon$p_value_wy_rmsd ## 0.8241758


###################################################
## perform decon on DESeq2 normalized bulk panNEN data

# normalize data
library(DESeq2)
library(edgeR)
library(limma)
DGE <- edgeR::DGEList(alvarez)
v <- limma::voom(DGE,design = NULL)
alvarez_norm <- v$E
DGE <- edgeR::DGEList(diedisheim)
v <- limma::voom(DGE,design = NULL)
diedisheim_norm <- v$E
DGE <- edgeR::DGEList(fadista)
v <- limma::voom(DGE,design = NULL)
fadista_norm <- v$E
DGE <- edgeR::DGEList(groetzinger)
v <- limma::voom(DGE,design = NULL)
groetzinger_norm <- v$E
DGE <- edgeR::DGEList(master)
v <- limma::voom(DGE,design = NULL)
master_norm <- v$E
DGE <- edgeR::DGEList(missiaglia)
v <- limma::voom(DGE,design = NULL)
missiaglia_norm <- v$E
DGE <- edgeR::DGEList(sadanandam)
v <- limma::voom(DGE,design = NULL)
sadanandam_norm <- v$E
DGE <- edgeR::DGEList(scarpa)
v <- limma::voom(DGE,design = NULL)
scarpa_norm <- v$E

start_time <- Sys.time()
alvarez_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = alvarez_norm, 
                                        bulk_meta = alvarez_meta, sc_data = qc_baron$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 55.03348 hours
alvarez_norm_baron_decon$p_value_wy_pearson ## 0.000999001
alvarez_norm_baron_decon$p_value_wy_spearman ## 0.000999001
alvarez_norm_baron_decon$p_value_wy_mad ## 0.2037962
alvarez_norm_baron_decon$p_value_wy_rmsd ## 0.8201798
start_time <- Sys.time()
diedisheim_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = diedisheim_norm, 
                                           bulk_meta = diedisheim_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 5.608439 Minuten
diedisheim_norm_baron_decon$p_value_wy_pearson ## 0.8711289
diedisheim_norm_baron_decon$p_value_wy_spearman ## 0.985015
diedisheim_norm_baron_decon$p_value_wy_mad ## 0.000999001
diedisheim_norm_baron_decon$p_value_wy_rmsd ## 0.987013
start_time <- Sys.time()
fadista_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = fadista_norm, 
                                        bulk_meta = fadista_meta, sc_data = qc_baron$sc.eset.qc, 
                                        cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 
fadista_norm_baron_decon$p_value_wy_pearson ## 
fadista_norm_baron_decon$p_value_wy_spearman ## 
fadista_norm_baron_decon$p_value_wy_mad ## 
fadista_norm_baron_decon$p_value_wy_rmsd ## 
start_time <- Sys.time()
groetzinger_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = groetzinger_norm, 
                                            bulk_meta = groetzinger_meta, sc_data = qc_baron$sc.eset.qc, 
                                            cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 7.573201 Minuten
groetzinger_norm_baron_decon$p_value_wy_pearson ## 0.000999001
groetzinger_norm_baron_decon$p_value_wy_spearman ## 0.000999001
groetzinger_norm_baron_decon$p_value_wy_mad ## 0.1828172
groetzinger_norm_baron_decon$p_value_wy_rmsd ## 1
start_time <- Sys.time()
master_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = master_norm, 
                                       bulk_meta = master_meta, sc_data = qc_baron$sc.eset.qc, 
                                       cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 4.431733 Minuten
master_norm_baron_decon$p_value_wy_pearson ## 0.08291708
master_norm_baron_decon$p_value_wy_spearman ## 0.5914086
master_norm_baron_decon$p_value_wy_mad ## 0.000999001
master_norm_baron_decon$p_value_wy_rmsd ## 0.994006
start_time <- Sys.time()
missiaglia_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = missiaglia_norm, 
                                           bulk_meta = missiaglia_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 11.23623 Minuten
missiaglia_norm_baron_decon$p_value_wy_pearson ## 0.000999001
missiaglia_norm_baron_decon$p_value_wy_spearman ## 0.6873127
missiaglia_norm_baron_decon$p_value_wy_mad ## 0.000999001
missiaglia_norm_baron_decon$p_value_wy_rmsd ## 0.974026
start_time <- Sys.time()
sadanandam_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = sadanandam_norm, 
                                           bulk_meta = sadanandam_meta, sc_data = qc_baron$sc.eset.qc, 
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 5.97341 Minuten
sadanandam_norm_baron_decon$p_value_wy_pearson ## 0.000999001
sadanandam_norm_baron_decon$p_value_wy_spearman ## 0.06793207
sadanandam_norm_baron_decon$p_value_wy_mad ## 0.000999001
sadanandam_norm_baron_decon$p_value_wy_rmsd ## 0.9490509
start_time <- Sys.time()
scarpa_norm_baron_decon <- Calculate_pvalue(nrep = reps, ncores = nc, bulk_data = scarpa_norm, 
                                       bulk_meta = scarpa_meta, sc_data = qc_baron$sc.eset.qc, 
                                       cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
end_time <- Sys.time()
end_time - start_time ## 5.360151 Minuten
scarpa_norm_baron_decon$p_value_wy_pearson ## 0.01598402
scarpa_norm_baron_decon$p_value_wy_spearman ## 0.000999001
scarpa_norm_baron_decon$p_value_wy_mad ## 0.03296703
scarpa_norm_baron_decon$p_value_wy_rmsd ## 1


###################################################
## compare results in boxplot
## per measure one box plot
## in each box plot you find one pair for each dataset (raw vs normalized)
library(reshape2)
library(ggplot2)

my_env <- ls()
my_env <- my_env[grep("decon", my_env)]
my_env <- my_env[grep("baron", my_env)]
my_env <- my_env[-grep("repset", my_env)]
my_env <- my_env[-grep("fadista", my_env)]

pearson_pval <- matrix(ncol = 3)
colnames(pearson_pval) <- c("dataset", "state", "value")
spearman_pval <- matrix(ncol = 3)
colnames(spearman_pval) <- c("dataset", "state", "value")
mad_pval <- matrix(ncol = 3)
colnames(mad_pval) <- c("dataset", "state", "value")
rmsd_pval <- matrix(ncol = 3)
colnames(rmsd_pval) <- c("dataset", "state", "value")
for (deconres in my_env) {
  name <- deconres
  decon_object <- get(deconres)
  if(grepl("norm", name)){
    state <- "norm"
  } else {
    state <- "raw"
  }
  
  dataset <- strsplit(name, "_")[[1]][1]
  
  tmp_matrix_pearson <- matrix(c(rep(dataset, nrow(decon_object$p_value_per_sample)), 
                         rep(state, nrow(decon_object$p_value_per_sample)), 
                         decon_object$p_value_per_sample$Pearson),
                       ncol = 3, nrow = nrow(decon_object$p_value_per_sample)) 
  tmp_matrix_spearman <- matrix(c(rep(dataset, nrow(decon_object$p_value_per_sample)), 
                                 rep(state, nrow(decon_object$p_value_per_sample)), 
                                 decon_object$p_value_per_sample$Spearman),
                               ncol = 3, nrow = nrow(decon_object$p_value_per_sample)) 
  tmp_matrix_mad <- matrix(c(rep(dataset, nrow(decon_object$p_value_per_sample)), 
                                 rep(state, nrow(decon_object$p_value_per_sample)), 
                                 decon_object$p_value_per_sample$mAD),
                               ncol = 3, nrow = nrow(decon_object$p_value_per_sample)) 
  tmp_matrix_rmsd <- matrix(c(rep(dataset, nrow(decon_object$p_value_per_sample)), 
                                 rep(state, nrow(decon_object$p_value_per_sample)), 
                                 decon_object$p_value_per_sample$RMSD),
                               ncol = 3, nrow = nrow(decon_object$p_value_per_sample)) 
  
  pearson_pval <- rbind(pearson_pval, tmp_matrix_pearson)
  spearman_pval <- rbind(spearman_pval, tmp_matrix_spearman)
  mad_pval <- rbind(mad_pval, tmp_matrix_mad)
  rmsd_pval <- rbind(rmsd_pval, tmp_matrix_rmsd)
}

pearson_pval <- as.data.frame(pearson_pval[-1,])
pearson_pval$value <- as.numeric(pearson_pval$value)
spearman_pval <- as.data.frame(spearman_pval[-1,])
spearman_pval$value <- as.numeric(spearman_pval$value)
rmsd_pval <- as.data.frame(rmsd_pval[-1,])
rmsd_pval$value <- as.numeric(rmsd_pval$value)
mad_pval <- as.data.frame(mad_pval[-1,])
mad_pval$value <- as.numeric(mad_pval$value)

pearson_pval$state <- factor(pearson_pval$state, levels = c("raw", "norm"))
spearman_pval$state <- factor(spearman_pval$state, levels = c("raw", "norm"))
mad_pval$state <- factor(mad_pval$state, levels = c("raw", "norm"))
rmsd_pval$state <- factor(rmsd_pval$state, levels = c("raw", "norm"))

# ggplot(pearson_pval, aes(x = dataset, y = value, fill = state)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(title = "p-value of pearson correlation")

ggplot(spearman_pval, aes(x = dataset, y = value, fill = state)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "p-value of spearman correlation")

ggplot(mad_pval, aes(x = dataset, y = value, fill = state)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "p-value of mad")

ggplot(rmsd_pval, aes(x = dataset, y = value, fill = state)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "p-value of rmsd")
