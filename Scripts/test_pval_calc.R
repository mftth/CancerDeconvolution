## Mastherthesis, Melanie Fattohi
## test Calculate_pvalue function
## 1) random data; expectation: no significant p-value
## 2) known/simulated data; expectation: p-value has expected value
## 3) real data

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")


#library(SCDC)
library(parallel)
library(robustbase)

set.seed(40)

## read in real bulk RNA-seq data and set some variables
# repset <- read.table(file = "~/Praktikum/Data/RepSet.S57.tsv", header = TRUE, 
#                      sep = "\t",)
# pannen_meta <- read.table(file = "~/Praktikum/Data/Clinial_data_Meta_information.tsv", 
#                           header = TRUE, sep = "\t")
# scarpa_meta <- pannen_meta[which(pannen_meta$Study == "Scarpa"),]
# riemer_meta <- pannen_meta[which(pannen_meta$Study == "Riemer"),]
# riemer_meta <- riemer_meta[which(riemer_meta$Subtype == "Cancer"),]
# repset_meta <- rbind(scarpa_meta, riemer_meta)
# repset_meta$Name <- as.character(repset_meta$Name)
# repset_meta$Name[1:55] <- paste("X", repset_meta$Name[1:55], sep = "")
# repset_meta <- repset_meta[match(colnames(repset), repset_meta$Name),]
# all(repset_meta$Name == colnames(repset))
# rownames(repset_meta) <- repset_meta$Name
repset <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")

cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")
#scpath1 <- "~/Praktikum/Deko_SCDC/Training_Data/Baron_qc_exo.RDS"
scpath2 <- "~/Praktikum/Deko_SCDC/Training_Data/Lawlor_qc_exo.RDS"
scpath3 <- "~/Praktikum/Deko_SCDC/Training_Data/Segerstolpe_qc_exo.RDS"
reps <- 50
ncores <- 15

## perform QC on Baron
#baron <- readRDS(file = "~/Praktikum/Data/Baron/Baron.RDS")
#baron_meta <- readRDS(file = "~/Praktikum/Data/Baron/Baron_meta.RDS")
#qc_baron <- Quality_control(sc_data = baron, sc_meta = baron_meta, sc_path = "~/Praktikum/Data/Baron/qc_baron_exo.RDS",
#                            multiple_donors = TRUE, ct.varname = "cluster", sample = "sample", ct.sub = cts)
qc_baron <- readRDS(file = "~/Praktikum/Data/Baron/qc_baron_exo.RDS")

## 1) random data; expectation: no significant p-value
## create matrix of same size as repset, but random values
random_bulk <- matrix(runif(nrow(repset)*ncol(repset)), nrow=nrow(repset))
rownames(random_bulk) <- rownames(repset)
colnames(random_bulk) <- colnames(repset)

decon_pval_random <- Calculate_pvalue(nrep = reps, ncores = ncores,  bulk_data = random_bulk, sc_data = qc_baron$sc.eset.qc,
                                      bulk_meta = repset_meta, cell_types = cts, ensemble = FALSE,
                                      multiple_donors = TRUE)
decon_pval_random$p_value_wy_pearson  
decon_pval_random$p_value_wy_spearman 
decon_pval_random$p_value_wy_mad 
decon_pval_random$p_value_wy_rmsd 
#View(decon_pval_random$decon_res$prop.est.mvw)


## 2) known/simulated data; expectation: p-value has expected value
#qc_baron <- readRDS(scpath1)
qc_segerstolpe <- readRDS(scpath3)
pseudo_bulk <- generateBulk_allcells(qc_segerstolpe$sc.eset.qc, ct.varname = "cluster", sample = "sample", 
                                     ct.sub = cts)
pseudo_bulk_rand <- generateBulk_norep(qc_segerstolpe$sc.eset.qc, ct.varname = "cluster", sample = "sample", 
                                       ct.sub = cts, nbulk = 50)

decon_pval_pseudo <- Calculate_pvalue(nrep = reps, ncores = ncores, bulk_data = pseudo_bulk$pseudo_eset@assayData$exprs, 
                                      sc_data = qc_baron$sc.eset.qc, bulk_meta = pseudo_bulk$pseudo_eset@phenoData@data,
                                      cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
decon_pval_pseudo$p_value_wy_pearson  
decon_pval_pseudo$p_value_wy_spearman 
decon_pval_pseudo$p_value_wy_mad 
decon_pval_pseudo$p_value_wy_rmsd 
SCDC_peval(ptrue = pseudo_bulk$truep, pest = decon_pval_pseudo$decon_res$prop.est.mvw, 
           pest.names = "pseudo_bulk")$evals.table
#               RMSD     mAD      R
# pseudo_bulk 0.0397 0.02933 0.9669
# 0.11334 0.08213 0.9292


decon_pval_pseudo_rand <- Calculate_pvalue(nrep = reps, ncores = ncores, bulk_data = pseudo_bulk_rand$pseudo_bulk, 
                                           sc_data = qc_baron$sc.eset.qc, bulk_meta = pseudo_bulk_rand$pseudo_eset@phenoData@data,
                                           cell_types = cts, ensemble = FALSE, multiple_donors = TRUE)
decon_pval_pseudo_rand$p_value_wy_pearson  
decon_pval_pseudo_rand$p_value_wy_spearman 
decon_pval_pseudo_rand$p_value_wy_mad 
decon_pval_pseudo_rand$p_value_wy_rmsd 
SCDC_peval(ptrue = pseudo_bulk_rand$true_p, pest = decon_pval_pseudo_rand$decon_res$prop.est.mvw, 
           pest.names = "pseudo_bulk_rand")$evals.table
#                     RMSD     mAD      R
# pseudo_bulk_rand 0.11276 0.08166 0.9341
# 0.117 0.08222 0.9025


## 3) real data
decon_pval <- Calculate_pvalue(nrep = reps, ncores = ncores, bulk_data = repset, bulk_meta = repset_meta,
                               cell_types =  cts, sc_data = qc_baron$sc.eset.qc,
                               ensemble = FALSE, multiple_donors = TRUE)
decon_pval$p_value_wy_pearson  
decon_pval$p_value_wy_spearman 
decon_pval$p_value_wy_mad 
decon_pval$p_value_wy_rmsd 

