## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## Deconvolve Flowers and Bailey
## nreps = 1000, ncores = 15
## survival analysis
## correlation analysis
## ML anaylsis

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/survival_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/correlation_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")

## bulk RNA-seq datasets
## 7 samples; survival; tumor grading; RNA-seq (mouse)
Flowers_bulk <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_bulk.RDS")
Flowers_meta <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_metadata.RDS")
rownames(Flowers_meta) <- colnames(Flowers_bulk)
Flowers_meta$COO <- sapply(rownames(Flowers_meta), function(x) strsplit(x, split = "\\.")[[1]][1])
rownames(Flowers_bulk) <- toupper(rownames(Flowers_bulk))

bulk_list <- list("Flowers" = Flowers_bulk)
bulk_meta_list <- list("Flowers" = Flowers_meta)

qc_baron_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_baron_exo.RDS")
qc_segerstolpe_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/Segerstolpe_qc_exo.RDS")
qc_lawlor_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/Lawlor_qc_exo.RDS")
qc_tosti_sc <- readRDS(file = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_tosti.RDS")

## perform deconvolution
reps <- 1000
ncores <- 15
cts <- c("alpha", "beta", "gamma", "delta", "acinar", "ductal")
cts_tosti <- c("sacinar", "racinar", "iacinar", "ductal", "mductal", "alpha", "beta", "gamma", "delta")

###
res_path_baron <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Baron"
decon_baron <- lapply(1:length(bulk_list), function(x) {
  decon_baron_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = qc_baron_sc$sc.eset.qc, cell_types = cts,
                                    ensemble = FALSE, multiple_donors = TRUE)
  saveRDS(decon_baron_x, file = paste(res_path_baron, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

decon_baron <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_baron, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_baron) <- names(bulk_list)
###
###
res_path_tosti <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti"
decon_tosti <- lapply(1:length(bulk_list), function(x) {
  decon_tosti_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = qc_tosti_sc$sc.eset.qc, cell_types = cts_tosti,
                                    ensemble = FALSE, multiple_donors = TRUE)
  saveRDS(decon_tosti_x, file = paste(res_path_tosti, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

decon_tosti <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_tosti, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_tosti) <- names(bulk_list)
###
###
baron_flowers_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Flowers,
                                                  clinical_characteristics = data.frame("COO" = Flowers_meta$COO,
                                                                                        row.names= rownames(Flowers_meta)))
baron_flowers_prop_bar <- barplot_proportions(decon_output = decon_baron$Flowers,
                                                  clinical_characteristics = data.frame("COO" = Flowers_meta$COO,
                                                                                        row.names= rownames(Flowers_meta)))

tosti_flowers_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Flowers,
                                                  clinical_characteristics = data.frame("COO" = Flowers_meta$COO,
                                                                                        row.names= rownames(Flowers_meta)))
tosti_flowers_prop_bar <- barplot_proportions(decon_output = decon_tosti$Flowers,
                                              clinical_characteristics = data.frame("COO" = Flowers_meta$COO,
                                                                                    row.names= rownames(Flowers_meta)))
