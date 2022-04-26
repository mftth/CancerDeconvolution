## Mastherthesis, Melanie Fattohi
## scRNA-seq data by SeneSys
## Deconvolve Chapuy, Schleich and Schmitz (maybe Reddy)
## nreps = 1000, ncores = 15
## survival analysis
## correlation analysis
## ML anaylsis


source("~/Masterthesis/CancerDeconvolution/Scripts/Quality_control.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/survival_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/correlation_analysis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_Deconvolution.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/Execute_MachineLearning.R")


## bulk data
## Chapuy: 137 samples;  OS, OS_stat; Cluster; COO
Chapuy_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_bulk.RDS")
Chapuy_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_metadata.RDS")
## Schleich: 109 samples; treatment response
Schleich_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_bulk.RDS")
Schleich_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_metadata.RDS")
rownames(Schleich_bulk) <- toupper(rownames(Schleich_bulk))
## Schmitz 481 samples; survival; COO; treatment
Schmitz_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_bulk.RDS")
Schmitz_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_metadata.RDS")

bulk_list <- list("Chapuy" = Chapuy_bulk,
                  "Schleich" = Schleich_bulk,
                  "Schmitz" = Schmitz_bulk)
bulk_meta_list <- list("Chapuy" = Chapuy_meta,
                       "Schleich" = Schleich_meta,
                       "Schmitz" = Schmitz_meta)

## SeneSys scRNA-seq data
mouse_senescence <- read.table("~/SeneSys_scRNA_mouse/senescence_mouse.tsv", sep = "\t", header = TRUE, row.names = 1)
mouse_senescence <- as.matrix(mouse_senescence)
rownames(mouse_senescence) <- toupper(rownames(mouse_senescence))
phenotypes <- sapply(colnames(mouse_senescence), function(x) strsplit(x, "_\\d+")[[1]][1])
phenotypes <- as.data.frame(phenotypes)
phenotypes$sample <- rep("donor", nrow(phenotypes))
colnames(phenotypes)[1] <- "cluster"
qc_mouse_senescence <- Quality_control(sc_data = mouse_senescence, sc_meta = phenotypes, 
                                       sc_path = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_mouse_senescence.RDS",
                                       cell_types = unique(phenotypes$cluster),
                                       multiple_donors = FALSE)

## deconvolution
reps <- 1000
ncores <- 15
cts <- c("ADR", "ADR_OHT", "PS") ## ueberpruefe nochmal

res_path <- "~/Masterthesis/CancerDeconvolution/Results/DLBCL_deconvolution2"
decon_dlbcl <- lapply(1:length(bulk_list), function(x) {
  decon_dlbcl_x <- Calculate_pvalue(nrep = reps, ncores = ncores, silent = FALSE, 
                                    bulk_data = bulk_list[[x]], bulk_meta = bulk_meta_list[[x]],
                                    sc_data = qc_mouse_senescence$sc.eset.qc, cell_types = cts,
                                    ensemble = FALSE, multiple_donors = FALSE)
  saveRDS(decon_dlbcl_x, file = paste(res_path_baron, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})

## plot p-values
technology = c("microarray", "microarray", "RNA-seq")
pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                        pvalue_type = "Spearman", technology = technology) 
pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                       pvalue_type = "Pearson", technology = technology) 
pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                   pvalue_type = "mAD", technology = technology) 
pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                    pvalue_type = "RMSD", technology = technology) 

## visualize ct props in heatmaps
## for chapuy: coo und cluster und treatment
## for schleich: treatment and treatment response
## for schmitz: coo, treatment and IPI group
baron_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Guo,
                                              clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, 
                                                                                    "MKI67" = Guo_mki67,
                                                                                    row.names = rownames(Guo_meta)),
                                              annotation_colors = guo_annot_colors, fontsize = 11)