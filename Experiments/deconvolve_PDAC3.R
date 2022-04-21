## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## scRNA-seq von PDAC werden noch gesucht
## Deconvolve Yang, PAAD, Moffitt, Guo
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
## 183 samples; survival; tumor grading; RNA-seq
PAAD_bulk <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_raw_bulk.RDS")
PAAD_meta <- readRDS("~/Masterthesis/Data/Bulk/PAAD/PAAD_metadata.RDS")
## 69 samples; survival; tumor grading; microarray
Yang_bulk <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_bulk.RDS")
Yang_meta <- readRDS("~/Masterthesis/Data/Bulk/Yang/Yang_metadata.RDS")
## 62 samples; survival; tumor subtype (basal, classical, hybrid); RNA-seq
Guo_bulk <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_bulk.RDS")
Guo_meta <- readRDS("~/Masterthesis/Data/Bulk/Guo/Guo_metadata.RDS")
## 357 samples; survival; tumor subtype; microarray
Moffitt_array_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_bulk.RDS")
Moffitt_array_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_metadata.RDS")
## 7 samples; survival; tumor grading; RNA-seq (mouse)
Flowers_bulk <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_bulk.RDS")
Flowers_meta <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_metadata.RDS")
rownames(Flowers_meta) <- colnames(Flowers_bulk)
Flowers_meta$COO <- sapply(rownames(Flowers_meta), function(x) strsplit(x, split = "\\.")[[1]][1])
rownames(Flowers_bulk) <- toupper(rownames(Flowers_bulk))

bulk_list <- list("PAAD" = PAAD_bulk,
                  "Yang" = Yang_bulk,
                  "Guo" = Guo_bulk,
                  #"Janky" = Janky_bulk,
                  #"Kirby" = Kirby_bulk,
                  #"Moffitt_seq" = Moffitt_seq_bulk,
                  "Moffitt_array" = Moffitt_array_bulk,
                  "Flowers" = Flowers_bulk)
bulk_meta_list <- list("PAAD" = PAAD_meta,
                       "Yang" = Yang_meta,
                       "Guo" = Guo_meta,
                       #"Janky" = Janky_meta,
                       #"Kirby" = Kirby_meta,
                       #"Moffitt_seq" = Moffitt_seq_meta,
                       "Moffitt_array" = Moffitt_array_meta,
                       "Flowers" = Flowers_meta)

res_path_baron <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Baron"
decon_baron <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_baron, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_baron) <- names(bulk_list)
res_path_tosti <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Tosti"
decon_tosti <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path_tosti, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_tosti) <- names(bulk_list)
# res_path_ensemble <- "~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Ensemble"
# decon_ensemble <- lapply(1:length(bulk_list), 
#                          function(x) readRDS(file = paste(res_path_ensemble, "/", names(bulk_list)[x], 
#                                                           "_decon.RDS", sep = "")))
# names(decon_ensemble) <- names(bulk_list)
guo_ensemble <- readRDS("~/Masterthesis/CancerDeconvolution/Results/PDAC_deconvolution/Ensemble/guo_ensemble.RDS")

## include only primary cancer samples in Moffitt
Moffitt_array_primary <- which(Moffitt_array_meta$source_name_ch2 == "Pancreas_Primary")
Moffitt_array_bulk <- Moffitt_array_bulk[,Moffitt_array_primary]
Moffitt_array_meta <- Moffitt_array_meta[Moffitt_array_primary,]
decon_baron$Moffitt_array$p_value_per_sample <- decon_baron$Moffitt_array$p_value_per_sample[Moffitt_array_primary,]
decon_baron$Moffitt_array$decon_res$prop.est.mvw <- decon_baron$Moffitt_array$decon_res$prop.est.mvw[Moffitt_array_primary,]
decon_baron$Moffitt_array$statistics_observed$pearson_vec <- decon_baron$Moffitt_array$statistics_observed$pearson_vec[Moffitt_array_primary]
decon_baron$Moffitt_array$statistics_observed$spearman_vec <- decon_baron$Moffitt_array$statistics_observed$spearman_vec[Moffitt_array_primary]
decon_baron$Moffitt_array$statistics_observed$mad_vec <- decon_baron$Moffitt_array$statistics_observed$mad_vec[Moffitt_array_primary]
decon_baron$Moffitt_array$statistics_observed$rmsd_vec <- decon_baron$Moffitt_array$statistics_observed$rmsd_vec[Moffitt_array_primary]
decon_tosti$Moffitt_array$p_value_per_sample <- decon_tosti$Moffitt_array$p_value_per_sample[Moffitt_array_primary,]
decon_tosti$Moffitt_array$decon_res$prop.est.mvw <- decon_tosti$Moffitt_array$decon_res$prop.est.mvw[Moffitt_array_primary,]
decon_tosti$Moffitt_array$statistics_observed$pearson_vec <- decon_tosti$Moffitt_array$statistics_observed$pearson_vec[Moffitt_array_primary]
decon_tosti$Moffitt_array$statistics_observed$spearman_vec <- decon_tosti$Moffitt_array$statistics_observed$spearman_vec[Moffitt_array_primary]
decon_tosti$Moffitt_array$statistics_observed$mad_vec <- decon_tosti$Moffitt_array$statistics_observed$mad_vec[Moffitt_array_primary]
decon_tosti$Moffitt_array$statistics_observed$rmsd_vec <- decon_tosti$Moffitt_array$statistics_observed$rmsd_vec[Moffitt_array_primary]
Moffitt_array_meta$tumor_subtype <- rep(NA, nrow(Moffitt_array_meta))
Moffitt_array_meta$tumor_subtype[Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2` == "1"] <- "Classical" 
Moffitt_array_meta$tumor_subtype[Moffitt_array_meta$`tumor_subtype_0na_1classical_2basal:ch2` == "2"] <- "Basal" 
write.table(Moffitt_array_bulk, "/mnt/home/melanie/Masterthesis/Data/Bulk/Moffitt/Moffitt_array_bulk2.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## include only G1,G2,G3 in PAAD and Yang
Yang_grading <- which(Yang_meta$`grading:ch1` %in% c("G2", "G3"))
Yang_bulk <- Yang_bulk[,Yang_grading]
Yang_meta <- Yang_meta[Yang_grading,]
decon_baron$Yang$p_value_per_sample <- decon_baron$Yang$p_value_per_sample[Yang_grading,]
decon_baron$Yang$decon_res$prop.est.mvw <- decon_baron$Yang$decon_res$prop.est.mvw[Yang_grading,]
decon_baron$Yang$statistics_observed$pearson_vec <- decon_baron$Yang$statistics_observed$pearson_vec[Yang_grading]
decon_baron$Yang$statistics_observed$spearman_vec <- decon_baron$Yang$statistics_observed$spearman_vec[Yang_grading]
decon_baron$Yang$statistics_observed$mad_vec <- decon_baron$Yang$statistics_observed$mad_vec[Yang_grading]
decon_baron$Yang$statistics_observed$rmsd_vec <- decon_baron$Yang$statistics_observed$rmsd_vec[Yang_grading]
decon_tosti$Yang$p_value_per_sample <- decon_tosti$Yang$p_value_per_sample[Yang_grading,]
decon_tosti$Yang$decon_res$prop.est.mvw <- decon_tosti$Yang$decon_res$prop.est.mvw[Yang_grading,]
decon_tosti$Yang$statistics_observed$pearson_vec <- decon_tosti$Yang$statistics_observed$pearson_vec[Yang_grading]
decon_tosti$Yang$statistics_observed$spearman_vec <- decon_tosti$Yang$statistics_observed$spearman_vec[Yang_grading]
decon_tosti$Yang$statistics_observed$mad_vec <- decon_tosti$Yang$statistics_observed$mad_vec[Yang_grading]
decon_tosti$Yang$statistics_observed$rmsd_vec <- decon_tosti$Yang$statistics_observed$rmsd_vec[Yang_grading]

PAAD_grading <- which(PAAD_meta$neoplasm_histologic_grade %in% c("g1", "g2", "g3"))
PAAD_bulk <- PAAD_bulk[,PAAD_grading]
PAAD_meta <- PAAD_meta[PAAD_grading,]
decon_baron$PAAD$p_value_per_sample <- decon_baron$PAAD$p_value_per_sample[PAAD_grading,]
decon_baron$PAAD$decon_res$prop.est.mvw <- decon_baron$PAAD$decon_res$prop.est.mvw[PAAD_grading,]
decon_baron$PAAD$statistics_observed$pearson_vec <- decon_baron$PAAD$statistics_observed$pearson_vec[PAAD_grading]
decon_baron$PAAD$statistics_observed$spearman_vec <- decon_baron$PAAD$statistics_observed$spearman_vec[PAAD_grading]
decon_baron$PAAD$statistics_observed$mad_vec <- decon_baron$PAAD$statistics_observed$mad_vec[PAAD_grading]
decon_baron$PAAD$statistics_observed$rmsd_vec <- decon_baron$PAAD$statistics_observed$rmsd_vec[PAAD_grading]
decon_tosti$PAAD$p_value_per_sample <- decon_tosti$PAAD$p_value_per_sample[PAAD_grading,]
decon_tosti$PAAD$decon_res$prop.est.mvw <- decon_tosti$PAAD$decon_res$prop.est.mvw[PAAD_grading,]
decon_tosti$PAAD$statistics_observed$pearson_vec <- decon_tosti$PAAD$statistics_observed$pearson_vec[PAAD_grading]
decon_tosti$PAAD$statistics_observed$spearman_vec <- decon_tosti$PAAD$statistics_observed$spearman_vec[PAAD_grading]
decon_tosti$PAAD$statistics_observed$mad_vec <- decon_tosti$PAAD$statistics_observed$mad_vec[PAAD_grading]
decon_tosti$PAAD$statistics_observed$rmsd_vec <- decon_tosti$PAAD$statistics_observed$rmsd_vec[PAAD_grading]
PAAD_meta$neoplasm_histologic_grade[which(PAAD_meta$neoplasm_histologic_grade == "g1")] <- "G1"
PAAD_meta$neoplasm_histologic_grade[which(PAAD_meta$neoplasm_histologic_grade == "g2")] <- "G2"
PAAD_meta$neoplasm_histologic_grade[which(PAAD_meta$neoplasm_histologic_grade == "g3")] <- "G3"

## import supplementary file from Hayashi et al. 2020
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

## visualize p-values
technology = c("RNA-seq", "microarray", "RNA-seq", "microarray", "RNA-seq")
baron_pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_baron,
                                              pvalue_type = "Spearman", technology = technology) + 
  ggtitle("Reference: Baron et al.") 
tosti_pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_tosti,
                                              pvalue_type = "Spearman", technology = technology) + 
  ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_spearman, tosti_pval_boxplot_spearman, common.legend = TRUE) 

baron_pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_baron,
                                             pvalue_type = "Pearson", technology = technology) + 
  ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_tosti,
                                             pvalue_type = "Pearson", technology = technology) + 
  ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_pearson, tosti_pval_boxplot_pearson, common.legend = TRUE)

baron_pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_baron,
                                         pvalue_type = "mAD", technology = technology) + 
  ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_tosti,
                                         pvalue_type = "mAD", technology = technology) + 
  ggtitle("Reference: Tosti et al.")
ggarrange(baron_pval_boxplot_mad, tosti_pval_boxplot_mad, common.legend = TRUE)

baron_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_baron,
                                          pvalue_type = "RMSD", technology = technology) + 
  ggtitle("Reference: Baron et al.")
tosti_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_tosti,
                                          pvalue_type = "RMSD", technology = technology) + 
  ggtitle("Reference: Tosti et al.") + theme(plot.title = element_text(size=12))
ggarrange(baron_pval_boxplot_rmsd, tosti_pval_boxplot_rmsd, common.legend = TRUE) 

guo_ensemble_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = list("Guo" = guo_ensemble),
                                                 pvalue_type = "RMSD", technology = "RNA-seq")
# guo_ensemble_pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = list("Guo" = guo_ensemble),
#                                                  pvalue_type = "Spearman", technology = "RNA-seq")
# guo_ensemble_pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = list("Guo" = guo_ensemble),
#                                                  pvalue_type = "Pearson", technology = "RNA-seq")
# guo_ensemble_pval_boxplot_mad <- boxplot_pvalue(decon_output_list = list("Guo" = guo_ensemble),
#                                                  pvalue_type = "mAD", technology = "RNA-seq")
decon_baron$Guo_ensemble <- guo_ensemble
baron_pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_baron,
                                          pvalue_type = "RMSD", technology = c(technology, "RNA-seq")) + 
  ggtitle("Reference: Baron et al.") + theme(plot.title = element_text(size=12))
ggarrange(baron_pval_boxplot_rmsd, tosti_pval_boxplot_rmsd, common.legend = TRUE) 


## visualize cell type props
Guo_mki67_thirds <- quantile(Guo_bulk["MKI67",], probs = seq(0, 1, 1/3))
Guo_mki67 <- as.numeric(Guo_bulk["MKI67",])
Guo_mki67[sapply(as.numeric(Guo_bulk["MKI67",]), function(x) x <= Guo_mki67_thirds[2])] <- "MKI67_low"
Guo_mki67[sapply(as.numeric(Guo_bulk["MKI67",]), function(x) x > Guo_mki67_thirds[2] && x <= Guo_mki67_thirds[3])] <- "MKI67_medium"
Guo_mki67[sapply(as.numeric(Guo_bulk["MKI67",]), function(x) x > Guo_mki67_thirds[3])] <- "MKI67_high"
Guo_mki67 <- factor(Guo_mki67, levels = c("MKI67_low", "MKI67_medium", "MKI67_high"))
guo_annot_colors <- list(tumor_subtype = c(Basal = "#f683ad", Classical = "#f8fc88", Hybrid ="#6eacf2"),
                         MKI67 = c(MKI67_low = "#77f387",  MKI67_medium = "#f19e5b", MKI67_high = "#d34545"))

baron_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Guo,
                                              clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description, 
                                                                                    "MKI67" = Guo_mki67,
                                                                                    row.names = rownames(Guo_meta)),
                                              annotation_colors = guo_annot_colors, fontsize = 11)
tosti_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Guo,
                                              clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                                    "MKI67" = Guo_mki67,
                                                                                    row.names = rownames(Guo_meta)),
                                              annotation_colors = guo_annot_colors, fontsize = 11)
ensemble_Guo_prop_heatmap <- heatmap_proportions(decon_output = guo_ensemble,
                                                 clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                                       "MKI67" = Guo_mki67,
                                                                                       row.names = rownames(Guo_meta)),
                                                 annotation_colors = guo_annot_colors, fontsize = 11,
                                                 clustering_method = "ward.D2")

# tosti_Guo_boxplot_prop <- boxplot_proportions(decon_output = decon_tosti$Guo,
#                                               clinical_characteristics_vec = Guo_meta$description,
#                                               cell_types = c("sacinar", "mductal", "racinar")) + 
#   geom_signif(comparisons = list(c("Hybrid", "Classical"), c("Hybrid", "Basal"), c("Classical", "Basal")), map_signif_level=TRUE)


PAAD_mki67_thirds <- quantile(PAAD_bulk["MKI67",], probs = seq(0, 1, 1/3))
PAAD_mki67 <- as.numeric(PAAD_bulk["MKI67",])
PAAD_mki67[sapply(as.numeric(PAAD_bulk["MKI67",]), function(x) x <= PAAD_mki67_thirds[2])] <- "MKI67_low"
PAAD_mki67[sapply(as.numeric(PAAD_bulk["MKI67",]), function(x) x > PAAD_mki67_thirds[2] && x <= PAAD_mki67_thirds[3])] <- "MKI67_medium"
PAAD_mki67[sapply(as.numeric(PAAD_bulk["MKI67",]), function(x) x > PAAD_mki67_thirds[3])] <- "MKI67_high"
PAAD_mki67 <- factor(PAAD_mki67, levels = c("MKI67_low", "MKI67_medium", "MKI67_high"))

PAAD_annot_colors <- list(Collisson = c(Exocrine = "#6eacf2", Classical = "#f8fc88", QM ="#f683ad"),
                          Bailey = c(ADEX = "#6eacf2", Immunogenic = "#f19e5b", Progenitor ="#f8fc88", Squamous = "#f683ad"),
                          Moffitt = c(Basal = "#f683ad", Classical = "#f8fc88"),
                          grading = c(G1 = "#77f387", G2 = "#f19e5b", G3 = "#d34545"),
                          MKI67 = c(MKI67_low = "#77f387", MKI67_medium ="#f19e5b",MKI67_high = "#d34545")) 

baron_PAAD_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                     "Moffitt" = PAAD_meta$tumor_moffitt,
                                                                                     "Collisson" = PAAD_meta$tumor_collisson,
                                                                                     "Bailey" = PAAD_meta$tumor_bailey,
                                                                                     "MKI67" = PAAD_mki67,
                                                                                     row.names = rownames(PAAD_meta)), 
                                               clustering_method = "ward.D2", annotation_colors = PAAD_annot_colors, fontsize = 11)

tosti_PAAD_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                     "Moffitt" = PAAD_meta$tumor_moffitt,
                                                                                     "Collisson" = PAAD_meta$tumor_collisson,
                                                                                     "Bailey" = PAAD_meta$tumor_bailey,
                                                                                     "MKI67" = PAAD_mki67,
                                                                                     row.names = rownames(PAAD_meta)), 
                                               clustering_method = "ward.D2", annotation_colors = PAAD_annot_colors, fontsize = 11)


baron_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)),
                                               fontsize = 11)
tosti_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)),
                                               fontsize = 11)


Moffitt_mki67_thirds <- quantile(Moffitt_array_bulk["MKI67",], probs = seq(0, 1, 1/3))
Moffitt_mki67 <- as.numeric(Moffitt_array_bulk["MKI67",])
Moffitt_mki67[sapply(as.numeric(Moffitt_array_bulk["MKI67",]), function(x) x <= Moffitt_mki67_thirds[2])] <- "MKI67_low"
Moffitt_mki67[sapply(as.numeric(Moffitt_array_bulk["MKI67",]), function(x) x > Moffitt_mki67_thirds[2] && x <= Moffitt_mki67_thirds[3])] <- "MKI67_medium"
Moffitt_mki67[sapply(as.numeric(Moffitt_array_bulk["MKI67",]), function(x) x > Moffitt_mki67_thirds[3])] <- "MKI67_high"
Moffitt_mki67 <- factor(Moffitt_mki67, levels = c("MKI67_low", "MKI67_medium", "MKI67_high"))
Moffitt_annot_colors <- list(tumor_subtype = c(Basal = "#f683ad", Classical = "#f8fc88"),
                         MKI67 = c(MKI67_low = "#77f387",  MKI67_medium = "#f19e5b", MKI67_high = "#d34545"))


baron_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Moffitt_array,
                                                        clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype, 
                                                                                              "MKI67" = Moffitt_mki67,
                                                                                              row.names = rownames(Moffitt_array_meta)),
                                                        fontsize = 11, annotation_colors = Moffitt_annot_colors)
tosti_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Moffitt_array,
                                                        clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype, 
                                                                                              "MKI67" = Moffitt_mki67,
                                                                                              row.names = rownames(Moffitt_array_meta)),
                                                        fontsize = 11, annotation_colors = Moffitt_annot_colors, clustering_method = "average")
## ANOVA: Guo, PAAD
tosti_guo_anova <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                        clinical_characteristic = Guo_meta$description)
ggarrange(tosti_guo_anova$aov_plots[[1]], tosti_guo_anova$aov_plots[[2]], tosti_guo_anova$aov_plots[[3]], nrow = 1, ncol = 3) 
tosti_guo_anova2 <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                         clinical_characteristic = as.character(Guo_mki67))

tosti_paad_anova1 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                         clinical_characteristic = PAAD_meta$tumor_moffitt)
ggarrange(tosti_paad_anova1$aov_plots$racinar, tosti_paad_anova1$aov_plots$ductal, nrow = 1, ncol = 2)
tosti_paad_anova2 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$tumor_bailey)
ggarrange(tosti_paad_anova2$aov_plots$sacinar, tosti_paad_anova2$aov_plots$racinar, tosti_paad_anova2$aov_plots$mductal, nrow = 2, ncol = 2)
tosti_paad_anova3 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$tumor_collisson)
ggarrange(tosti_paad_anova3$aov_plots$sacinar, tosti_paad_anova3$aov_plots$racinar, tosti_paad_anova3$aov_plots$mductal, nrow = 1, ncol = 3)
tosti_paad_anova4 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$neoplasm_histologic_grade)
ggarrange(tosti_paad_anova4$aov_plots$sacinar, tosti_paad_anova4$aov_plots$racinar, tosti_paad_anova4$aov_plots$mductal, nrow = 1, ncol = 3)
tosti_paad_anova5 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = as.character(PAAD_mki67))
ggarrange(tosti_paad_anova5$aov_plots$sacinar, tosti_paad_anova5$aov_plots$mductal, nrow = 1, ncol = 2)


tosti_guo_anova <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                        clinical_characteristic = Guo_meta$description)
ggarrange(tosti_guo_anova$aov_plots[[1]], tosti_guo_anova$aov_plots[[2]], tosti_guo_anova$aov_plots[[3]], nrow = 1, ncol = 3) 
tosti_guo_anova2 <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                         clinical_characteristic = as.character(Guo_mki67))

## survival analysis
Guo_OS <- Guo_meta$Days
Guo_Zensur <- Guo_meta$status.1
Guo_Zensur[Guo_Zensur == 1] <- 0
Guo_Zensur[Guo_Zensur == 2] <- 1
tosti_guo_survival <- survival_analysis(decon_output = decon_tosti$Guo, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                              "MKI67" = as.character(Guo_mki67),
                                                                              row.names = rownames(Guo_meta)))
ggpar(tosti_guo_survival$single_kp$tumor_subtype, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months") # + guides(colour = guide_legend(nrow = 3))

PAAD_OS <- rep(NA, nrow(PAAD_meta))
PAAD_OS[which(is.na(PAAD_meta$days_to_death))] <- PAAD_meta$days_to_last_followup[which(is.na(PAAD_meta$days_to_death))]
PAAD_OS[which(is.na(PAAD_meta$days_to_last_followup))] <- PAAD_meta$days_to_death[which(is.na(PAAD_meta$days_to_last_followup))]
PAAD_OS <- as.numeric(PAAD_OS)
PAAD_censor <- rep(NA, nrow(PAAD_meta))
PAAD_censor[which(PAAD_meta$vital_status == "alive")] <- 0
PAAD_censor[which(PAAD_meta$vital_status == "dead")] <- 1
tosti_PAAD_survival <- survival_analysis(decon_output = decon_tosti$PAAD, 
                                         OS = PAAD_OS, censor = PAAD_censor, 
                                         clinical_characteristics =  data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                "Moffitt" = PAAD_meta$tumor_moffitt,
                                                                                "Bailey" = PAAD_meta$tumor_bailey,
                                                                                "Collisson" = PAAD_meta$tumor_collisson,
                                                                                "MKI67" = as.character(PAAD_mki67),
                                                                                row.names = rownames(PAAD_meta)))
ggpar(tosti_PAAD_survival$single_kp$Moffitt, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")

