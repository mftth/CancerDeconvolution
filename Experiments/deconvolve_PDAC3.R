## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
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
#Flowers_bulk <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_bulk.RDS")
#Flowers_meta <- readRDS("~/Masterthesis/Data/Bulk/Flowers/Flowers_metadata.RDS")
#rownames(Flowers_meta) <- colnames(Flowers_bulk)
#Flowers_meta$COO <- sapply(rownames(Flowers_meta), function(x) strsplit(x, split = "\\.")[[1]][1])
#rownames(Flowers_bulk) <- toupper(rownames(Flowers_bulk))

bulk_list <- list("PAAD" = PAAD_bulk,
                  "Yang" = Yang_bulk,
                  "Guo" = Guo_bulk,
                  #"Janky" = Janky_bulk,
                  #"Kirby" = Kirby_bulk,
                  #"Moffitt_seq" = Moffitt_seq_bulk,
                  "Moffitt_array" = Moffitt_array_bulk)
                  #"Flowers" = Flowers_bulk)
bulk_meta_list <- list("PAAD" = PAAD_meta,
                       "Yang" = Yang_meta,
                       "Guo" = Guo_meta,
                       #"Janky" = Janky_meta,
                       #"Kirby" = Kirby_meta,
                       #"Moffitt_seq" = Moffitt_seq_meta,
                       "Moffitt_array" = Moffitt_array_meta)
                       #"Flowers" = Flowers_meta)

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
technology = c("RNA-seq", "microarray", "RNA-seq", "microarray")
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

pdf(onefile = FALSE, width = 7, height = 5)
ggarrange(baron_pval_boxplot_rmsd, tosti_pval_boxplot_rmsd, common.legend = TRUE) 
dev.off()


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
pdf(width = 4, height = 4)
barplot_proportions(decon_output = decon_baron$Guo,
                    clinical_characteristics_vec =  Guo_meta$description)
dev.off()
tosti_Guo_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Guo,
                                              clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                                    "MKI67" = Guo_mki67,
                                                                                    row.names = rownames(Guo_meta)),
                                              annotation_colors = guo_annot_colors, fontsize = 11)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_tosti$Guo, clinical_characteristic_vec = Guo_meta$description) 
dev.off()
ensemble_Guo_prop_heatmap <- heatmap_proportions(decon_output = guo_ensemble,
                                                 clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                                       "MKI67" = Guo_mki67,
                                                                                       row.names = rownames(Guo_meta)),
                                                 annotation_colors = guo_annot_colors, fontsize = 11,
                                                 clustering_method = "ward.D2")
pdf(width = 5, height = 5)
umap_plot(decon_output = guo_ensemble, clinical_characteristic_vec = Guo_meta$description) 
dev.off()
# pdf(width = 5, height = 5)
# umap_plot(decon_output = guo_ensemble, clinical_characteristic_vec = Guo_mki67) 
# dev.off()
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
pdf(width = 4, height = 4)
barplot_proportions(decon_output = decon_tosti$PAAD,
                    clinical_characteristics_vec = PAAD_meta$tumor_moffitt)
dev.off()

tosti_PAAD_prop_heatmap <-   heatmap_proportions(decon_output = decon_tosti$PAAD,
                                               clinical_characteristics = data.frame("grading" = PAAD_meta$neoplasm_histologic_grade,
                                                                                     "Moffitt" = PAAD_meta$tumor_moffitt,
                                                                                     "Collisson" = PAAD_meta$tumor_collisson,
                                                                                     "Bailey" = PAAD_meta$tumor_bailey,
                                                                                     "MKI67" = PAAD_mki67,
                                                                                     row.names = rownames(PAAD_meta)), 
                                               clustering_method = "ward.D2", annotation_colors = PAAD_annot_colors, fontsize = 11)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_tosti$PAAD, clinical_characteristic_vec = PAAD_meta$tumor_moffitt) 
dev.off()
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_tosti$PAAD, clinical_characteristic_vec = PAAD_meta$tumor_collisson) #+ stat_ellipse()
dev.off()
pdf(width = 6, height = 6)
umap_plot(decon_output = decon_tosti$PAAD, clinical_characteristic_vec = PAAD_meta$tumor_bailey) 
dev.off()


baron_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_baron$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)),
                                               fontsize = 11)
pdf(width = 3, height = 3)
barplot_proportions(decon_output = decon_baron$Yang,
                    clinical_characteristics_vec = Yang_meta$`grading:ch1`)
dev.off()

tosti_yang_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Yang,
                                               clinical_characteristics = data.frame("grading" = Yang_meta$`grading:ch1`, 
                                                                                     #stage = Yang_meta$`Stage:ch1`, 
                                                                                     row.names = rownames(Yang_meta)),
                                               fontsize = 11)
pdf(width = 3, height = 3)
barplot_proportions(decon_output = decon_tosti$Yang,
                    clinical_characteristics_vec = Yang_meta$`grading:ch1`)
dev.off()


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
pdf(width = 3, height = 3)
barplot_proportions(decon_output = decon_baron$Moffitt_array,
                    clinical_characteristics_vec = Moffitt_array_meta$tumor_subtype)
dev.off()

tosti_Moffitt_array_prop_heatmap <- heatmap_proportions(decon_output = decon_tosti$Moffitt_array,
                                                        clinical_characteristics = data.frame("tumor_subtype" = Moffitt_array_meta$tumor_subtype, 
                                                                                              "MKI67" = Moffitt_mki67,
                                                                                              row.names = rownames(Moffitt_array_meta)),
                                                        fontsize = 11, annotation_colors = Moffitt_annot_colors, clustering_method = "average")
pdf(width = 3, height = 3)
barplot_proportions(decon_output = decon_tosti$Moffitt_array,
                    clinical_characteristics_vec = Moffitt_array_meta$tumor_subtype)
dev.off()


## ANOVA: Guo, PAAD
tosti_guo_anova <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                        clinical_characteristic = Guo_meta$description)
#ggarrange(tosti_guo_anova$aov_plots[[1]], tosti_guo_anova$aov_plots[[2]], tosti_guo_anova$aov_plots[[3]], nrow = 1, ncol = 3) 
pdf(width = 3.5, height = 5)
tosti_guo_anova$aov_plots$sacinar + 
  geom_signif(comparisons = tosti_guo_anova$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_guo_anova$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 3.5, height = 5)
tosti_guo_anova$aov_plots$racinar + 
  geom_signif(comparisons = tosti_guo_anova$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_guo_anova$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 3.5, height = 5)
tosti_guo_anova$aov_plots$mductal +  
  geom_signif(comparisons = tosti_guo_anova$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_guo_anova$comparison_list)*0.05),0.05))
dev.off()

tosti_guo_anova2 <- correlation_analysis(decon_output = decon_tosti$Guo, 
                                         clinical_characteristic = as.character(Guo_mki67))


tosti_paad_anova1 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                         clinical_characteristic = PAAD_meta$tumor_moffitt)
ggarrange(tosti_paad_anova1$aov_plots$racinar, tosti_paad_anova1$aov_plots$ductal, nrow = 1, ncol = 2)
pdf(width = 5, height = 5)
tosti_paad_anova1$aov_plots$racinar + 
  geom_signif(comparisons = tosti_paad_anova1$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova1$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
tosti_paad_anova1$aov_plots$ductal + 
  geom_signif(comparisons = tosti_paad_anova1$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova1$comparison_list)*0.05),0.05)) 
dev.off()


tosti_paad_anova2 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$tumor_bailey)
ggarrange(tosti_paad_anova2$aov_plots$sacinar, tosti_paad_anova2$aov_plots$racinar, tosti_paad_anova2$aov_plots$mductal, nrow = 2, ncol = 2)
pdf(width = 5, height = 5)
tosti_paad_anova2$aov_plots$sacinar +
  geom_signif(comparisons = tosti_paad_anova2$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova2$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
tosti_paad_anova2$aov_plots$racinar + 
  geom_signif(comparisons = tosti_paad_anova2$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova2$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
tosti_paad_anova2$aov_plots$mductal + 
  geom_signif(comparisons = tosti_paad_anova2$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova2$comparison_list)*0.05),0.05)) 
dev.off()


tosti_paad_anova3 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$tumor_collisson)
ggarrange(tosti_paad_anova3$aov_plots$sacinar, tosti_paad_anova3$aov_plots$racinar, tosti_paad_anova3$aov_plots$mductal, nrow = 1, ncol = 3)
pdf(width = 3.5, height = 5)
tosti_paad_anova3$aov_plots$sacinar + 
  geom_signif(comparisons = tosti_paad_anova3$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova3$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 3.5, height = 5)
tosti_paad_anova3$aov_plots$racinar + 
  geom_signif(comparisons = tosti_paad_anova3$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova3$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 3.5, height = 5)
tosti_paad_anova3$aov_plots$mductal +  
  geom_signif(comparisons = tosti_paad_anova3$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(tosti_paad_anova3$comparison_list)*0.05),0.05))
dev.off()


tosti_paad_anova4 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = PAAD_meta$neoplasm_histologic_grade)
ggarrange(tosti_paad_anova4$aov_plots$sacinar, tosti_paad_anova4$aov_plots$racinar, tosti_paad_anova4$aov_plots$mductal, nrow = 1, ncol = 3)
tosti_paad_anova5 <- correlation_analysis(decon_output = decon_tosti$PAAD, 
                                          clinical_characteristic = as.character(PAAD_mki67))
ggarrange(tosti_paad_anova5$aov_plots$sacinar, tosti_paad_anova5$aov_plots$mductal, nrow = 1, ncol = 2)


ensemble_guo_anova <- correlation_analysis(decon_output = guo_ensemble, 
                                        clinical_characteristic = Guo_meta$description)
ggarrange(ensemble_guo_anova$aov_plots$acinar, ensemble_guo_anova$aov_plots$ductal) 
pdf(width = 5, height = 5)
ensemble_guo_anova$aov_plots$acinar + 
  geom_signif(comparisons = ensemble_guo_anova$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(ensemble_guo_anova$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
ensemble_guo_anova$aov_plots$ductal + 
  geom_signif(comparisons = ensemble_guo_anova$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(ensemble_guo_anova$comparison_list)*0.05),0.05)) 
dev.off()

ensemble_guo_anova2 <- correlation_analysis(decon_output = guo_ensemble, 
                                         clinical_characteristic = as.character(Guo_mki67))
ggarrange(ensemble_guo_anova2$aov_plots$acinar, ensemble_guo_anova2$aov_plots$ductal) 



## survival analysis
Guo_OS <- Guo_meta$Days
Guo_Zensur <- Guo_meta$status.1
Guo_Zensur[Guo_Zensur == 1] <- 0
Guo_Zensur[Guo_Zensur == 2] <- 1
tosti_guo_survival <- survival_analysis(decon_output = decon_tosti$Guo, OS = Guo_OS, censor = Guo_Zensur, 
                                        clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                              "MKI67" = as.character(Guo_mki67),
                                                                              row.names = rownames(Guo_meta)))
pdf(width = 8, height = 6)
ggpar(tosti_guo_survival$single_kp$tumor_subtype, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months") # + guides(colour = guide_legend(nrow = 3))
dev.off()
ensemble_guo_survival <- survival_analysis(decon_output = guo_ensemble, OS = Guo_OS, censor = Guo_Zensur, 
                                           clinical_characteristics = data.frame("tumor_subtype" = Guo_meta$description,
                                                                              "MKI67" = as.character(Guo_mki67),
                                                                              row.names = rownames(Guo_meta)))

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
pdf(width = 8, height = 6)
ggpar(tosti_PAAD_survival$single_kp$Moffitt, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in days")
dev.off()
pdf(width = 8, height = 6)
ggpar(tosti_PAAD_survival$single_kp$grading, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in days")
dev.off()
pdf(width = 8, height = 6)
ggpar(tosti_PAAD_survival$single_kp$MKI67, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in days")
dev.off()


## ML analysis (Guo-ENSEMBLE, Guo-Tosti, PAAD-Tosti)
## Guo (ensemble + tosti) train on whole dataset (predict PDAC subtype)
## PAAD train on part of it (predict PDAC subtype, and grading)

tosti_guo_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$Guo, 
                                       clinical_char = Guo_meta$description)
tosti_guo_ml_model <- train_ML_model(trainData = tosti_guo_prepped)
ensemble_guo_prepped <- prepare_decon_res(p_value = TRUE, decon_res = guo_ensemble, 
                                          clinical_char = Guo_meta$description)
ensemble_guo_ml_model <- train_ML_model(trainData = ensemble_guo_prepped)
guo_signature_genes <- data.frame("KRAS" = as.numeric(Guo_bulk["KRAS",]), "GATA6" = as.numeric(Guo_bulk["GATA6",]), 
                                   "response" = Guo_meta$description, row.names = rownames(Guo_meta))
guo_signature_genes$response <- factor(guo_signature_genes$response)
guo_baseline_model <- train_ML_model(trainData = guo_signature_genes, feature_selection = FALSE)
guo_baseline_comparison <- boxplot_ML_sd(ml_model_list = 
                                           list("Guo_ENSEMBLE" = ensemble_guo_ml_model$rf_model_whole,
                                           "Guo_Tosti" = tosti_guo_ml_model$rf_model_reduced,
                                           "Guo_baseline" = guo_baseline_model$rf_model_whole),
                                         levels = c("Classical", "Basal", "Hybrid"))
guo_baseline_comparison$boxplots + theme(legend.position="top")
mean(guo_baseline_comparison$boxplots$data$value[guo_baseline_comparison$boxplots$data$Model == "Guo_ENSEMBLE" &
                                        guo_baseline_comparison$boxplots$data$variable == "Sensitivity"])
plot(tosti_guo_ml_model$varimp_whole)
plot(ensemble_guo_ml_model$varimp_whole)
plot.roc(guo_baseline_comparison$ROCcurves$Guo_Tosti, print.auc = TRUE, col = "blue")
plot.roc(guo_baseline_comparison$ROCcurves$Guo_ENSEMBLE, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.4, col = "green", add = TRUE)
plot.roc(guo_baseline_comparison$ROCcurves$Guo_baseline, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.3,col = "red", add = TRUE)
legend(1, 1, legend=c("Guo-Tosti", "Guo-ENSEMBLE", "Guo-baseline"),
       col=c("blue", "green", "red"), lt = 1,cex=0.8)


tosti_PAAD_prepped_grading <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, 
                                                clinical_char = PAAD_meta$neoplasm_histologic_grade)
tosti_PAAD_model_grading <- train_ML_model(trainData = tosti_PAAD_prepped_grading)
PAAD_mki67_prepped <- data.frame("MKI67" = as.numeric(PAAD_bulk["MKI67",]),
                                 "response" = PAAD_meta$neoplasm_histologic_grade, row.names = rownames(PAAD_meta))
PAAD_mki67_prepped$response <- factor(PAAD_mki67_prepped$response)
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "down",
  savePred = TRUE
)
PAAD_mki67_ml_model <- train(x = data.frame("MKI67" = PAAD_mki67_prepped$MKI67, 
                                            row.names = rownames(PAAD_mki67_prepped)), 
                             y = PAAD_mki67_prepped$response, 
                             method = "rf", metric = "Accuracy", trControl = fitControl,
                             type = "Classification", ntree = 500)
PAAD_grading_ml_evaluation <- boxplot_ML_sd(list("PAAD_Tosti" = tosti_PAAD_ml_model_grading$rf_model_reduced ,
                                                 "MKI67_baseline_PAAD" = PAAD_mki67_ml_model),
                                            levels = c("G1", "G2", "G3"))
PAAD_grading_ml_evaluation$boxplots + theme(legend.position="top")
mean(PAAD_grading_ml_evaluation$boxplots$data$value[PAAD_grading_ml_evaluation$boxplots$data$Model == "PAAD_Tosti" &
                                                      PAAD_grading_ml_evaluation$boxplots$data$variable == "Accuracy"])
plot(tosti_PAAD_model_grading$varimp_whole)
plot.roc(PAAD_grading_ml_evaluation$ROCcurves$PAAD_Tosti, print.auc = TRUE, col = "blue")
plot.roc(PAAD_grading_ml_evaluation$ROCcurves$MKI67_baseline_PAAD, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.4,col = "red", add = TRUE)
legend(1, 1, legend=c("PAAD-Tosti", "PAAD-baseline"),
       col=c("blue", "red"), lt = 1,cex=0.8)


tosti_PAAD_prepped_moffitt <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, 
                                                clinical_char = PAAD_meta$tumor_moffitt)
empty_hayashi2 <- which(is.na(tosti_PAAD_prepped_moffitt$response))
tosti_PAAD_prepped_moffitt <- tosti_PAAD_prepped_moffitt[-empty_hayashi2,]
tosti_PAAD_model_moffitt <- train_ML_model(trainData = tosti_PAAD_prepped_moffitt)
PAAD_signature_genes2 <- data.frame("KRAS" = as.numeric(PAAD_bulk["KRAS",]), "GATA6" = as.numeric(PAAD_bulk["GATA6",]), 
                                    "response" = PAAD_meta$tumor_moffitt, row.names = rownames(PAAD_meta))
PAAD_signature_genes2 <- PAAD_signature_genes2[-empty_hayashi2,]
PAAD_signature_genes2$response <- factor(PAAD_signature_genes2$response)
PAAD_signature_genes_model2 <- train_ML_model(trainData = PAAD_signature_genes2, feature_selection = FALSE)
PAAD_baseline_comparison <- boxplot_ML_sd(ml_model_list = 
                                            list("PAAD_Tosti" = tosti_PAAD_model_moffitt$rf_model_reduced,
                                                 "PAAD_baseline" = PAAD_signature_genes_model2$rf_model_whole),
                                          levels = c("Classical", "Basal"))
PAAD_baseline_comparison$boxplots + theme(legend.position="top")
mean(PAAD_baseline_comparison$boxplots$data$value[PAAD_baseline_comparison$boxplots$data$Model == "PAAD_Tosti" &
                                                   PAAD_baseline_comparison$boxplots$data$variable == "Accuracy"])
plot(tosti_PAAD_model_moffitt$varimp_whole)
plot.roc(PAAD_baseline_comparison$ROCcurves$PAAD_Tosti, print.auc = TRUE, col = "blue")
plot.roc(PAAD_baseline_comparison$ROCcurves$PAAD_baseline, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.4,col = "red", add = TRUE)
legend(1, 1, legend=c("PAAD-Tosti", "PAAD-baseline"),
       col=c("blue", "red"), lt = 1,cex=0.8)


tosti_PAAD_prepped_collisson <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, 
                                                clinical_char = PAAD_meta$tumor_collisson)
empty_hayashi3 <- which(is.na(tosti_PAAD_prepped_collisson$response))
tosti_PAAD_prepped_collisson <- tosti_PAAD_prepped_collisson[-empty_hayashi3,]
tosti_PAAD_model_collisson <- train_ML_model(trainData = tosti_PAAD_prepped_collisson)
PAAD_signature_genes3 <- data.frame("KRAS" = as.numeric(PAAD_bulk["KRAS",]), "GATA6" = as.numeric(PAAD_bulk["GATA6",]), 
                                    "response" = PAAD_meta$tumor_collisson, row.names = rownames(PAAD_meta))
PAAD_signature_genes3 <- PAAD_signature_genes2[-empty_hayashi3,]
PAAD_signature_genes3$response <- factor(PAAD_signature_genes3$response)
PAAD_signature_genes_model3 <- train_ML_model(trainData = PAAD_signature_genes3, feature_selection = FALSE)
PAAD_baseline_comparison <- boxplot_ML_sd(ml_model_list = 
                                            list("PAAD_Tosti" = tosti_PAAD_model_collisson$rf_model_whole,
                                                 "PAAD_baseline" = PAAD_signature_genes_model3$rf_model_whole),
                                          levels = c("Exocrine", "Classical", "QM"))
PAAD_baseline_comparison$boxplots + theme(legend.position="top")
mean(PAAD_baseline_comparison$boxplots$data$value[PAAD_baseline_comparison$boxplots$data$Model == "PAAD_Tosti" &
                                                    PAAD_baseline_comparison$boxplots$data$variable == "Accuracy"])
plot(tosti_PAAD_model_collisson$varimp_whole)
plot.roc(PAAD_baseline_comparison$ROCcurves$PAAD_Tosti, print.auc = TRUE, col = "blue")
plot.roc(PAAD_baseline_comparison$ROCcurves$PAAD_baseline, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.4,col = "red", add = TRUE)
legend(1, 1, legend=c("PAAD-Tosti", "PAAD-baseline"),
       col=c("blue", "red"), lt = 1,cex=0.8)


tosti_PAAD_prepped_bailey <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, 
                                                clinical_char = PAAD_meta$tumor_bailey)
empty_hayashi4 <- which(is.na(tosti_PAAD_prepped_bailey$response))
tosti_PAAD_prepped_bailey <- tosti_PAAD_prepped_bailey[-empty_hayashi2,]
tosti_PAAD_model_bailey <- train_ML_model(trainData = tosti_PAAD_prepped_bailey)
PAAD_signature_genes4 <- data.frame("KRAS" = as.numeric(PAAD_bulk["KRAS",]), "GATA6" = as.numeric(PAAD_bulk["GATA6",]), 
                                    "response" = PAAD_meta$tumor_bailey, row.names = rownames(PAAD_meta))
PAAD_signature_genes4 <- PAAD_signature_genes4[-empty_hayashi4,]
PAAD_signature_genes4$response <- factor(PAAD_signature_genes4$response)
PAAD_signature_genes_model4 <- train_ML_model(trainData = PAAD_signature_genes4, feature_selection = FALSE)
PAAD_baseline_comparison <- boxplot_ML_sd(ml_model_list = 
                                            list("Tosti_PAAD" = tosti_PAAD_model_bailey$rf_model_reduced,
                                                 "baseline_PAAD" = PAAD_signature_genes_model4$rf_model_whole),
                                          levels = c("ADEX", "Immunogenic", "Squamous", "Progenitor"))
plot.roc(PAAD_baseline_comparison$ROCcurves$Tosti_PAAD, print.auc = TRUE)
plot.roc(PAAD_baseline_comparison$ROCcurves$baseline_PAAD, print.auc = TRUE)



# tosti_PAAD_prepped_collisson <- prepare_decon_res(p_value = TRUE, decon_res = decon_tosti$PAAD, 
#                                                   clinical_char = PAAD_meta$tumor_collisson)
# empty_hayashi <- which(is.na(tosti_PAAD_prepped_collisson$response))
# tosti_PAAD_prepped_collisson <- tosti_PAAD_prepped_collisson[-empty_hayashi,]
# tosti_PAAD_collisson_trainRowNumbers <- createDataPartition(tosti_PAAD_prepped_collisson$response, p = 0.8, list = FALSE)
# tosti_PAAD_collisson_train <- tosti_PAAD_prepped_collisson[tosti_PAAD_collisson_trainRowNumbers,]
# tosti_PAAD_collisson_test <- tosti_PAAD_prepped_collisson[-tosti_PAAD_collisson_trainRowNumbers,]
# tosti_PAAD_ml_model_collisson <- train_ML_model(trainData = tosti_PAAD_collisson_train)
# tosti_PAAD_ml_pred_collisson <- test_ML_model(train_output = tosti_PAAD_ml_model_collisson, 
#                                               testData = tosti_PAAD_collisson_test[,-ncol(tosti_PAAD_collisson_test)], 
#                                               truth_vec = tosti_PAAD_collisson_test$response)
# tosti_PAAD_collisson_roc_curve <- roc_curve(labels = tosti_PAAD_collisson_test$response, 
#                                             predictions = tosti_PAAD_ml_pred_collisson$predicted_whole, 
#                                             levels = c("Classical", "Exocrine", "QM"))
# PAAD_signature_genes <- data.frame("KRAS" = as.numeric(PAAD_bulk["KRAS",]), "GATA6" = as.numeric(PAAD_bulk["GATA6",]), 
#                                   "response" = PAAD_meta$tumor_collisson, row.names = rownames(PAAD_meta))
# PAAD_signature_genes <- PAAD_signature_genes[-empty_hayashi,]
# PAAD_signature_genes$response <- factor(PAAD_signature_genes$response)
# PAAD_signature_genes_trainRowNumbers <- createDataPartition(PAAD_signature_genes$response, p = 0.8, list = FALSE)
# PAAD_signature_genes_train <- PAAD_signature_genes[PAAD_signature_genes_trainRowNumbers,]
# PAAD_signature_genes_test <- PAAD_signature_genes[-PAAD_signature_genes_trainRowNumbers,]
# PAAD_signature_genes_model <- train_ML_model(trainData = PAAD_signature_genes_train, feature_selection = FALSE)
# PAAD_signature_genes_pred <- test_ML_model(train_output = PAAD_signature_genes_model, 
#                                            testData = PAAD_signature_genes_test[,-ncol(PAAD_signature_genes_test)], 
#                                            truth_vec = PAAD_signature_genes_test$response,
#                                            feature_selection = FALSE)
# PAAD_signature_genes_roc_curve <- roc_curve(labels = PAAD_signature_genes_test$response, 
#                                       predictions = PAAD_signature_genes_pred$predicted_whole, 
#                                       levels = c("Classical", "Exocrine", "QM"))
# PAAD_collisson_ml_evaluation <- boxplot_ML_sd(list("Tosti_PAAD" = tosti_PAAD_ml_model_collisson$rf_model_whole,
#                                                    "baseline_PAAD" = PAAD_signature_genes_model$rf_model_whole))
# PAAD_collisson_ml_evaluation2 <- barplot_ML_evaluation(list("Tosti_PAAD" = tosti_PAAD_ml_pred_collisson$evaluation_whole,
#                                                             "baseline_PAAD" = PAAD_signature_genes_pred$evaluation_whole))