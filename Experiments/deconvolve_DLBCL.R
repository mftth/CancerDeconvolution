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
## Chapuy: 137 samples;  OS, OS_stat; Cluster; COO; IPI
Chapuy_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_bulk.RDS")
Chapuy_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Chapuy/Chapuy_metadata.RDS")
## Schleich: 109 samples; treatment response
Schleich_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_bulk.RDS")
Schleich_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schleich/Schleich_metadata.RDS")
rownames(Schleich_bulk) <- toupper(rownames(Schleich_bulk))
## Schmitz 481 samples; survival; COO; treatment, IPI
Schmitz_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_bulk.RDS")
Schmitz_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Schmitz_counts/Schmitz_metadata.RDS")
## Reddy samples 773; survival, COO, treatment response, IPI
Reddy_bulk <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Reddy_counts/Reddy_bulk.RDS")
Reddy_meta <- readRDS("~/Masterthesis/Data/Bulk/DLBCL/Reddy_counts/Reddy_metadata.RDS")


bulk_list <- list("Chapuy" = Chapuy_bulk,
                  "Schleich" = Schleich_bulk,
                  "Schmitz" = Schmitz_bulk,
                  "Reddy" = Reddy_bulk)
bulk_meta_list <- list("Chapuy" = Chapuy_meta,
                       "Schleich" = Schleich_meta,
                       "Schmitz" = Schmitz_meta,
                       "Reddy" = Reddy_meta)

## SeneSys scRNA-seq data
# mouse_senescence <- read.table("~/SeneSys_scRNA_mouse/senescence_mouse.tsv", sep = "\t", header = TRUE, row.names = 1)
# mouse_senescence <- as.matrix(mouse_senescence)
# rownames(mouse_senescence) <- toupper(rownames(mouse_senescence))
# phenotypes <- sapply(colnames(mouse_senescence), function(x) strsplit(x, "_\\d+")[[1]][1])
# phenotypes <- as.data.frame(phenotypes)
# phenotypes$sample <- rep("donor", nrow(phenotypes))
# colnames(phenotypes)[1] <- "cluster"
# qc_mouse_senescence <- Quality_control(sc_data = mouse_senescence, sc_meta = phenotypes, 
#                                        sc_path = "~/Masterthesis/CancerDeconvolution/Data/SingleCell/qc_mouse_senescence.RDS",
#                                        cell_types = unique(phenotypes$cluster),
#                                        multiple_donors = FALSE)

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
  saveRDS(decon_dlbcl_x, file = paste(res_path, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})
decon_dlbcl <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_dlbcl) <- names(bulk_list)

schmitz_null_statistics <- sapply(decon_dlbcl$Schmitz$statistics_sampled, function(x) is.null(x))
decon_dlbcl$Schmitz$statistics_sampled <- decon_dlbcl$Schmitz$statistics_sampled[!schmitz_null_statistics]

pearson_matrix_sampled <- sapply(decon_dlbcl$Schmitz$statistics_sampled, function(x) x$pearson_vec)
spearman_matrix_sampled <- sapply(decon_dlbcl$Schmitz$statistics_sampled, function(x) x$spearman_vec)
mad_matrix_sampled <- sapply(decon_dlbcl$Schmitz$statistics_sampled, function(x) x$mad_vec)
rmsd_matrix_sampled <- sapply(decon_dlbcl$Schmitz$statistics_sampled, function(x) x$rmsd_vec)

statistics_obs <- decon_dlbcl$Schmitz$statistics_observed
p_value_wy_pearson <- (sum(colMedians(abs(pearson_matrix_sampled)) >= median(abs(statistics_obs$pearson_vec)))+1)/(ncol(pearson_matrix_sampled)+1)
p_value_wy_spearman <- (sum(colMedians(abs(spearman_matrix_sampled)) >= median(abs(statistics_obs$spearman_vec)))+1)/(ncol(spearman_matrix_sampled)+1)
p_value_wy_mad <- (sum(colMedians(abs(mad_matrix_sampled)) <= median(abs(statistics_obs$mad_vec)))+1)/(ncol(mad_matrix_sampled)+1)
p_value_wy_rmsd <- (sum(colMedians(abs(rmsd_matrix_sampled)) <= median(abs(statistics_obs$rmsd_vec)))+1)/(ncol(rmsd_matrix_sampled)+1)

p_value_wy_pearson_per_sample <- sapply(1:nrow(pearson_matrix_sampled), 
                                        function(x) (sum(abs(pearson_matrix_sampled[x,]) >= abs(statistics_obs$pearson_vec[x]))+1)/(ncol(pearson_matrix_sampled)+1))
p_value_wy_spearman_per_sample <- sapply(1:nrow(spearman_matrix_sampled), 
                                         function(x) (sum(abs(spearman_matrix_sampled[x,]) >= abs(statistics_obs$spearman_vec[x]))+1)/(ncol(spearman_matrix_sampled)+1))
p_value_wy_mad_per_sample <- sapply(1:nrow(mad_matrix_sampled), 
                                    function(x) (sum(abs(mad_matrix_sampled[x,]) <= abs(statistics_obs$mad_vec[x]))+1)/(ncol(mad_matrix_sampled)+1))
p_value_wy_rmsd_per_sample <- sapply(1:nrow(rmsd_matrix_sampled), 
                                     function(x) (sum(abs(rmsd_matrix_sampled[x,]) <= abs(statistics_obs$rmsd_vec[x]))+1)/(ncol(rmsd_matrix_sampled)+1))

p_value_per_sample <- data.frame(Pearson = p_value_wy_pearson_per_sample,
                                 Spearman = p_value_wy_spearman_per_sample,
                                 mAD = p_value_wy_mad_per_sample,
                                 RMSD = p_value_wy_rmsd_per_sample,
                                 row.names = rownames(decon_dlbcl$Schmitz$decon_res$prop.est.mvw))
decon_dlbcl$Schmitz$p_value_per_sample <- p_value_per_sample
decon_dlbcl <- decon_dlbcl[1:3]

## plot p-values
technology = c("microarray", "microarray", "RNA-seq")#, "RNA-seq")
pval_boxplot_spearman <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                        pvalue_type = "Spearman", technology = technology) 
pval_boxplot_pearson <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                       pvalue_type = "Pearson", technology = technology) + 
  theme(plot.title = element_text(size=12)) +  theme(legend.position="top")
pval_boxplot_mad <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                   pvalue_type = "mAD", technology = technology) 
pval_boxplot_rmsd <- boxplot_pvalue(decon_output_list = decon_dlbcl,
                                    pvalue_type = "RMSD", technology = technology) 


## import ecotyper results
ecotyper_results_path <- list.files("~/Ecotyper/Results2", full.names = TRUE)
names(ecotyper_results_path) <- c("Chapuy", "Reddy", "Schleich", "Schmitz")
ecotyper <- lapply(ecotyper_results_path, function(x) {
  ecotyper_result <- unzip(x)
  ecotyper_bcell_assignment_file <- ecotyper_result[grep("B.cells_Cell_State_Assignment", ecotyper_result)]
  return(read.table(ecotyper_bcell_assignment_file, header = TRUE))
})
Chapuy_meta$ecotyper_bcell_state <- rep(NA, nrow(Chapuy_meta))
Chapuy_meta$ecotyper_bcell_state <- ecotyper$Chapuy$Cell.State[match(Chapuy_meta$geo_accession, ecotyper$Chapuy$ID)]
Schleich_meta$ecotyper_bcell_state <- rep(NA, nrow(Schleich_meta))
Schleich_meta$ecotyper_bcell_state <- ecotyper$Schleich$Cell.State[match(Schleich_meta$geo_accession, ecotyper$Schleich$ID)]
Schmitz_meta$ecotyper_bcell_state <- rep(NA, nrow(Schmitz_meta))
Schmitz_meta$ecotyper_bcell_state <- ecotyper$Schmitz$Cell.State[match(Schmitz_meta$sample, ecotyper$Schmitz$ID)]
#Reddy_meta$ecotyper_bcell_state <- rep(NA, nrow(Reddy_meta))
#Reddy_meta$ecotyper_bcell_state <- ecotyper$Reddy$Cell.State[match(Reddy_meta$X, ecotyper$Reddy$ID)]


## visualize ct props in heatmaps
## muss dafür na und ? entfernen
## plotte mit den gaps, kb die ueberall rauszunehmen
## for chapuy: coo und cluster, ecotyper bcell states
## for schleich: treatment and treatment response, ecotyper bcell states
## for schmitz: coo, treatment and IPI group, ecotyper bcell states
## for reddy: coo, treatment response, ipi, ecotyper bcell states

Chapuy_meta$`R-CHOP-like Chemo`[Chapuy_meta$`R-CHOP-like Chemo` == "na"] <- NA
#Chapuy_meta$Cluster <- sapply(Chapuy_meta$Cluster, function(x) paste0("C", x, collapse = ""))
Chapuy_annot_colors <- list(COO = c(ABC = "#6eacf2", GCB = "#f8fc88", Unclassified ="#f683ad"),
                            #cluster = c(C0 = "#ffed69", C1 = "#69ff6b", C2 = "#69f0ff",
                            #            C3 = "#8269ff", C4 = "#ff69dc", C5 = "#ff6969"),
                            treatment = c(yes = "#c4bcbc", no = "#000000"),
                            IPI = c(Low = "#77f387", Intermediate ="#f19e5b", High = "#d34545"),
                            Bcell.state = c(S01 ="#8269ff",  S02 = "#ff69dc", S03 =  "#69f0ff",
                                            S04 = "#69ff6b", S05 = "#ff6969"))
Chapuy_meta$IPI[Chapuy_meta$IPI == "na"] <- NA
Chapuy_meta$IPI[Chapuy_meta$IPI == "0" | Chapuy_meta$IPI == "1"] <- "Low"
Chapuy_meta$IPI[Chapuy_meta$IPI == "2" | Chapuy_meta$IPI == "3"] <- "Intermediate"
Chapuy_meta$IPI[Chapuy_meta$IPI == "4" | Chapuy_meta$IPI == "5"] <- "High"
Chapuy_meta$IPI <- factor(Chapuy_meta$IPI, levels = c("Low", "Intermediate", "High"))
Chapuy_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Chapuy,
                                           clinical_characteristics = data.frame("COO" = Chapuy_meta$COO_byGEP, 
                                                                                 #"cluster" = Chapuy_meta$Cluster,
                                                                                 "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                                 "IPI" = Chapuy_meta$IPI,
                                                                                 "Bcell state" = Chapuy_meta$ecotyper_bcell_state,
                                                                                 row.names = rownames(Chapuy_meta)),
                                           fontsize_rows = 11, angle_col = 45, annotation_colors = Chapuy_annot_colors)


Schleich_annot_colors <- list(treatment_response = c(NR = "#6eacf2", RES = "#f8fc88", RP ="#f683ad"),
                              treatment = c(CTX = "#c4bcbc", native = "#000000"),
                              Bcell.state = c(S01 ="#8269ff",  S02 = "#ff69dc", S03 =  "#69f0ff",
                                              S04 = "#69ff6b", S05 = "#ff6969"))
Schleich_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Schleich,
                                             clinical_characteristics = data.frame("t.response" = Schleich_meta$treatment_response,
                                                                                   "treatment" = Schleich_meta$treatment,
                                                                                   "Bcell state" = Schleich_meta$ecotyper_bcell_state,
                                                                                   row.names = rownames(Schleich_meta)),
                                             fontsize = 11, angle_col = 45, annotation_colors = Schleich_annot_colors, clustering_method = "average")


Schmitz_annot_colors <- list(COO = c(ABC = "#6eacf2", GCB = "#f8fc88", Unclass ="#f683ad"),
                             IPI = c(Low = "#77f387", Intermediate ="#f19e5b", High = "#d34545"),
                             Bcell.state = c(S01 ="#8269ff",  S02 = "#ff69dc", S03 =  "#69f0ff",
                                             S04 = "#69ff6b", S05 = "#ff6969"))
Schmitz_meta$IPI.Group[Schmitz_meta$IPI.Group == ""] <- NA
Schmitz_meta$IPI.Group <- factor(Schmitz_meta$IPI.Group, levels = c("Low", "Intermediate", "High"))
Schmitz_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Schmitz,
                                            clinical_characteristics = data.frame("COO" = Schmitz_meta$Gene.Expression.Subgroup ,
                                                                                  "IPI" = Schmitz_meta$IPI.Group,
                                                                                  "Bcell state" = Schmitz_meta$ecotyper_bcell_state,
                                                                                   row.names = rownames(Schmitz_meta)),
                                             fontsize = 11, angle_col = 45, annotation_colors = Schmitz_annot_colors, clustering_method = "ward.D2")


# Reddy_meta$IPI_group <- rep(NA, nrow(Reddy_meta))
# Reddy_meta$IPI_group[Reddy_meta$IPI == "0" | Reddy_meta$IPI == "1"] <- "Low"
# Reddy_meta$IPI_group[Reddy_meta$IPI == "2" | Reddy_meta$IPI == "3"] <- "Intermediate"
# Reddy_meta$IPI_group[Reddy_meta$IPI == "4" | Reddy_meta$IPI == "5"] <- "High"
# Reddy_meta$IPI_group <- factor(Reddy_meta$IPI_group, levels = c("Low", "Intermediate", "High"))
# Reddy_meta$Response.to.initial.therapy[Reddy_meta$Response.to.initial.therapy == ""] <- NA
# #Reddy_meta$treatment_response <- rep(NA, nrow(Reddy_meta))
# Reddy_meta$Overall.Survival.years <- as.numeric(gsub(",", ".", Reddy_meta$Overall.Survival.years))
# Reddy_annot_colors <- list(COO = c(ABC = "#6eacf2", GCB = "#f8fc88", Unclassified ="#f683ad"),
#                            Treatment.response = c("Complete response" = "#6eacf2", "No response" = "#f8fc88", 
#                                                   "Partial response" ="#f683ad"),
#                            IPI = c(Low = "#77f387", Intermediate ="#f19e5b", High = "#d34545"),
#                            Bcell.state = c(S01 ="#8269ff",  S02 = "#ff69dc", S03 =  "#69f0ff",
#                                            S04 = "#69ff6b", S05 = "#ff6969"))
# Reddy_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Reddy,
#                                             clinical_characteristics = data.frame("COO" = Reddy_meta$ABC.GCB..RNAseq.,
#                                                                                   "IPI" = Reddy_meta$IPI_group,
#                                                                                   "Treatment response" = Reddy_meta$Response.to.initial.therapy,
#                                                                                   "Bcell state" = Reddy_meta$ecotyper_bcell_state,
#                                                                                   row.names = rownames(Reddy_meta)),
#                                             fontsize = 11, annotation_colors = Reddy_annot_colors)#, clustering_method = "ward.D2")



## ANOVA
chapuy_anova_coo <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                         clinical_characteristic = Chapuy_meta$COO_byGEP)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Chapuy,
                             clinical_characteristic_vec = Chapuy_meta$COO_byGEP)
dev.off()
#chapuy_anova_cluster <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
#                                             clinical_characteristic = Chapuy_meta$Cluster)
#chapuy_cluster_umap <- umap_plot(decon_output = decon_dlbcl$Chapuy,
#                             clinical_characteristic_vec = Chapuy_meta$Cluster)
chapuy_anova_treatment <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                               clinical_characteristic = Chapuy_meta$`R-CHOP-like Chemo`)
chapuy_treatment_umap <- umap_plot(decon_output = decon_dlbcl$Chapuy,
                                   clinical_characteristic_vec = Chapuy_meta$`R-CHOP-like Chemo`)
chapuy_anova_ipi <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                         clinical_characteristic = as.character(Chapuy_meta$IPI))
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Chapuy,
                             clinical_characteristic_vec = Chapuy_meta$IPI)
dev.off()
chapuy_anova_bcellstate <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                                clinical_characteristic = as.character(Chapuy_meta$ecotyper_bcell_state))
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Chapuy,
                                    clinical_characteristic_vec = Chapuy_meta$ecotyper_bcell_state)
dev.off()


schleich_anova_treatment <- correlation_analysis(decon_output = decon_dlbcl$Schleich, 
                                                 clinical_characteristic = Schleich_meta$treatment)
schleich_treatment_umap <- umap_plot(decon_output = decon_dlbcl$Schleich,
                                     clinical_characteristic_vec = Schleich_meta$treatment) 
schleich_anova_treatmentoutcome <- correlation_analysis(decon_output = decon_dlbcl$Schleich, 
                                                        clinical_characteristic = Schleich_meta$treatment_response)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Schleich,
                                            clinical_characteristic_vec = Schleich_meta$treatment_response) 
dev.off()
schleich_anova_bcellstate <- correlation_analysis(decon_output = decon_dlbcl$Schleich, 
                                                  clinical_characteristic = Schleich_meta$ecotyper_bcell_state)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Schleich,
                                      clinical_characteristic_vec = Schleich_meta$ecotyper_bcell_state) 
dev.off()
pdf(width = 5, height = 5)
schleich_anova_treatmentoutcome$aov_plots[[1]] + 
  geom_signif(comparisons = schleich_anova_treatmentoutcome$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schleich_anova_treatmentoutcome$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
schleich_anova_treatmentoutcome$aov_plots[[2]] +
  geom_signif(comparisons = schleich_anova_treatmentoutcome$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schleich_anova_treatmentoutcome$comparison_list)*0.05),0.05)) 
dev.off()


schmitz_anova_coo <- correlation_analysis(decon_output = decon_dlbcl$Schmitz, 
                                          clinical_characteristic = Schmitz_meta$Gene.Expression.Subgroup)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Schmitz,
                              clinical_characteristic_vec = Schmitz_meta$Gene.Expression.Subgroup) 
dev.off()
schmitz_anova_ipi <- correlation_analysis(decon_output = decon_dlbcl$Schmitz, 
                                          clinical_characteristic = as.character(Schmitz_meta$IPI.Group))
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Schmitz,
                              clinical_characteristic_vec = Schmitz_meta$IPI.Group)
dev.off()
schmitz_anova_bcellstate <- correlation_analysis(decon_output = decon_dlbcl$Schmitz, 
                                                 clinical_characteristic = Schmitz_meta$ecotyper_bcell_state)
pdf(width = 5, height = 5)
umap_plot(decon_output = decon_dlbcl$Schmitz,
                                     clinical_characteristic_vec = Schmitz_meta$ecotyper_bcell_state)
dev.off()
pdf(width = 5, height = 5)
schmitz_anova_ipi$aov_plots$ADR_OHT + 
  geom_signif(comparisons = schmitz_anova_ipi$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schmitz_anova_ipi$comparison_list)*0.05),0.05)) 
dev.off()
pdf(width = 5, height = 5)
schmitz_anova_bcellstate$aov_plots$ADR_OHT + 
  geom_signif(comparisons = schmitz_anova_bcellstate$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schmitz_anova_bcellstate$comparison_list)*0.08),0.08)) 
dev.off()
pdf(width = 5, height = 5)
schmitz_anova_bcellstate$aov_plots$ADR + 
  geom_signif(comparisons = schmitz_anova_bcellstate$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schmitz_anova_bcellstate$comparison_list)*0.08),0.08)) 
dev.off()
pdf(width = 5, height = 5)
schmitz_anova_bcellstate$aov_plots$PS + 
  geom_signif(comparisons = schmitz_anova_bcellstate$comparison_list, 
              map_signif_level=TRUE, tip_length = 0,
              y_position = seq(1,(1+length(schmitz_anova_bcellstate$comparison_list)*0.08),0.08)) 
dev.off()
#ggarrange(schmitz_anova_bcellstate$aov_plots[[1]], #so2 mit jedem und s1-s3, s3-s4, s3-s5
#          schmitz_anova_bcellstate$aov_plots[[2]], #so2 mit jedem und s1-s3, s3-s4, s3-s5
#          schmitz_anova_bcellstate$aov_plots[[3]], nrow =  3) # s2-s3, s2-s4, s1-s3, s3-s5


## survival analysis
chapuy_OS <- Chapuy_meta$OS
chapuy_OS[chapuy_OS == "na"] <- NA
chapuy_OS <- as.numeric(chapuy_OS)
chapuy_os_zensur <- Chapuy_meta$OS_STAT
chapuy_os_zensur[chapuy_os_zensur == "na"] <- NA
chapuy_os_zensur <- as.numeric(chapuy_os_zensur)
chapuy_PFS <- Chapuy_meta$PFS
chapuy_PFS[chapuy_PFS == "na"] <- NA
chapuy_PFS <- as.numeric(chapuy_PFS)
chapuy_pfs_zensur <- Chapuy_meta$PFS_STAT
chapuy_pfs_zensur[chapuy_pfs_zensur == "na"] <- NA
chapuy_pfs_zensur <- as.numeric(chapuy_pfs_zensur)
chapuy_os_survival <- survival_analysis(decon_output = decon_dlbcl$Chapuy, OS = chapuy_OS, censor = chapuy_os_zensur, 
                                        clinical_characteristics = data.frame("COO" = Chapuy_meta$COO_byGEP, 
                                                                              #"cluster" = Chapuy_meta$Cluster,
                                                                              "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                              "IPI" = Chapuy_meta$IPI,
                                                                              "Bcell state" = Chapuy_meta$ecotyper_bcell_state,
                                                                              row.names = rownames(Chapuy_meta)))
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_os_survival$single_kp$ADR_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_os_survival$single_kp$ADR_OHT_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_os_survival$single_kp$IPI, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()

chapuy_pfs_survival <- survival_analysis(decon_output = decon_dlbcl$Chapuy, OS = chapuy_PFS, censor = chapuy_pfs_zensur, 
                                         clinical_characteristics = data.frame("COO" = Chapuy_meta$COO_byGEP, 
                                                                               #"cluster" = Chapuy_meta$Cluster,
                                                                               "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                               "IPI" = Chapuy_meta$IPI,
                                                                               "Bcell state" = Chapuy_meta$ecotyper_bcell_state,
                                                                               row.names = rownames(Chapuy_meta)))
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_pfs_survival$single_kp$ADR_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_pfs_survival$single_kp$ADR_OHT_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(chapuy_pfs_survival$single_kp$IPI, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()


schmitz_PFS <- Schmitz_meta$`Progression_Free Survival _PFS_ Time _yrs`
schmitz_pfs_zensur <- Schmitz_meta$`Progression_Free Survival _PFS_ Status_ 0 No Progressoin_ 1 Progression`
schmitz_pfs_survival <- survival_analysis(decon_output = decon_dlbcl$Schmitz, OS = schmitz_PFS, censor = schmitz_pfs_zensur,
                                          clinical_characteristics = data.frame("COO" = Schmitz_meta$Gene.Expression.Subgroup ,
                                                                                "IPI" = Schmitz_meta$IPI.Group,
                                                                                "Bcell state" = Schmitz_meta$ecotyper_bcell_state,
                                                                                row.names = rownames(Schmitz_meta)))
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(schmitz_pfs_survival$single_kp$ADR_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(schmitz_pfs_survival$single_kp$ADR_OHT_high_low , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(schmitz_pfs_survival$single_kp$IPI, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(schmitz_pfs_survival$single_kp$COO , 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()
pdf(onefile = FALSE, width = 8, height = 6.08)
ggpar(schmitz_pfs_survival$single_kp$Bcell.state, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(11),font.tickslab = c(12), 
      xlab = "Time in months")
dev.off()


# reddy_OS <- Reddy_meta$Overall.Survival.years
# reddy_censored <- Reddy_meta$Censored
# reddy_os_survival <- survival_analysis(decon_output = decon_dlbcl$Reddy, OS = reddy_OS, censor = reddy_censored,
#                                        clinical_characteristics = data.frame("COO" = Reddy_meta$ABC.GCB..RNAseq.,
#                                                                              "IPI" = Reddy_meta$IPI_group,
#                                                                              "Treatment response" = Reddy_meta$Response.to.initial.therapy,
#                                                                              "Bcell state" = Reddy_meta$ecotyper_bcell_state,
#                                                                              row.names = rownames(Reddy_meta)))
## 


### ML
## train a model on schleich for prediction of treatment response
## baseline model: survarness genes
## apply to schmitz and chapuy and see if survival was correctly predicted (1 nr, 2 rp, 3 res)

## load suvarness gene set
library(GSA)
senesys_gmt <- GSA.read.gmt(filename = "~/SeneSys/Chapuy_Enrichment/SeneSys_gene_sets.gmt.tsv")
names(senesys_gmt$genesets) <- senesys_gmt$geneset.names
senesys_gmt$genesets <- lapply(senesys_gmt$genesets, function(x) {
  x <- x[which(x != "")]
})
suvarness <- senesys_gmt$genesets$SUVARNESS
suvarness_reduced <- intersect(suvarness, rownames(Schleich_bulk))
suvarness_reduced <- intersect(suvarness_reduced, rownames(Schmitz_bulk))
suvarness_reduced <- intersect(suvarness_reduced, rownames(Chapuy_bulk))


schleich_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_dlbcl$Schleich, 
                                      clinical_char = Schleich_meta$treatment_response)
schleich_ml_model <- train_ML_model(trainData = schleich_prepped)
suvarness_prepped <- data.frame(t(Schleich_bulk[suvarness_reduced,]), "response" = Schleich_meta$treatment_response, 
                                row.names = rownames(Schleich_meta))
suvarness_prepped$response <- factor(suvarness_prepped$response)
suvarness_baseline_model <- train_ML_model(trainData = suvarness_prepped)
schleich_baseline_comparison <- boxplot_ML_sd(ml_model_list = 
                                              list("Schleich" = schleich_ml_model$rf_model_whole,
                                                   "SUVARness_baseline" = suvarness_baseline_model$rf_model_whole),
                                         levels = c("NR", "RP", "RES"))
schleich_baseline_comparison$boxplots + theme(legend.position="top")

mean(schleich_baseline_comparison$boxplots$data$value[schleich_baseline_comparison$boxplots$data$Model == "Schleich" &
                                                      schleich_baseline_comparison$boxplots$data$variable == "Sensitivity"])
plot(schleich_ml_model$varimp_whole)
plot(suvarness_baseline_model$varimp_whole)
plot.roc(schleich_baseline_comparison$ROCcurves$Schleich, print.auc = TRUE, col = "red")
plot.roc(schleich_baseline_comparison$ROCcurves$SUVARness_baseline, print.auc = TRUE,
         print.auc.x = 0.5, print.auc.y = 0.4,col = "blue", add = TRUE)
legend(0.5, 0.2, legend=c("Schleich", "SUVARness-baseline"),
       col=c("red", "blue"), lt = 1,cex=0.8)

## apply both models to chapuy and schmitz
chapuy_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_dlbcl$Chapuy, 
                                    clinical_char = Chapuy_meta$IPI)
chapuy_treatment_response_prediction1 <- predict(schleich_ml_model$rf_model_whole, 
                                                 chapuy_prepped[,-ncol(chapuy_prepped)])
chapuy_prepped2 <- data.frame(t(Chapuy_bulk[suvarness,]), row.names = rownames(Chapuy_meta))
chapuy_treatment_response_prediction2 <- predict(suvarness_baseline_model$rf_model_whole, 
                                                  chapuy_prepped2)

schmitz_prepped <- prepare_decon_res(p_value = TRUE, decon_res = decon_dlbcl$Schmitz, 
                                    clinical_char = Schmitz_meta$IPI.Group)
schmitz_treatment_response_prediction1 <- predict(schleich_ml_model$rf_model_whole, 
                                                 schmitz_prepped[,-ncol(schmitz_prepped)])
schmitz_prepped2 <- data.frame(t(Schmitz_bulk[suvarness,]), row.names = rownames(Schmitz_meta))
schmitz_treatment_response_prediction2 <- predict(suvarness_baseline_model$rf_model_whole, 
                                                  schmitz_prepped2)

Chapuy_meta$pred_treatment_response <- chapuy_treatment_response_prediction1
Chapuy_meta$pred_treatment_response2 <- chapuy_treatment_response_prediction2
Schmitz_meta$pred_treatment_response <- schmitz_treatment_response_prediction1


## repeat survival analysis
chapuy_os_survival2 <- survival_analysis(decon_output = decon_dlbcl$Chapuy, OS = chapuy_OS, censor = chapuy_os_zensur, 
                                        clinical_characteristics = data.frame("response_pred_schleich" = Chapuy_meta$pred_treatment_response, 
                                                                              "response_pred_baseline" = Chapuy_meta$pred_treatment_response2, 
                                                                              row.names = rownames(Chapuy_meta)))
ggpar(chapuy_os_survival2$single_kp$response_pred_schleich, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
ggpar(chapuy_os_survival2$single_kp$response_pred_baseline, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")

chapuy_pfs_survival2 <- survival_analysis(decon_output = decon_dlbcl$Chapuy, OS = chapuy_PFS, censor = chapuy_pfs_zensur, 
                                         clinical_characteristics = data.frame("response_pred_schleich" = Chapuy_meta$pred_treatment_response, 
                                                                               "response_pred_baseline" = Chapuy_meta$pred_treatment_response2,
                                                                               row.names = rownames(Chapuy_meta)))
ggpar(chapuy_pfs_survival2$single_kp$response_pred_schleich, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")
ggpar(chapuy_pfs_survival2$single_kp$response_pred_baseline, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months")

schmitz_pfs_survival2 <- survival_analysis(decon_output = decon_dlbcl$Schmitz, OS = schmitz_PFS, censor = schmitz_pfs_zensur,
                                          clinical_characteristics = data.frame("response_pred_schleich" = Schmitz_meta$pred_treatment_response,
                                                                                row.names = rownames(Schmitz_meta)))
ggpar(schmitz_pfs_survival2$single_kp$response_pred_schleich, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in years")


umap_plot(decon_output = decon_dlbcl$Schmitz, clinical_characteristic_vec = Schmitz_meta$pred_treatment_response)# + 
  stat_ellipse()
umap_plot(decon_output = decon_dlbcl$Chapuy, clinical_characteristic_vec = Chapuy_meta$pred_treatment_response) + 
  stat_ellipse()
umap_plot(decon_output = decon_dlbcl$Chapuy, clinical_characteristic_vec = Chapuy_meta$pred_treatment_response2) 