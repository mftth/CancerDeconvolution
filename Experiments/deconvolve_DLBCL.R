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
  saveRDS(decon_dlbcl_x, file = paste(res_path, "/", names(bulk_list)[x], "_decon.RDS", sep = ""))
})
decon_dlbcl <- lapply(1:length(bulk_list), 
                      function(x) readRDS(file = paste(res_path, "/", names(bulk_list)[x], 
                                                       "_decon.RDS", sep = "")))
names(decon_dlbcl) <- names(bulk_list)

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
## muss dafür na und ? entfernen
## plotte mit den gaps, kb die ueberall rauszunehmen
## for chapuy: coo und cluster
## for schleich: treatment and treatment response
## for schmitz: coo, treatment and IPI group

Chapuy_meta$`R-CHOP-like Chemo`[Chapuy_meta$`R-CHOP-like Chemo` == "na"] <- NA
Chapuy_meta$Cluster <- sapply(Chapuy_meta$Cluster, function(x) paste0("C", x, collapse = ""))
Chapuy_annot_colors <- list(COO = c(ABC = "#6eacf2", GCB = "#f8fc88", Unclassified ="#f683ad"),
                            cluster = c(C0 = "#ffed69", C1 = "#69ff6b", C2 = "#69f0ff",
                                        C3 = "#8269ff", C4 = "#ff69dc", C5 = "#ff6969"),
                            treatment = c(yes = "#c4bcbc", no = "#000000"),
                            IPI = c(Low = "#77f387", Intermediate ="#f19e5b", High = "#d34545"))
Chapuy_meta$IPI[Chapuy_meta$IPI == "na"] <- NA
Chapuy_meta$IPI[Chapuy_meta$IPI == "0" | Chapuy_meta$IPI == "1"] <- "Low"
Chapuy_meta$IPI[Chapuy_meta$IPI == "2" | Chapuy_meta$IPI == "3"] <- "Intermediate"
Chapuy_meta$IPI[Chapuy_meta$IPI == "4" | Chapuy_meta$IPI == "5"] <- "High"
Chapuy_meta$IPI <- factor(Chapuy_meta$IPI, levels = c("Low", "Intermediate", "High"))
Chapuy_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Chapuy,
                                           clinical_characteristics = data.frame("COO" = Chapuy_meta$COO_byGEP, 
                                                                                 "cluster" = Chapuy_meta$Cluster,
                                                                                 "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                                 "IPI" = Chapuy_meta$IPI,
                                                                                 row.names = rownames(Chapuy_meta)),
                                           fontsize = 11, annotation_colors = Chapuy_annot_colors)


Schleich_annot_colors <- list(treatment_response = c(NR = "#6eacf2", RES = "#f8fc88", RP ="#f683ad"),
                              treatment = c(CTX = "#c4bcbc", native = "#000000"))
Schleich_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Schleich,
                                             clinical_characteristics = data.frame("treatment_response" = Schleich_meta$treatment_response,
                                                                                   "treatment" = Schleich_meta$treatment,
                                                                                   row.names = rownames(Schleich_meta)),
                                             fontsize = 11, annotation_colors = Schleich_annot_colors, clustering_method = "average")


Schmitz_annot_colors <- list(COO = c(ABC = "#6eacf2", GCB = "#f8fc88", Unclass ="#f683ad"),
                             IPI = c(Low = "#77f387", Intermediate ="#f19e5b", High = "#d34545"))
Schmitz_meta$IPI.Group[Schmitz_meta$IPI.Group == ""] <- NA
Schmitz_meta$IPI.Group <- factor(Schmitz_meta$IPI.Group, levels = c("Low", "Intermediate", "High"))
Schmitz_prop_heatmap <- heatmap_proportions(decon_output = decon_dlbcl$Schmitz,
                                            clinical_characteristics = data.frame("COO" = Schmitz_meta$Gene.Expression.Subgroup ,
                                                                                  "IPI" = Schmitz_meta$IPI.Group,
                                                                                   row.names = rownames(Schmitz_meta)),
                                             fontsize = 11, annotation_colors = Schmitz_annot_colors)


## ANOVA
chapuy_anova_coo <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                         clinical_characteristic = Chapuy_meta$COO_byGEP)
chapuy_anova_cluster <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                             clinical_characteristic = Chapuy_meta$Cluster)
chapuy_anova_treatment <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                               clinical_characteristic = Chapuy_meta$`R-CHOP-like Chemo`)
chapuy_anova_ipi <- correlation_analysis(decon_output = decon_dlbcl$Chapuy, 
                                         clinical_characteristic = as.character(Chapuy_meta$IPI))

schleich_anova_treatment <- correlation_analysis(decon_output = decon_dlbcl$Schleich, 
                                                 clinical_characteristic = Schleich_meta$treatment)
schleich_anova_treatmentoutcome <- correlation_analysis(decon_output = decon_dlbcl$Schleich, 
                                                        clinical_characteristic = Schleich_meta$treatment_response)
ggarrange(schleich_anova_treatmentoutcome$aov_plots[[1]], 
          schleich_anova_treatmentoutcome$aov_plots[[2]]) 

schmitz_anova_coo <- correlation_analysis(decon_output = decon_dlbcl$Schmitz, 
                                          clinical_characteristic = Schmitz_meta$Gene.Expression.Subgroup)
schmitz_anova_ipi <- correlation_analysis(decon_output = decon_dlbcl$Schmitz, 
                                          clinical_characteristic = as.character(Schmitz_meta$IPI.Group))
schmitz_anova_ipi$aov_plots$ADR_OHT


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
                                                                              "cluster" = Chapuy_meta$Cluster,
                                                                              "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                              "IPI" = Chapuy_meta$IPI,
                                                                              row.names = rownames(Chapuy_meta)))
# IPI, adr-oht, adr, coo
chapuy_pfs_survival <- survival_analysis(decon_output = decon_dlbcl$Chapuy, OS = chapuy_PFS, censor = chapuy_pfs_zensur, 
                                         clinical_characteristics = data.frame("COO" = Chapuy_meta$COO_byGEP, 
                                                                               "cluster" = Chapuy_meta$Cluster,
                                                                               "treatment" = Chapuy_meta$`R-CHOP-like Chemo`,
                                                                               "IPI" = Chapuy_meta$IPI,
                                                                               row.names = rownames(Chapuy_meta)))
## IPI, adr, adroht
ggpar(chapuy_os_survival$single_kp$IPI, 
      font.main = c(12), font.x = c(14), font.y = c(14),
      font.caption = c(12), font.legend = c(12),font.tickslab = c(12), 
      xlab = "Time in months") # + guides(colour = guide_legend(nrow = 3))


schmitz_PFS <- Schmitz_meta$`Progression_Free Survival _PFS_ Time _yrs`
schmitz_pfs_zensur <- Schmitz_meta$`Progression_Free Survival _PFS_ Status_ 0 No Progressoin_ 1 Progression`
schmitz_pfs_survival <- survival_analysis(decon_output = decon_dlbcl$Schmitz, OS = schmitz_PFS, censor = schmitz_pfs_zensur,
                                          clinical_characteristics = data.frame("COO" = Schmitz_meta$Gene.Expression.Subgroup ,
                                                                                "IPI" = Schmitz_meta$IPI.Group,
                                                                                row.names = rownames(Schmitz_meta)))
## COO, IPI