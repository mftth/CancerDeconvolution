## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## Deconvolve repset (Scarpa + Riemer), Fadista, Missiaglia, Alvarez, Sadanandam
## alle bulks reinladen
## alle sc reinladen (entweder als table oder schon qc'ed by SCDC)
## alle bulks in eine liste speichern
## mit lapply compute_pvalue ausfuehren
## nreps = 1000, ncores = 15
## ergebnisse pro dataset darstellen
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
bulk_meta <- read.table("~/Masterthesis/Data/Bulk/PanNEN/Meta_information.tsv", header = TRUE, sep = "\t")

alvarez_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Alvarez.S104.HGNC.tsv", header = TRUE, sep = "\t")
alvarez_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(alvarez_bulk)),]
rownames(alvarez_meta) <- alvarez_meta$Sample
alvarez_meta <- alvarez_meta[match(colnames(alvarez_bulk), rownames(alvarez_meta)),]
alvarez_idx <- which(alvarez_meta$Histology_Primary == "Pancreatic")
alvarez_bulk <- alvarez_bulk[,alvarez_idx]
alvarez_meta <- alvarez_meta[alvarez_idx,]

fadista_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Fadista.S89.tsv", header = TRUE, sep = "\t", row.names = 1)
fadista_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(fadista_bulk)),]
rownames(fadista_meta) <- fadista_meta$Sample

master_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Master.S34.HGNC.tsv", header = TRUE, sep = "\t", check.names = FALSE)
master_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(master_bulk)),]
rownames(master_meta) <- master_meta$Sample
master_meta <- master_meta[match(colnames(master_bulk), rownames(master_meta)),]

missiaglia_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Missaglia.S75.tsv", header = TRUE, sep = "\t", check.names = FALSE)
missiaglia_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(missiaglia_bulk)),]
rownames(missiaglia_meta) <- missiaglia_meta$Sample

sadanandam_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Sadanandam.S29.tsv", header = TRUE, sep = "\t", check.names = FALSE)
sadanandam_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(sadanandam_bulk)),]
rownames(sadanandam_meta) <- sadanandam_meta$Sample

scarpa_bulk <- read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Scarpa.New.HGNC.tsv", header = TRUE, sep = "\t", check.names = FALSE)
scarpa_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(scarpa_bulk)),]
rownames(scarpa_meta) <- scarpa_meta$Sample

diedisheim_bulk <-  read.table(file = "~/Masterthesis/Data/Bulk/PanNEN/Diedisheim.S66.HGNC.tsv", header = TRUE, sep = "\t", row.names = 1)
diedisheim_meta <- bulk_meta[which(bulk_meta$Sample %in% colnames(diedisheim_bulk)),]
rownames(diedisheim_meta) <- diedisheim_meta$Sample
diedisheim_meta <- diedisheim_meta[match(colnames(diedisheim_bulk), rownames(diedisheim_meta)),]

repset_bulk <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset.RDS")
repset_meta <- readRDS("~/Masterthesis/Data/Bulk/RepSet/repset_meta.RDS")

bulk_list <- list("Alvarez" = alvarez_bulk,
                  "Diedisheim" = diedisheim_bulk,
                  "Fadista" = fadista_bulk,
                  "Master" = master_bulk,
                  "Missiaglia" = missiaglia_bulk,
                  "RepSet" = repset_bulk,
                  "Sadanandam" = sadanandam_bulk,
                  "Scarpa" = scarpa_bulk)
bulk_meta_list <- list("Alvarez" = alvarez_meta,
                       "Diedisheim" = diedisheim_meta,
                       "Fadista" = fadista_meta,
                       "Master" = master_meta,
                       "Missiaglia" = missiaglia_meta,
                       "RepSet" = repset_meta,
                       "Sadanandam" = sadanandam_meta,
                       "Scarpa" = scarpa_meta)


## single-cell RNA-seq datasets
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
res_path_baron <- "~/Masterthesis/CancerDeconvolution/Results/PanNEN_deconvolution/Baron"
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
res_path_tosti <- "~/Masterthesis/CancerDeconvolution/Results/PanNEN_deconvolution/Tosti"
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