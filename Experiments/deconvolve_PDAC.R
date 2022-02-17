## Mastherthesis, Melanie Fattohi
## scRNA-seq data by Tosti, Baron
## scRNA-seq von PDAC werden noch gesucht
## Deconvolve Yang, PAAD, Moffitt (Seq + array), Guo, Kirby, Janky
## alle bulks reinladen
## alle sc reinladen (entweder als table oder schon qc'ed by SCDC)
## alle bulks in eine liste speichern
## mit lapply compute_pvalue ausfuehren
## nreps = 1000, ncores = 15
## ergebnisse pro dataset darstellen
## survival analysis machen
## vllt iwas mit basal, classical und hybrid type
## ecotyper zeugs?

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/visualization.R")
source("~/Masterthesis/CancerDeconvolution/Scripts/survival_analysis.R")

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
## 131 samples; microarray
Janky_bulk <- readRDS("~/Masterthesis/Data/Bulk/Janky/Janky_bulk.RDS")
Janky_meta <- readRDS("~/Masterthesis/Data/Bulk/Janky/Janky_metadata.RDS")
## 51 samples; survival; RNA-seq
Kirby_bulk <- readRDS("~/Masterthesis/Data/Bulk/Kirby/Kirby_bulk.RDS")
Kirby_meta <- readRDS("~/Masterthesis/Data/Bulk/Kirby/Kirby_metadata.RDS")
## 61 samples; RNA-seq
Moffitt_seq_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_seq_bulk.RDS")
Moffitt_seq_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_seq_metadata.RDS")
## 357 samples; survival; tumor subtype; microarray
Moffitt_array_bulk <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_bulk.RDS")
Moffitt_array_meta <- readRDS("~/Masterthesis/Data/Bulk//Moffitt/Moffitt_array_metadata.RDS")

bulk_list <- list("PAAD" = PAAD_bulk,
                  "Yang" = Yang_bulk,
                  "Guo" = Guo_bulk,
                  "Janky" = Janky_bulk,
                  "Kirby" = Kirby_bulk,
                  "Moffitt_seq" = Moffitt_seq_bulk,
                  "Moffitt_array" = Moffitt_array_bulk)
bulk_meta_list <- list("PAAD" = PAAD_meta,
                       "Yang" = Yang_meta,
                       "Guo" = Guo_meta,
                       "Janky" = Janky_meta,
                       "Kirby" = Kirby_meta,
                       "Moffitt_seq" = Moffitt_seq_meta,
                       "Moffitt_array" = Moffitt_array_meta)

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