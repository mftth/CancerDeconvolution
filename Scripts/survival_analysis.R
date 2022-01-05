## Mastherthesis, Melanie Fattohi
## survival analysis 
## test if relevant cell type proportions are predictive in regards to disease-related survival
## relevant cell types are determined through Machine Learning (feature importance) or otherwise
## just testing cell types of interest
## also: visualization using Kaplan-Meier plots

source("~/Masterthesis/CancerDeconvolution/Scripts/Permute_basis.R")
library(survival)
library(survminer)
library(ggplot2)

## need input meta data 
# angabe über welche spalten der meta data in survival analysis getestet werden soll
#+ cell type(s) of interest + additional clinical characteristic of interest as a whole 

repset_meta$OS_Tissue <- gsub(",", ".", repset_meta$OS_Tissue)
repset_meta$OS_Tissue <- as.numeric(repset_meta$OS_Tissue)
repset_meta$ductal <-  Repset_scdc_baron$decon_res$prop.est.mvw[,"ductal"]
repset_meta$ductal_high_low <- NA
repset_meta$ductal_high_low[which(repset_meta$ductal <= mean(repset_meta$ductal))] <- "ductal_low"
repset_meta$ductal_high_low[which(repset_meta$ductal > mean(repset_meta$ductal))] <- "ductal_high"
repset_meta$alpha <-  Repset_scdc_baron$decon_res$prop.est.mvw[,"alpha"]
repset_meta$alpha_high_low <- NA
repset_meta$alpha_high_low[which(repset_meta$alpha <= mean(repset_meta$alpha))] <- "alpha_low"
repset_meta$alpha_high_low[which(repset_meta$alpha > mean(repset_meta$alpha))] <- "alpha_high"


#km <- with(repset_meta, Surv(OS_Tissue, Zensur))
#km_fit <- survfit(Surv(OS_Tissue, Zensur) ~ 1, data = repset_meta)
#plot(km_fit)
km_fit_ductal <- survfit(Surv(OS_Tissue, Zensur) ~ ductal_high_low, data = repset_meta)
km_fit_alpha <- survfit(Surv(OS_Tissue, Zensur) ~ alpha_high_low, data = repset_meta)
#plot(km_fit_ductal)
surv_ductal <- survminer::surv_pvalue(km_fit_ductal, data = repset_meta)
surv_ductal$pval
surv_alpha <- survminer::surv_pvalue(km_fit_alpha, data = repset_meta)
surv_alpha$pval

survminer::ggsurvplot(km_fit_ductal, data = repset_meta, pval = TRUE)#, 
#                      xlab = "month", ylab = "survival rate")

survminer::ggsurvplot(list(ductal = km_fit_ductal, alpha = km_fit_alpha), data = repset_meta,
                      combine = TRUE)
