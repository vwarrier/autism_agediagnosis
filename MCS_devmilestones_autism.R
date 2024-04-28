#MCS scripts - Autism + dev milestones

#1: Library
library("tidyr")
library("reshape2")
library("dplyr")
library("ggplot2")
library("data.table")
library(purrr)
setwd("/mnt/beegfs/home4/arc/vw260/MCS")


#2: Read in the data
link = fread("Phenotypes/GDAC_2023_05_WARRIER_mcs_genoID_projectID_link.txt")
setnames(link, "geno_sid", "IID")

demographics = fread("EGA_mcs_basic_demographics_v0003.csv")
setnames(demographics, "geno_sid", "IID")

over10 = fread("./PGS/MCS_autismover10IPSYCH_finalscore.profile")
over10 = over10[,c("FID","IID", "SCORE")]
setnames(over10, "SCORE", "over10_PGS")

under11 = fread("./PGS/MCS_autismunder11IPSYCH_finalscore.profile")
under11 = under11[,c("IID", "SCORE")]
setnames(under11, "SCORE", "under11_PGS")

PCs = fread("MCS_PCsforGWAS_european.txt")
setnames(PCs, "Sample_name", "IID")
PCs = PCs[,c(1:16)]

df_list <- list(demographics, over10, under11, PCs, link)

merged = df_list %>% reduce(full_join, by='IID') # 7841 children, 12416 adults/[parents]


#3: Check and clean phenotypes for the dev milestones, run regression
pheno = fread("./Phenotypes/GDAC_2023_05_WARRIER_mcs_cm_structure_pheno_data_2023-08-16_16-22-56.txt", fill = TRUE)
pheno2 = pheno[pheno$WARRIER_SID %in% merged$project_sid,]
pheno2$ID_merge = paste0(pheno2$WARRIER_FID, "_", pheno2$PNUM)


pheno_raw = fread("./Phenotypes/GDAC_2023_05_WARRIER_mcs_parent_cm_structure_pheno_data_2023-08-16_16-22-56.txt", fill= TRUE)

pheno_raw = pheno_raw %>%  filter(!PNUM=='mcs6_family_with_saliva')
pheno_raw$PNUM = as.numeric(as.character(pheno_raw$PNUM))
pheno_raw$PNUM2 = pheno_raw$PNUM*100
pheno_raw$ID_merge = paste0(pheno_raw$WARRIER_FID, "_", pheno_raw$PNUM2)
setnames(pheno_raw, "WARRIER_SID", "WARRIER_SID_parents")

pheno_merge = pheno2[,c("WARRIER_SID", "WARRIER_FID", "ID_merge")]
pheno_raw_merged = merge(pheno_merge, pheno_raw, by = "ID_merge")


development1 = pheno_raw_merged[,c("WARRIER_SID","ACSMIL00", "ACSITU00", "ACSTAN00","ACHAND00","ACGRAB00",
                             "ACPICK00","ACPTOY00", "ACWALK00")]

development1 = development1 %>% mutate(across(ACSMIL00:ACNODS00, recode, '1'= 0, "2" = 0, "3"=1))
development1$fine = development1$ACHAND00 + development1$ACGRAB00 + development1$ACPICK00 + development1$ACPTOY00
development1$gross = development1$ACSITU00 + development1$ACSTAN00 + development1$ACWALK00 
development1$cdi = development1$ACGIVE00 + development1$ACWAVE00 + development1$ACARMS00  + development1$ACNODS00


merged_development1 = merge(merged, development1, by.x = "project_sid", by.y = "WARRIER_SID")


list1 = c("fine", "gross")
results_over10 = NULL
results_under11 = NULL
Samplesize = NULL

for(i in list1){
  results_all = summary(lm(scale(merged_development1[[i]]) ~ scale(over10_PGS) + scale(under11_PGS) + sex + age + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged_development1))
  results_over10 = rbind(results_over10, cbind(i, t(results_all$coefficients[2,])))
  Samplesize = rbind(Samplesize, results_all$df.null)
  results_under11 = rbind(results_under11, cbind(i, t(results_all$coefficients[3,])))
}



#4: Check and clean phenotypes for the autism, run regression
autism = pheno_raw_merged[,c("WARRIER_SID","CPAUTS00", "DPAUTS00", "EPAUTS00","FPAUTS00")]
autism = autism %>% mutate(across(CPAUTS00:FPAUTS00, recode, '1'= 1, "2" = 0))

autism$under8 = ifelse(autism$CPAUTS00 ==1 | autism$DPAUTS00 ==1, 1, 0)
autism$under8 = ifelse(autism$under8 ==0 & autism$FPAUTS00 ==1, NA, autism$under8)
autism$under8 = ifelse(autism$under8 ==0 & autism$EPAUTS00 ==1, NA, autism$under8)

autism$under11 = ifelse(autism$CPAUTS00 ==1 | autism$DPAUTS00 == 1 | autism$EPAUTS00 == 1, 1, 0)
autism$under11 = ifelse(autism$under11 == 0 & autism$FPAUTS00 ==1, NA, autism$under11)

autism$over12 = ifelse(autism$FPAUTS00 == 1, 1, 0)
autism$over12 = ifelse(autism$over12 == 1 & autism$CPAUTS00 ==1, NA, autism$over12)
autism$over12 = ifelse(autism$over12 == 1 & autism$DPAUTS00 ==1, NA, autism$over12)
autism$over12 = ifelse(autism$over12 == 1 & autism$EPAUTS00 ==1, NA, autism$over12)


merged_autism = merge(merged, autism, by.x = "project_sid", by.y = "WARRIER_SID")
library(logistf)


summary(logistf(under8 ~ scale(over10_PGS) + scale(under11_PGS) + sex + age + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,family = "binomial", data = merged_autism))
summary(logistf(under11 ~ scale(over10_PGS) + scale(under11_PGS) + sex + age + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,family = "binomial", data = merged_autism))
summary(logistf(over12 ~ scale(over10_PGS) + scale(under11_PGS) + sex + age + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,family = "binomial", data = merged_autism))


