#Polygenic scoring in SPARK"

library(data.table)
library(MASS)
library(lme4)
library(tidyverse)
library(boot)

Autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_autism_finalscore.sscore") #47170
Autism = Autism[,c("#FID", "IID", "SCORE1_AVG")]
setnames(Autism, old = c("SCORE1_AVG", "#FID"), new = c("Autism_PGS", "FID"))

ADHD = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_ADHD_finalscore.profile")
ADHD = ADHD[,c("IID", "SCORE")]
setnames(ADHD, "SCORE", "ADHD_PGS")

scz = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_scz_finalscore.profile")
scz = scz[,c("IID", "SCORE")]
setnames(scz, "SCORE", "scz_PGS")


bipolar = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_bipolar_finalscore.profile")
bipolar = bipolar[,c("IID", "SCORE")]
setnames(bipolar, "SCORE", "bipolar_PGS")


edu = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_edu_finalscore.profile")
edu = edu[,c("IID", "SCORE")]
setnames(edu, "SCORE", "edu_PGS")


IQ = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_IQ_finalscore.profile")
IQ = IQ[,c("IID", "SCORE")]
setnames(IQ, "SCORE", "IQ_PGS")

depression = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_dep_finalscore.profile")
depression = depression[,c("IID", "SCORE")]
setnames(depression, "SCORE", "depression_PGS")

over10 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHover10_finalscore.profile")
over10 = over10[,c("IID", "SCORE")]
setnames(over10, "SCORE", "over10_PGS")

under11 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHunder11_finalscore.profile")
under11 = under11[,c("IID", "SCORE")]
setnames(under11, "SCORE", "under11_PGS")



df_list <- list(Autism, ADHD, scz, IQ, edu, bipolar,  depression, over10, under11)

merged = df_list %>% reduce(full_join, by='IID')

PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt") # 47170
setnames(PCs, "Sample_name", "IID")
merged = merge(merged, PCs, by = "IID")


###Read the pheno files###

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total = merge(merged, autism_eligible, by.x = "IID", by.y = "subject_sp_id")
merged_total$diagnosis_age3 = scale(merged_total$diagnosis_age/12)



#Run regressions
#Model 1: ID included as a covariate

list1 = c("ADHD_PGS", "Autism_PGS", "bipolar_PGS", "depression_PGS", "scz_PGS", "edu_PGS", "IQ_PGS", "over10_PGS", "under11_PGS")

results_model1 = NULL
results_model2 = NULL
results_model3 = NULL
results_model4 = NULL
results_model5 = NULL

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
  results_model1 = rbind(results_model1, cbind(i, t(results_all$coefficients[2,])))
}




#Model 2: Sensitivity analyses - developmental milestones
merged_total$walked_age_mos = ifelse(merged_total$walked_age_mos == "888", NA, merged_total$walked_age_mos)
merged_total$used_words_age_mos = ifelse(merged_total$used_words_age_mos == "888", NA, merged_total$used_words_age_mos)

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + walked_age_mos + used_words_age_mos, data = merged_total))
  results_model2 = rbind(results_model2, cbind(i, t(results_all$coefficients[2,])))
}


#Model 3: With SES

SES =  fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/background_history_child_2022-12-12.csv")
SES2 = SES[,c("subject_sp_id", "mother_highest_education", "father_highest_education", "mother_occupation", "father_occupation", "annual_household_income")]

setnames(SES2, "subject_sp_id", "IID")
merged_total_ses= merge(merged_total, SES2, by = "IID")

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total_ses[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + walked_age_mos + used_words_age_mos + mother_occupation + father_occupation + annual_household_income + father_highest_education + mother_highest_education, data = merged_total_ses))
  results_model3 = rbind(results_model3, cbind(i, t(results_all$coefficients[2,])))
}



#Model 4: With SES and deprivation

deprivation = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/area_deprivation_index_2022-12-12.csv")
setnames(deprivation, "subject_sp_id", "IID")

merged_total_dep = merge(merged_total_ses, deprivation, by = "IID")

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total_dep[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + walked_age_mos + used_words_age_mos + mother_occupation + father_occupation + annual_household_income + adi_national_rank_percentile, data = merged_total_dep))
  results_model4 = rbind(results_model4, cbind(i, t(results_all$coefficients[2,])))
}


##Sex-diff

summary(lm(diagnosis_age3 ~ scale(ADHD_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(depression_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(edu_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(scz_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(over10_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(under11_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(under9_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
summary(lm(diagnosis_age3 ~ scale(over11_PGS)*sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + age_at_registration_years + cognitive_impairment_latest, data = merged_total))


# Accounting for ADHD diagnosis
health = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/basic_medical_screening_2022-12-12.csv")
health2 = health[,c("subject_sp_id", "attn_behav")]
health2[is.na(health2)] <- 0

setnames(health2, "subject_sp_id", "IID")

merged_total_health = merge(merged_total, health2, by = "IID")

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total_health[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + attn_behav, data = merged_total_health))
  results_model5 = rbind(results_model5, cbind(i, t(results_all$coefficients[2,])))
}


#Sensitivity analyses - restricting to individuals without ID and who can speak in sentences
merged_noID = subset(merged_total, cognitive_impairment_latest  == FALSE)
merged_noID_fullsentence = subset(merged_noID, language_level_latest  == "Uses longer sentences of his/her own and is able to tell you something that happened")


results_model6 = NULL

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_noID_fullsentence[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years, data = merged_noID_fullsentence))
  results_model6 = rbind(results_model6, cbind(i, t(results_all$coefficients[2,])))
}



