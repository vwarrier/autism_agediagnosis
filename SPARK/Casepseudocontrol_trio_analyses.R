#R, case-pseudocontrol regression}
rm(list = ls())

library(data.table)
library(MASS)
library(lme4)
library(tidyverse)
library(boot)
set.seed(1029384756)

Autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_autism_pseudocons_finalscore.profile")
Autism = Autism[,c("FID", "IID", "SCORE")]
setnames(Autism, "SCORE", "Autism_PGS")

ADHD = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_ADHD_pseudocons_finalscore.profile")
ADHD = ADHD[,c("IID", "SCORE")]
setnames(ADHD, "SCORE", "ADHD_PGS")

scz = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_scz_pseudocons_finalscore.profile")
scz = scz[,c("IID", "SCORE")]
setnames(scz, "SCORE", "scz_PGS")


bipolar = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_bipolar_pseudocons_finalscore.profile")
bipolar = bipolar[,c("IID", "SCORE")]
setnames(bipolar, "SCORE", "bipolar_PGS")


edu = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_edu_pseudocons_finalscore.profile")
edu = edu[,c("IID", "SCORE")]
setnames(edu, "SCORE", "edu_PGS")


IQ = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_IQ_pseudocons_finalscore.profile")
IQ = IQ[,c("IID", "SCORE")]
setnames(IQ, "SCORE", "IQ_PGS")

depression = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_dep_pseudocons_finalscore.profile")
depression = depression[,c("IID", "SCORE")]
setnames(depression, "SCORE", "depression_PGS")

df_list <- list(Autism, ADHD, scz, IQ, edu, bipolar,  depression)

merged = df_list %>% reduce(full_join, by='IID')

check1 = data.frame(do.call('rbind', strsplit(as.character(merged$IID),'_',fixed=TRUE)))

merged = cbind(merged, check1)
merged = merged[,c("FID", "X1", "Autism_PGS", "ADHD_PGS", "bipolar_PGS", "depression_PGS", "scz_PGS", "edu_PGS", "IQ_PGS", "X2" )]
setnames(merged, old = c("X1", "X2"), new = c("IID", "Transmission_status"))

transmitted = subset(merged, Transmission_status == "T")
untransmitted = subset(merged, Transmission_status == "U")

merged = merge(transmitted, untransmitted, by = "IID")
setnames(merged, old = c("FID.x", "Autism_PGS.x", "ADHD_PGS.x", "bipolar_PGS.x", "depression_PGS.x", "scz_PGS.x","edu_PGS.x", "IQ_PGS.x", "Autism_PGS.y", "ADHD_PGS.y", "bipolar_PGS.y", "depression_PGS.y", "scz_PGS.y","edu_PGS.y", "IQ_PGS.y"), new = c("FID", "Autism_PGS.T", "ADHD_PGS.T", "bipolar_PGS.T", "depression_PGS.T", "scz_PGS.T", "edu_PGS.T","IQ_PGS.T", "Autism_PGS.U", "ADHD_PGS.U", "bipolar_PGS.U", "depression_PGS.U", "scz_PGS.U", "edu_PGS.U", "IQ_PGS.U"))

merged = merged[,-c("Transmission_status.x", "FID.y", "Transmission_status.y")]

PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt")
setnames(PCs, "Sample_name", "IID")
merged = merge(merged, PCs, by = "IID")



###Read the pheno files###

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442


merged_total = merge(merged, autism_eligible, by.x = "IID", by.y = "subject_sp_id")
merged_total$diagnosis_age3 = scale(merged_total$diagnosis_age/12)


summary(lmer(diagnosis_age3 ~ scale(edu_PGS.T) + scale(edu_PGS.U) +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + (1|FID), data = merged_total))

summary(lmer(diagnosis_age3 ~ scale(ADHD_PGS.T) + scale(ADHD_PGS.U) +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest + (1|FID), data = merged_total))

#Cor between T and U:  ADHD: 0.022 (-0.002 - 0.04), Edu: 0.19 (0.17 - 0.21)

###Bootstrapping - scripts modified from Perline Demange (see: https://github.com/PerlineDemange/GeneticNurtureNonCog)

# Bootstrapping
nboot <- 10000

trio_lm<-function(data,index){
  datx<-data[index,]
  mod <- lm(diagnosis_age3 ~ scale(ADHD_PGS.T) + scale(ADHD_PGS.U) +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = datx)
  coef(mod) #get fixed effects
}

summary( lm(diagnosis_age3 ~ scale(ADHD_PGS.T) + scale(ADHD_PGS.U) +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = merged_total))

trio_lm_boot<-boot(merged_total, trio_lm, nboot, parallel = "multicore", ncpus=20) 



# Plot to check bootstrapping
png("SPARK.trios.EA.bootstrap_lm_2023jul05.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(trio_EA_lm_boot)
dev.off()

# Save t output of boot
bootoutput <- as.data.frame(trio_lm_boot$t)
colnames(bootoutput) <- rownames(as.data.frame(trio_lm_boot$t0))

# Get values out of boot.out for all estimates + create direct and ratio estimates
original <- as.data.frame(t(trio_lm_boot$t0)) # estimates of the original sample #best estimates of the effects

original$total_effect <- original$`scale(ADHD_PGS.T)`
original$indirect_effect <- original$`scale(ADHD_PGS.U)`
original$direct_effect <- original$total_effect - original$indirect_effect
original$ratio_effect <- original$indirect_effect / original$direct_effect
original$ratio_tot_effect <- original$indirect_effect / original$total_effect

bootoutput$total_effect <- bootoutput$`scale(ADHD_PGS.T)`
bootoutput$indirect_effect <- bootoutput$`scale(ADHD_PGS.U)`
bootoutput$direct_effect <- bootoutput$total_effect - bootoutput$indirect_effect
bootoutput$ratio_effect <- bootoutput$indirect_effect / bootoutput$direct_effect
bootoutput$ratio_tot_effect <- bootoutput$indirect_effect / bootoutput$total_effect

mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
tot <- statsoutput[,(ncol(statsoutput)-10):ncol(statsoutput)] # get only summary statistics 
tot <- tot[,c(ncol(tot), 1:(ncol(tot)-1))]
tot


# Comparing estimates

diffeffect <- original$direct_effect - original$indirect_effect
SD_sampling_diffeffect <- sd(bootoutput$direct_effect - bootoutput$indirect_effect)
Z_diffeffect <- diffeffect/SD_sampling_diffeffect
P_diffeffect <- 2*pnorm(-abs(Z_diffeffect))

```


```{r Cleaning data, trio vs non-trio}
#Here we check if there are differences in PGS between trios and non-trios

library(data.table)
library(MASS)
library(lme4)
library(tidyverse)



Autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_autism_finalscore.profile") #47170
Autism = Autism[,c("FID", "IID", "SCORE")]
setnames(Autism, "SCORE", "Autism_PGS")

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

df_list <- list(Autism, ADHD, scz, IQ, edu, bipolar,  depression)

merged = df_list %>% reduce(full_join, by='IID')

PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt") # 38280 
setnames(PCs, "Sample_name", "IID")
merged = merge(merged, PCs, by = "IID")




###Read the pheno files###

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total = merge(merged, autism_eligible, by.x = "IID", by.y = "subject_sp_id")
merged_total$diagnosis_age2 = log(merged_total$diagnosis_age/12 + 1)
merged_total$diagnosis_age3 = scale(merged_total$diagnosis_age/12)
merged_total$diagnosis_age4 = merged_total$diagnosis_age/12

registration = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")

trios =  registration[!(registration$biofather_id =="" | registration$biomother_id==""), ] # 65143
asd_trios = subset(trios, asd == "TRUE") #44236


merged_pgs_trio = merge(merged_total, asd_trios, by.x = "IID", by.y = "subject_sp_id") #11079

merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biomother_id.x %in% merged$IID,] #8312, 2980
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biofather_id.x %in% merged$IID,] #6554

merged_total$trio = ifelse(merged_total$IID %in% merged_pgs_trio$IID, 1, 0)


list1 = c("ADHD_PGS", "Autism_PGS", "bipolar_PGS", "depression_PGS", "scz_PGS", "edu_PGS", "IQ_PGS")

results_model1 = NULL

for(i in list1){
  results_all = summary(glm(trio ~ scale(merged_total[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = merged_total,  family = "binomial"))
  
  results_model1 = rbind(results_model1, cbind(i, t(results_all$coefficients[2,])))
}



results_model2 = NULL

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total[[i]]) + trio + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
  
  results_model2 = rbind(results_model2, cbind(i, t(results_all$coefficients[2,])))
}

results_model3 = NULL

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(merged_total[[i]])*trio + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = merged_total))
  
  results_model3 = rbind(results_model3, cbind(i, t(results_all$coefficients[17,])))
}


trios = subset(merged_total, trio == "1")

results_model4 = NULL

for(i in list1){
  results_all = summary(lm(diagnosis_age3 ~ scale(trios[[i]]) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + sex + age_at_registration_years + cognitive_impairment_latest, data = trios))
  
  results_model4 = rbind(results_model4, cbind(i, t(results_all$coefficients[2,])))
}

