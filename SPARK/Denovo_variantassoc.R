```{R, Association with de novo SNVs}
phenos = fread("~/Autism_heterogeneity/Mahmoud_denovos/SPARK_iWES2_trios_phen_mm60_20240131-2.csv")
autism = subset(phenos, ASD == 2)

denovo = fread("~/Autism_heterogeneity/Mahmoud_denovos/SPARK_iWES2_trios_dnms_mm60_20240131.csv")
PTV = subset(denovo, CSQ == "frameshift_variant" | CSQ == "start_lost" | CSQ == "stop_gained" | CSQ == "splice_donor_variant" | CSQ == "splice_acceptor_variant" | CSQ == "stop_lost")
PTV = subset(PTV, oe_lof_upper_bin == 0)


missense = subset(denovo, CSQ == "missense_variant")
missense = subset(missense, MPC > 1.9999)
missense = subset(missense, oe_lof_upper_bin == 0)

denovo_constraint = rbind(PTV, missense)

autism$denovo = ifelse(autism$IID %in% denovo_constraint$IID, 1, 0)
autism$PTV = ifelse(autism$IID %in% PTV$IID, 1, 0)
autism$missense = ifelse(autism$IID %in% missense$IID, 1, 0)


autism2 = autism[,c("IID", "denovo", "PTV", "missense")]

Autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_autism_finalscore.sscore") #47170
Autism = Autism[,c("#FID", "IID", "SCORE1_AVG")]
setnames(Autism, "SCORE1_AVG", "Autism_PGS")
setnames(Autism, "#FID", "FID")

over10 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHover10_finalscore.profile")
over10 = over10[,c("IID", "SCORE")]
setnames(over10, "SCORE", "over10_PGS")

under11 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHunder11_finalscore.profile")
under11 = under11[,c("IID", "SCORE")]
setnames(under11, "SCORE", "under11_PGS")

merged = merge(Autism, under11, by = "IID") #Hardcoding it as reduce isn't happening
merged = merge(merged, over10, by = "IID")


PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt") # 38280 
setnames(PCs, "Sample_name", "IID")
merged = merge(merged, PCs, by = "IID")

merged = merge(merged, autism2, by = "IID") #6634


summary(glm(denovo ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))
summary(glm(PTV ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))
summary(glm(missense ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total = merge(merged, autism_eligible, by.x = "IID", by.y = "subject_sp_id")

merged_total$diagnosis_age2 = log(merged_total$diagnosis_age/12 + 1)
merged_total$diagnosis_age3 = scale(merged_total$diagnosis_age/12) #6206

summary(lm(diagnosis_age3 ~ denovo + sex , data = merged_total))
summary(lm(diagnosis_age3 ~ PTV + sex, data = merged_total))
summary(lm(diagnosis_age3 ~ missense + sex, data = merged_total))


#inherited
inherited = fread("~/Autism_heterogeneity/Mahmoud_denovos/SPARK_iWES2_trios_inherited_mm60_20240131.csv")

autism = subset(inherited, ASD == "Probands")


PTV = subset(inherited, CSQ == "frameshift_variant" | CSQ == "start_lost" | CSQ == "stop_gained" | CSQ == "splice_donor_variant" | CSQ == "splice_acceptor_variant" | CSQ == "stop_lost")
PTV = subset(PTV, oe_lof_upper_bin == 0)


missense = subset(inherited, CSQ == "missense_variant")
missense = subset(missense, MPC > 1.9999)
missense = subset(missense, oe_lof_upper_bin == 0)

inherited_constraint = rbind(PTV, missense)

autism$inherited= ifelse(autism$IID %in% inherited_constraint$IID, 1, 0)
autism$PTV = ifelse(autism$IID %in% PTV$IID, 1, 0)
autism$missense = ifelse(autism$IID %in% missense$IID, 1, 0)


autism2 = autism[,c("IID", "inherited", "PTV", "missense")]

Autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_autism_finalscore.sscore") #47170
Autism = Autism[,c("#FID", "IID", "SCORE1_AVG")]
setnames(Autism, "SCORE1_AVG", "Autism_PGS")
setnames(Autism, "#FID", "FID")

over10 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHover10_finalscore.profile")
over10 = over10[,c("IID", "SCORE")]
setnames(over10, "SCORE", "over10_PGS")

under11 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHunder11_finalscore.profile")
under11 = under11[,c("IID", "SCORE")]
setnames(under11, "SCORE", "under11_PGS")

merged = merge(Autism, under11, by = "IID") #Hardcoding it as reduce isn't happening
merged = merge(merged, over10, by = "IID")


PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt") # 38280 
setnames(PCs, "Sample_name", "IID")
merged = merge(merged, PCs, by = "IID")

merged = merge(merged, autism2, by = "IID") #6634

summary(glm(inherited ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))
summary(glm(PTV ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))
summary(glm(missense ~ scale(over10_PGS) + scale(under11_PGS) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged))

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total = merge(merged, autism_eligible, by.x = "IID", by.y = "subject_sp_id")

merged_total$diagnosis_age2 = log(merged_total$diagnosis_age/12 + 1)
merged_total$diagnosis_age3 = scale(merged_total$diagnosis_age/12) #6206

summary(lm(diagnosis_age3 ~ inherited + sex +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged_total))
summary(lm(diagnosis_age3 ~ PTV + sex +  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged_total))
summary(lm(diagnosis_age3 ~ missense + sex + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged_total))

