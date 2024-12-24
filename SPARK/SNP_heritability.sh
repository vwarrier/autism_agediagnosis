#SNP-based heritability scripts
#Discovery_SPARK
qcovar = fread("/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_noage.txt")
covar = fread("/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_covar.txt")

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

merged_total = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total$walked_age_mos[merged_total$walked_age_mos == 888] <- NA
merged_total$used_words_age_mos[merged_total$used_words_age_mos == 888] <- NA
merged_total$sex = ifelse(merged_total$sex == "Male", 1, 0)
merged_total$cognitive_impairment_latest = ifelse(merged_total$cognitive_impairment_latest == "TRUE", 1, 0)
merged_total$reported_cog_test_score_numeric <- factor(merged_total$reported_cog_test_score,
                                                       levels = c("24 - below", "25 - 39", "40 - 54", 
                                                                  "55 - 69", "70 - 79", "80 - 89",
                                                                  "90 - 109", "110 - 119", "120 - 129",
                                                                  "130 - above"),
                                                       labels = 1:10)

# Convert to numeric
merged_total$reported_cog_test_score_numeric <- as.numeric(as.character(merged_total$reported_cog_test_score_numeric))

data3 = merged_total[,c("subject_sp_id", "reported_cog_test_score_numeric",  "scq_total_final_score",  "rbsr_total_final_score", "walked_age_mos",  "used_words_age_mos")]

qcovar2 = merge(qcovar, data3, by.x = "IID", by.y = "subject_sp_id")
qcovar2 = qcovar2[,c(2,1,3:17)]

data4 = merged_total[,c("subject_sp_id", "regress_lang_y_n",  "regress_other_y_n")]
covar2 = merge(covar, data4, by.x = "IID", by.y = "subject_sp_id")
covar2 = covar2[,c(2,1,3:6)]

write.table(covar2, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_covar_fullclinical.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_fullclinical.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "walked_age_mos", "used_words_age_mos")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_devonly.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "walked_age_mos", "used_words_age_mos", "reported_cog_test_score_numeric")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_devonlywithIQ.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "scq_total_final_score", "rbsr_total_final_score")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_autismseverityonly", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "scq_total_final_score", "rbsr_total_final_score", "walked_age_mos", "used_words_age_mos")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_autismseverity_devmilestones.txt", row.names = F, col.names = T, quote = F)


./gcta-1.94.1 \
--threads 15 \
--reml \
--grm /mnt/home4/arc/vw260/SPARK/SPARK_v3/Postimputation_genotyping/omega/alpha_omega_hg19_autismonly_chr \
--pheno /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_pheno.txt \
--grm-cutoff 0.05 \
--qcovar /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_autismseverity_devmilestones.txt \
--covar /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_covar.txt \
--out /mnt/home4/arc/vw260/Autism_heterogeneity/Age_diagnosis_analyses/GCTA_results/h2_autismage_autismseverity_devmilestones


#SNP-based heritability scripts
#Replication


qcovar = fread("/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar.txt")
covar = fread("/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_covar.txt")

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

merged_total = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total$walked_age_mos[merged_total$walked_age_mos == 888] <- NA
merged_total$used_words_age_mos[merged_total$used_words_age_mos == 888] <- NA
merged_total$sex = ifelse(merged_total$sex == "Male", 1, 0)
merged_total$cognitive_impairment_latest = ifelse(merged_total$cognitive_impairment_latest == "TRUE", 1, 0)
merged_total$reported_cog_test_score_numeric <- factor(merged_total$reported_cog_test_score,
                                                       levels = c("24 - below", "25 - 39", "40 - 54", 
                                                                  "55 - 69", "70 - 79", "80 - 89",
                                                                  "90 - 109", "110 - 119", "120 - 129",
                                                                  "130 - above"),
                                                       labels = 1:10)

# Convert to numeric
merged_total$reported_cog_test_score_numeric <- as.numeric(as.character(merged_total$reported_cog_test_score_numeric))

data3 = merged_total[,c("subject_sp_id", "reported_cog_test_score_numeric",  "scq_total_final_score",  "rbsr_total_final_score", "walked_age_mos",  "used_words_age_mos")]

qcovar2 = merge(qcovar, data3, by.x = "IID", by.y = "subject_sp_id")
qcovar2 = qcovar2[,c(2,1,3:17)]

data4 = merged_total[,c("subject_sp_id", "regress_lang_y_n",  "regress_other_y_n")]
covar2 = merge(covar, data4, by.x = "IID", by.y = "subject_sp_id")
covar2 = covar2[,c(2,1,3:7)]

write.table(covar2, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_covar_fullclinical.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_fullclinical.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "walked_age_mos", "used_words_age_mos")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_devonly.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "walked_age_mos", "used_words_age_mos", "reported_cog_test_score_numeric")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_devonlywithIQ.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "scq_total_final_score", "rbsr_total_final_score")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_autismseverityonly", row.names = F, col.names = T, quote = F)
write.table(qcovar2[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "scq_total_final_score", "rbsr_total_final_score", "walked_age_mos", "used_words_age_mos")], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_autismseverity_devmilestones.txt", row.names = F, col.names = T, quote = F)

SES =  fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/background_history_child_2022-12-12.csv")
SES2 = SES[,c("subject_sp_id", "mother_highest_education", "father_highest_education", "mother_occupation", "father_occupation", "annual_household_income")]


covar3 = merge(covar, SES2, by.x = "IID", by.y = "subject_sp_id")
covar3[covar3 == ""] <- NA
write.table(covar3[,c(2,1,3:10)], file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_covar_withSES.txt", row.names = F, col.names = T, quote = F)

deprivation = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/area_deprivation_index_2022-12-12.csv")
setnames(deprivation, "subject_sp_id", "IID")
deprivation = deprivation[,c("IID", "adi_national_rank_percentile")]

qcovar3 = merge(qcovar2, deprivation, by = "IID")

qcovar3 = qcovar3[,c(2,1,3:12,16:18)]
qcovar3[qcovar3 == ""] <- NA

write.table(qcovar3, file= "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_withdeprivation.txt", row.names = F, col.names = T, quote = F)


./gcta-1.94.1 \
--threads 15 \
--reml \
--grm ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_autismonly \
--pheno ~/SPARK/SPARK_v4/GWAS/age_diagnosis_pheno.txt \
--grm-cutoff 0.05 \
--qcovar ~/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar_withdeprivation.txt \
--covar ~/SPARK/SPARK_v4/GWAS/age_diagnosis_covar_withSES.txt \
--out /mnt/home4/arc/vw260/Autism_heterogeneity/Age_diagnosis_analyses/GCTA_results/h2_replication_autism_agediagnosis_SES
