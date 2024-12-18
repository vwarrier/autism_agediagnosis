#########Notes###########
## GWAS - age at diagnosis
###Read the imputed GWAS files###

data1 = fread("~/SPARK/SPARK_v4/postimputation/alpha/sequencing_alpha_omega_hg19_allchrs.fam")

###Read the pheno files###

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

data1_autismeligible_autism = data1[data1$V2 %in% core_autism$subject_sp_id,] #12084

data1_autismeligible = data1[data1$V2 %in% autism_eligible$subject_sp_id,] #9701
write.table(data1_autismeligible[,c(1:2)], file = "~/SPARK/SPARK_v4/postimputation/sequencing_alphaomega_autismonly.txt", row.names = F, col.names = T, quote = F )

### Lets make a GRM

for i in {1..22}; do ./gcta-1.94.1  --bfile ~/SPARK/SPARK_v4/postimputation/alpha/sequencing_alpha_omega_hg19_allchrs --keep ~/SPARK/SPARK_v4/postimputation/sequencing_alphaomega_autismonly.txt --autosome  --make-grm  --maf 0.01 --thread-num 20 --chr ${i} --out ~/SPARK/SPARK_v4/postimputation/sequencing_alphaomega_hg19_autismonly_chr${i}; done
./gcta-1.94.1 --mgrm ~/SPARK/SPARK_v4/postimputation/mbfile.txt --make-grm  --out ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_autismonly --thread-num 20
./gcta-1.94.1 --grm ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_autismonly --thread-num 20 --make-bK-sparse 0.05 --out ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_autismonly_sparse


###Create phenotype and covariate file for GWAS

PCs = fread("~/SPARK/SPARK_v4/SPARK_twistsequencing_PCsforGWAS.txt")
setnames(PCs, "Sample_name", "IID")
core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442
autism_eligible = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total = merge(PCs, autism_eligible, by.x = "IID", by.y = "subject_sp_id")
merged_total$diagnosis_age = scale(merged_total$diagnosis_age/12)

setnames(merged_total, "family_sf_id", "FID")

pheno = merged_total[,c("FID", "IID", "diagnosis_age")]
qcovar = merged_total[,c("FID", "IID", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")]
covar = merged_total[,c("FID", "IID", "sex", "cognitive_impairment_latest")]

fam1 = fread("~/SPARK/SPARK_v4/postimputation/alpha/alpha_chr1_hg19_cleaned_v2.fam")
fam2 = fread("~/SPARK/SPARK_v4/postimputation/omega/omega_chr1_hg19_cleaned_v2.fam")

covar$batch = NA
covar$batch = ifelse(covar$IID %in% fam1$V2, "alpha", "omega")


write.table(pheno, file = "~/SPARK/SPARK_v4/GWAS/age_diagnosis_pheno.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar, file = "~/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar.txt", row.names = F, col.names = T, quote = F)
write.table(covar, file = "~/SPARK/SPARK_v4/GWAS/age_diagnosis_covar.txt", row.names = F, col.names = T, quote = F)


###Run GWAS

./gcta-1.94.1 --mbfile /mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v4/GWAS/mbfile_agediagnosis_GWAS.txt --grm-sparse ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse --fastGWA-mlm --pheno ~/SPARK/SPARK_v4/GWAS/age_diagnosis_pheno.txt --qcovar ~/SPARK/SPARK_v4/GWAS/age_diagnosis_qcovar.txt --covar ~/SPARK/SPARK_v4/GWAS/age_diagnosis_covar.txt --thread-num 20 --maf 0.01 --out ~/SPARK/SPARK_v4/GWAS/Agediagnosis_linear_noage --exclude ~/SPARK/SPARK_v4/GWAS/batch_SNPs.txt


##Read GWAS output

data1 = fread("/mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v3/GWAS/Age_of_diagnosis/geno_assoc.fastGWA")




###################################################################
#### GWAS for different age at diagnosis groups, using fastGWA#####
###################################################################


### Lets make a GRM

for i in {1..22}; do ./gcta-1.94.1  --bfile ~/SPARK/SPARK_v4/postimputation/alpha/sequencing_alpha_omega_hg19_allchrs --autosome  --make-grm  --maf 0.01 --thread-num 20 --chr ${i} --out ~/SPARK/SPARK_v4/postimputation/sequencing_alphaomega_hg19_chr${i}; done
./gcta-1.94.1 --mgrm ~/SPARK/SPARK_v4/postimputation/mbfile_allpeople.txt --make-grm  --out ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs --thread-num 20
./gcta-1.94.1 --grm ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs --thread-num 20 --make-bK-sparse 0.05 --out ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse


### Create phenotype files
fam = fread("~/SPARK/SPARK_v4/postimputation/alpha/sequencing_alpha_omega_hg19_allchrs.fam")

core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/individuals_registration_2022-12-12.csv")
covar = core[,c("subject_sp_id", "sex", "age_at_registration_years", "diagnosis_age", "asd")]

merged = merge(fam, covar, by.x = "V2", by.y = "subject_sp_id")
merged$autism = ifelse(merged$asd  == "TRUE", 1, 0)

over10 = subset(merged, diagnosis_age > 120) #3358, 2885 repli
under10 = subset(merged, diagnosis_age < 121) #18719, 9162 repli
under6 = subset(merged, diagnosis_age < 72) #14587, 6857 repli
controls = subset(merged, autism == 0) #24965, 33302


merged$over10 = ifelse(merged$V2 %in% over10$V2, 1, NA)
merged$under10 = ifelse(merged$V2 %in% under10$V2, 1, NA)
merged$under6 = ifelse(merged$V2 %in% under6$V2, 1, NA)

merged$over10 = ifelse(merged$V2 %in% controls$V2, 0, merged$over10)
merged$under10 = ifelse(merged$V2 %in% controls$V2, 0, merged$under10)
merged$under6 = ifelse(merged$V2 %in% controls$V2, 0, merged$under6)

pcs = fread("~/SPARK/SPARK_v4/SPARK_twistsequencing_PCsforGWAS.txt")

pheno_all = merged[,c("V1", "V2", "autism")]
pheno_under10 = merged[,c("V1", "V2", "under10")]
pheno_over10 = merged[,c("V1", "V2", "over10")]
pheno_under6 = merged[,c("V1", "V2", "under6")]


write.table(pheno_all, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_all.txt", row.names = F, col.names = F, quote = F)
write.table(pheno_under10, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_under10.txt", row.names = F, col.names = F, quote = F)
write.table(pheno_over10, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_over10.txt", row.names = F, col.names = F, quote = F)
write.table(pheno_under6, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_under6.txt", row.names = F, col.names = F, quote = F)

merged = merge(merged, pcs, by.x = "V2", by.y = "Sample_name")

covar = merged[,c("V1", "V2", "sex")]
qcovar = merged[,c("V1", "V2", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "age_at_registration_years")]

fam1 = fread("~/SPARK/SPARK_v4/postimputation/alpha/alpha_chr1_hg19_cleaned_v2.fam")
fam2 = fread("~/SPARK/SPARK_v4/postimputation/omega/omega_chr1_hg19_cleaned_v2.fam")

covar$batch = NA
covar$batch = ifelse(covar$V2 %in% fam1$V2, "alpha", "omega")

write.table(covar, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/covar.txt", row.names = F, col.names = F, quote = F)
write.table(qcovar, file = "/mnt/home4/arc/vw260/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/qcovar.txt", row.names = F, col.names = F, quote = F)


### Run GWAS
./gcta-1.94.1 --fastGWA-mlm-binary --mbfile ~/SPARK/SPARK_v4/GWAS/mbfile_agediagnosis_GWAS.txt --grm-sparse ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse --pheno ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_all.txt --qcovar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/qcovar.txt --covar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/covar.txt  --thread-num 20 --maf 0.01 --out ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/All_autism_SPARKreplication --exclude ~/SPARK/SPARK_v4/GWAS/batch_SNPs.txt
./gcta-1.94.1 --fastGWA-mlm-binary --mbfile ~/SPARK/SPARK_v4/GWAS/mbfile_agediagnosis_GWAS.txt --grm-sparse ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse --pheno ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_under10.txt --qcovar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/qcovar.txt --covar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/covar.txt  --thread-num 20 --maf 0.01 --out ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/under10_autism_SPARKreplication --exclude ~/SPARK/SPARK_v4/GWAS/batch_SNPs.txt
./gcta-1.94.1 --fastGWA-mlm-binary --mbfile ~/SPARK/SPARK_v4/GWAS/mbfile_agediagnosis_GWAS.txt --grm-sparse ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse --pheno ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_over10.txt --qcovar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/qcovar.txt --covar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/covar.txt  --thread-num 20 --maf 0.01 --out ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/over10_autism_SPARKreplication --exclude ~/SPARK/SPARK_v4/GWAS/batch_SNPs.txt
./gcta-1.94.1 --fastGWA-mlm-binary --mbfile ~/SPARK/SPARK_v4/GWAS/mbfile_agediagnosis_GWAS.txt --grm-sparse ~/SPARK/SPARK_v4/postimputation/sequencing_alpha_omega_hg19_allchrs_sparse --pheno ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/pheno_under6.txt --qcovar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/qcovar.txt --covar ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/covar.txt  --thread-num 20 --maf 0.01 --out ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/under6_autism_SPARKreplication --exclude ~/SPARK/SPARK_v4/GWAS/batch_SNPs.txt


### Meta-analyse the GWAS
./plink --meta-analysis ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/All_autism_SPARKreplication.fastGWA ~/SPARK/SPARK_v3/GWAS/Autism_SPARK_June2024/All_autism.fastGWA + qt report-all --out ./metaanalysis/SPARK_allautism --meta-analysis-bp-field POS
./plink --meta-analysis ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/under10_autism_SPARKreplication.fastGWA ~/SPARK/SPARK_v3/GWAS/Autism_SPARK_June2024/under10_autism.fastGWA + qt report-all --out ./metaanalysis/SPARK_under10 --meta-analysis-bp-field POS
./plink --meta-analysis ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/under6_autism_SPARKreplication.fastGWA ~/SPARK/SPARK_v3/GWAS/Autism_SPARK_June2024/under6_autism.fastGWA + qt report-all --out ./metaanalysis/SPARK_under6 --meta-analysis-bp-field POS
./plink --meta-analysis ~/SPARK/SPARK_v4/GWAS/Autism_SPARK_December2024/over10_autism_SPARKreplication.fastGWA ~/SPARK/SPARK_v3/GWAS/Autism_SPARK_June2024/over10_autism.fastGWA + qt report-all --out ./metaanalysis/SPARK_over10 --meta-analysis-bp-field POS


