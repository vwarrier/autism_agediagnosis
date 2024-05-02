###pTDT analyses 
library("tidyr")
library("reshape2")
library("dplyr")
library("ggplot2")
library("data.table")

prs_autism = fread("~/Autism_heterogeneity/PGS/SPARKv3_alphaomega_ADHD_finalscore.profile", header = TRUE)
setnames(prs_autism, "SCORE", "autism_prs")
#setnames(prs_autism, old = c("#FID", "SCORE1_AVG"), new = c("FID", "autism_prs"))
merged = prs_autism[,c("IID", "autism_prs", "FID")]
merged = unique(merged)

registration = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")

trios =  registration[!(registration$biofather_id =="" | registration$biomother_id==""), ] # 65143
asd_trios = subset(trios, asd == "TRUE") #44236


merged_pgs_trio = merge(merged, asd_trios, by.x = "IID", by.y = "subject_sp_id") #11763, 3638
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biomother_id %in% merged$IID,] #8312, 2980
merged_pgs_trio = merged_pgs_trio[merged_pgs_trio$biofather_id %in% merged$IID,] #6956, 2477


fatherpgs = merged[merged$IID %in% asd_trios$biofather_id,] # 7763
fatherpgs = fatherpgs[,c("IID", "FID", "autism_prs")]
motherpgs = merged[merged$IID %in% asd_trios$biomother_id,] # 8430
motherpgs = motherpgs[,c("IID", "FID","autism_prs")]

casespgs = merged_pgs_trio #6956, 2477

parentpgs = merge(motherpgs, fatherpgs, by = "FID") #6473
parentpgs$midparent_autism = (parentpgs$autism_prs.x + parentpgs$autism_prs.y)/2

triopgs = merge(parentpgs, casespgs, by = "FID") # 6981
triopgs = triopgs[!duplicated(triopgs[,c("IID")]),] # 6513 (alpha_omega), 2477 (theta_tau), total = 6513 + 2477 = 8,990


Sd_autism = sd(triopgs$midparent_autism) 

triopgs$autism_diff = (triopgs$autism_prs - triopgs$midparent_autism)/Sd_autism

#PTDT overtransmission analyses

N = sqrt(nrow(triopgs))
Z = mean(triopgs$autism_diff)/(sd(triopgs$autism_diff)/N)
mean(triopgs$autism_diff)
(sd(triopgs$autism_diff)/N)
2*pnorm(-abs(Z))

males = subset(triopgs, sex == "Male")
N = sqrt(nrow(males))
Z = mean(males$autism_diff)/(sd(males$autism_diff)/N)
mean(males$autism_diff)
(sd(males$autism_diff)/N)
2*pnorm(-abs(Z))

females = subset(triopgs, sex == "Female")
N = sqrt(nrow(females))
Z = mean(females$autism_diff)/(sd(females$autism_diff)/N)
mean(females$autism_diff)
(sd(females$autism_diff)/N)
2*pnorm(-abs(Z))


###subset to no_ID/ID
noID = subset(triopgs, reported_cog_test_score == "110 - 119" | reported_cog_test_score == "120 - 129" | reported_cog_test_score == "130 - above" | reported_cog_test_score == "70 - 79" | reported_cog_test_score == "80 - 89"  | reported_cog_test_score == "90 - 109")
t.test(noID$autism_diff ~ noID$sex)

ID = subset(triopgs, reported_cog_test_score == "24 - below" | reported_cog_test_score == "25 - 39" | reported_cog_test_score == "40 - 54" | reported_cog_test_score == "55 - 69")
t.test(ID$autism_diff ~ ID$sex)

N = sqrt(nrow(ID))
Z = mean(ID$autism_diff)/(sd(ID$autism_diff)/N)
mean(ID$autism_diff)
(sd(ID$autism_diff)/N)
2*pnorm(-abs(Z))

N = sqrt(nrow(noID))
Z = mean(noID$autism_diff)/(sd(noID$autism_diff)/N)
mean(noID$autism_diff)
(sd(noID$autism_diff)/N)
2*pnorm(-abs(Z))


males = subset(noID, sex == "Male")
females = subset(noID, sex == "Female")

N = sqrt(nrow(males))
Z = mean(males$autism_diff)/(sd(males$autism_diff)/N)
mean(males$autism_diff)
(sd(males$autism_diff)/N)
2*pnorm(-abs(Z))


N = sqrt(nrow(females))
Z = mean(females$autism_diff)/(sd(females$autism_diff)/N)
mean(females$autism_diff)
(sd(females$autism_diff)/N)
2*pnorm(-abs(Z))

males = subset(ID, sex == "Male")
females = subset(ID, sex == "Female")

N = sqrt(nrow(males))
Z = mean(males$autism_diff)/(sd(males$autism_diff)/N)
mean(males$autism_diff)
(sd(males$autism_diff)/N)
2*pnorm(-abs(Z))


N = sqrt(nrow(females))
Z = mean(females$autism_diff)/(sd(females$autism_diff)/N)
mean(females$autism_diff)
(sd(females$autism_diff)/N)
2*pnorm(-abs(Z))


#Analyses with age at diagnosis

updated_trio = subset(triopgs, age_at_registration_years < 22)

summary(lm(scale(diagnosis_age/12) ~ scale(autism_diff) + scale(midparent_autism) + sex + age_at_registration_years + cognitive_impairment_latest, data = updated_trio))




