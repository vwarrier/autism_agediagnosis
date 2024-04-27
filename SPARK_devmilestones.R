#Developmental milestones SPARK -  PGS script

#Install packages 
library("tidyr")
library("reshape2")
library("dplyr")
library("ggplot2")
library("data.table")
library(purrr)

#1 PGS in Autistic individuals
data1 = fread("~/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/background_history_child_2022-12-12.csv")

setnames(data1, 1, "IID")
setnames(data1, 2, "FID")

data1$walked_age_mos = ifelse(data1$walked_age_mos == "888", NA, data1$walked_age_mos)
data1$used_words_age_mos = ifelse(data1$used_words_age_mos == "888", NA, data1$used_words_age_mos)
data1$bowel_trained_age_mos = ifelse(data1$bowel_trained_age_mos == "888", NA, data1$bowel_trained_age_mos)
data1$bladder_trained_age_mos = ifelse(data1$bladder_trained_age_mos == "888", NA, data1$bladder_trained_age_mos)
data1$smiled_age_mos = ifelse(data1$bladder_trained_age_mos == "888", NA, data1$smiled_age_mos)
data1$sat_wo_support_age_mos = ifelse(data1$sat_wo_support_age_mos == "888", NA, data1$sat_wo_support_age_mos)
data1$crawled_age_mos = ifelse(data1$crawled_age_mos == "888", NA, data1$crawled_age_mos)
data1$combined_phrases_age_mos = ifelse(data1$combined_phrases_age_mos == "888", NA, data1$combined_phrases_age_mos)
data1$fed_self_spoon_age_mos = ifelse(data1$fed_self_spoon_age_mos == "888", NA, data1$fed_self_spoon_age_mos)

over10 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHover10_finalscore.profile")
over10 = over10[,c("FID", "IID", "SCORE")]
setnames(over10, "SCORE", "over10_PGS")

under11 = fread("~/Autism_heterogeneity/PGS/SPARK_alphaomega_iPSYCHunder11_finalscore.profile")
under11 = under11[,c("IID", "SCORE")]
setnames(under11, "SCORE", "under11_PGS")


df_list <- list( over10, under11, data1)
merged = df_list %>% reduce(full_join, by='IID')

autism = subset(merged, asd == "TRUE")

PCs = fread("~/SPARK/SPARK_v3/Preimputation_genotype/SPARKalphaomega_PCsforGWAS_v2.txt") # 47170
setnames(PCs, "Sample_name", "IID")
merged2 = merge(autism, PCs, by = "IID") 


list1 = c("walked_age_mos", "used_words_age_mos", "bowel_trained_age_mos", "bladder_trained_age_mos", "smiled_age_mos", "sat_wo_support_age_mos", "crawled_age_mos", "combined_phrases_age_mos", "fed_self_spoon_age_mos")

replace_outliers <- function(x) {
  median_x <- median(x, na.rm = TRUE)
  mad_x <- mad(x, na.rm = TRUE)
  threshold <- median_x + 5 * mad_x
  x[x > threshold] <- NA
  return(x)
}

merged2 <- merged2 %>%
  mutate(across(all_of(list1), replace_outliers))

# Loop through each data frame in the list

results_over10 = NULL
results_under11 = NULL
Samplesize = NULL


for(i in list1){
  results_all = summary(lm(scale(merged2[[i]]) ~ scale(over10_PGS) + scale(under11_PGS) + sex + age_at_eval_years + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = merged2))
  results_over10 = rbind(results_over10, cbind(i, t(results_all$coefficients[2,])))
  results_under11 = rbind(results_under11, cbind(i, t(results_all$coefficients[3,])))
  dim =  nrow(merged2) - length(results_all$na.action)
  Samplesize = rbind(Samplesize, dim)
}
