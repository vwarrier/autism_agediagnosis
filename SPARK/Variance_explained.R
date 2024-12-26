#Assessing sources of variance and visualisation

library(data.table)
library(relaimpo)
core = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/core_descriptive_variables_2022-12-12.csv")
core_autism = subset(core, asd == "TRUE") # 123442

merged_total = subset(core_autism, age_at_registration_years < 22) # 109657, persons born in 1994 or afterwards after DSM -4 came out, and subsequently, the use of Asperger Syndrome

merged_total$diagnosis_age3 = scale(merged_total$age_onset_mos/12)

ggplot(merged_total, aes(x=diagnosis_age3, fill=sex)) +
  geom_density(alpha=0.2) +
  scale_fill_manual(values=c("green", "darkred")) +
  theme_minimal() +
  labs(
    title="Distribution of Age at Diagnosis by Sex",
    x="Age at Diagnosis",
    y="Density"
  ) +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position="top"
  )

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

data3 = merged_total[,c("subject_sp_id", "diagnosis_age3", "sex", "reported_cog_test_score_numeric", "cognitive_impairment_latest",  "scq_total_final_score",  "rbsr_total_final_score", "walked_age_mos",  "used_words_age_mos", "regress_lang_y_n",  "regress_other_y_n")]

clean_data <- na.omit(data3)

model <- lm(diagnosis_age3 ~ sex + reported_cog_test_score_numeric + cognitive_impairment_latest+ scq_total_final_score + rbsr_total_final_score + walked_age_mos + used_words_age_mos + regress_lang_y_n + regress_other_y_n, data = clean_data)

calc.relimp(model, type = "lmg", rela = TRUE)

#Now add other metrics

SES =  fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/background_history_child_2022-12-12.csv")
ethnicity = fread("/mnt/beegfs/home4/arc/vw260/SPARK/Phenotypes/V9/SPARK_collection_v9_2022-12-12/individuals_registration_2022-12-12.csv")
ethnicity$race <- NA  # Create new column with NA as default

# Fill in values based on conditions
ethnicity$race[ethnicity$race_asian == 1] <- 1
ethnicity$race[ethnicity$race_white == 1] <- 2
ethnicity$race[ethnicity$race_african_amer == 1] <- 3
ethnicity$race[ethnicity$race_native_amer == 1] <- 4
ethnicity$race[ethnicity$race_native_hawaiian == 1] <- 5
ethnicity$race[ethnicity$race_other == 1] <- 6
ethnicity$race[ethnicity$race_more_than_one_calc == 1] <- 7



SES2 = SES[,c("subject_sp_id", "mother_highest_education", "father_highest_education", "mother_occupation", "father_occupation", "annual_household_income")]

SES2$mother_education_numeric <- factor(SES2$mother_highest_education,
                                        levels = c("did_not_attend_high_school",
                                                   "some_high_school",
                                                   "ged_diploma",
                                                   "high_school_graduate",
                                                   "trade_school",
                                                   "some_college",
                                                   "associate_degree",
                                                   "baccalaureate_degree",
                                                   "graduate_or_professional_degree"),
                                        labels = 1:9)

# Convert to numeric
SES2$mother_education_numeric <- as.numeric(as.character(SES2$mother_education_numeric))

SES2$father_education_numeric <- factor(SES2$father_highest_education,
                                        levels = c("did_not_attend_high_school",
                                                   "some_high_school",
                                                   "ged_diploma",
                                                   "high_school_graduate",
                                                   "trade_school",
                                                   "some_college",
                                                   "associate_degree",
                                                   "baccalaureate_degree",
                                                   "graduate_or_professional_degree"),
                                        labels = 1:9)

# Convert to numeric
SES2$father_education_numeric <- as.numeric(as.character(SES2$father_education_numeric))

SES2$income_numeric <- factor(SES2$annual_household_income,
                              levels = c("less_than_20000",
                                         "21000_35000",
                                         "36000_50000",
                                         "51000_65000",
                                         "66000_80000",
                                         "81000_100000",
                                         "101000_130000",
                                         "131000_160000",
                                         "over_161000"),
                              labels = 1:9)

# Convert to numeric
SES2$income_numeric <- as.numeric(as.character(SES2$income_numeric))

SES3 = SES2[,c("subject_sp_id", "income_numeric", "mother_education_numeric", "father_education_numeric")]
ethnicity2 = ethnicity[,c("subject_sp_id", "race")]

SES_ethnicity = merge(SES3, ethnicity2, by = "subject_sp_id")
clean_data2 = merge(clean_data, SES_ethnicity, by = "subject_sp_id")


model <- lm(diagnosis_age3 ~ race + income_numeric + mother_education_numeric + father_education_numeric + sex + reported_cog_test_score_numeric + cognitive_impairment_latest+ scq_total_final_score + rbsr_total_final_score + walked_age_mos + used_words_age_mos + regress_lang_y_n + regress_other_y_n, data = clean_data2)

calc.relimp(model, type = "lmg", rela = TRUE)
