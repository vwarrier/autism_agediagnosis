
library(GenomicSEM)

setwd("/mnt/beegfs6/home4/arc/vw260/ldsc/sumstats")
list1 = c("ADHD_demontis_2022.sumstats.gz", "anorexia2019.sumstats.gz","MDD_als.sumstats.gz", "scz3.sumstats.gz", "BIP_2021.sumstats", "PTSD_2024_natgen.sumstats.gz", "anxietycc.sumstats", "childhoodtrauma_retro_prosp_meta.sumstats", "selfharmscoreUKB.sumstats.gz",  "savage_intelligence.sumstats.gz", "leeedu2018.sumstats.gz", "Autismagediagnosis_linear.sumstats.gz", "asd_agediagnosis_iPSYCH.sumstats.gz", "ASDf.sumstats", "ASDm.sumstats", "MoBa_impulsivity_age3.sumstats.gz", "MoBa_inattention_age3.sumstats.gz", "MoBa_language_age3.sumstats.gz", "MoBa_motor_age3.sumstats.gz", "MoBa_NVcom_age3.sumstats.gz", "MoBa_play_age3.sumstats.gz", "MoBa_prosocial_age3.sumstats.gz", "MoBa_RepBehavior_age3.sumstats.gz", "MoBa_RepSpeech_age3.sumstats.gz", "MoBa_SocialAtt_age3.sumstats.gz", "MoBa_waiting_age3.sumstats.gz", "WALK_gwas.sumstats.gz", "Latephase_EV_2024.sumstats.gz", "Latephase_RV_2024.sumstats.gz", "Earlyphase_EV_2024.sumstats.gz", "NDD_meta_HCM2024.sumstats.gz")

ED = NULL
LD = NULL
for (i in list1){
  traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "Autism_iPSYCH_over12.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under6.sumstats", "SPARK_over10.sumstats", paste0(i))
  sample.prev <- c(NA, NA, NA, NA, NA, NA, NA)
  population.prev <- c(NA, NA, NA, NA, NA, NA, NA)
  ld <- "~/ldsc/eur_w_ld_chr/"
  wld <- "~/ldsc/eur_w_ld_chr/"
  trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over12", "FinnGen",  "SPARK_under6",  "SPARK_over10", "T")
  LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
  
  modelrg<-'ED=~NA*PGC + SPARK_under6 + iPSYCH_under9
  LD=~NA*iPSYCH_over12 + iPSYCH_under9 + FinnGen + SPARK_over10
  LT=~NA*T
  LD~~1*LD
  ED~~1*ED
  LT~~1*LT
  T~~0*T
  ED~~LD
  LD~~0*T
  ED~~0*T
  PGC ~~ a*PGC
  a > .001
  iPSYCH_over12 ~~ b*iPSYCH_over12
  b > .001'
  
  Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = modelrg, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
  
  
  ED = rbind(ED, Autism_CFA$results[16,])
  LD = rbind(LD, Autism_CFA$results[17,])
  
}




setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("dep_anx1.sumstats.gz", "dep_anx2.sumstats.gz", "dep_anx3.sumstats.gz")
sample.prev <- c(NA, NA, NA)
population.prev <- c(NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("depression",  "anxiety", "comorbid")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)


GWISmodel <- 'comorbid ~ depression + anxiety
depression~~anxiety'

#run the model using the usermodel function
GWISoutput<-usermodel(LDSCoutput, estimation = "DWLS", model = GWISmodel, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

#print the output
GWISoutput
