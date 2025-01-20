library(GenomicSEM)
library(lavaan)

# Set working directory
setwd("/mnt/beegfs6/home4/arc/vw260/ldsc/")
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"


traits <- c("./sumstats/Autism_PGC2_2017.sumstats.gz", "./sumstats/Autism_iPSYCHonly_Warrier2022.sumstats.gz", "./sumstats/ADHD_demontis_2022.sumstats.gz", "./sumstats/scz3.sumstats", "./sumstats/wraymdd.sumstats", "./sumstats/PTSD_2024_natgen.sumstats.gz", "./sumstats/anorexia2019.sumstats.gz", "./sumstats/BIP_2021.sumstats", "./sumstats/Finngen_r10_autism.sumstats.gz", "./sumstats/SPARK_over10_meta.sumstats", "./sumstats/SPARK_under6_meta.sumstats", "./sumstats/autism_ipsych_over10.sumstats.gz", "./sumstats/Autism_iPSYCH_under9.sumstats.gz")
sample.prev <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
population.prev <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
trait.names<-c("pgc", "ipsych", "adhd", "scz", "mdd", "ptsd", "anorexia", "bipolar", "finngen", "SPARK_over10", "SPARK_under6", "iPSYCH_over10", "iPSYCH_under9")
GWIS<-ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
save(GWIS, file="Autism_MH_LDSC_2024.RData")


GWIS <- LDSCoutput

#specify a model in which Educational Achievement is regressed on Bipolar Disorder and Schizophrenia, and the correlation is modeled #between Bipolar and Schizophrenia


baseline_pgc <- 'ipsych ~ pgc + adhd + mdd + ptsd + scz + anorexia + bipolar '
usermodel(GWIS, estimation = "DWLS", model = baseline_pgc)


baseline_spark <- 'ipsych ~ adhd + mdd + ptsd + scz + anorexia + bipolar + SPARK_under6'
usermodel(GWIS, estimation = "DWLS", model = baseline_spark)


baseline_adultspark <- 'SPARK_over10 ~ SPARK_under6 + adhd + mdd + ptsd + scz + anorexia + bipolar'
usermodel(GWIS, estimation = "DWLS", model = baseline_adultspark)


full_pgc <- 'ipsych ~ pgc + adhd + mdd + ptsd + scz + anorexia + bipolar + SPARK_over10 '
usermodel(GWIS, estimation = "DWLS", model = full_pgc)


full_spark <- 'ipsych ~ adhd + mdd + ptsd + scz + anorexia + bipolar + SPARK_under6 + SPARK_over10'
usermodel(GWIS, estimation = "DWLS", model = full_spark)


full_adultspark <- 'SPARK_over10 ~ SPARK_under6 + adhd + mdd + ptsd + scz + anorexia + bipolar + finngen'
usermodel(GWIS, estimation = "DWLS", model = full_adultspark)


full_adultiPSYCH <- 'iPSYCH_over10 ~ SPARK_under6 + adhd + mdd + ptsd + scz + anorexia + bipolar + SPARK_over10'
usermodel(GWIS, estimation = "DWLS", model = full_adultiPSYCH)

full_iPSYCH<- 'ipsych ~ SPARK_under6 + SPARK_over10'
usermodel(GWIS, estimation = "DWLS", model = full_iPSYCH)

