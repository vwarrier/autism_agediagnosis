#GenomicSEM_basedgeneticcor

library(data.table)
library(GenomicSEM)

#Genetic correlation - genomicSEM

traits <- c("Autism_PGC2_2017.sumstats.gz",  "SPARK_matoba_multiancestry.sumstats.gz", "autism_ipsych_under11.sumstats.gz", "autism_ipsych_over10.sumstats.gz", "Finngen_r10_autism.sumstats.gz", "leeea2018.sumstats.gz")
sample.prev <- c(0.47, 0.5, 0.2, 0.18, 0.002, NA)
population.prev <- c(0.027, 0.027, 0.027, 0.027, 0.027, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC", "SPARK", "Under11", "Over10", "FinnGen", "T")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
##optional command to save the ldsc output in case you want to use it in a later R session. 


modelrg<-'ED=~NA*PGC + SPARK + Under11
  LD=~NA*Over10 + Under11 + FinnGen
  LT=~NA*T
  LD~~1*LD
  ED~~1*ED
  LT~~1*LT
  T~~0*T
  ED~~LD
  LD~~0*T
  ED~~0*T'

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = modelrg, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA

setwd("/mnt/beegfs/home4/arc/vw260/ldsc")

#Genetic correlation - GWAS by subtraction

traits <- c("./sumstats/autism_ipsych_over10.sumstats.gz","./sumstats/PTSD_eur_2019.sumstats","./sumstats/anorexia2019.sumstats.gz")
sample.prev <- c(NA,NA,NA)
population.prev <- c(NA,NA,NA)
ld<-"eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
trait.names<-c("Autism", "ADHD", "T")

LDSCoutput <- ldsc(traits, 
                   sample.prev, 
                   population.prev, 
                   ld, 
                   wld, 
                   trait.names)



modelrg<-'add=~NA*Autism + start(0.4)*ADHD
  aut=~NA*Autism
  LT=~NA*T
  add~~1*add
  aut~~1*aut
  LT~~1*LT
  Autism~~0*Autism
  ADHD~~0*ADHD
  T~~0*T
  add~~0*aut
  Autism~~0*ADHD
  Autism~~0*T
  ADHD~~0*T
  '
output<-usermodel(LDSCoutput,estimation="DWLS",model=modelrg)
output
