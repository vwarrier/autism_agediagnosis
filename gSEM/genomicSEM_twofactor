#Scripts for genomicSEM


###Scripts for global EFA and CFA using genomic SEM

library(data.table)
require(GenomicSEM)
library(dplyr)


setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz",  "SPARK_matoba_multiancestry.sumstats.gz", "autism_ipsych_under11.sumstats.gz", "autism_ipsych_over10.sumstats.gz", "Finngen_r10")
sample.prev <- c(0.47, 0.5, 0.2, 0.18, 0.002)
population.prev <- c(0.027, 0.027, 0.027, 0.027, 0.027)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC", "SPARK", "Under11", "Over10", "FinnGen")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
##optional command to save the ldsc output in case you want to use it in a later R session. 
save(LDSCoutput, file="Autism_all.RData")

CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")


CFAofEFA <- "F1 =~ NA*PGC + SPARK + Under11
             F2 =~ NA*Over10 + Under11 + FinnGen
            F1~~F2"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)





