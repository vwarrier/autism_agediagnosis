#Scripts for genomicSEM


###Scripts for global EFA and CFA using genomic SEM

library(data.table)
require(GenomicSEM)
library(dplyr)


setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "Autism_iPSYCH_over12.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under6.sumstats", "SPARK_over10.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over12", "FinnGen",  "SPARK_under6",  "SPARK_over10")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
##optional command to save the ldsc output in case you want to use it in a later R session. 
save(LDSCoutput, file="Autism_all_nonoverlapping_GWAS_withSPARK_extremevalue.RData")



#Two correlated factor
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + iPSYCH_under9
             F2 =~ NA*iPSYCH_over12  + FinnGen + SPARK_over10 + iPSYCH_under9
            F1~~F2
            PGC ~~ a*PGC
            a > .001
            iPSYCH_over12 ~~ b*iPSYCH_over12
            b > .001"

Autism_CFA_cor<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)




