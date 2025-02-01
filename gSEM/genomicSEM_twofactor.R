#Scripts for genomicSEM
###Scripts for global EFA and CFA using genomic SEM - preprint

setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "autism_ipsych_over10.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under6_meta.sumstats", "SPARK_over10_meta.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over10", "FinnGen",  "SPARK_under6",  "SPARK_over10")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
##optional command to save the ldsc output in case you want to use it in a later R session. 
save(LDSCoutput, file="Autism_all_nonoverlapping_GWAS_withSPARK_meta_final.RData")


CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS

#Two correlated factor
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + iPSYCH_under9
             F2 =~ NA*iPSYCH_over10  + FinnGen + SPARK_over10 + iPSYCH_under9
            F1~~F2
            iPSYCH_over10 ~~ a*iPSYCH_over10
            a > .001"

Autism_CFA_cor<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


Autism_CFA_cor

#Two correlated factor based on geography
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + SPARK_over10
             F2 =~ NA*iPSYCH_over10  + FinnGen + iPSYCH_under9
            F1~~F2
            iPSYCH_under9 ~~ b*iPSYCH_under9
            b > 0.001"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA


#Bifactor - geography - Model not identified
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + SPARK_over10
             F2 =~ NA*iPSYCH_over10  + FinnGen + iPSYCH_under9
             F3 =~ NA*PGC + SPARK_under6 + SPARK_over10 + iPSYCH_over10  + FinnGen + iPSYCH_under9
            F1~~0*F2
            F3 ~~ 0*F1
            F3 ~~ 0*F2
            iPSYCH_under9 ~~ a*iPSYCH_under9
            a > 0.001"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA 

#Bifactor - age
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 +  iPSYCH_under9
             F2 =~ NA*iPSYCH_over10  + FinnGen + SPARK_over10
             F3 =~ NA*PGC + SPARK_under6 + SPARK_over10 + iPSYCH_over10  + FinnGen + iPSYCH_under9
            F3 ~~ 0*F1
            F3 ~~ 0*F2
            F1~~0*F2
            iPSYCH_under9 ~~ b*iPSYCH_under9
            b > 0.001
            PGC ~~ a*PGC
            a > 0.001
            iPSYCH_over10 ~~ c*iPSYCH_over10
            c > 0.001
            SPARK_over10 ~~ d*SPARK_over10
            d > 0.001"

Autism_CFA_bifacotr<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA_bifacotr 




###Sensitivity analyses 1 - without FinnGen

setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "autism_ipsych_over10.sumstats.gz",  "SPARK_under6_meta.sumstats", "SPARK_over10_meta.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over10", "SPARK_under6",  "SPARK_over10")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)


CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS

#Two correlated factor
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + iPSYCH_under9
             F2 =~ NA*iPSYCH_over10  + SPARK_over10 + iPSYCH_under9
            F1~~F2
            iPSYCH_over10 ~~ a*iPSYCH_over10
            a > .001"

Autism_CFA_cor<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


Autism_CFA_cor


###Sensitivity analyses 1 - less stringent difference

 setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "autism_ipsych_under11.sumstats.gz", "autism_ipsych_over10.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under10_meta.sumstats", "SPARK_over10_meta.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under11", "iPSYCH_over10", "FinnGen",  "SPARK_under11",  "SPARK_over10")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)



CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS

#Two correlated factor
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under11 + iPSYCH_under11 
             F2 =~ NA*iPSYCH_over10  + FinnGen + SPARK_over10 + iPSYCH_under11 + SPARK_under11
            F1~~F2
          iPSYCH_over10 ~~ a*iPSYCH_over10
            a > .001"

Autism_CFA_cor<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Autism_CFA_cor


