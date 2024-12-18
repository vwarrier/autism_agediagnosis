#Scripts for genomicSEM
###Scripts for global EFA and CFA using genomic SEM - preprint

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

CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

#Two correlated factor
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + iPSYCH_under9
             F2 =~ NA*iPSYCH_over12  + FinnGen + SPARK_over10 + iPSYCH_under9
            F1~~F2
            PGC ~~ a*PGC
            a > .001
            iPSYCH_over12 ~~ b*iPSYCH_over12
            b > .001"

Autism_CFA_cor<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)



#Two correlated factor based on geography
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + SPARK_over10
             F2 =~ NA*iPSYCH_over12  + FinnGen + iPSYCH_under9
            F1~~F2
            iPSYCH_under9 ~~ b*iPSYCH_under9
            b > 0.001"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)




#Bifactor - geography - Model not identified
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 + SPARK_over10
             F2 =~ NA*iPSYCH_over12  + FinnGen + iPSYCH_under9
             F3 =~ NA*PGC + SPARK_under6 + SPARK_over10 + iPSYCH_over12  + FinnGen + iPSYCH_under9
            F1~~0*F2
            F3 ~~ 0*F1
            F3 ~~ 0*F2"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA 

#Bifactor - age
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 +  iPSYCH_under9
             F2 =~ NA*iPSYCH_over12  + FinnGen + SPARK_over10
             F3 =~ NA*PGC + SPARK_under6 + SPARK_over10 + iPSYCH_over12  + FinnGen + iPSYCH_under9
            F3 ~~ 0*F1
            F3 ~~ 0*F2
            F1~~0*F2
            iPSYCH_under9 ~~ b*iPSYCH_under9
            b > 0.001
            PGC ~~ a*PGC
            a > 0.001
            iPSYCH_over12 ~~ c*iPSYCH_over12
            c > 0.001"

Autism_CFA_bifacotr<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

Autism_CFA 



#Heirarchical factor - age
CFAofEFA <- "F1 =~ NA*PGC + SPARK_under6 +  iPSYCH_under9
             F2 =~ NA*iPSYCH_over12  + FinnGen + SPARK_over10
             F3 =~ F1 + F2
            iPSYCH_under9 ~~ b*iPSYCH_under9
            b > 0.001
            PGC ~~ a*PGC
            a > 0.001
            iPSYCH_over12 ~~ c*iPSYCH_over12
            c > 0.001"

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Autism_CFA 


#####New - with the SPARK meta-analysed GWAS - just change the sumstats as follows and rerun everything above

setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "Autism_iPSYCH_over12.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under6_meta.sumstats", "SPARK_over10_meta.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over12", "FinnGen",  "SPARK_under6",  "SPARK_over10")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
##optional command to save the ldsc output in case you want to use it in a later R session. 
save(LDSCoutput, file="Autism_all_nonoverlapping_GWAS_withSPARK_extremevalue_meta.RData")

