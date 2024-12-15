library(genomicSEM)
library(data.table)



setwd("/mnt/beegfs/home4/arc/vw260/ldsc/sumstats")
traits <- c("Autism_PGC2_2017.sumstats.gz", "Autism_iPSYCH_under9.sumstats.gz", "autism_ipsych_over10.sumstats.gz", "Finngen_r10_autism.sumstats.gz",  "SPARK_under6.sumstats", "SPARK_over10.sumstats", "Menarche.sumstats")
sample.prev <- c(NA, NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA, NA)
ld <- "~/ldsc/eur_w_ld_chr/"
wld <- "~/ldsc/eur_w_ld_chr/"
trait.names<-c("PGC",  "iPSYCH_under9", "iPSYCH_over10", "FinnGen",  "SPARK_under6",  "SPARK_over10", "T")
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)



modelrg<-'ED =~ NA*PGC + SPARK_under6 + iPSYCH_under9 
            LD =~ NA*iPSYCH_over10  + FinnGen + SPARK_over10 + iPSYCH_under9
            ED~~LD
            iPSYCH_over10 ~~ b*iPSYCH_over10
            b > .001
            PGC ~~ a*PGC
            a > .001
  LT=~NA*T
  LD~~1*LD
  ED~~1*ED
  LT~~1*LT
  T~~0*T
  LD~~0*T
  ED~~0*T'

Autism_CFA<-usermodel(LDSCoutput, estimation = "DWLS", model = modelrg, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


ED = Autism_CFA$results[16,] 

LD = Autism_CFA$results[17,] 



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
