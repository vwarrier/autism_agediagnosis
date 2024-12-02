#Example script. Modify as needed

./gcta64
--threads 15
--reml
--grm /mnt/home4/arc/vw260/SPARK/SPARK_v3/Postimputation_genotyping/omega/alpha_omega_hg19_autismonly_chr
--pheno /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_pheno.txt
--grm-cutoff0.05
--qcovar /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_qcovar_devmilestones.txt
--covar /mnt/home4/arc/vw260/SPARK/SPARK_v3/GWAS/age_diagnosis_covar.txt
--out /mnt/home4/arc/vw260/Autism_heterogeneity/Age_diagnosis_analyses/GCTA_results/h2_autismage_devmilestones
