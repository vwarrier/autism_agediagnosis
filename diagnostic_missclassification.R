###Script for checking the impact of diagnostic missclassification

#Scripts generated from doi: 10.1038/ejhg.2011.257.(Wray et al., 2012, EJHG).

# Input parameters
N_caseTA =  # Number of true cases for disorder A
N_controlA =  # Number of controls for disorder A
N_caseTB =  # Number of true cases for disorder B
N_controlB =  # Number of controls for disorder B
M_TA =  # Proportion of disorder A cases misclassified as disorder B
M_TB =  # Proportion of disorder B cases misclassified as disorder A
SNP_h2_A =  # Observed scale true SNP heritability of disorder A
SNP_h2_B =  # Observed scale true SNP heritability of disorder B

#Calculate the genetic variance for binary GWAS when SNP h2 is available (observed scale)
calculate_genetic_variance <- function(h2_snp, N_cases, N_controls) {
  P <- N_cases / (N_cases + N_controls)
  s2_u <- h2_snp * P * (1 - P)
  return(s2_u)
}



s2_uTA <- calculate_genetic_variance(SNP_h2_A, N_caseTA, N_caseTA) # True genetic variance for disorder A
s2_uTB <- calculate_genetic_variance(SNP_h2_B, N_caseTB, N_caseTB) # True genetic variance for disorder B

# Calculate genetic covariance when genetic correlation and SNP heritabilities (in observed scale) are available
rg = 
s_uTA_uTB <- rg*(sqrt(s2_uTA) * sqrt(s2_uTB))


# Calculate parameters
N_caseDA <- (1 - M_TA) * N_caseTA + M_TB * N_caseTB
N_caseDB <- (1 - M_TB) * N_caseTB + M_TA * N_caseTA
P_DA <- N_caseDA / (N_caseDA + N_controlA)
P_DB <- N_caseDB / (N_caseDB + N_controlB)

s2_uDA <- ((1 - M_TA)^2 * N_caseTA^2 * s2_uTA + M_TB^2 * N_caseTB^2 * s2_uTB + 2 * (1 - M_TA) * M_TB * N_caseTA * N_caseTB * s_uTA_uTB) / N_caseDA^2
s2_uDB <- ((1 - M_TB)^2 * N_caseTB^2 * s2_uTB + M_TA^2 * N_caseTA^2 * s2_uTA + 2 * (1 - M_TB) * M_TA * N_caseTA * N_caseTB * s_uTA_uTB) / N_caseDB^2
s_uDA_uDB <- ((((1 - M_TA) * (1 - M_TB) + M_TA * M_TB) * N_caseTA * N_caseTB * s_uTA_uTB + (1 - M_TA) * M_TA * N_caseTA^2 * s2_uTA + (1 - M_TB) * M_TB * N_caseTB^2 * s2_uTB) / (N_caseDA * N_caseDB))

# Heritability estimates
h2_DA <- s2_uDA / (P_DA * (1 - P_DA))
h2_DB <- s2_uDB / (P_DB * (1 - P_DB))

# Genetic correlation estimate
rg_D <- s_uDA_uDB / (sqrt(s2_uDA) * sqrt(s2_uDB))

# Print results
cat("Heritability for disorder A (diagnosed):", h2_DA, "\n")
cat("Heritability for disorder B (diagnosed):", h2_DB, "\n")
cat("Genetic correlation between disorders A and B (diagnosed):", rg_D, "\n")
