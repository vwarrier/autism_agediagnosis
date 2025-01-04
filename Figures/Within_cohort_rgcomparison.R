# Create separate dataframes for iPSYCH and SPARK
iPSYCH_data <- data.frame(
  Phenotype = factor(rep(c("ADHD", "Bipolar disorder", "Educational attainment", "PTSD", 
                           "Schizophrenia", "Depression", "Anorexia", "Cognitive aptitude",
                           "Childhood maltreatment", "Self harm score", "Anxiety"), each=2),
                     levels = c("ADHD", "Bipolar disorder", "Educational attainment", "PTSD", 
                                "Schizophrenia", "Depression", "Anorexia", "Cognitive aptitude",
                                "Childhood maltreatment", "Self harm score", "Anxiety")),
  Timing = rep(c("early", "late"), times=11),
  rg = c(0.37, 0.52, 0.25, 0.28, 0.24, 0, 0.27, 0.54, 0.20, 0.30, 
         0.26, 0.52, 0.19, 0.21, 0.22, 0.15, 0.41, 0.56, 0.37, 0.59, 0.15, 0.39),
  SE = c(0.05, 0.04, 0.05, 0.05, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04,
         0.04, 0.04, 0.06, 0.06, 0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.15, 0.17),
  zdiff = rep(c(-2.34, -0.42, 5.66, -4.77, -1.77, -4.60, -0.24, 1.24, -2.12, -2.59, -1.06), each=2),
  Cohort = "iPSYCH before 9/iPSYCH after 10"
)

SPARK_data <- data.frame(
  Phenotype = factor(rep(c("ADHD", "Bipolar disorder", "Educational attainment", "PTSD", 
                           "Schizophrenia", "Depression", "Anorexia", "Cognitive aptitude",
                           "Childhood maltreatment", "Self harm score", "Anxiety"), each=2),
                     levels = c("ADHD", "Bipolar disorder", "Educational attainment", "PTSD", 
                                "Schizophrenia", "Depression", "Anorexia", "Cognitive aptitude",
                                "Childhood maltreatment", "Self harm score", "Anxiety")),
  Timing = rep(c("early", "late"), times=11),
  rg = c(0.2218, 0.4593, 0.147, 0.3154, 0.1366, 0.2929, 0.1755, 0.4983, 
         0.2632, 0.1885, 0.1645, 0.5569, 0.143, 0.2128, 0.0606, 0.4519, 
         0.2467, 0.629, 0.1676, 0.7481, 0.2326, 0.1261),
  SE = c(0.0594, 0.0838, 0.053, 0.0672, 0.0459, 0.0547, 0.0476, 0.0837,
         0.0494, 0.0616, 0.0439, 0.0877, 0.0658, 0.0825, 0.0455, 0.0768,
         0.0595, 0.1056, 0.0702, 0.1253, 0.1885, 0.2174),
  zdiff = rep(c(-2.31, -2.01, -2.19, -3.35, 0.95, -4.00, -0.66, -4.38, -3.15, -4.04, 0.37), each=2),
  Cohort = "SPARK before 6/SPARK after 10"
)

# Combine the data
plot_data <- rbind(iPSYCH_data, SPARK_data)

# Create the plot
ggplot(plot_data, aes(x = rg, y = factor(Phenotype, levels = rev(levels(Phenotype))), color = Timing)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = rg - 1.96*SE, xmax = rg + 1.96*SE),
                 position = position_dodge(width = 0.5),
                 height = 0.2) +
  facet_wrap(~Cohort) +
  scale_color_manual(values = c("early" = "#36A047", "late" = "#6B3F98"),
                     labels = c("Early diagnosed", "Late diagnosed")) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(y = "Phenotype",
       x = "Genetic Correlation",
       color = "Diagnosis Timing")




# Create dataframe for late diagnosed comparison
late_comparison <- data.frame(
  Phenotype = c("ADHD", "Bipolar disorder", "Educational attainment", "PTSD", 
                "Schizophrenia", "Depression", "Anorexia", "Cognitive aptitude",
                "Childhood maltreatment", "Self harm score", "Anxiety"),
  SPARK_rg = c(0.4593, 0.3154, 0.2929, 0.4983, 0.1885, 0.5569, 0.2128, 0.4519, 
               0.629, 0.7481, 0.1261),
  SPARK_SE = c(0.0838, 0.0672, 0.0547, 0.0837, 0.0616, 0.0877, 0.0825, 0.0768, 
               0.1056, 0.1253, 0.2174),
  iPSYCH_rg = c(0.52, 0.28, 0, 0.54, 0.30, 0.52, 0.21, 0.15, 0.56, 0.59, 0.39),
  iPSYCH_SE = c(0.04, 0.05, 0.03, 0.04, 0.04, 0.04, 0.06, 0.04, 0.05, 0.06, 0.17)
)

# Create the plot
ggplot(late_comparison, aes(x = SPARK_rg, y = iPSYCH_rg)) +
  geom_point() +
  geom_errorbar(aes(ymin = iPSYCH_rg - 1.96*iPSYCH_SE, 
                    ymax = iPSYCH_rg + 1.96*iPSYCH_SE)) +
  geom_errorbarh(aes(xmin = SPARK_rg - 1.96*SPARK_SE, 
                     xmax = SPARK_rg + 1.96*SPARK_SE)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = Phenotype), hjust = -0.1, vjust = -0.5, size = 3) +
  theme_classic() +
  labs(x = "SPARK after 10",
       y = "iPSYCH after 10") +
  coord_equal(xlim = c(-0.2, 1), ylim = c(-0.2, 1))
