library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

# Create phenotype levels with wrapped text
phenotype_levels <- rev(c(
  "Sex and 10\ngenetic PCs",
  "Sex, ID and 10\ngenetic PCs (baseline)",
  "+ Dev.\nmilestones",
  "+ Dev. milestones &\ndev. regression",
  "+ Dev. milestones, IQ scores,\n& dev. regression",
  "+ SCQ and\nRBS-R scores",
  "+ SCQ and RBS-R scores,\ndev. milestones, \n& dev. regression",
  "+ Dev. milestones\n& SES",
  "+ Dev.milestones, SES,\n& area deprivation"
))

data <- data.frame(
  Phenotype = factor(rep(phenotype_levels, 3), levels = phenotype_levels),
  h2 = c(0.11, 0.11, 0.12, 0.10, 0.09, 0.09, 0.08, 0.10, 0.07,  # Discovery
         0.11, 0.10, 0.13, 0.10, 0.17, 0.13, 0.16, 0.14, 0.09,  # Replication
         0.11, 0.11, 0.13, 0.10, 0.10, 0.09, 0.09, 0.11, 0.08), # Meta
  SE = c(0.02, 0.02, 0.03, 0.03, 0.07, 0.03, 0.03, 0.03, 0.04,  # Discovery
         0.04, 0.04, 0.07, 0.07, 0.16, 0.07, 0.07, 0.07, 0.08,  # Replication
         0.02, 0.02, 0.03, 0.03, 0.06, 0.03, 0.03, 0.03, 0.03), # Meta
  Cohort = factor(rep(c("Discovery", "Replication", "Meta-analysis"), each = 9),
                  levels = c("Discovery", "Replication", "Meta-analysis"))
)

# Calculate original CIs without constraints
data$CI_lower_orig <- data$h2 - 1.96*data$SE
data$CI_upper_orig <- data$h2 + 1.96*data$SE

# Constrain CIs to plot limits
x_min <- 0
x_max <- 0.25
data$CI_lower <- pmax(x_min, data$CI_lower_orig)
data$CI_upper <- pmin(x_max, data$CI_upper_orig)
data$cohorttype <-ifelse(data$Cohort == "Meta-analysis", 1, 0)

# Create dataset for arrows
arrows_data <- data %>%
  mutate(
    needs_right_arrow = CI_upper_orig > x_max,
    needs_left_arrow = CI_lower_orig < x_min,
    y_position = as.numeric(Phenotype)
  ) %>%
  group_by(Phenotype, Cohort) %>%
  mutate(
    y_position = case_when(
      Cohort == "Discovery" ~ y_position - 0.25,
      Cohort == "Replication" ~ y_position,
      Cohort == "Meta-analysis" ~ y_position + 0.25
    )
  )

B = ggplot(data, aes(y = Phenotype, x = h2, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray5") +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper,
                      linetype = Cohort),
                  position = position_dodge(width = 0.65),
                  size = 0.8) +
  # Add right arrows
  geom_segment(data = subset(arrows_data, needs_right_arrow),
               aes(x = x_max - 0.01, xend = x_max,
                   y = y_position, yend = y_position),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 0.8) +
  scale_color_manual(values = c("#fc8d62","#8da0cb", "#56B4E9")) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid")) +
  theme_classic() +
  labs(x = expression(paste("SNP-based ", h^2)), y = "") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10, lineheight = 0.8),
        legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 0.25))


###Second figure
data <- data.frame(
  Cohort = factor(c("iPSYCH", "SPARK discovery", "SPARK replication"),
                  levels = rev(c("SPARK discovery", "SPARK replication", "iPSYCH"))),
  h2 = c(0.1, 0.11, 0.11),
  SE = c(0.03, 0.02, 0.04)
)

data$CI_lower <- pmax(0, data$h2 - 1.96*data$SE)
data$CI_upper <- data$h2 + 1.96*data$SE

set2_colors <- brewer.pal(3, "Set2")
custom_colors <- c(set2_colors[1], set2_colors[3], set2_colors[2])

A = ggplot(data, aes(y = Cohort, x = h2, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray5") +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper),
                  size = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(x = expression(paste("SNP-based ", h^2)), y = "") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(limits = c(0, max(data$CI_upper) * 1.1))



#Third figure


# Create initial dataframe
data <- data.frame(
  Phenotype = c("PGC 2017", "SPARK (Matoba 2020)", "Autism (Grove 2019)", 
                "Autism (FinnGen)", "Autism (iPSYCH, before_11)",
                "Autism (iPSYCH, males only)", "Autism (iPSYCH only)",
                "Autism (iPSYCH, after_10)", "Autism (iPSYCH, females only)"))

# Add correlation and SE columns
data$SPARK_rg <- c(-0.70, -0.55, -0.27, -0.02, -0.06, 0.07, 0.12, 0.27, 0.23)
data$SPARK_SE <- c(0.11, 0.13, 0.09, 0.16, 0.10, 0.11, 0.09, 0.10, 0.13)
data$iPSYCH_rg <- c(-0.57, -0.59, -0.50, -0.03, -0.56, -0.33, -0.27, 0.09, -0.01)
data$iPSYCH_SE <- c(0.14, 0.15, 0.11, 0.21, 0.08, 0.12, 0.10, 0.13, 0.14)

# Order by SPARK rg
data <- data[order(data$SPARK_rg), ]
data$Phenotype <- factor(data$Phenotype, levels = data$Phenotype)

# Create long format correctly preserving SEs
long_data <- data.frame(
  Phenotype = rep(data$Phenotype, 2),
  Cohort = rep(c("SPARK", "iPSYCH"), each = nrow(data)),
  rg = c(data$SPARK_rg, data$iPSYCH_rg),
  SE = c(data$SPARK_SE, data$iPSYCH_SE)
)

long_data$CI_lower <- long_data$rg - 1.96 * long_data$SE
long_data$CI_upper <- long_data$rg + 1.96 * long_data$SE

C = ggplot(long_data, aes(y = Phenotype, x = rg, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray5") +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper),
                  position = position_dodge(0.5),
                  size = 0.8) +
  scale_color_manual(values = c("iPSYCH" = brewer.pal(3, "Set2")[1],
                                "SPARK" = "#56B4E9")) +
  theme_classic() +
  labs(x = expression(r[g]), y = "", color = "") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom") +
  scale_x_continuous(limits = c(min(long_data$CI_lower) * 1.1,
                                max(long_data$CI_upper) * 1.1))

leftcol <- plot_grid(A, C, labels = c('A', 'C'), label_size = 12, ncol = 1, rel_heights = c(1, 2.5))
plot_grid(leftcol, B, labels = c('', 'B'), label_size = 12, ncol = 2, rel_widths = c(1, 1.2))
