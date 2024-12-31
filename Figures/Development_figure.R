library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(patchwork)


palette1 = c("#6B3F98", "#36A047")

palette2 = c("#95D840FF", "#440154FF")


pd <- position_dodge(width = 0.5)
pd2 <- position_dodge(width = 1)
pd3 <- position_dodge(width = -0.5)

data_devmilestone1 = read.delim(pipe("pbpaste"))
data_devmilestone1$Phenotype <- factor(data_devmilestone1$Phenotype, levels = unique(data_devmilestone1$Phenotype))


A5 = ggplot(data_devmilestone1, aes(x = Phenotype, y = rg, colour = Factor)) + geom_point(position = pd, aes(shape = BY.adjusted.P)) + scale_shape_manual(values = c(16,8)) +
  geom_errorbar(width=.3, aes(ymin = rg-1.96*SE, ymax = rg+1.96*SE), position=pd) +
  theme_classic() + xlab("") + ylab("Genetic correlation") + geom_hline(yintercept = 0) + theme(axis.text.x=element_text(size=rel(1))) + coord_flip()  +
  theme(strip.text.y = element_blank()) + scale_colour_manual(values=palette2) + theme(legend.position="bottom") +
  theme(strip.text.y = element_blank()) + theme(text = element_text(size = 12)) + guides(shape = "none")


data_devmilestone2 = read.delim(pipe("pbpaste"))
data_devmilestone2$Phenotype <- factor(data_devmilestone2$Phenotype, levels = unique(data_devmilestone2$Phenotype))
data_devmilestone2$adjusted = ifelse(data_devmilestone2$BY.adjusted.P < 0.05, "*", "")

B5 = ggplot(data_devmilestone2, aes(x = Phenotype, y = Beta, colour = PGS)) + 
  geom_point(position = pd, aes(shape = adjusted)) + scale_shape_manual(values = c(16,8)) + geom_errorbar(width=.3, aes(ymin = Beta-1.96*SE, ymax = Beta+1.96*SE), position=pd) +
  theme_classic() + xlab("") + ylab("PGS association") + geom_hline(yintercept = 0) + theme(axis.text.x=element_text(size=rel(1))) + coord_flip() +
  theme(strip.text.y = element_blank()) + scale_colour_manual(values=palette1) + theme(legend.position="bottom") +
  theme(strip.text.y = element_blank()) + theme(text = element_text(size = 12)) + guides(size = "none") + guides(shape = "none")




MCS_data = read.delim(pipe("pbpaste"))
MCS_data = subset(MCS_data, Trait == "Total")
ALSPAC_data = read.delim(pipe("pbpaste"))
ALSPAC_data = subset(ALSPAC_data, Trait == "Total")



# Function to process data
process_data <- function(data, ages) {
  # Calculate age mean and sd 
  age_mean <- mean(ages)
  age_sd <- sd(ages)
  
  # Standardise ages
  ages_std <- (ages - age_mean) / age_sd
  
  # Create processed dataset
  processed <- data %>%
    mutate(
      version = ifelse(grepl("after10", PGS), "After 2010", "Before 2011"),
      trait = gsub("Autism.*", "", Trait)
    ) %>%
    select(trait, version, Beta, SE) %>%
    # Create data for each age point
    crossing(age = ages) %>%
    mutate(
      age_std = (age - age_mean) / age_sd,
      effect = Beta * age_std,
      lower = effect - 1.96 * SE * abs(age_std),
      upper = effect + 1.96 * SE * abs(age_std)
    )
  
  return(processed)
}



# Process both datasets
MCS_processed <- process_data(MCS_data, c(3, 5, 7, 11, 14, 17))
ALSPAC_processed <- process_data(ALSPAC_data, c(7, 10, 12, 13, 17))

# Function to create plot for each cohort
create_cohort_plot <- function(data, title_prefix, ages) {
  ggplot(data, aes(x = age, y = effect, 
                   color = version, fill = version)) +
    facet_wrap(~trait, scales = "free_y") +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha = 0.2, color = NA) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = palette1) +
    scale_fill_manual(values = palette1) +
    theme_minimal() +
    labs(x = "Age (years)",
         y = "Standardized Effect",
         color = "none",
         fill = "none") +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = ages)
}

# Create individual plots
p1 <- create_cohort_plot(MCS_processed, "A. MCS", c(3, 5, 7, 11, 14, 17))
p2 <- create_cohort_plot(ALSPAC_processed, "B. ALSPAC", c(7, 10, 12, 13, 17))

# Combine plots
C5 <- p1 / p2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "none") 



rightcol <- plot_grid(B5, C5, labels = c('B', 'C'), label_size = 12, ncol = 1, rel_heights = c(1, 1.5))
plot_grid(A5, rightcol, labels = c('A', ''), label_size = 12, ncol = 2, rel_widths = c(1.5, 1))



