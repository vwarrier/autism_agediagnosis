library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


MCS_data = read.delim(pipe("pbpaste"))
ALSPAC_data = read.delim(pipe("pbpaste"))


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

# Define color palette
palette1 <- c("#6B3F98", "#36A047")

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
         color = "PGS Version",
         fill = "PGS Version",
         title = paste(title_prefix)) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11)) +
    scale_x_continuous(breaks = ages)
}

# Create individual plots
p1 <- create_cohort_plot(MCS_processed, "A. MCS", c(3, 5, 7, 11, 14, 17))
p2 <- create_cohort_plot(ALSPAC_processed, "B. ALSPAC", c(7, 10, 12, 13, 17))

# Combine plots
combined_plot <- p1 / p2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)


