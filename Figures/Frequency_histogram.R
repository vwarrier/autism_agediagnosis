library(tidyverse)

#Data column headers - Age, Count, Cohort, long format. "SPARK Discovery", "SPARK replication", "iPSYCH"

data = read.delim(pipe("pbpaste"))


weighted_median <- function(x, w) {
  order_idx <- order(x)
  x_sorted <- x[order_idx]
  w_sorted <- w[order_idx]
  w_cum <- cumsum(w_sorted)
  w_total <- sum(w)
  median_idx <- which(w_cum >= w_total/2)[1]
  return(x_sorted[median_idx])
}

weighted_mad <- function(x, w, median) {
  deviations <- abs(x - median)
  return(weighted_median(deviations, w))
}

stats_by_cohort <- data %>%
  group_by(Cohort) %>%
  summarise(
    median = weighted_median(Age, Count),
    mad = weighted_mad(Age, Count, weighted_median(Age, Count))
  )

ggplot(data, aes(x = Age, y = Count, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_vline(data = stats_by_cohort, 
             aes(xintercept = median, color = Cohort),
             linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(
    title = "Age at Diagnosis Distribution with Medians",
    x = "Age at Diagnosis",
    y = "Count"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  geom_text(data = stats_by_cohort,
            aes(x = median, 
                y = max(data$Count),
                color = Cohort,
                label = sprintf("Median: %.1f\nMAD: %.1f", median, mad)),
            vjust = 1,
            show.legend = FALSE)
