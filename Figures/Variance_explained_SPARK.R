# Create the data frame with corrected values and explicit ordering
data <- data.frame(
  Subcategory = c("Race/ethnicity", "Family income", "Mother's highest education", 
                  "Father's highest education", "Sex",
                  "Reported cognitive test score", "Reported cognitive impairment",
                  "SCQ total score", "RBS-R total score", "Age at walking independently",
                  "Age at using words", "Language regression", "Other regression",
                  "SNP-based heritability", "Unexplained"),
  Category = c(rep("Sociodemographic factors", 5),
               rep("Clinical Factors", 8),
               "Genetics", "Unexplained"),
  Variance = c(0.011883012, 0.15059013, 0.07300062, 0.40860654, 1.550329078,
               2.9524485, 0.666184706, 8.309129944, 0.601664745, 0.121051424,
               2.880973253, 5.36291336, 0.181224689, 11, 65.73)
)

# Order within categories
data <- data %>%
  group_by(Category) %>%
  arrange(Category, desc(Variance)) %>%
  ungroup()

# Add unexplained variance

# Create the explicit order we want (from bottom to top)
ordered_levels <- c(
  # Sociodemographic (ordered by variance)
  "Race/ethnicity", "Mother's highest education", "Family income", 
  "Father's highest education", "Sex",
  # Clinical (ordered by variance)
  "Age at walking independently", "Other regression", "RBS-R total score",
  "Reported cognitive impairment", "Age at using words", 
  "Reported cognitive test score", "Language regression", "SCQ total score",
  # Genetics
  "SNP-based heritability",
  # Unexplained at top
  "Unexplained"
)

# Set the factor levels
data$Subcategory <- factor(data$Subcategory, levels = ordered_levels)

# Define colors explicitly matched to categories
color_mapping <- c(
  # Sociodemographic - blues
  "Race/ethnicity" = "#003f5c",
  "Mother's highest education" = "#2f4b7c",
  "Family income" = "#465c7c",
  "Father's highest education" = "#665191",
  "Sex" = "#7b6195",
  # Clinical - greens/yellows
  "Age at walking independently" = "#1F968BFF",
  "Other regression" = "#98c9b0",
  "RBS-R total score" = "#73D055FF",
  "Reported cognitive impairment" = "lightgreen",
  "Age at using words" = "#B8DE29FF",
  "Reported cognitive test score" = "#FDE725FF" ,
  "Language regression" = "#ffbc42",
  "SCQ total score" = "#ffa600", 
  # Genetics - purple
  "SNP-based heritability" = "darkred",
  # Unexplained - grey
  "Unexplained" = "#cccccc"
)


# Create a layout with two plots side by side
# First, create the data for the second plot (without unexplained)
data_explained <- data %>%
  filter(Subcategory != "Unexplained")


# Create two plots
p1 <- ggplot(data, aes(x = "Total Variance", y = Variance, fill = Subcategory)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = color_mapping) +
  theme_minimal() +
  labs(y = "Variance (%)", 
       x = "",
       title = "Total Variance") +
  theme(legend.position = "right",
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0, 100))

p2 <- ggplot(data_explained, aes(x = "Explained Variance", y = Variance, fill = Subcategory)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), 
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_manual(values = color_mapping) +
  theme_minimal() +
  labs(y = "Explained Variance (%)", 
       x = "",
       title = "Explained Variance Only") +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0, 34.27))

# Combine the plots
library(gridExtra)
grid.arrange(p1, p2, ncol = 2, widths = c(1.3, 1))

