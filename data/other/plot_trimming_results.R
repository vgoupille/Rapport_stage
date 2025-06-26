# Script to create plots for trimming results
# Author: Bioinformatics analysis
# Date: 2024

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(grid)

# Create the data frame with sequencing results
trimming_data <- data.frame(
  Sample = c("BC_0076", "BC_0077", "BC_0079", "BC_0080"),
  Before_trimming = c(631.4, 325.5, 379.1, 397.7),
  After_trimming = c(450.8, 248.4, 285.2, 300.1),
  Change_percent = c(-28.6, -23.7, -24.8, -24.5)
)

# Convert to long format for plotting
trimming_long <- trimming_data %>%
  pivot_longer(cols = c(Before_trimming, After_trimming), 
               names_to = "Condition", 
               values_to = "Sequences_M") %>%
  mutate(Condition = factor(Condition, levels = c("Before_trimming", "After_trimming")))

# Create bar plot comparing before and after trimming
p1 <- ggplot(trimming_long, aes(x = Sample, y = Sequences_M, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Before_trimming" = "#3498db", "After_trimming" = "#e74c3c"),
                    labels = c("Before trimming", "After trimming")) +
  labs(title = "Sequencing Depth Before and After Trimming",
       x = "Sample",
       y = "Sequences (M)",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom") +
  scale_y_continuous(labels = comma)

# Create percentage change plot
p2 <- ggplot(trimming_data, aes(x = Sample, y = Change_percent)) +
  geom_bar(stat = "identity", fill = "#f39c12", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Change_percent)), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  labs(title = "Percentage Change in Sequences After Trimming",
       x = "Sample",
       y = "Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylim(min(trimming_data$Change_percent) - 2, 2)

# Create a combined plot with both metrics

# Save plots
ggsave("figures/trimming_comparison.png", p1, width = 10, height = 6, dpi = 300)
ggsave("figures/trimming_percentage_change.png", p2, width = 8, height = 6, dpi = 300)

# Create a comprehensive plot
combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                             top = textGrob("Sequencing Data Quality Assessment", 
                                           gp = gpar(fontsize = 16, fontface = "bold")))

ggsave("figures/trimming_comprehensive_analysis.png", combined_plot, 
       width = 16, height = 8, dpi = 300)

# Print summary statistics
cat("Summary of trimming results:\n")
cat("Total sequences before trimming:", sum(trimming_data$Before_trimming), "M\n")
cat("Total sequences after trimming:", sum(trimming_data$After_trimming), "M\n")
cat("Overall reduction:", round((sum(trimming_data$After_trimming) / sum(trimming_data$Before_trimming) - 1) * 100, 1), "%\n")
cat("Mean reduction per sample:", round(mean(trimming_data$Change_percent), 1), "%\n")

# Create additional visualization: sequence length comparison
length_data <- data.frame(
  Sample = rep(c("BC_0076", "BC_0077", "BC_0079", "BC_0080"), 2),
  Read = rep(c("R1", "R2"), each = 4),
  Before_length = c(241, 241, 241, 241, 91, 91, 91, 91),
  After_length = c(127, 157, 152, 132, 90, 90, 90, 90)
)

length_long <- length_data %>%
  pivot_longer(cols = c(Before_length, After_length),
               names_to = "Condition",
               values_to = "Length_bp") %>%
  mutate(Condition = factor(Condition, levels = c("Before_length", "After_length")))

p3 <- ggplot(length_long, aes(x = Sample, y = Length_bp, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~Read, scales = "free_y") +
  scale_fill_manual(values = c("Before_length" = "#2ecc71", "After_length" = "#e67e22"),
                    labels = c("Before trimming", "After trimming")) +
  labs(title = "Read Length Before and After Trimming",
       x = "Sample",
       y = "Length (bp)",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("figures/read_length_comparison.png", p3, width = 12, height = 6, dpi = 300)

cat("\nPlots have been saved to the figures/ directory.\n") 