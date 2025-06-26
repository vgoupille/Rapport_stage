# Simple trimming plot for Quarto integration
# This script creates a clean visualization of trimming results

library(ggplot2)
library(dplyr)
library(tidyr)

# Data from the results
trimming_data <- data.frame(
  Sample = c("BC_0076", "BC_0077", "BC_0079", "BC_0080"),
  Before_trimming = c(631.4, 325.5, 379.1, 397.7),
  After_trimming = c(450.8, 248.4, 285.2, 300.1),
  Change_percent = c(-28.6, -23.7, -24.8, -24.5)
)

# Create comparison plot
trimming_long <- trimming_data %>%
  pivot_longer(cols = c(Before_trimming, After_trimming), 
               names_to = "Condition", 
               values_to = "Sequences_M") %>%
  mutate(Condition = factor(Condition, levels = c("Before_trimming", "After_trimming")))

# Main comparison plot
trimming_plot <- ggplot(trimming_long, aes(x = Sample, y = Sequences_M, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Before_trimming" = "#3498db", "After_trimming" = "#e74c3c"),
                    labels = c("Before trimming", "After trimming")) +
  labs(title = "Sequencing Depth Before and After Trimming",
       subtitle = "Comparison of sequence counts across samples",
       x = "Sample",
       y = "Sequences (M)",
       fill = "Condition",
       caption = "Data shows consistent reduction across all samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 9, color = "gray50")) +
  scale_y_continuous(labels = scales::comma) +
  geom_text(aes(label = sprintf("%.1fM", Sequences_M)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3, fontface = "bold")

# Percentage change plot
percent_plot <- ggplot(trimming_data, aes(x = Sample, y = Change_percent)) +
  geom_bar(stat = "identity", fill = "#f39c12", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Change_percent)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Percentage Reduction After Trimming",
       x = "Sample",
       y = "Reduction (%)",
       caption = "All samples show similar reduction rates (~24-29%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.caption = element_text(hjust = 0, size = 9, color = "gray50")) +
  ylim(min(trimming_data$Change_percent) - 2, 2) +
  geom_hline(yintercept = mean(trimming_data$Change_percent), 
             linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 2.5, y = mean(trimming_data$Change_percent) + 1, 
           label = paste("Mean:", round(mean(trimming_data$Change_percent), 1), "%"), 
           color = "red", fontface = "bold")

# Print summary statistics
cat("=== TRIMMING RESULTS SUMMARY ===\n")
cat("Total sequences before trimming:", sum(trimming_data$Before_trimming), "M\n")
cat("Total sequences after trimming:", sum(trimming_data$After_trimming), "M\n")
cat("Overall reduction:", round((sum(trimming_data$After_trimming) / sum(trimming_data$Before_trimming) - 1) * 100, 1), "%\n")
cat("Mean reduction per sample:", round(mean(trimming_data$Change_percent), 1), "%\n")
cat("Range of reduction:", round(min(trimming_data$Change_percent), 1), "% to", round(max(trimming_data$Change_percent), 1), "%\n")
cat("================================\n")

# Display plots
print(trimming_plot)
print(percent_plot) 