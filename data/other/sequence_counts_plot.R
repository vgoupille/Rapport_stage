# Plot for sequence counts table visualization
# Based on the data from the results table

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data from the table in the results
sequence_data <- data.frame(
  Sample = c("BC_0076", "BC_0077", "BC_0079", "BC_0080", "Total"),
  Before_trimming = c(631.4, 325.5, 379.1, 397.7, 1733.7),
  After_trimming = c(450.8, 248.4, 285.2, 300.1, 1284.5),
  Change_percent = c(-28.6, -23.7, -24.8, -24.5, -25.4)
)

# Create long format for plotting (excluding Total for main plot)
plot_data <- sequence_data %>%
  filter(Sample != "Total") %>%
  pivot_longer(cols = c(Before_trimming, After_trimming), 
               names_to = "Condition", 
               values_to = "Sequences_M") %>%
  mutate(Condition = factor(Condition, levels = c("Before_trimming", "After_trimming")))

# Main comparison plot
p1 <- ggplot(plot_data, aes(x = Sample, y = Sequences_M, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Before_trimming" = "#3498db", "After_trimming" = "#e74c3c"),
                    labels = c("Before trimming", "After trimming")) +
  labs(title = "Sequencing Depth Before and After Trimming",
       subtitle = "Number of sequences (M) per sample",
       x = "Sample",
       y = "Sequences (M)",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  scale_y_continuous(labels = comma) +
  geom_text(aes(label = sprintf("%.1fM", Sequences_M)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5, fontface = "bold")

# Percentage change plot
p2 <- ggplot(sequence_data %>% filter(Sample != "Total"), 
             aes(x = Sample, y = Change_percent)) +
  geom_bar(stat = "identity", fill = "#f39c12", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Change_percent)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Percentage Change in Sequences",
       x = "Sample",
       y = "Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ylim(min(sequence_data$Change_percent[sequence_data$Sample != "Total"]) - 2, 2) +
  geom_hline(yintercept = -25.4, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 2.5, y = -23, 
           label = "Total: -25.4%", color = "red", fontface = "bold")

# Total comparison plot
total_data <- sequence_data %>%
  filter(Sample == "Total") %>%
  pivot_longer(cols = c(Before_trimming, After_trimming), 
               names_to = "Condition", 
               values_to = "Sequences_M") %>%
  mutate(Condition = factor(Condition, levels = c("Before_trimming", "After_trimming")))

p3 <- ggplot(total_data, aes(x = "Total", y = Sequences_M, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Before_trimming" = "#3498db", "After_trimming" = "#e74c3c"),
                    labels = c("Before trimming", "After trimming")) +
  labs(title = "Total Sequences Across All Samples",
       x = "",
       y = "Sequences (M)",
       fill = "Condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom") +
  scale_y_continuous(labels = comma) +
  geom_text(aes(label = sprintf("%.1fM", Sequences_M)), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5, size = 4, fontface = "bold")

# Print summary
cat("=== SEQUENCE COUNTS SUMMARY ===\n")
cat("Sample\tBefore\tAfter\tChange\n")
for(i in 1:nrow(sequence_data)) {
  cat(sprintf("%s\t%.1fM\t%.1fM\t%.1f%%\n", 
              sequence_data$Sample[i], 
              sequence_data$Before_trimming[i], 
              sequence_data$After_trimming[i], 
              sequence_data$Change_percent[i]))
}
cat("==============================\n")

# Display plots
print(p1)
print(p2)
print(p3)

# Save plots
ggsave("figures/sequence_counts_comparison.png", p1, width = 10, height = 6, dpi = 300)
ggsave("figures/sequence_percentage_change.png", p2, width = 8, height = 6, dpi = 300)
ggsave("figures/total_sequences.png", p3, width = 6, height = 6, dpi = 300) 