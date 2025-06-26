# R code chunk for Quarto integration
# Copy this code directly into your .qmd file

```{r}
#| label: fig-sequence-counts
#| fig-cap: "Comparison of sequencing depth before and after trimming across all samples"
#| fig-width: 10
#| fig-height: 6

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data from the results table
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
ggplot(plot_data, aes(x = Sample, y = Sequences_M, fill = Condition)) +
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
```

```{r}
#| label: fig-percentage-change
#| fig-cap: "Percentage reduction in sequences after trimming"
#| fig-width: 8
#| fig-height: 6

# Percentage change plot
ggplot(sequence_data %>% filter(Sample != "Total"), 
       aes(x = Sample, y = Change_percent)) +
  geom_bar(stat = "identity", fill = "#f39c12", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Change_percent)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Percentage Change in Sequences After Trimming",
       x = "Sample",
       y = "Change (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ylim(min(sequence_data$Change_percent[sequence_data$Sample != "Total"]) - 2, 2) +
  geom_hline(yintercept = -25.4, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 2.5, y = -23, 
           label = "Total: -25.4%", color = "red", fontface = "bold")
```

```{r}
#| label: fig-total-sequences
#| fig-cap: "Total sequences across all samples"
#| fig-width: 6
#| fig-height: 6

# Total comparison plot
total_data <- sequence_data %>%
  filter(Sample == "Total") %>%
  pivot_longer(cols = c(Before_trimming, After_trimming), 
               names_to = "Condition", 
               values_to = "Sequences_M") %>%
  mutate(Condition = factor(Condition, levels = c("Before_trimming", "After_trimming")))

ggplot(total_data, aes(x = "Total", y = Sequences_M, fill = Condition)) +
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
``` 