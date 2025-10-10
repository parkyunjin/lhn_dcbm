library(ggplot2)
library(readr)
library(dplyr)
library(rstatix)
library(ggpubr)
library(broom)


dat <- read_csv("/result/permloss_100.csv")
dat <- dat %>%
  mutate(algorithm = if_else(algorithm == "scoreplus", "SCORE+", algorithm))%>%
  mutate(metric = if_else(metric == "jac_coef", "Jaccard", metric))%>%
  mutate(metric = if_else(metric == "nrc", "DCRC", metric))%>%
  mutate(metric = toupper(metric))%>%
  mutate(metric = if_else(metric == "JACCARD", "Jaccard", metric))%>%
  mutate(graph = toupper(graph))%>%
  mutate(graph = str_replace(graph, "^N", "Graph "))


dat <- dat %>%
  group_by(seed, graph, algorithm) %>%
  mutate(difference_from_dcrc = difference[metric == 'DCRC'] -  difference) #bigger the difference, it is worse than DCRC


# exp1
dat1 <- dat[dat['graph'] == "Graph 1",]
dat12<- dat[dat['graph'] == "Graph 1" & dat['metric'] != "DCRC",]


# exp3
dat3 <- dat[dat['metric'] == "DCRC",]


### EXP1 Figure ###

significance_results <- dat12 %>%
  group_by(metric, algorithm) %>%
  do(tidy(t.test(.$difference_from_dcrc, mu = 0, alternative = "greater")))


print(significance_results)

label_data <- significance_results %>%
  mutate(
    signif_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "" 
    )
  ) %>%
  left_join(dat12 %>% group_by(metric) %>% summarize(max_val = max(difference_from_dcrc)), by = "metric")

print(label_data)

exp1plot <- ggplot(dat12, aes(x = metric, y = difference_from_dcrc, fill = metric)) +
  geom_boxplot() +
  facet_grid(graph ~ algorithm, labeller = labeller(
    algorithm = toupper
  )) + 
  labs(
    x = "Curvature",
    y = "Relative performance of DCRC"
  ) +
  theme_light() + 
  theme(legend.position = "none",  axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_blank(),
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  geom_text(
    data = label_data,
    aes(x = metric, label = signif_label),
    y= max(dat1$difference_from_dcrc) + 1,
    size = 6
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1)

### EXP3 Figure ###

significance_results3 <- dat3 %>%
  group_by(graph, algorithm) %>%
  do(tidy(t.test(.$difference, mu = 0, alternative = "greater")))

print(significance_results)

label_data3 <- significance_results3 %>%
  mutate(
    signif_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "" 
    )
  ) %>%
  left_join(dat3 %>% group_by(graph) %>% summarize(max_val = max(difference)), by = "graph")

print(label_data)



exp3plot <- ggplot(dat3, aes(x = graph, y = difference, fill = graph)) +
  geom_boxplot() +
  facet_grid(metric~ algorithm, labeller = labeller(
    algorithm = toupper
  )) + 
  labs(
    x = "Algorithm",
    y = "Improvement of DCRC"
  ) +
  theme_light() + # A clean theme is nice for faceted plots
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_blank(),
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))+
  geom_text(
    data = label_data3,
    aes(x = graph, label = signif_label),
    y= max(dat3$difference) + 1,
    size = 6
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1)



# Save as PDF
setwd("/work/users/y/j/yjinpark/Tracy/neurips_2025/experiment_3/fig")
ggsave(
  "exp1plot.pdf",
  plot = exp1plot,
  width = 16,
  height = 5
)

ggsave(
  "exp3plot.pdf",
  plot = exp3plot,
  width = 16,
  height = 5
)
