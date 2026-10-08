# UXM plots
library(ggplot2)
library(readxl)
library(ggpubr)
library(tidyverse)

uxm_prop<-read_excel("V:/NeuroEpigen_031114/39. EM-seq/BIOP-AD/3. Major study/4. Bioinformática/2. Deconvolution/UXM/results/output_all_groups.xlsx")
uxm_prop <- read_xlsx("C:/Users/alvarez.166446/Desktop/Fragmentomic/uxm/output_all_groups.xlsx")
samples$Sample <- gsub("_AD$|_CT$", "", samples$Sample)

colnames(uxm_prop)
# Step 1: Pivot uxm_prop to long format
uxm_long<- uxm_prop %>%
  pivot_longer(
    cols = -CellType,
    names_to = "Sample",
    values_to = "Proportion"
  )

# Step 2: Clean Sample names (remove '.deduplicated' suffix)

uxm_long <- uxm_long %>%
  mutate(Sample = str_remove(Sample, "\\.deduplicated.*"))

# Step 3: Join with group info from `samples`
colnames(samples) <- c("Sample", "Group")  # rename for joining
uxm_long <- uxm_long %>%
  left_join(samples, by = "Sample")

# Set 4: select the groups

uxm_ad_ct <- uxm_long %>%
  filter(Group %in% c("AD", "CT")) %>%
  mutate(
    Group = factor(Group, levels = c("AD", "CT"))
  )

my_comparisons <- list(c("AD", "CT"))

palette_colors <- c(
  "AD" = "cornflowerblue",
  "CT" = "yellow2")

# Guardar la matriz estdístuica de celltype 

summary_celltypes <- uxm_ad_ct %>%
  group_by(CellType, Group) %>%
  summarise(
    N = sum(!is.na(Proportion)),
    Min = min(Proportion, na.rm = TRUE),
    Q1 = quantile(Proportion, 0.25, na.rm = TRUE),
    Median = median(Proportion, na.rm = TRUE),
    Mean = mean(Proportion, na.rm = TRUE),
    Q3 = quantile(Proportion, 0.75, na.rm = TRUE),
    Max = max(Proportion, na.rm = TRUE),
    SD = sd(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_celltypes)

write.csv2(summary_celltypes, "ADvsCT_cell_proportions_summary.csv")
summary_celltypes<- read.csv2("ADvsCT_cell_proportions_summary.tsv")
#test estadístico wilcoxon (u de mann) 

wilcox_results <- uxm_ad_ct %>%
  group_by(CellType) %>%
  summarise(
    p_value = wilcox.test(Proportion ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_FDR = p.adjust(p_value, method = "fdr"),
    significance = case_when(
      p_adj_FDR <= 0.0001 ~ "****",
      p_adj_FDR < 0.001 ~ "***",
      p_adj_FDR < 0.01 ~ "**",
      p_adj_FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
wilcox_results2 <- uxm_ad_ct %>%
  group_by(CellType) %>%
  summarise(
    n_CT = sum(Group == "CT" & !is.na(Proportion)),
    n_AD = sum(Group == "AD" & !is.na(Proportion)),
    mean_CT = mean(Proportion[Group == "CT"], na.rm = TRUE),
    mean_AD = mean(Proportion[Group == "AD"], na.rm = TRUE),
    median_CT = median(Proportion[Group == "CT"], na.rm = TRUE),
    median_AD = median(Proportion[Group == "AD"], na.rm = TRUE),
    p_value = ifelse(
      n_CT > 1 & n_AD > 1,
      wilcox.test(Proportion ~ Group, exact = FALSE)$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_FDR = p.adjust(p_value, method = "fdr"),
    significance = case_when(
      is.na(p_adj_FDR) ~ NA_character_,
      p_adj_FDR <= 0.0001 ~ "****",
      p_adj_FDR < 0.001 ~ "***",
      p_adj_FDR < 0.01 ~ "**",
      p_adj_FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    direction = case_when(
      mean_AD > mean_CT ~ "Higher in AD",
      mean_AD < mean_CT ~ "Lower in AD",
      TRUE ~ "No change"
    )
  )

write.csv2(wilcox_results2, "test_wilcox_ADvsCT_UXM")
celltype_interest <- "Thyroid-Ep" # poner el tipo celular que nos interesa graficar.
plot_df <- uxm_ad_ct %>%
  filter(CellType == celltype_interest)

stats_cell <- wilcox_results %>%
  filter(CellType == celltype_interest)
p<- ggboxplot(
  plot_df,
  x = "Group",
  y = "Proportion",
  fill = "Group",
  palette = palette_colors,
  alpha = 0.7,
  outlier.shape = NA
) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.7,
    color = "black"
  ) +
  stat_compare_means(
    comparisons = list(c("AD", "CT")),
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.signif"
  ) +
  theme_minimal() +
  labs(
    title = paste0(celltype_interest, " proportion"),
    subtitle = paste0(
      "p = ", signif(stats_cell$p_value, 3),
      " | FDR = ", signif(stats_cell$p_adj_FDR, 3)
    ),
    x = "",
    y = paste0(celltype_interest," proportion")
  ) +
  theme(
    legend.position = "none"
  )
print(p)
ggsave(
  filename = file.path(output_dir, paste0("AD_vs_CT_", celltype_interest, "_boxplot.pdf")),
  plot = p,
  width = 4,
  height = 4,
  units = "in"
)
output_dir <- "C:/Users/alvarez.166446/Desktop/Fragmentomic/uxm"
