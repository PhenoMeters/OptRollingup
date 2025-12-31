library(tidyverse)
library(ggpubr)

text_size = 16


#LIP Figures

## Methods

## Simulation Results

lip_res = readRDS(here::here("Out", "sim_res", "lip_rollup_sims.RDS"))
lip_res = lip_res %>%
  mutate(aggregation_method = ifelse(Analysis_Method == "TrypticPeptides", "", aggregation_method)) %>%
  mutate(Analysis_Method = paste(Analysis_Method, aggregation_method, sep = ' '))


lip_res %>%
  filter(
    trypsin_missed_prob == 0,
    proteinase_k_missed_prob == 0) %>%
  ggplot(
    aes(
      x = masked,
      y = MSE,
      color = Analysis_Method
    )
  ) +
  geom_boxplot() +
  facet_grid(n_per_group ~ n_protein_replicates) +
  theme_bw()

lip_res %>%
  filter(
    n_protein_replicates == 100,
    n_per_group == 25) %>%
  ggplot(
    aes(
      x = masked,
      y = MSE,
      color = Analysis_Method
    )
  ) +
  geom_boxplot() +
  facet_grid(trypsin_missed_prob ~ proteinase_k_missed_prob) +
  theme_bw()


p_overall = lip_res %>%
  mutate(`Region Type` = ifelse(masked, "Differential Masking", "Non-differential Masking")) %>%
  mutate(RMSE = sqrt(MSE)) %>%
  ggplot(
    aes(
      x = `Region Type`,
      y = RMSE,
      color = Analysis_Method
    )
  ) +
  geom_boxplot() +
  # facet_grid(trypsin_missed_prob ~ proteinase_k_missed_prob) +
  theme_bw() +
  guides(color = guide_legend(title = "Analysis Method")) +
  scale_color_viridis_d(labels = c("Site Rollup (Mean)", "Site Rollup (Median)", "Site Rollup (Sum)", "Tryptic Peptides")) +
  # scale_color_manual(
  #   values = c("#7AD151FF", '#414487FF'),
  #   labels = c("Site Rollup", "Tryptic Peptides")
  # ) +
  theme(
    legend.position = 'bottom',
    text = element_text(size = text_size)
    )


ggsave(here::here("Out", "plots", "lip_overall_plot_new.png"), p_overall, width = 6, height = 4)
p_overall


## Supplementary

res = lip_res

this_n_per_group = 5
# this_n_per_group = 10
# this_n_per_group = 25
# this_n_per_group = 50

df = res %>%
  filter(
    n_per_group == this_n_per_group
  ) %>%
  mutate(masked = ifelse(masked, "Differential Masking", 'Non-differential Masking')) %>%
  rename(
    `Method` = Analysis_Method,
    `Region Type` = masked,
    `Samples per Group` = n_per_group,
    `Protein Replicates` = n_protein_replicates,
    `p1` = trypsin_missed_prob,
    `p2` = proteinase_k_missed_prob
  )
df

p = df %>%
  ggplot(
    aes(
      x = `Region Type`,
      y = MSE,
      color = `Method`
    )
  ) +
  geom_boxplot() +
  scale_color_viridis_d(labels = c("Site Rollup (Mean)", "Site Rollup (Median)", "Site Rollup (Sum)", "Tryptic Peptides")) +
  facet_grid(p1 + p2 ~ `Protein Replicates`) +
  ylim(0, 1.5) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotates x-axis ticks by 45 degrees
    text = element_text(size = text_size)
  ) +
  ggtitle(paste0(n_per_group, ' Samples per Group'))


ggsave(here::here('Out', 'plots', paste0('lip_supp_', this_n_per_group, '_new.png')), p, width = 8, height = 12)



### Gap figure

lip_masked_sims <- readRDS("~/Documents/project_files/ppi_systems/workspace/Code/proteosim/Out/sim_res/lip_masked_sims.RDS")

gap_plot = lip_masked_sims %>%
  mutate(`Region Type` = ifelse(masked, "Differential Masking", "Non-differential Masking")) %>%
  mutate(RMSE = sqrt(MSE)) %>%
  ggplot(
    aes(
      x = `Region Type`,
      y = RMSE,
      color = Analysis_Method
    )
  ) +
  geom_boxplot() +
  # facet_grid(trypsin_missed_prob ~ proteinase_k_missed_prob) +
  theme_bw() +
  facet_grid(gap_parameter ~ max_masked_region_size) +
  guides(color = guide_legend(title = "Analysis Method")) +
  scale_color_viridis_d(labels = c("Site Rollup (Median)", "Tryptic Peptides")) +
  # scale_color_manual(
  #   values = c("#7AD151FF", '#414487FF'),
  #   labels = c("Site Rollup", "Tryptic Peptides")
  # ) +
  theme(
    legend.position = 'bottom',
    text = element_text(size = text_size),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotates x-axis ticks by 45 degrees
    )
ggsave(here::here("Out", "plots", "lip_gap_plot_new.png"), gap_plot, width = 6, height = 6)


#PTM Figures

## Methods

## Simulation Results

res <- readRDS("~/Documents/project_files/ppi_systems/workspace/Code/proteosim/Out/sim_res/init_sim_res_ptm3.RDS")


df = res %>% group_by(
  Method,
  Combine_Function,
  total_abundance,
  protein_sequence_length,
  missingness_prop,
  prop_to_miss
) %>%
  rename(
    `Combine Function` = Combine_Function,
    `Protein Abundance` = total_abundance,
    `Protein Length` = protein_sequence_length,
    `p1` = prop_to_miss
  ) %>%
  mutate(
    `Combine Function` = case_when(
      `Combine Function` == "mean" ~ "Mean",
      `Combine Function` == "median" ~ "Median",
      `Combine Function` == "sum" ~ "Sum"
    )
  )


p = df %>%
  ggplot(
    aes(
      x = `Method`,
      y = rmse,
      color = `Combine Function`
    )
  ) +
  ylim(0, 5) +
  geom_boxplot() +
  scale_color_viridis_d(begin = 0, end = 0.8, direction = -1) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotates x-axis ticks by 45 degrees
    text = element_text(size = text_size)
  ) +
  ylab("RMSE")# +
  # ggtitle(paste0("Proportion of Inducted Missingness: ", missingness_prop))

ggsave(here::here('Out', 'plots', paste0('ptm_overall_paper.png')), p, width = 6, height = 6)



this_missingness_prop = 0
# this_missingness_prop = 0.25
# this_missingness_prop = 0.5

df = res %>% group_by(
  Method,
  Combine_Function,
  total_abundance,
  protein_sequence_length,
  # missingness_prop,
  prop_to_miss
) %>%
  filter(
    missingness_prop == this_missingness_prop
  ) %>%
  rename(
    `Combine Function` = Combine_Function,
    `Protein Abundance` = total_abundance,
    `Protein Length` = protein_sequence_length,
    `p1` = prop_to_miss
  ) %>%
  mutate(
    `Combine Function` = case_when(
      `Combine Function` == "mean" ~ "Mean",
      `Combine Function` == "median" ~ "Median",
      `Combine Function` == "sum" ~ "Sum"
    )
  )


p = df %>%
  ggplot(
    aes(
      x = `Method`,
      y = rmse,
      color = `Combine Function`
    )
  ) +
  ylim(0, 5) +
  geom_boxplot() +
  scale_color_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(p1 ~ `Protein Abundance` + `Protein Length`) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotates x-axis ticks by 45 degrees
    text = element_text(size = text_size)
  ) +
  ylab("RMSE") +
  ggtitle(paste0("Proportion of Inducted Missingness: ", this_missingness_prop))

ggsave(here::here('Out', 'plots', paste0('ptm_supp_', this_missingness_prop, '.png')), p, width = 8, height = 10)


