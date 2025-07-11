### LIP


n_per_group = 25

df = res %>% group_by(
  Analysis_Method,
  masked,
  n_per_group,
  n_protein_replicates,
  trypsin_missed_prob,
  proteinase_k_missed_prob
) %>%
  filter(
    n_per_group == n_per_group
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


p = df %>%
  ggplot(
    aes(
      x = `Region Type`,
      y = MSE,
      color = `Method`
    )
  ) +
  geom_boxplot() +
  scale_color_viridis_d(begin = 0, end = 0.8, direction = -1) +
  facet_grid(p1 + p2 ~ `Protein Replicates`) +
  ylim(0, 1.5) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotates x-axis ticks by 45 degrees
  ) +
  ggtitle(paste0(n_per_group, ' Samples per Group'))

ggsave(here::here('Out', 'plots', paste0('lip_supp_', n_per_group, '.png')), p, width = 8, height = 10)



### PTM

missingness_prop = 0

df = res %>% group_by(
  Method,
  Combine_Function,
  total_abundance,
  protein_sequence_length,
  missingness_prop,
  prop_to_miss
) %>%
  filter(
    missingness_prop == missingness_prop
  ) %>%
  rename(
    `Combine Function` = Combine_Function,
    `Protein Abundance` = total_abundance,
    `Protein Length` = protein_sequence_length,
    `p1` = prop_to_miss
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
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotates x-axis ticks by 45 degrees
  ) +
  ggtitle(paste0("Proportion of Inducted Missingness: ", missingness_prop))

ggsave(here::here('Out', 'plots', paste0('lip_supp_', missingness_prop, '.png')), p, width = 8, height = 10)



