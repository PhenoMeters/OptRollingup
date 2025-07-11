library(ggplot2)
library(dplyr)
library(OrgMassSpecR)
library(pmartR)
library(tidyverse)


source(here::here("code_to_migrate", "ptm_utils.R"))
source(here::here("code_to_migrate", "LiP_utils.R"))
source(here::here("code_to_migrate", "missing_data_utils.R"))

######################################################
### Simulation
######################################################

# Create amino acid distribution
lip_data = readr::read_tsv(here::here("Data", "double_pept_lytic_sites.tsv"))
all_peptides = unlist(lapply(lip_data$Sequence, \(x) strsplit(x, split = '; ')[[1]]))
unique_peptides = unique(all_peptides)
amino_acids = unlist(lapply(unique_peptides, \(x) strsplit(x, split = "")[[1]]))
amino_acid_distribution = create_amino_acid_distribution(amino_acids)

# Tuning parameters for protein masking
set.seed(70)
sample_delay = function() rpois(1, 10)
downsample_fn = function(vec, prop) sample(vec, floor(prop * length(vec)))
min_num_trypsin_sites = 5
n_masking_patterns = 2
sequence_length = 1000
pk_mu = 8
tr_mu = 3

#stuff I want to change
n_per_group_vec = c(50)
n_protein_replicates_vec = c(1000)
trypsin_missed_prob_vec = c(0.25)
proteinase_k_missed_prob_vec = c(0.25)
max_masked_region_size_vec = c(25)
gap_parameter_vec = c(50)
functions = list(
  'mean' = mean,
  'median' = median,
  'sum' = sum
  )

sim_dat = expand.grid(
  n_per_group = n_per_group_vec,
  trypsin_missed_prob = trypsin_missed_prob_vec,
  proteinase_k_missed_prob = proteinase_k_missed_prob_vec,
  n_protein_replicates = n_protein_replicates_vec,
  max_masked_region_size = max_masked_region_size_vec,
  gap_parameter = gap_parameter_vec,
  function_idx = seq_along(functions)
  )

n_reps = 30

#iterate though replicates
res = c()
for(r in seq_len(n_reps)){
  #iterate though scenarios
  for(i in seq_len(nrow(sim_dat))){
    #get scenario parameters
    n_g1 = sim_dat$n_per_group[i]
    n_g2 = sim_dat$n_per_group[i]
    n_protein_replicates = sim_dat$n_protein_replicates[i]
    trypsin_missed_prob = sim_dat$trypsin_missed_prob[i]
    proteinase_k_missed_prob = sim_dat$proteinase_k_missed_prob[i]
    max_masked_region_size = sim_dat$max_masked_region_size[i]
    gap_parameter = sim_dat$gap_parameter[i]
    function_idx = sim_dat$function_idx[i]

    sample_gap = function() rpois(1, gap_parameter)

    #Generate one synthetic protein
    synthetic_protein = generate_protein(sequence_length, amino_acid_distribution)

    # Perform full masking
    full_masking = mask_sequence(
      sequence = synthetic_protein,
      delay_fn = sample_delay,
      gap_fn = sample_gap,
      max_masked_region_size = max_masked_region_size,
      min_num_trypsin_sites = min_num_trypsin_sites
    )

    get_masking_pattern_proportions = function(masking_pattern_set, mu = 0){
      unifs = runif(length(masking_pattern_set))
      unifs = rnorm(length(masking_pattern_set), mu)
      proportions = unifs / sum(unifs)
      return(proportions)
    }

    masking_pattern_set = create_masking_set(full_masking, n_masking_patterns)
    g1_masking_pattern_proportions = n_masking_patterns:1 / sum(1:n_masking_patterns)
    g2_masking_pattern_proportions = 1:n_masking_patterns / sum(1:n_masking_patterns)
    expected_fc = g1_masking_pattern_proportions / g2_masking_pattern_proportions

    g1_dat = simulate_proteins_for_group(
      n_subjects = n_g1,
      group_id = 1,
      n_reps = n_protein_replicates,
      masking_set = masking_pattern_set,
      masking_proportions = g1_masking_pattern_proportions,
      full_masking = full_masking,
      pk_missed_prop = proteinase_k_missed_prob,
      tr_missed_prop = trypsin_missed_prob,
      sequences = synthetic_protein,
      pk_prob = NULL,
      tr_prob = NULL,
      pk_mu = pk_mu,
      tr_mu = tr_mu)

    g2_dat = simulate_proteins_for_group(
      n_subjects = n_g2,
      group_id = 2,
      n_reps = n_protein_replicates,
      masking_set = masking_pattern_set,
      masking_proportions = g2_masking_pattern_proportions,
      full_masking = full_masking,
      pk_missed_prop = proteinase_k_missed_prob,
      tr_missed_prop = trypsin_missed_prob,
      sequences = synthetic_protein,
      pk_prob = NULL,
      tr_prob = NULL,
      pk_mu = pk_mu,
      tr_mu = tr_mu)

    full_pep_data = get_full_pep_data(g1_dat, g2_dat)
    full_cleavage_mapping = get_full_cleavage_mapping(g1_dat, g2_dat)
    f_data = generate_f_data(full_pep_data)

    #add missingness capabilities here

    site_res = get_site_res(full_pep_data, full_cleavage_mapping, functions[[function_idx]])
    site_plot_res = create_site_plot(site_res, masking_pattern_set)
    site_p = site_plot_res$p +
      theme(
        text = element_text(size = 24)
      )
    site_line_data = site_plot_res$line_data

    tryptic_res = get_tryptic_res(full_pep_data, full_cleavage_mapping)
    tryptic_plot_res = create_tryptic_plot(tryptic_res, masking_pattern_set)
    tryptic_p = tryptic_plot_res$p +
      theme(
        text = element_text(size = 24)
      )
    tryptic_line_data = tryptic_plot_res$line_data
    plot_data = tryptic_plot_res$plot_data
    tmp_res = get_comparison_results(site_line_data, tryptic_line_data, plot_data, full_masking)


    compare_p = tmp_res$p +
      theme(
        text = element_text(size = 24)
      ) +
      scale_linetype_discrete(labels = c("Estimated FC", "log2 FC")) +
      scale_color_viridis_d(labels = c("Estimated FC", "log2 FC"), begin = 0.0, end = 0.7)

    res = tmp_res$res %>%
      mutate(
        replicate = r,
        n_per_group = sim_dat$n_per_group[i],
        n_protein_replicates = n_protein_replicates,
        trypsin_missed_prob = trypsin_missed_prob,
        proteinase_k_missed_prob = proteinase_k_missed_prob,
        overall_missingness = full_pep_data %>% dplyr::select(-peptide, -protein, -start_site) %>% apply(., 1, \(x) mean(is.na(x))) %>% mean(),
        sd_across_peptides = full_pep_data %>% dplyr::select(-peptide, -protein, -start_site) %>% apply(., 1, \(x) mean(is.na(x))) %>% sd(),
        max_masked_region_size = max_masked_region_size,
        gap_parameter = gap_parameter,
        aggregation_method = names(functions)[function_idx]
      ) %>%
      rbind(res, .)
    saveRDS(res, here::here("Out", "sim_res", "lip_rollup_sims.RDS"))
    ggsave(here::here("Out", "plots", paste0("comparison_",r, ".png")), compare_p, width = 14, height = 6)
    ggsave(here::here("Out", "plots", paste0("tryptic_lineplot_",r, ".png")), tryptic_p, width = 14, height = 6)
    ggsave(here::here("Out", "plots", paste0("site_lineplot_",r, ".png")), site_p, width = 14, height = 6)
  }
}

res = readRDS(here::here("Out", "sim_res", "lip_rollup_sims.RDS"))


res %>% ggplot(aes(x = Analysis_Method, y = MSE)) + geom_boxplot() + facet_grid(max_masked_region_size~gap_parameter)



ggdat = res %>% mutate(masked = ifelse(masked, "Masked Region(s)", "Unmasked Region(s)"))
compare_p = ggplot(ggdat, aes(x = masked, y = MSE, color = Analysis_Method)) +
  geom_boxplot() +
  facet_grid(proteinase_k_missed_prob~trypsin_missed_prob) +
  guides(color = guide_legend(title = "Analysis Method")) +
  scale_color_manual(
    values = c("#7AD151FF", '#414487FF'),
    labels = c("Site Rollup", "Tryptic Peptides")
    ) +
  xlab("") +
  theme_bw()
compare_p


ggplot(ggdat, aes(y = overall_missingness, x = n_protein_replicates)) +
  geom_point() + facet_wrap( ~n_per_group) +
  theme_bw()



ggsave(here::here("Out", "plots", paste0("method_compare.png")), compare_p, width = 14, height = 6)
compare_p


# Immediate:
# [x] potential masking patterns and proportions of each within each group
# [x] lower bound on the peptide length
# [x] fix analysis and rollup code to match new simulation code
# - [x] site level
# - [x] peptide level
# [x] Fix PK bug with site level
# [x] Get code in format so its easy to tweak some nobs
# [x] Code to compare between methods (wrt fold change/true fold change)
# [x] Create plots/metrics comparing methods for LiP
# [x] Create plots/metrics comparing methods for PTM
# [x] Run some initial LiP sims
# [ ] Run some initial PTM sims



# Stretch:
# [ ] synthetic label free missingness
# [ ] more sophisticated imperfect digestion




#lineplot code

ggdat = tryptic_res %>%
  select(Peptide_Site,
         Fold_change_1_vs_2) %>%
  rename(FC = Fold_change_1_vs_2) %>%
  mutate(
    Peptide = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][2]),
    Site = as.numeric(sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][1])),
    # protein = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][3]),
    Enzyme = "Trypsin"
  ) %>%
  select(-Peptide_Site)

masking_data = c()
for(idx in seq_along(masking_pattern_set)){
  curr_masking_pattern = list(get_masking_pattern(full_masking, masking_pattern_set[[idx]]))
  tmp_masking_data = data.frame(
    Site = seq_len(nchar(curr_masking_pattern[[1]]$sequence)),
    AA = strsplit(curr_masking_pattern[[1]]$masked_sequence, '')[[1]],
    pattern_id = idx
  ) %>%
    mutate(Site_Type = ifelse(AA == "#", "Masked", "Unmasked"))
  masking_data = rbind(masking_data, tmp_masking_data)
}

all_sites <- data.frame(
  Site = rep(seq(min(masking_data$Site), max(masking_data$Site))), length(masking_pattern_set),
  pattern_id = rep(seq_along(masking_pattern_set), each = max(masking_data$Site)))

plot_data_1 <- all_sites %>%
  left_join(masking_data, by = c("Site", "pattern_id")) %>%
  mutate(Enzyme = "Trypsin")
plot_data_2 <- all_sites %>%
  left_join(masking_data, by = c("Site", "pattern_id")) %>%
  mutate(Enzyme = "Proteinase K")
plot_data = rbind(plot_data_1, plot_data_2) %>%
  left_join(., ggdat, by = c("Site", "Enzyme"))

tryptic_line_data <- plot_data %>%
  filter(!is.na(FC))

plot_data = plot_data %>% filter(Enzyme == "Trypsin")



plot_data = plot_data |>
  mutate(
    Site_Type = ifelse(Site_Type == "Unmasked", Site_Type, ifelse(pattern_id == 1, "Group 1 Masking Pattern", "Group 2 Masking Pattern"))
  )

p1 = ggplot() +
  geom_tile(data = plot_data, aes(x = Site, y = (max(site_line_data$FC, na.rm = TRUE) + min(site_line_data$FC, na.rm = TRUE)) / 2, fill = Site_Type), height = max(site_line_data$FC, na.rm = TRUE) - min(site_line_data$FC, na.rm = TRUE), alpha = 0.2) +
  geom_line(data = tryptic_line_data, aes(x = Site, y = FC)) +
  theme_bw() +
  scale_fill_manual(values = c("Unmasked" = "gray90", "Group 1 Masking Pattern" = "#FDE725FF", "Group 2 Masking Pattern" = "#21918c")) +
  labs(fill = "Region Masked?") +
  ylab("log2 Fold Change")
p1


ggdat_site = site_res %>%
  select(
    Site_Type,
    Fold_change_1_vs_2
  ) %>%
  rename(FC = Fold_change_1_vs_2) %>%
  mutate(
    Site = as.numeric(sapply(Site_Type, \(x) strsplit(x, split = '_')[[1]][1])),
    Enzyme = sapply(Site_Type, \(x) strsplit(x, split = '_')[[1]][2])
  ) %>%
  select(-Site_Type)

masking_data = c()
for(idx in seq_along(masking_pattern_set)){
  curr_masking_pattern = list(get_masking_pattern(full_masking, masking_pattern_set[[idx]]))
  tmp_masking_data = data.frame(
    Site = seq_len(nchar(curr_masking_pattern[[1]]$sequence)),
    AA = strsplit(curr_masking_pattern[[1]]$masked_sequence, '')[[1]],
    pattern_id = idx
  ) %>%
    mutate(Site_Type = ifelse(AA == "#", "Masked", "Unmasked"))
  masking_data = rbind(masking_data, tmp_masking_data)
}

all_sites <- data.frame(
  Site = rep(seq(min(masking_data$Site), max(masking_data$Site))), length(masking_pattern_set),
  pattern_id = rep(seq_along(masking_pattern_set), each = max(masking_data$Site)))

plot_data_1_site <- all_sites %>%
  left_join(masking_data, by = c("Site", "pattern_id")) %>%
  mutate(Enzyme = "Trypsin")
plot_data_2_site <- all_sites %>%
  left_join(masking_data, by = c("Site", "pattern_id")) %>%
  mutate(Enzyme = "Proteinase K")
plot_data_site = rbind(plot_data_1_site, plot_data_2_site) %>%
  left_join(., ggdat_site, by = c("Site", "Enzyme"))

plot_data_site = plot_data_site %>% filter(Enzyme != "Trypsin")
site_line_data <- plot_data_site %>%
  filter(!is.na(FC))

plot_data_site = plot_data_site |>
  mutate(
    Site_Type = ifelse(Site_Type == "Unmasked", Site_Type, ifelse(pattern_id == 1, "Group 1 Masking Pattern", "Group 2 Masking Pattern"))
  )

p2 = ggplot() +
  geom_tile(data = plot_data_site, aes(x = Site, y = (max(site_line_data$FC, na.rm = TRUE) + min(site_line_data$FC, na.rm = TRUE)) / 2, fill = Site_Type), height = max(site_line_data$FC, na.rm = TRUE) - min(site_line_data$FC, na.rm = TRUE), alpha = 0.2) +
  geom_line(data = site_line_data, aes(x = Site, y = FC)) +
  theme_bw() +
  scale_fill_manual(values = c("Unmasked" = "gray90", "Group 1 Masking Pattern" = "#FDE725FF", "Group 2 Masking Pattern" = "#21918c")) +
  labs(fill = "Region Masked?") +
  ylab("log2 Fold Change")
p2


plot_data = plot_data |>
  mutate(
    Site_Type = ifelse(Site_Type == "Unmasked", Site_Type, ifelse(pattern_id == 1, "Group 1 Masking Pattern", "Group 2 Masking Pattern"))
  )

p1 = ggplot() +
  geom_tile(data = plot_data_site, aes(x = Site, y = (max(site_line_data$FC, na.rm = TRUE) + min(site_line_data$FC, na.rm = TRUE)) / 2, fill = Site_Type), height = max(site_line_data$FC, na.rm = TRUE) - min(site_line_data$FC, na.rm = TRUE), alpha = 0.2) +
  geom_line(data = tryptic_line_data, aes(x = Site, y = FC)) +
  theme_bw() +
  scale_fill_manual(values = c("Unmasked" = "gray90", "Group 1 Masking Pattern" = "#FDE725FF", "Group 2 Masking Pattern" = "#21918c")) +
  labs(fill = "Region Masked?") +
  ylab("log2 Fold Change")
p1



library(ggpubr)

plegend = ggpubr::get_legend(p1 + theme(text = element_text(size = 22)))

pmain = ggarrange(
  p1 + theme(legend.position = "none", text = element_text(size = 22)) + ggtitle("Tryptic Peptides"),
  p2 + theme(legend.position = "none", text = element_text(size = 22)) + ggtitle("Proteinase K Sites"),
  nrow = 2
)

pfull = ggarrange(
  pmain,
  plegend,
  widths = c(0.75, 0.25)
)
pfull


ggsave(here::here("Out", "plots", "lineplotfull.png"), pfull, width = 14, height = 7)


