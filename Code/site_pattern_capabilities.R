library(OrgMassSpecR)
library(tidyverse)
library(pmartR)

here::i_am("Code/ptm_simulation_and_dials.Rmd")

source(here::here("Code", "protein_quant2.R"))
source(here::here("Code", "ptm_utils.R"))
source(here::here("Code", "TMT_simulation", "Code", "utils.R"))

samples_per_group = 20
samples_per_plex = 4
n_plex = (2 * samples_per_group) / samples_per_plex
site_sd = 0.1
subject_sd = 0.1
prop_splits_to_miss = 0.1
missingness_prop = 0.2

n_reps = 30
for(rr in seq_len(n_reps)){
  cat("Replicate", rr, "\n")
  
  sequences = c(generate_random_protein(1000),
                generate_random_protein(1000),
                generate_random_protein(1000))#,
                # generate_random_protein(1000),
                # generate_random_protein(1000))
  
  
  sequence_abundances = rbind(c(1000, 1000),
                              c(750, 750),
                              c(1000, 1000))#,
                              # c(1250, 1250),
                              # c(1000, 1000))
  
  
  out = get_samples_for_sequences(sequences = sequences, 
                                  sequence_abundances = sequence_abundances, 
                                  samps_per_group = samples_per_group, 
                                  site_sd = site_sd,
                                  subject_sd = subject_sd,
                                  prop_to_miss = prop_splits_to_miss,
                                  amino_acid = "S", 
                                  identifier = "@")
  
  pepdat = out$pepdat
  ground_truth = out$ground_truth_site
  ground_truth_sp = out$ground_truth_site_pattern
  
  e_meta_site_pattern = out$pepdat$e_meta %>% 
    separate(., 'site', into = c('site', 'protein'), sep = ';') %>%
    group_by(peptide, protein) %>%
    summarise(site_pattern = paste(site, collapse = ';')) %>%
    mutate(site_pattern = paste(site_pattern, protein, sep = ';'))
  
  pepdat_sp = as.pepData(e_data = out$pepdat$e_data,
                         f_data = out$pepdat$f_data,
                         e_meta = e_meta_site_pattern,
                         edata_cname = 'peptide',
                         fdata_cname = 'sample',
                         emeta_cname = 'site_pattern')
  
  rm(out)
  gc()
  
  # Randomly assign plexes and induce TMT missingness
  pepdat = assign_plexes(pepdat, samples_per_plex, n_plex)
  missingness_bins = get_missingness_bins(pepdat)
  pepdat = simulate_missingness(pepdat, 
                                missingness_prop, 
                                missingness_bins$binned_peps, 
                                missingness_bins$bin_probs,
                                verbose = F)
  
  
  # Pre-process data:
  pepdat = pmartR::edata_transform(pepdat, "log2")
  pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
  pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
  pepdat = pmartR::group_designation(pepdat, main_effects = 'group')
  pepdat = pmartR::applyFilt(pmartR::imdanova_filter(pepdat), pepdat, min_nonmiss_anova = 2)
  
  #Site Rollup:
  out_dat_site = c()
  methods = c('rollup', 'rrollup', 'zrollup')
  combine_fns = c("mean", "median", "sum")
  combinations = expand.grid(method = methods, combine_fn = combine_fns)
  for(i in seq_len(nrow(combinations))){
    method = combinations$method[i]
    combine_fn = combinations$combine_fn[i]
    single_pep = method == "zrollup"
    
    sitedat = protein_quant2(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
    tmp_res = imd_anova(sitedat, test_method = 'anova')
    res = data.frame(tmp_res %>% 
                       select(site, starts_with("Fold_change")) %>%
                       left_join(., ground_truth %>% mutate(site = paste(site, protein, sep = ';')), by = 'site'), 
                     Method = method, 
                     Combine_Function = combine_fn,
                     replicate = rr)
    
    out_dat_site = rbind(out_dat_site, res)
  }
  out_dat_site = as_tibble(out_dat_site)
  
  #Site Pattern Rollup:
  # pepdat$e_meta %>% mutate(site = strsplit(site, split = ";")[[1]][1])
  
  rm(pepdat)
  gc()
  
  pepdat_sp = pmartR::edata_transform(pepdat_sp, "log2")
  pepdat_sp = pmartR::normalize_global(pepdat_sp, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
  pepdat_sp = pmartR::applyFilt(pmartR::molecule_filter(pepdat_sp), pepdat_sp, min_num = 2)
  pepdat_sp = pmartR::group_designation(pepdat_sp, main_effects = 'group')
  pepdat_sp = pmartR::applyFilt(pmartR::imdanova_filter(pepdat_sp), pepdat_sp, min_nonmiss_anova = 2)
  
  out_dat_site_pattern = c()
  methods = c('rollup', 'rrollup', 'zrollup')
  combine_fns = c("mean", "median", "sum")
  combinations = expand.grid(method = methods, combine_fn = combine_fns)
  for(i in seq_len(nrow(combinations))){
    method = combinations$method[i]
    combine_fn = combinations$combine_fn[i]
    single_pep = method == "zrollup"
    
    sitedat = protein_quant2(pepdat_sp, method, combine_fn = combine_fn, single_pep = single_pep)
    tmp_res = imd_anova(sitedat, test_method = 'anova')
    res = data.frame(tmp_res %>% 
                       select(site_pattern, starts_with("Fold_change")) %>%
                       left_join(., ground_truth_sp %>% mutate(site_pattern = paste(site_pattern, protein, sep = ';')), by = 'site_pattern'), 
                     Method = method, 
                     Combine_Function = combine_fn,
                     replicate = rr)
    
    out_dat_site_pattern = rbind(out_dat_site_pattern, res)
  }
  out_dat_site_pattern = as_tibble(out_dat_site_pattern)
  
  rm(pepdat_sp)
  gc()
  
  # Performance for each method compared to ground truth:
  saveRDS(list(out_dat_site = out_dat_site, out_dat_site_pattern = out_dat_site_pattern), 
          here::here("Results", "ptm_test_sim.RDS"))
}

ggdat = out_dat_site %>% mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2) %>% group_by(replicate, Method, Combine_Function) %>% summarize(SEL = mean(SEL))
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  ggtitle("Site-wise rollup")
theme_bw()

ggdat = out_dat_site_pattern %>% mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2) %>% group_by(replicate, Method, Combine_Function) %>% summarize(SEL = mean(SEL))
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  ggtitle("Site-pattern rollup") + 
  theme_bw()
