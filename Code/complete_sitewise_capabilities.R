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

sequences = c(generate_random_protein(1000),
              generate_random_protein(1000),
              generate_random_protein(1000),
              generate_random_protein(1000),
              generate_random_protein(1000))


sequence_abundances = rbind(c(1000, 1000),
                            c(750, 750),
                            c(1000, 1000),
                            c(1250, 1250),
                            c(1000, 1000))

n_reps = 30

out_dat = c()
for(rr in seq_len(n_reps)){
  out = get_samples_for_sequences(sequences = sequences, 
                                  sequence_abundances = sequence_abundances, 
                                  samps_per_group = samples_per_group, 
                                  site_sd = 0.1,
                                  subject_sd = 0.1,
                                  amino_acid = "S", 
                                  identifier = "@")
  
  
  pepdat = out$pepdat
  ground_truth = out$ground_truth
  
  # Randomly assign plexes and induce TMT missingness
  missingness_prop = 0.5
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
  
  #Rollup:
  
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
                       left_join(., ground_truth, by = 'site'), 
                     Method = method, 
                     Combine_Function = combine_fn,
                     replicate = rr)
    
    out_dat = rbind(out_dat, res)
  }
  out_dat = as_tibble(out_dat)
}


# Performance for each method compared to ground truth:

ggdat = out_dat %>% mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2) %>% group_by(replicate, Method, Combine_Function) %>% summarize(SEL = mean(SEL))
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  theme_bw()

