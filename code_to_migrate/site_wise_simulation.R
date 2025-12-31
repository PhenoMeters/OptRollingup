library(OrgMassSpecR)
library(tidyverse)

here::i_am("Code/site_wise_simulation.R")

add_char_at_n = function(original_string, new_char, n) {
  # Split the string into two parts
  part1 = substr(original_string, 1, n - 1)
  part2 = substr(original_string, n, nchar(original_string))
  
  # Concatenate the parts with the new character
  modified_string = paste0(part1, new_char, part2)
  
  return(modified_string)
}

generate_random_protein = function(n) {
  amino_acids = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  random_sequence = sample(amino_acids, size = n, replace = TRUE)
  random_sequence_string = paste(random_sequence, collapse = "")
  return(random_sequence_string)
}

add_ptms = function(sequence, amino_acid, identifier) {
  modified_string = str_replace_all(sequence, amino_acid, paste0(amino_acid, identifier))
  return(modified_string)
}

get_ptm_locations = function(sequence, amino_acid){
  return(stringr::str_locate_all(sequence, amino_acid)[[1]][,1])
}

get_ptm_site_mapping = function(sequence, amino_acid, drop_na = TRUE){
  digest_results = Digest(sequence)
  ptm_locations = get_ptm_locations(sequence, amino_acid)
  tmp_df = c()
  for(i in seq_along(ptm_locations)){
    peptide_id = which(ptm_locations[i] >= digest_results$start & ptm_locations[i] <= digest_results$stop)
    tmp_df = rbind(tmp_df, c(peptide = digest_results$peptide[peptide_id], site = paste0("@", amino_acid, ptm_locations[i])))
  }
  out = digest_results %>% left_join(., as_tibble(tmp_df), by = 'peptide') %>% mutate(protein = sequence)  %>% select(site, peptide, protein, start, stop)
  if(drop_na) out = out %>% filter(!is.na(site))
  return(out)
}

get_modify_indices = function(e_meta, sequence_abundance){
  unique_sites = unique(e_meta$site)
  num_sites = length(unique_sites)
  ptm_abundances = sample(seq_len(sequence_abundance), num_sites)
  modify_indices = sapply(ptm_abundances, \(x) seq_len(sequence_abundance) %in% sample(seq_len(sequence_abundance), x))
  colnames(modify_indices) = unique_sites
  return(modify_indices)
}

ptm_at_site = function(peptide, site, identifier, start){ 
  location = as.numeric(substr(site, 3, nchar(site))) - start + 2
  return(add_char_at_n(peptide, identifier, location))
}

ptm_at_sites = function(peptide, sites, identifier, start){
  for(site in sites){
    peptide = ptm_at_site(peptide, site, identifier = identifier, start = start)
    start = start - 1
  }
  return(peptide)
}

modify_and_count = function(e_meta, modify_indices){
  peps = unique(e_meta$peptide)
  out_dat = c()
  for(i in seq_along(peps)){
    pep_e_meta = e_meta %>% filter(peptide == peps[i])
    tmp_site = pep_e_meta %>% pull(site)
    site_modifications = matrix(modify_indices[,colnames(modify_indices) %in% tmp_site], ncol = length(tmp_site))
    colnames(site_modifications) = tmp_site
    site_modification_counts = table(apply(site_modifications, 1, \(x) paste0(colnames(site_modifications)[x], collapse = ";")))
    #drop no modification case here
    site_modification_counts = site_modification_counts[names(site_modification_counts) != '']
    
    for(j in seq_along(site_modification_counts)){
      curr_ptm = site_modification_counts[j]
      sites = strsplit(names(curr_ptm), split = ';')[[1]]
      ptm = ptm_at_sites(peps[i], sites, "#", pep_e_meta %>% filter(site %in% sites) %>% pull(start) %>% unique)
      out_dat = rbind(out_dat, tibble(peptide = ptm, site = sites, peptide_count = curr_ptm))
    }
  }
  return(out_dat)
}

#random subject effects
#random group effects
#random site effects

#function to generate data for one protein per subject per group given

generate_ptm_samples_for_sequence = function(ptm_site_mapping,
                                             sequence_abundances, 
                                             subject_coefs, 
                                             site_coefs,
                                             group_ids){
  #assert length(sequence_abundances) = length(subject_coef)
  
  e_data_full = c()
  e_meta_full = c()
  ground_truth = c()
  for(j in seq_along(sequence_abundances)){
    modify_indices = get_modify_indices(ptm_site_mapping, sequence_abundances[j]) #sequence should depend on group
    ground_truth = cbind(ground_truth, colSums(modify_indices))
    
    for(i in seq_len(sum(group_ids == 1))){
      e_data = modify_and_count(ptm_site_mapping, modify_indices)
      e_meta = e_data %>% rename(true_peptide_count = peptide_count)
      
      e_data_subj = e_data %>% mutate(peptide_count = peptide_count * exp(subject_coefs[i]))
      names(e_data_subj)[3] = paste0("Subject_", i, "_G", j)
      
      if(i == 1 & j == 1){
        e_data_full = e_data_subj
      } else {
        e_data_full = e_data_full %>% 
          full_join(e_data_subj, by = c("peptide", "site")) %>%
          mutate_if(is.numeric,coalesce, 0)
      }
      
      if(i == 1) e_meta_full = rbind(e_meta_full, e_meta %>% mutate(group = j))
    }
  }
  
  coef_dat = tibble(Mult_Factor = exp(site_coefs), site = unique(e_meta_full$site))
  e_data_full = e_data_full %>% 
    left_join(., coef_dat) %>% 
    select(-site) %>%
    group_by(peptide) %>% 
    summarise(Mult_Factor = prod(Mult_Factor),
              across(starts_with("Subject"), ~ mean(.x))) %>% 
    mutate(across(starts_with("Subject"), ~ . * Mult_Factor)) %>% 
    select(-Mult_Factor)
  
  return(list(e_data = e_data_full, e_meta = e_meta_full, ground_truth = ground_truth))
}

# set.seed(25)
amino_acid = "S"
identifier = "@"

sequence = generate_random_protein(300)
sequence_abundances = c(100, 100)

ptm_site_mapping = get_ptm_site_mapping(sequence, amino_acid)
site_sd = 0.1
site_coefs = rep(0, length(unique(ptm_site_mapping$site)))


n_per_group = 2
subject_effect_sd = 0.1
# subject_coefs = rnorm(n_per_group * 2, sd = subject_effect_sd)
subject_coefs = rep(0, n_per_group * 2)


out = generate_ptm_samples_for_sequence(ptm_site_mapping,
                                        sequence_abundances,
                                        subject_coefs,
                                        site_coefs,
                                        group_ids = rep(c(1, 2), each = n_per_group))

ground_truth = data.frame(site = rownames(out$ground_truth), 
                          G1 = out$ground_truth[,1], 
                          G2 = out$ground_truth[,2]) %>%
  mutate(true_fold_change = log2(G1 / G2)) %>%
  select(site, true_fold_change)

f_data_full = data.frame(sample = out$e_data %>% select(-peptide) %>% names,
                         group = rep(c(1, 2), each = n_per_group))


pepdat = as.pepData(e_data = out$e_data,
                    f_data = f_data_full,
                    e_meta = out$e_meta %>% select(-group, -true_peptide_count) %>% distinct(),
                    edata_cname = "peptide",
                    fdata_cname = "sample",
                    emeta_cname = "site")

pepdat = pmartR::edata_transform(pepdat, "log2")
pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
pepdat = pmartR::group_designation(pepdat, main_effects = 'group')

sitedat = protein_quant(pepdat, method = 'rrollup')
sitedat


methods = c('rollup', 'rrollup', 'zrollup')
combine_fns = c("mean", "median")
combinations = expand.grid(method = methods, combine_fn = combine_fns)
out_dat = c()
for(i in seq_len(nrow(combinations))){
  method = combinations$method[i]
  combine_fn = combinations$combine_fn[i]
  single_pep = method == "zrollup"
  
  sitedat = pmartR::protein_quant(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
  tmp_res = imd_anova(sitedat, test_method = 'anova')
  res = data.frame(tmp_res %>% 
                     select(site, starts_with("Fold_change")) %>%
                     left_join(., ground_truth, by = 'site'), 
                   Method = method, 
                   Combine_Function = combine_fn)
  
  out_dat = rbind(out_dat, res)
}
out_dat = as_tibble(out_dat)

ggdat = out_dat %>% mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2)
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  theme_bw()



####################################################################################
####################################################################################
# Old Code
####################################################################################
####################################################################################












# set.seed(25)
amino_acid = "S"
identifier = "@"

test_seq = generate_random_protein(300)
ptm_labels = get_ptm_site_mapping(test_seq, amino_acid)

#simple case for getting peptide abunances
sequence_abundance = 100
e_meta = ptm_labels# %>% mutate(peptide = add_ptms(peptide, "S", "@"))


# set.seed(1)
modify_indices_G1 = get_modify_indices(e_meta, sequence_abundance)
e_data_truth_G1 = modify_and_count(e_meta, modify_indices_G1)
e_meta_ptm_G1 = e_data_truth_G1 %>% select(-peptide_count)

# set.seed(1)
modify_indices_G2 = get_modify_indices(e_meta, sequence_abundance)
e_data_truth_G2 = modify_and_count(e_meta, modify_indices_G2)
e_meta_ptm_G2 = e_data_truth_G2 %>% select(-peptide_count)

e_data_truth = e_data_truth_G1 %>% 
  full_join(., e_data_truth_G1, by = c("peptide", "site")) %>% 
  full_join(., e_data_truth_G2, by = c("peptide", "site"))  %>% 
  full_join(., e_data_truth_G2, by = c("peptide", "site"))  %>%
  mutate_if(is.numeric,coalesce, 0)
e_meta_full = rbind(e_meta_ptm_G1, e_meta_ptm_G2) %>% distinct()

names(e_data_truth) = c("peptide", "site", "Samp1_1", "Samp2_1", "Samp1_2", "Samp2_2")

coef_sd = 0.1
coefs = rnorm(unique(e_meta_full$site), sd = coef_sd)
coef_dat = tibble(Mult_Factor = exp(coefs), site = unique(e_meta_full$site))

e_data = e_data_truth %>% 
  left_join(., coef_dat) %>% 
  select(-site) %>%
  group_by(peptide) %>% 
  summarise(Mult_Factor = prod(Mult_Factor),
            across(starts_with("Samp"), ~ mean(.x))) %>% 
  mutate(across(starts_with("Samp"), ~ . * Mult_Factor)) %>% 
  select(-Mult_Factor)


f_data_full = data.frame(sample = c("Samp1_1", "Samp2_1", "Samp1_2", "Samp2_2"), group = c("G1", "G1", "G2", "G2"))

library(pmartR)
pepdat = as.pepData(e_data = e_data,
                    f_data = f_data_full,
                    e_meta = e_meta_full,
                    edata_cname = "peptide",
                    fdata_cname = "sample",
                    emeta_cname = "site")

# Standard transformation and median centering
pepdat = pmartR::edata_transform(pepdat, "log2")
pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
pepdat = pmartR::group_designation(pepdat, main_effects = 'group')

sitedat = protein_quant(pepdat, method = 'rollup')
sitedat

ground_truth = log2(colSums(modify_indices_G1) / colSums(modify_indices_G2))
ground_truth_df = data.frame(site = names(ground_truth), Expected_fold_change = ground_truth)

res = imd_anova(sitedat, test_method = "anova")
res %>% 
  select(site, Fold_change_G1_vs_G2) %>% 
  left_join(., ground_truth_df) %>% view

