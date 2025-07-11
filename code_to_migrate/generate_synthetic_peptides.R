# Packages
library(OrgMassSpecR)
library(foreach)
library(pmartR)
library(seqinr)
library(tidyverse)

# Set diectory and load helper functions
here::i_am("Code/generate_synthetic_peptides.R")
source(here::here("Code", "pmart_rollup_debuggers.R"))
source(here::here("Code", "missing_data_utils.R"))

# Functions to simulate mass spec proteomics
digest_proteins = function(sequence, enzyme = "trypsin") Digest(sequence, enzyme = enzyme)$peptide

digest_n_proteins = function(sequence, n, noise_fn, enzyme = "trypsin"){
  peptide_counts = table(digest_proteins(sequence, enzyme = enzyme)) * n
  peptide_counts_adj = noise_fn(peptide_counts)
  return(data.frame(Count = c(peptide_counts_adj), Peptide = names(peptide_counts)))
} 

digest_sequences = function(sequences, mean_abundances, noise_fn, enzyme = "trypsin"){
  digest_results = c()
  for(i in seq_along(sequences)){
    if(i %% 100 == 0) cat("Sequence ", i, " out of ", length(sequences), '\r')
    sequence = sequences[i]
    mean_abundance = mean_abundances[i]
    digest_results = rbind(digest_results, digest_n_proteins(sequence, mean_abundance, enzyme = enzyme, noise_fn = noise_fn))
    digest_results = digest_results %>% group_by(Peptide) %>% summarise(Count = sum(Count)) %>% ungroup
  }
  return(digest_results)
}

generate_n_samples = function(sequences, mean_abundances, n_samples, noise_fn, samp_id_base = NULL, enzyme = "trypsin"){
  for(samp in seq_len(n_samples)){
    if(samp == 1) {
      peptide_abundances = digest_sequences(sequences, mean_abundances, noise_fn, enzyme = enzyme)
    } else {
      peptide_abundances = full_join(peptide_abundances, digest_sequences(sequences, mean_abundances, noise_fn, enzyme = enzyme), by = 'Peptide')
    }
  }
  if(!is.null(samp_id_base)){
    names(peptide_abundances)[which(names(peptide_abundances) != 'Peptide')] = paste0(samp_id_base, seq_len(n_samples))
  } else {
    names(peptide_abundances)[which(names(peptide_abundances) != 'Peptide')] = paste0("Samp_", seq_len(n_samples))
  }
  return(peptide_abundances)
}

get_peptide_mapping = function(sequences){
  peptide_mapping = c()
  for(i in seq_along(sequences)){
    sequence = sequences[i]
    digest_results = digest_proteins(sequence)
    peptide_mapping = rbind(peptide_mapping, cbind(Peptide = unique(digest_results), Protein = sequences[i]))
    
  }
  return(data.frame(peptide_mapping))
}

abundances_to_probs = function(abundance_dat){
  mdns = log2(abundance_dat %>% apply(., 1, median))
  mdns_cent = mdns - median(mdns)
  probs = exp(-mdns_cent) / sum(exp(-mdns_cent))
  return(probs)
}

generate_label_free_missingness = function(e_data, missingness_prop, e_data_cname){
  abundance_dat = e_data %>% select(-all_of(e_data_cname))
  peptides = e_data %>% pluck(e_data_cname)
  samples = names(e_data)[which(names(e_data) != e_data_cname)]
  probs = abundances_to_probs(abundance_dat = abundance_dat)
  while(mean(is.na(abundance_dat)) < missingness_prop){
    peptide_to_go_missing = sample(x = peptides, size = 1, prob = probs)
    pep_id = which(peptides == peptide_to_go_missing)
    
    tmp_samples = samples[!is.na(abundance_dat[pep_id, ])]
    
    if(identical(tmp_samples, character(0))){
      probs[pep_id] = 0
      probs = probs / sum(probs)
      next
    }
    
    samp_to_go_missing = sample(x = tmp_samples, size = 1)
    abundance_dat[pep_id, which(samples == samp_to_go_missing)] = NA
  }
  e_data[, which(names(e_data) != e_data_cname)] = abundance_dat
  return(e_data)
}

# Pick some sequences
sequence1 = OrgMassSpecR::example.sequence # Human Serum Albumin
sequence2 = "MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"
sequences = c(sequence1, sequence2)

G1_mean_abundances = c(1000, 2000)
G2_mean_abundances = c(1000, 1000)

expected_FC = data.frame(Protein = sequences, Expected_FC = log2(G1_mean_abundances / G2_mean_abundances))

X1 = generate_n_samples(sequences = sequences, 
                   mean_abundances = G1_mean_abundances, 
                   n_samples = 2, 
                   noise_fn = function(peptide_counts) peptide_counts, 
                   samp_id_base = NULL, 
                   enzyme = "trypsin")


X2 = generate_n_samples(sequences = sequences, 
                        mean_abundances = G2_mean_abundances, 
                        n_samples = 2, 
                        noise_fn = function(peptide_counts) peptide_counts, 
                        samp_id_base = NULL, 
                        enzyme = "trypsin")

X = full_join(X1, X2, by = 'Peptide')

# Generate peptide data
pep_e_data = X
pep_f_data = data.frame(Samp = names(pep_e_data)[names(pep_e_data) != "Peptide"], group = rep(c("1", "2"), each = 2))
pep_e_meta = get_peptide_mapping(sequences)

# Make into pmart object
pepdat = pmartR::as.pepData(e_data = pep_e_data,
                            f_data = pep_f_data,
                            e_meta = pep_e_meta,
                            edata_cname = "Peptide",
                            emeta_cname = "Protein",
                            fdata_cname = "Samp")

# Standard transformation and median centering
pepdat = pmartR::edata_transform(pepdat, "log2")
pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = F)
pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
pepdat = pmartR::group_designation(pepdat, main_effects = 'group')

# Rollup
methods = c('rollup', 'rrollup', 'zrollup')
combine_fns = c("mean", "median")
combinations = expand.grid(method = methods, combine_fn = combine_fns)
out_dat = c()
for(i in seq_len(nrow(combinations))){
  method = combinations$method[i]
  combine_fn = combinations$combine_fn[i]
  single_pep = method == "zrollup"
  
  protdat = pmartR::protein_quant(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
  tmp_res = imd_anova(protdat, test_method = 'anova')
  res = data.frame(tmp_res %>% 
                           select(Protein, starts_with("Fold_change")) %>%
                           left_join(., expected_FC, by = 'Protein'), 
                   Method = method, 
                   Combine_Function = combine_fn)
  
  out_dat = rbind(out_dat, res)
}
out_dat = as_tibble(out_dat)

ggdat = out_dat %>% mutate(SEL = (Fold_change_1_vs_2 - Expected_FC)^2)
ggplot(ggdat, aes(x = Method, y = SEL, color = Combine_Function)) + 
  geom_boxplot() + 
  ylab("Squared Error Loss of log2FC Recovery") + 
  labs(color = 'Combine Function') + 
  theme_bw()











# # Load real data
# test_fasta = read.fasta(here::here("Data", "Bt_97-27_GCF_000008505.1_ASM850v1_cont_protein.faa"), seqtype = "AA", as.string = T)
# 
# 
# # Set simulatin parameters
# n_proteins = 200
# mean_abundances = round(seq(50, 20000, length.out = n_proteins))
# n_samples = 10
# noise_fn = function(peptide_counts, retention_rate = 0.8) rbinom(n = length(peptide_counts), size = peptide_counts, prob = retention_rate)
# missingness_vector = seq(0, 0.5, by = 0.25)
# 
# 
# set.seed(42)
# 
# # Subset real data down to desired number of proteins
# sequences_full = sapply(test_fasta, \(x) x[1])
# sequences = sequences_full[sample(seq_along(sequences_full), n_proteins, replace = FALSE)]
# names(sequences) = NULL
# 
# 
# # Loop over missing data proportions
# out_dat = c()
# truth_dat = data.frame(Protein = sequences, Truth = mean_abundances)
# for(jj in seq_along(missingness_vector)){
#   prop_missingness = missingness_vector[jj]
#   cat("Missingness value", prop_missingness, '\n')
#   
#   # Generate peptide data
#   pep_e_data = generate_n_samples(sequences = sequences, 
#                                   mean_abundances = mean_abundances, 
#                                   n_samples = n_samples, 
#                                   noise_fn = noise_fn)
#   pep_f_data = data.frame(Samp = names(pep_e_data)[names(pep_e_data) != "Peptide"], group = "1")
#   pep_e_meta = get_peptide_mapping(sequences)
#   
#   # Induce label free missingness
#   pep_e_data_mi = generate_label_free_missingness(pep_e_data, prop_missingness, e_data_cname = "Peptide")
#   
#   # Make into pmart object
#   pepdat = pmartR::as.pepData(e_data = pep_e_data_mi,
#                               f_data = pep_f_data,
#                               e_meta = pep_e_meta,
#                               edata_cname = "Peptide",
#                               emeta_cname = "Protein",
#                               fdata_cname = "Samp")
#   
#   # Standard transformation and median centering
#   pepdat = pmartR::edata_transform(pepdat, "log2")
#   pepdat = pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = T, backtransform = T)
#   pepdat = pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
#   
#   
#   # Rollup
#   methods = c('rollup', 'rrollup', 'zrollup')
#   combine_fns = c("mean", "median")
#   combinations = expand.grid(method = methods, combine_fn = combine_fns)
#   for(i in seq_len(nrow(combinations))){
#     method = combinations$method[i]
#     combine_fn = combinations$combine_fn[i]
#     single_pep = method == "zrollup"
#     
#     protdat = pmartR::protein_quant(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
#     res = apply(protdat$e_data[-which(names(protdat$e_data) == "Protein")], 1, \(x) median(x, na.rm = T))
#     tmp_res = data.frame(Protein = protdat$e_data$Protein, 
#                          Results = res, 
#                          Method = method, 
#                          Combine_Function = combine_fn,
#                          Missingness = prop_missingness)
#     
#     out_dat = rbind(out_dat, left_join(tmp_res, truth_dat, by = 'Protein'))
#   }
# }
# 
# 
# # Plot results
# ggplot(data = out_dat, aes(x = Truth, y = Results, color = Method, pch = Combine_Function, linetype = Combine_Function)) + 
#   geom_point() + 
#   geom_line() + 
#   stat_function(fun = function(x) log2(x), color = "black") + 
#   xlab("True Abundance") + 
#   ylab("Estimated log2 Abundance") + 
#   facet_grid(~Missingness) + 
#   # ylim(0, 1) + 
#   theme_bw()



