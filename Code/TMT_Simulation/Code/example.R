library(pmartR)
library(tidyverse)

here::i_am("Code/example.R")
source(here::here("Code", "utils.R"))

# Load Data
e_data = readRDS(here::here("Data", "e_data.RDS"))
f_data = readRDS(here::here("Data", "f_data.RDS"))
e_meta = readRDS(here::here("Data", "e_meta.RDS"))

# Define some required objects
edata_cname = "Sequence"    
fdata_cname = "Subject_ID"  
emeta_cname = "Proteins"    

groupvar_cname = "Ecoli_FR" 
samps_per_plex = 4          
plex_cname = "Plex"         

# Create pmart object
pep_dat = as.isobaricpepData(e_data, 
                             f_data, 
                             e_meta, 
                             edata_cname = edata_cname,
                             fdata_cname = fdata_cname,
                             emeta_cname = emeta_cname,
                             isobaric_norm = T)

# Simulate plex
set.seed(349582)
pep_dat$f_data$Plex = c(sample(seq_len(5)), sample(seq_len(5)), sample(seq_len(5)), sample(seq_len(5)))

#### Log2 scaling and median normalization
pep_dat = edata_transform(pep_dat, 'log2')
pep_dat = normalize_global(omicsData = pep_dat,
                           subset_fn = "all",
                           norm_fn = "median",
                           apply_norm = T)

#### Get missingness bins
MI_bins = get_missingness_bins(pep_dat)
bin_medians_QE_B5 = MI_bins$bin_medians
binned_peps_QE_B5 = MI_bins$binned_peps
bin_probs_QE_B5 = MI_bins$bin_probs

#### Filter down to zero missingness
pep_dat = applyFilt(molecule_filter(pep_dat), pep_dat, min_num = 2)
pep_dat = applyFilt(proteomics_filter(omicsData = pep_dat), pep_dat, redundancy = T)
pep_dat = group_designation(pep_dat, main_effects = "Ecoli_FR")
pep_dat = applyFilt(imdanova_filter(omicsData = pep_dat), pep_dat, min_nonmiss_anova = 2)
pep_dat_complete <- applyFilt(molecule_filter(pep_dat), pep_dat, min_num = get_data_info(pep_dat)$num_samps)

#Update bins to only include those peptides that have no missingness
binned_peps_QE_B5_nomiss = binned_peps_QE_B5
for(i in seq_len(10)){
  binned_peps_QE_B5_nomiss[[i]] = binned_peps_QE_B5_nomiss[[i]][binned_peps_QE_B5_nomiss[[i]]  %in% pep_dat_complete$e_data[,1]]
}

# Specify missingness mechanism and total proportion of data to go missing
TMT = T
target_missingness = 0.2

# Simulate missingness
pep_dat_w_missingness = simulate_missingness(pepObj = pep_dat_complete, 
                                             prop_missingness = target_missingness, 
                                             binned_peps = binned_peps_QE_B5_nomiss, 
                                             bin_probs = bin_probs_QE_B5,
                                             verbose = F,
                                             p_adjust = T,
                                             TMT = TMT)

# Plot for sanity check
pd = get_plex_dist(pep_dat_w_missingness, plex_cname = plex_cname, samps_per_plex = samps_per_plex)
plot(pd)





