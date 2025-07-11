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

sim_dat = expand.grid(
  n_per_group = n_per_group_vec,
  trypsin_missed_prob = trypsin_missed_prob_vec,
  proteinase_k_missed_prob = proteinase_k_missed_prob_vec,
  n_protein_replicates = n_protein_replicates_vec,
  max_masked_region_size = max_masked_region_size_vec,
  gap_parameter = gap_parameter_vec
)

n_reps = 30

sim_dat |>
  mutate(
    scenario_id = 1
  )

#iterate though replicates
res = c()
r = 1
while(r <=n_reps){
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

    out = list(
      pep_data = full_pep_data,
      cleavage_mapping = full_cleavage_mapping,
      metadata = f_data,
      sim_data = tibble(scenario_id = i, replicate = r),
      ground_truth = expected_fc
      )

    saveRDS(out, here::here("synthetic_data", paste0("lip_", i, "_", r, ".RDS")))
  }
  r = r + 1
}



