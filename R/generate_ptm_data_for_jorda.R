library(OrgMassSpecR)
library(tidyverse)
library(pmartR)

here::i_am("code_to_migrate/ptm_simulation_and_dials.Rmd")

source(here::here("code_to_migrate", "protein_quant2.R"))
source(here::here("code_to_migrate", "ptm_utils.R"))
source(here::here("code_to_migrate", "TMT_simulation", "Code", "utils.R"))

samples_per_group_vec = c(50)
n_synth_proteins_vec = c(10)
total_abundance_vec = c(100)
protein_sequence_length_vec = c(1000)
site_sd_vec = c(0.1)
subject_sd_vec = c(0.1)
prop_to_miss_vec = c(0)
missingness_prop_vec = c(0)

samples_per_plex = 4

sim_scenarios = expand.grid(
  n_synth_proteins = n_synth_proteins_vec,
  total_abundance = total_abundance_vec,
  protein_sequence_length = protein_sequence_length_vec,
  site_sd = site_sd_vec,
  subject_sd = subject_sd_vec,
  prop_to_miss = prop_to_miss_vec,
  missingness_prop = missingness_prop_vec,
  samples_per_group = samples_per_group_vec
)

n_reps = 70
out_dat = c()

set.seed(87)

for (rr in seq_len(n_reps)) {
  for (ss in seq_len(nrow(sim_scenarios))) {
    tryCatch({
      # Set scenario parameters
      n_synth_proteins <- sim_scenarios$n_synth_proteins[ss]
      total_abundance <- sim_scenarios$total_abundance[ss]
      site_sd <- sim_scenarios$site_sd[ss]
      protein_sequence_length <- sim_scenarios$protein_sequence_length[ss]
      subject_sd <- sim_scenarios$subject_sd[ss]
      prop_to_miss <- sim_scenarios$prop_to_miss[ss]
      missingness_prop <- sim_scenarios$missingness_prop[ss]
      samples_per_group <- sim_scenarios$samples_per_group[ss]
      sequences <- sapply(seq_len(n_synth_proteins), \(x) generate_random_protein(protein_sequence_length))
      n_plex = (2 * samples_per_group) / samples_per_plex

      # Get protein data
      out <- tryCatch({
        get_samples_for_sequences(
          sequences = sequences,
          sequence_abundances = rbind(
            round(rbeta(n_synth_proteins, 3, 3) * total_abundance),
            round(rbeta(n_synth_proteins, 3, 3) * total_abundance)
          ) %>% t(),
          samps_per_group = samples_per_group,
          site_sd = site_sd,
          subject_sd = subject_sd,
          prop_to_miss = prop_to_miss,
          amino_acid = "S",
          identifier = "@"
        )
      }, error = function(e) {
        message(sprintf("Error in get_samples_for_sequences at rr = %d, ss = %d: %s", rr, ss, e$message))
        NULL
      })

      if (is.null(out)) next

      pepdat <- out$pepdat
      ground_truth_site <- out$ground_truth_site
      ground_truth_site_pattern <- out$ground_truth_site_pattern

      pepdat <- assign_plexes(pepdat, samples_per_plex, n_plex)
      pepdat <- pmartR::edata_transform(pepdat, "log2")
      pepdat <- pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = TRUE, backtransform = FALSE)

      tryCatch({
        pepdat <- pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
        pepdat <- pmartR::applyFilt(pmartR::proteomics_filter(pepdat), pepdat, min_num_peps = 2)
      }, error = function(e) {
        message(sprintf("Error in filtering at rr = %d, ss = %d: %s", rr, ss, e$message))
      })



      saveRDS(list(
        pep_data = pepdat,
        ground_truth = ground_truth_site
      ),
      here::here("synthetic_data", paste0("ptm_", ss, "_", rr, ".rds"))
      )
      #saveRDS(ground_truth_site_pattern, here::here("ground_truth_pattern", paste0("ptm_", ss, "_", rr, ".rds")))

    }, error = function(e) {
      message(sprintf("Error in outer loop rr = %d, ss = %d: %s", rr, ss, e$message))
    })
  }
}
