library(OrgMassSpecR)
library(tidyverse)
library(pmartR)

here::i_am("code_to_migrate/ptm_simulation_and_dials.Rmd")

source(here::here("code_to_migrate", "protein_quant2.R"))
source(here::here("code_to_migrate", "ptm_utils.R"))
source(here::here("code_to_migrate", "TMT_simulation", "Code", "utils.R"))

#things to change
samples_per_group_vec = c(10, 30, 50)
n_synth_proteins_vec = c(5, 10)
total_abundance_vec = c(100, 250)
protein_sequence_length_vec = c(100, 200)
site_sd_vec = c(0, 1) # sd on generating normal coefficients for site effects
subject_sd_vec = c(0, 1) # sd on generating normal coefficients for subject effects
prop_to_miss_vec = c(0, 0.25, 0.5)
missingness_prop_vec = c(0, 0.25, 0.5)

# So do we have DA proteins across groups and DA site modifications?
# Can probably confirm this by just simulating some data and checking

#misc
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

n_reps = 2
out_dat = c()

set.seed(42)

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

      # Pre-process data
      pepdat <- assign_plexes(pepdat, samples_per_plex, n_plex)
      pepdat <- pmartR::edata_transform(pepdat, "log2")
      pepdat <- pmartR::normalize_global(pepdat, subset_fn = 'all', norm_fn = "median", apply_norm = TRUE, backtransform = FALSE)

      tryCatch({
        pepdat <- pmartR::applyFilt(pmartR::molecule_filter(pepdat), pepdat, min_num = 2)
        pepdat <- pmartR::applyFilt(pmartR::proteomics_filter(pepdat), pepdat, min_num_peps = 2)
      }, error = function(e) {
        message(sprintf("Error in filtering at rr = %d, ss = %d: %s", rr, ss, e$message))
      })

      pepdat <- pmartR::group_designation(pepdat, main_effects = 'group')

      # Rollup
      methods <- c('rollup', 'rrollup', 'zrollup')
      combine_fns <- c("mean", "median", "sum")
      combinations <- expand.grid(method = methods, combine_fn = combine_fns)

      for (i in seq_len(nrow(combinations))) {
        tryCatch({
          method <- combinations$method[i]
          combine_fn <- combinations$combine_fn[i]
          single_pep <- method == "zrollup"

          sitedat <- protein_quant2(pepdat, method, combine_fn = combine_fn, single_pep = single_pep)
          sitedat <- pmartR::applyFilt(pmartR::imdanova_filter(sitedat), sitedat, min_nonmiss_anova = 2)

          tmp_res <- imd_anova(sitedat, test_method = 'anova') %>%
            mutate(
              protein = sapply(site, \(x) strsplit(x, split = ';')[[1]][2]),
              site = sapply(site, \(x) strsplit(x, split = ';')[[1]][1])
            )

          res <- data.frame(
            tmp_res %>%
              select(site, protein, starts_with("Fold_change")) %>%
              left_join(., ground_truth_site, by = c('site', "protein")),
            Method = method,
            Combine_Function = combine_fn
          )

          out_dat <- rbind(
            out_dat,
            res %>%
              mutate(se = (Fold_change_1_vs_2 - true_fold_change)^2) %>%
              select(-site, -protein, -Fold_change_1_vs_2, -true_fold_change) %>%
              group_by(Method, Combine_Function) %>%
              summarize(rmse = sqrt(mean(se))) %>%
              mutate(
                replicate = rr,
                scenario = ss,
                n_synth_proteins = n_synth_proteins,
                total_abundance = total_abundance,
                site_sd = site_sd,
                protein_sequence_length = protein_sequence_length,
                subject_sd = subject_sd,
                prop_to_miss = prop_to_miss,
                missingness_prop = missingness_prop
              )
          )

          saveRDS(out_dat, here::here("Out", "sim_res", "init_sim_res_ptm3.RDS"))

        }, error = function(e) {
          message(sprintf("Error in combination i = %d, rr = %d, ss = %d: %s", i, rr, ss, e$message))
        })
      }
    }, error = function(e) {
      message(sprintf("Error in outer loop rr = %d, ss = %d: %s", rr, ss, e$message))
    })
  }
}

out_dat <- as_tibble(out_dat)

# Overall Boxplot

# ggdat = out_dat %>%
#   mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2) %>%
#   group_by(replicate, Method, Combine_Function) %>%
#   summarize(RMSE = sqrt(mean(SEL)))
# p = ggplot(ggdat, aes(x = Method, y = RMSE, color = Combine_Function)) +
#   geom_boxplot() +
#   ylab("RMSE of log2FC Recovery") +
#   labs(color = 'Combine Function') +
#   ylim(0, NA) +
#   theme_bw()
# p

# Sitewise Lineplot

# ggdat = out_dat %>%
#   mutate(SEL = (Fold_change_1_vs_2 - true_fold_change)^2) %>%
#   group_by(site, Method, Combine_Function) %>%
#   summarize(RMSE = sqrt(mean(SEL))) %>%
#   mutate(site_num = as.numeric(substr(site, 3, nchar(site))))
#
# ggplot(ggdat, aes(x = site_num, y = RMSE, color = Combine_Function)) +
#   geom_line() +
#   facet_wrap(~Method) +
#   theme_bw()


# Sitewise Boxplot

# ggdat = out_dat %>%
#   mutate(AE = abs(Fold_change_1_vs_2 - true_fold_change)) %>%
#   group_by(site, replicate, Method, Combine_Function) %>%
#   mutate(site_num = as.numeric(substr(site, 3, nchar(site))))
#
# ggplot(ggdat, aes(x = as.factor(site_num), y = AE, color = Combine_Function)) +
#   geom_boxplot() +
#   facet_wrap(~Method) +
#   theme_bw()

