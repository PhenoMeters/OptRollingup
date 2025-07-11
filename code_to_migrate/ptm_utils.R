
#################################
### Sequence utils
#################################

add_char_at_n = function(original_string, new_char, n) {
  # Split the string into two parts
  part1 = substr(original_string, 1, n - 1)
  part2 = substr(original_string, n, nchar(original_string))

  # Concatenate the parts with the new character
  modified_string = paste0(part1, new_char, part2)

  return(modified_string)
}

create_amino_acid_distribution = function(amino_acids){
  return(table(amino_acids) / length(amino_acids))
}

generate_protein = function(n, amino_acid_distribution = NULL){
  if(is.null(amino_acid_distribution)) return(generate_random_protein(n))
  return(generate_protein_from_distribution(n, amino_acid_distribution))
}

generate_protein_from_distribution = function(n, amino_acid_distribution){
  paste(sample(names(amino_acid_distribution), size = n, replace = TRUE, prob = amino_acid_distribution), collapse = '')
}

generate_random_protein = function(n) {
  amino_acids = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  random_sequence = sample(amino_acids, size = n, replace = TRUE)
  random_sequence_string = paste(random_sequence, collapse = "")
  return(random_sequence_string)
}

#################################
### Site mapping utils
#################################

get_ptm_site_mapping = function(sequence, digest_results, amino_acid, drop_na = TRUE){
  ptm_locations = get_ptm_locations(sequence, amino_acid)
  tmp_df = c()
  for(i in seq_along(ptm_locations)){
    peptide_id = which(ptm_locations[i] >= as.numeric(digest_results$start) & ptm_locations[i] <= as.numeric(digest_results$stop))
    tmp_df = rbind(tmp_df, c(peptide = digest_results$peptide[peptide_id], site = paste0("@", amino_acid, ptm_locations[i])))
  }
  out = digest_results %>%
    full_join(., as_tibble(tmp_df), by = 'peptide') %>%
    mutate(protein = sequence)  %>%
    select(site, peptide, protein, start, stop) %>%
    mutate(site = replace_na(site, ""))
  # if(drop_na) out = out %>% filter(!is.na(site))
  return(out)
}

get_ptm_locations = function(sequence, amino_acid){
  return(stringr::str_locate_all(sequence, amino_acid)[[1]][,1])
}

add_ptms = function(sequence, amino_acid, identifier) {
  modified_string = str_replace_all(sequence, amino_acid, paste0(amino_acid, identifier))
  return(modified_string)
}

#################################
### Imperfect Digestion
#################################

add_peptide_effect = function(ptm_site_mapping, peptide_coef = NULL) {
  if(is.null(peptide_coef)){
    unique_peptides = unique(ptm_site_mapping$peptide)
    peptide_coef = rep(0, length(unique_peptides))
    names(peptide_coef) = unique_peptides
  }
  coef_dat = data.frame(peptide = names(peptide_coef), Mult_Factor = exp(peptide_coef))
  ptm_site_mapping = ptm_site_mapping %>% left_join(., coef_dat, by = "peptide")
  return(ptm_site_mapping)
}

choose_peptides_to_merge = function(ptm_site_mapping, prop_to_miss) {
  ptm_site_mapping = ptm_site_mapping %>% mutate(
    inverted_counts = (1 / nchar(peptide)),
    probs = inverted_counts / sum(inverted_counts)
  )
  n_peps = nrow(ptm_site_mapping)
  ptm_site_mapping = ptm_site_mapping %>%
    mutate(merge_with_next = (seq_len(n_peps) %in% sample(seq_len(n_peps)[-n_peps], round(prop_to_miss * (n_peps - 1)), prob = probs[-n_peps]))) %>%
    select(-inverted_counts, -probs)
  return(ptm_site_mapping)
}

merge_peptides = function(ptm_site_mapping) {
  remove_vec = rep(FALSE, nrow(ptm_site_mapping))
  ptm_site_mapping = ptm_site_mapping %>% as_tibble()

  for(i in seq_len(nrow(ptm_site_mapping) - 1)){
    if(ptm_site_mapping$merge_with_next[i]){

      #merge current values with next
      peptide = paste0(ptm_site_mapping$peptide[i], ptm_site_mapping$peptide[i + 1])
      start = ptm_site_mapping$start[i]
      stop = ptm_site_mapping$stop[i + 1]
      remove_vec[i] = TRUE
      merge_with_next = ptm_site_mapping$merge_with_next[i+1]
      site = paste(ptm_site_mapping$site[i:(i + 1)][ptm_site_mapping$site[i:(i + 1)] != ""], collapse = ";")

      #set current values
      ptm_site_mapping$peptide[i + 1] = peptide
      ptm_site_mapping$merge_with_next[i + 1] = merge_with_next
      ptm_site_mapping$start[i + 1] = start
      ptm_site_mapping$stop[i + 1] = stop
      ptm_site_mapping$site[i + 1] = site

      #merge peptide effects if introduced into the perfect ptm site mapping
      if ("Mult_Factor" %in% names(ptm_site_mapping)){
        mult_factor = ptm_site_mapping$Mult_Factor[i] * ptm_site_mapping$Mult_Factor[i + 1]
        ptm_site_mapping$Mult_Factor[i + 1] = mult_factor
      }
    }
  }
  return(ptm_site_mapping[!remove_vec,])
}

imperfect_digest = function(ptm_site_mapping, prop_to_miss = 0){
  ptm_site_mapping = choose_peptides_to_merge(ptm_site_mapping, prop_to_miss)
  ptm_site_mapping = merge_peptides(ptm_site_mapping)
  return(ptm_site_mapping)
}

#################################
### Functions to get site and peptide counts
#################################

get_modify_indices = function(ptm_site_mapping, sequence_abundance){
  # browser()
  unique_sites = unique(ptm_site_mapping %>% filter(site != "") %>% pull(site))
  num_sites = length(unique_sites)
  ptm_abundances = sample(seq_len(sequence_abundance), num_sites, replace = T)
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

get_peptide_counts = function(ptm_site_mapping, modify_indices){
  if (!("Mult_Factor" %in% names(ptm_site_mapping))) ptm_site_mapping = add_peptide_effect(ptm_site_mapping)

  peptides = unique(ptm_site_mapping$peptide)
  out_dat = c()
  for(i in seq_along(peptides)){
    pep_e_meta = ptm_site_mapping %>% filter(peptide == peptides[i]) %>% distinct()
    tmp_site_pattern = pep_e_meta %>% pull(site)
    tmp_sites = pep_e_meta %>% pull(site) %>% strsplit(., split = ';') %>% unlist()


    site_modifications = matrix(modify_indices[,colnames(modify_indices) %in% tmp_sites], ncol = length(tmp_sites))
    # site_modifications = matrix(modify_indices[,sapply(colnames(modify_indices), \(x) grepl(x, tmp_site))], ncol = length(tmp_site))

    colnames(site_modifications) = tmp_sites
    site_modification_counts = table(apply(site_modifications, 1, \(x) paste0(colnames(site_modifications)[x], collapse = ";")))
    #drop no modification case here
    site_modification_counts = site_modification_counts[names(site_modification_counts) != '']

    for(j in seq_along(site_modification_counts)){
      curr_ptm = site_modification_counts[j]
      sites = strsplit(names(curr_ptm), split = ';')[[1]]
      site_pattern =
      ptm = ptm_at_sites(peptides[i], sites, "#", pep_e_meta %>% filter(site %in% tmp_site_pattern) %>% pull(start) %>% unique %>% as.numeric())
      out_dat = rbind(
        out_dat,
        tibble(
          peptide = ptm,
          site = sites,
          site_pattern = paste(tmp_site_pattern, collapse = ";"),
          peptide_count = curr_ptm * ptm_site_mapping$Mult_Factor[i]
          )
        )
    }
  }
  return(out_dat)
}

#################################
### Uncategorized
#################################



generate_ptm_samples_for_sequence = function(ptm_site_mapping,
                                             sequence_abundances,
                                             subject_coefs,
                                             site_coefs,
                                             peptide_coef,
                                             group_ids,
                                             prop_to_miss){
  # browser()
  #assert length(sequence_abundances) = length(subject_coef)


  e_data_full = c()
  e_meta_full = c()
  ground_truth_site = c()
  # ground_truth_site_pattern = tibble(site_pattern = matrix(character(0), ncol = 1, nrow = 0))
  ground_truth_site_pattern = c()
  for(j in seq_along(sequence_abundances)){
    # Site modification is identical across subjects within a given group. While this may not exactly reflect reality, we can leave it for now and argue that the site modifcations within each group should be somewhat correlated, and this just assumes the correlation is 1. Potentially relax this later
    # I think its hard to argue that digestion should be the same for all subjects in a group. We can assume this is a completely random process, then modify the below objects based on the imperfect digestion
    modify_indices = get_modify_indices(ptm_site_mapping, sequence_abundances[j]) #sequence should depend on group
    ground_truth_site = cbind(ground_truth_site, colSums(modify_indices))

    for(i in seq_len(sum(group_ids == 1))){

      #modify ptm_site_mapping before passing into modify_and_count
      ptm_site_mapping_imperfect = imperfect_digest(ptm_site_mapping, prop_to_miss)

      e_data = get_peptide_counts(ptm_site_mapping_imperfect, modify_indices)
      e_meta = e_data %>% rename(true_peptide_count = peptide_count)

      #this maybe needs to be altered
      ground_truth_site_pattern = e_data  %>%
        group_by(peptide, site_pattern) %>%
        summarise(true_count = sum(peptide_count)) %>%
        ungroup() %>%
        select(-peptide) %>%
        mutate(Subject = i, Group = j)  %>%
        rbind(ground_truth_site_pattern, .)

      # ground_truth_site_pattern = e_data %>% group_by(peptide) %>% summarise(site_pattern = paste(site, collapse = ";"), true_count = sum(peptide_count)) %>% select(-peptide) %>% full_join(ground_truth_site_pattern, ., by = 'site_pattern')

      e_data_subj = e_data %>% mutate(peptide_count = peptide_count * exp(subject_coefs[i]))
      names(e_data_subj)[4] = paste0("Subject_", i, "_G", j)

      if(i == 1 & j == 1){
        e_data_full = e_data_subj %>% select(-site_pattern)
      } else {
        e_data_full = e_data_full %>%
          full_join(e_data_subj %>% select(-site_pattern), by = c("peptide", "site")) %>%
          mutate_if(is.numeric,coalesce, 0)
      }

      # if(i == 1) e_meta_full = rbind(e_meta_full, e_meta %>% mutate(group = j))
      e_meta_full = rbind(e_meta_full, e_meta %>% mutate(group = j))
    }
  }
  ground_truth_site_pattern = ground_truth_site_pattern %>%
    group_by(site_pattern, Group) %>%
    summarise(true_count = sum(true_count)) %>%
    pivot_wider(., values_from = true_count, names_from = Group)

  coef_dat = tibble(Mult_Factor = exp(site_coefs), site = names(site_coefs))

  e_data_full = e_data_full %>%
    left_join(., coef_dat) %>%
    select(-site) %>%
    group_by(peptide) %>%
    summarise(Mult_Factor = prod(Mult_Factor),
              across(starts_with("Subject"), ~ mean(.x))) %>%
    mutate(across(starts_with("Subject"), ~ . * Mult_Factor)) %>%
    select(-Mult_Factor)

  return(list(e_data = e_data_full, e_meta = e_meta_full, ground_truth_site = ground_truth_site, ground_truth_site_pattern = ground_truth_site_pattern))
}


get_samples_for_sequence = function(sequence,
                                    sequence_abundances,
                                    samps_per_group,
                                    site_sd,
                                    subject_sd,
                                    prop_to_miss,
                                    amino_acid = "S",
                                    identifier = "@"){
  # browser()
  digest_results = Digest(sequence)

  peptide_sd = 10
  unique_peptides = unique(digest_results$peptide)
  peptide_effects = rnorm(length(unique_peptides), 0, peptide_sd)
  peptide_effects = peptide_effects - (sum(peptide_effects) / length(unique_peptides)) #ensure they sum to 0
  names(peptide_effects) = unique_peptides

  ptm_site_mapping = get_ptm_site_mapping(sequence, digest_results, amino_acid)
  sites = unique(ptm_site_mapping$site) %>% stringi::stri_remove_empty()
  site_coefs = rnorm(length(sites), sd = site_sd)
  names(site_coefs) = sites
  subject_coefs = rnorm(samps_per_group * 2, sd = subject_sd)

  group_ids = rep(seq_len(2), each = samps_per_group)

  out = generate_ptm_samples_for_sequence(ptm_site_mapping,
                                          sequence_abundances,
                                          subject_coefs,
                                          site_coefs,
                                          peptide_coef = peptide_effects,
                                          group_ids = group_ids,
                                          prop_to_miss = prop_to_miss)

  #generate ground truth
  ground_truth_site = data.frame(site = rownames(out$ground_truth_site),
                                 G1 = out$ground_truth_site[,1],
                                 G2 = out$ground_truth_site[,2]) %>%
    mutate(true_fold_change = log2(G1 / G2)) %>%
    select(site, true_fold_change)

  ground_truth_site_pattern = out$ground_truth_site_pattern %>%
                              replace_na(., list(0, 0, 0)) %>%
                              rename(G1 = `1`, G2 = `2`) %>%
                              mutate(true_fold_change = log2(G1 / G2)) %>%
                              select(site_pattern, true_fold_change)

  f_data_full = data.frame(sample = out$e_data %>% select(-peptide) %>% names,
                           group = group_ids)

  pepdat = as.pepData(e_data = out$e_data,
                      f_data = f_data_full,
                      e_meta = out$e_meta %>% select(-group, -true_peptide_count, -site_pattern) %>% distinct(),
                      edata_cname = "peptide",
                      fdata_cname = "sample",
                      emeta_cname = "site")

  return(list(pepdat = pepdat, ground_truth_site = ground_truth_site, ground_truth_site_pattern = ground_truth_site_pattern))
}


get_samples_for_sequences = function(sequences, sequence_abundances, samps_per_group, site_sd, subject_sd, prop_to_miss, amino_acid = "S", identifier = "@"){
  e_data_full = c()
  f_data_full = c()
  e_meta_full = c()
  ground_truth_site_full = c()
  ground_truth_site_pattern_full = c()

  for(i in seq_along(sequences)){
    tmp_out = get_samples_for_sequence(sequence = sequences[i],
                                       sequence_abundances = sequence_abundances[i,],
                                       samps_per_group = samps_per_group,
                                       site_sd = site_sd,
                                       subject_sd = subject_sd,
                                       prop_to_miss = prop_to_miss,
                                       amino_acid = amino_acid,
                                       identifier = identifier)
    e_data_full = rbind(e_data_full, tmp_out$pepdat$e_data)
    e_meta_full = rbind(e_meta_full, tmp_out$pepdat$e_meta %>% mutate(site = paste(site, sequences[i], sep = ';')))
    ground_truth_site_full = rbind(ground_truth_site_full, tmp_out$ground_truth_site %>% mutate(protein = sequences[i]))
    ground_truth_site_pattern_full = rbind(ground_truth_site_pattern_full, tmp_out$ground_truth_site_pattern %>% mutate(protein = sequences[i]))
  }
  f_data_full = tmp_out$pepdat$f_data
  e_data_full = e_data_full %>% group_by(peptide) %>% summarise(across(everything(), ~sum(.x, na.rm = TRUE)))

  pepdat_full = as.pepData(e_data = e_data_full,
                           f_data = f_data_full,
                           e_meta = e_meta_full %>% distinct(),
                           edata_cname = "peptide",
                           fdata_cname = "sample",
                           emeta_cname = "site")
  return(list(pepdat = pepdat_full, ground_truth_site = ground_truth_site_full, ground_truth_site_pattern = ground_truth_site_pattern_full))
}

assign_plexes = function(omicsData, samples_per_plex, n_plex){
  #randomly assign plexes
  omicsData$f_data = omicsData$f_data %>% mutate(Plex = sample(rep(seq_len(n_plex), each = samples_per_plex)))
  return(omicsData)
}

# imperfect_digest = function(digest_data, prop_to_miss){
#   n_peps = nrow(digest_data)
#
#   # TODO does the probability of a sucessful split depend on the peptide length?
#   digest_data = digest_data %>% mutate(merge_with_next = (seq_len(n_peps) %in% sample(seq_len(n_peps)[-n_peps], round(prop_to_miss * (n_peps - 1)))))
#
#   remove_vec = rep(FALSE, nrow(digest_data))
#   for(i in seq_len(nrow(digest_data) - 1)){
#     if(digest_data$merge_with_next[i]){
#       peptide = paste0(digest_data$peptide[i], digest_data$peptide[i + 1])
#       start = digest_data$start[i]
#       stop = digest_data$stop[i + 1]
#       remove_vec[i] = TRUE
#       merge_wth_next = digest_data$merge_with_next[i+1]
#
#       digest_data[i + 1,] = c(peptide, start, stop, 0, 0, 0, 0, merge_wth_next)
#     }
#   }
#   return(digest_data[!remove_vec,])
# }





