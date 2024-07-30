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

get_modify_indices = function(e_meta, sequence_abundance){
  unique_sites = unique(e_meta %>% drop_na() %>% pull(site))
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
  # browser()
  peps = unique(e_meta$peptide)
  out_dat = c()
  for(i in seq_along(peps)){
    pep_e_meta = e_meta %>% filter(peptide == peps[i]) %>% distinct()
    tmp_site_pattern = pep_e_meta %>% pull(site)
    tmp_sites = pep_e_meta %>% pull(site) %>% strsplit(., split = ';') %>% unlist() %>% stringi::stri_remove_empty()
    
    
    site_modifications = matrix(modify_indices[,colnames(modify_indices) %in% tmp_sites], ncol = length(tmp_sites))
    # site_modifications = matrix(modify_indices[,sapply(colnames(modify_indices), \(x) grepl(x, tmp_site))], ncol = length(tmp_site))
    
    colnames(site_modifications) = tmp_sites
    site_modification_counts = table(apply(site_modifications, 1, \(x) paste0(colnames(site_modifications)[x], collapse = ";")))
    #drop no modification case here
    site_modification_counts = site_modification_counts[names(site_modification_counts) != '']
    
    for(j in seq_along(site_modification_counts)){
      curr_ptm = site_modification_counts[j]
      sites = strsplit(names(curr_ptm), split = ';')[[1]]
      ptm = ptm_at_sites(peps[i], sites, "#", pep_e_meta %>% filter(site %in% tmp_site_pattern) %>% pull(start) %>% unique %>% as.numeric())
      out_dat = rbind(out_dat, tibble(peptide = ptm, site = sites, peptide_count = curr_ptm))
    }
  }
  return(out_dat)
}

generate_ptm_samples_for_sequence = function(ptm_site_mapping,
                                             sequence_abundances, 
                                             subject_coefs, 
                                             site_coefs,
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
    #Site modification is identical across subjects within a given group. While this may not exactly reflect reality, we can leave it for now and argue that the site modifcations within each group should be somewhat correlated, and this just assumes the correlation is 1. Potentially relax this later
    # I think its hard to argue that digestion should be the same for all subjects in a group. We can assume this is a completely random process, then modify the below objects based on the imperfect digestion
    modify_indices = get_modify_indices(ptm_site_mapping, sequence_abundances[j]) #sequence should depend on group
    ground_truth_site = cbind(ground_truth_site, colSums(modify_indices))
    
    for(i in seq_len(sum(group_ids == 1))){
      # TODO: add imperfect digestion here
      
      #modify ptm_site_mapping before passing into modify_and_count
      ptm_site_mapping_imperfect = imperfect_digest(ptm_site_mapping, prop_to_miss)
      
      e_data = modify_and_count(ptm_site_mapping_imperfect, modify_indices)
      e_meta = e_data %>% rename(true_peptide_count = peptide_count)
      
      #this maybe needs to be altered
      ground_truth_site_pattern = e_data %>% 
        group_by(peptide) %>% 
        summarise(site_pattern = paste(site, collapse = ";"), true_count = sum(peptide_count)) %>% 
        select(-peptide) %>% 
        mutate(Subject = i, Group = j) %>%
        rbind(ground_truth_site_pattern, .)
      
      # ground_truth_site_pattern = e_data %>% group_by(peptide) %>% summarise(site_pattern = paste(site, collapse = ";"), true_count = sum(peptide_count)) %>% select(-peptide) %>% full_join(ground_truth_site_pattern, ., by = 'site_pattern')

      e_data_subj = e_data %>% mutate(peptide_count = peptide_count * exp(subject_coefs[i]))
      names(e_data_subj)[3] = paste0("Subject_", i, "_G", j)
      
      if(i == 1 & j == 1){
        e_data_full = e_data_subj
      } else {
        e_data_full = e_data_full %>% 
          full_join(e_data_subj, by = c("peptide", "site")) %>%
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
  ptm_site_mapping = get_ptm_site_mapping(sequence, amino_acid)
  sites = unique(ptm_site_mapping$site) %>% stringi::stri_remove_empty()
  site_coefs = rnorm(length(sites), sd = site_sd)
  names(site_coefs) = sites
  subject_coefs = rnorm(samps_per_group * 2, sd = subject_sd)
  
  group_ids = rep(seq_len(2), each = samps_per_group)
  
  out = generate_ptm_samples_for_sequence(ptm_site_mapping,
                                          sequence_abundances,
                                          subject_coefs,
                                          site_coefs,
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
                      e_meta = out$e_meta %>% select(-group, -true_peptide_count) %>% distinct(),
                      edata_cname = "peptide",
                      fdata_cname = "sample",
                      emeta_cname = "site")
  
  return(list(pepdat = pepdat, ground_truth_site = ground_truth_site, ground_truth_site_pattern = ground_truth_site_pattern))
}


get_samples_for_sequences = function(sequences, sequence_abundances, samps_per_group, site_sd, subject_sd, prop_to_miss, amino_acid = "S", identifier = "@"){
  # browser()
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

imperfect_digest = function(ptm_site_mapping, prop_to_miss){
  n_peps = nrow(ptm_site_mapping)
  
  # TODO does the probability of a sucessful split depend on the peptide length?
  ptm_site_mapping = ptm_site_mapping %>% mutate(merge_with_next = (seq_len(n_peps) %in% sample(seq_len(n_peps)[-n_peps], round(prop_to_miss * (n_peps - 1)))))
  protein = ptm_site_mapping$protein %>% unique()
  
  remove_vec = rep(FALSE, nrow(ptm_site_mapping))
  for(i in seq_len(nrow(ptm_site_mapping) - 1)){
    if(ptm_site_mapping$merge_with_next[i]){
      peptide = paste0(ptm_site_mapping$peptide[i], ptm_site_mapping$peptide[i + 1])
      start = ptm_site_mapping$start[i]
      stop = ptm_site_mapping$stop[i + 1]
      remove_vec[i] = TRUE
      merge_wth_next = ptm_site_mapping$merge_with_next[i+1]
      site =  paste(ptm_site_mapping$site[i], ptm_site_mapping$site[i + 1], sep = ';')
      
      ptm_site_mapping[i + 1,] = c(site, peptide, protein, start, stop, merge_wth_next)
    }
  }
  return(ptm_site_mapping[!remove_vec,])
}



