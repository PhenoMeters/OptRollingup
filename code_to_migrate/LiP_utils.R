get_mask_indices = function(
    sequence,
    start,
    delay_fn,
    max_masked_region_size,
    min_num_trypsin_sites
){
  if(start > length(amino_acids)) stop("Start position is larger than sequence length")
  adjusted_start = start + delay_fn()
  curr_ind = adjusted_start
  amino_acids = unlist(lapply(sequence, \(x) strsplit(x, split = "")[[1]]))
  n_trypsin_sites = 0
  region_size = 0

  # Check if the adjusted start is valid
  if(adjusted_start > length(amino_acids)){
    return(list(start = start, stop = Inf)) #stop at inf to flag finished procedure (this one will be tossed out)
  }

  #check if one of the conditions can be satisfied
  remaining_amino_acids = amino_acids[adjusted_start:length(amino_acids)]
  remaining_spit_sites = sum(grepl("R", remaining_amino_acids) | grepl("K", remaining_amino_acids))
  if(length(remaining_amino_acids) < max_masked_region_size & remaining_spit_sites < min_num_trypsin_sites){
    return(list(start = adjusted_start, stop = Inf)) #stop at inf to flag finished procedure (this one will be tossed out)
  }


  while(n_trypsin_sites < min_num_trypsin_sites & region_size < max_masked_region_size){
    curr_amino_acid = amino_acids[curr_ind]
    if(curr_amino_acid == "K" | curr_amino_acid == "R"){
      n_trypsin_sites = n_trypsin_sites + 1
    }
    region_size = region_size + 1
    curr_ind = curr_ind + 1
  }
  end_ind = curr_ind + delay_fn()
  return(list(start = adjusted_start, stop = end_ind))
}

mask_sequence = function(
    sequence,
    delay_fn,
    gap_fn,
    max_masked_region_size,
    min_num_trypsin_sites
){
  start_vec = c()
  stop_vec = c()
  curr_start = 1 + gap_fn()
  last_stop = 0
  masked_sequence = sequence
  while(last_stop < nchar(sequence)){
    mask_indices = get_mask_indices(
      sequence =sequence,
      start = curr_start,
      delay_fn = delay_fn,
      max_masked_region_size = max_masked_region_size,
      min_num_trypsin_sites = min_num_trypsin_sites
    )
    last_stop = mask_indices$stop
    curr_start = mask_indices$stop + gap_fn()
    if(!(last_stop > nchar(sequence) | curr_start > nchar(sequence))){
      start_vec = c(start_vec, mask_indices$start)
      stop_vec = c(stop_vec, mask_indices$stop)
      masked_sequence = replace_between_indices(masked_sequence, mask_indices$start, mask_indices$stop, '#')
    }
  }
  return(list(sequence = sequence, masked_sequence = masked_sequence, start_vec = start_vec, stop_vec = stop_vec))
}

replace_between_indices = function(x, start, stop, replacement) {
  before = substr(x, 1, start - 1)
  after = substr(x, stop + 1, nchar(x))
  replacement_string = strrep(replacement, stop - start + 1)
  result = paste0(before, replacement_string, after)
  return(result)
}

mask_all_proteins = function(
    sequences,
    delay_fn,
    gap_fn,
    max_masked_region_size,
    min_num_trypsin_sites){
  out_list = vector('list', length = length(sequences))
  for(i in seq_along(sequences)){
    out_list[[i]] = mask_sequence(
      sequence = sequences[[i]],
      delay_fn = sample_delay,
      gap_fn = sample_gap,
      max_masked_region_size = max_masked_region_size,
      min_num_trypsin_sites = min_num_trypsin_sites
    )
  }
  return(out_list)
}

downsample_masking = function(
    masking_output,
    downsample_fn,
    ...
){
  indices = downsample_fn(seq_along(masking_output$start_vec), ...)

  masking_output$start_vec = masking_output$start_vec[indices]
  masking_output$stop_vec = masking_output$stop_vec[indices]
  masking_output$masked_sequence = remask_sequence(
    masking_output$sequence,
    masking_output$start_vec,
    masking_output$stop_vec
  )
  return(masking_output)
}

remask_sequence = function(sequence, start_vec, stop_vec){
  masked_sequence = sequence
  for(i in seq_along(start_vec)){
    masked_sequence = replace_between_indices(masked_sequence, start_vec[i], stop_vec[i], '#')
  }
  return(masked_sequence)
}

downsample_maskings = function(
    masking_list,
    downsample_fn,
    ...){
  out_list = vector('list', length = length(masking_list))
  for(i in seq_along(masking_list)){
    out_list[[i]] = downsample_masking(
      masking_output = masking_list[[i]],
      downsample_fn = downsample_fn,
      ...
    )
  }
  return(out_list)
}

simulate_proteinase_k_cleavage = function(sequence) {
  # Define amino acids considered for cleavage (aliphatic, aromatic, hydrophobic)
  cleavage_aa = c("A", "V", "L", "I", "F", "Y", "W", "M", "P")  # Aliphatic, aromatic, and hydrophobic residues

  # Convert the sequence to a character vector
  sequence_vector = unlist(strsplit(sequence, ""))

  # Initialize vectors to store peptide information
  peptides = c()
  start_indices = c()
  end_indices = c()

  # Initialize start index
  start = 1

  # Loop through the sequence to find potential cleavage sites
  for (i in 1:(length(sequence_vector) - 1)) {
    if (sequence_vector[i] %in% cleavage_aa) {
      peptides = c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices = c(start_indices, start)
      end_indices = c(end_indices, i)
      start = i + 1
    }
  }

  # Add the final peptide after the last cleavage site
  if (start <= length(sequence_vector)) {
    peptides = c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices = c(start_indices, start)
    end_indices = c(end_indices, length(sequence_vector))
  }

  # Create a dataframe
  cleavage_df = data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}

insert_char_at_indices = function(string, indices, char = '|') {
  indices = sort(indices)
  indices = indices[indices >= 1 & indices <= nchar(string)]

  chars = unlist(strsplit(string, ""))

  offset = 0
  for (index in indices) {
    chars = append(chars, char, after = index + offset - 1)
    offset = offset + 1
  }

  new_string = paste(chars, collapse = "")

  return(new_string)
}

add_cleavage_characters = function(masked_sequence_dat, digest_results, char = '|'){
  masked_sequence_dat$masked_sequence = insert_char_at_indices(masked_sequence_dat$masked_sequence, digest_results$start, char = char)
  masked_sequence_dat$sequence = insert_char_at_indices(masked_sequence_dat$sequence, digest_results$start, char = char)

  return(masked_sequence_dat)
}

simulate_trypsin_cleavage = function(sequence) {
  cleavage_aa = c("K", "R")
  skip_aa = "P"

  cleavage_rule = function(prev_aa, curr_aa) {
    return(curr_aa %in% cleavage_aa && prev_aa != skip_aa)
  }

  sequence_vector = unlist(strsplit(sequence, ""))

  peptides = c()
  start_indices = c()
  end_indices = c()

  start = 1

  for (i in 1:(length(sequence_vector) - 1)) {
    if (cleavage_rule(sequence_vector[i], sequence_vector[i + 1])) {
      peptides = c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices = c(start_indices, start)
      end_indices = c(end_indices, i)
      start = i + 1
    }
  }

  if (start <= length(sequence_vector)) {
    peptides = c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices = c(start_indices, start)
    end_indices = c(end_indices, length(sequence_vector))
  }

  cleavage_df = data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}


get_peptides_before_cleavage = function(sequence, cleavage_site_mapping){
  curr_start = 1
  peptides = c()
  for(i in seq_len(nrow(cleavage_site_mapping))){
    curr_stop = cleavage_site_mapping$Site[i] - 1
    peptides = c(peptides, substr(sequence, start = curr_start, stop = curr_stop))
    curr_start = curr_stop + 1
  }
  return(peptides)
}

get_peptides_after_cleavage = function(sequence, cleavage_site_mapping){
  peptides = c()
  for(i in seq_len(nrow(cleavage_site_mapping) - 1)){
    curr_start = cleavage_site_mapping$Site[i]
    curr_stop = cleavage_site_mapping$Site[i + 1] - 1
    peptides = c(peptides, substr(sequence, start = curr_start, stop = curr_stop))
  }
  curr_start = cleavage_site_mapping$Site[i + 1]
  curr_stop = nchar(sequence)
  peptides = c(peptides, substr(sequence, start = curr_start, stop = curr_stop))

  return(peptides)
}

get_peptide_abundances = function(sequences, cleavage_site_mapping, sequence_abundances = 1){
  prior_peptide_vec = c()
  out = vector("list", length(sequences))
  for(i in seq_along(sequences)){
    peptides_before_cleavage = get_peptides_before_cleavage(
      sequences[[i]],
      cleavage_site_mapping[[i]]
    )
    peptides_after_cleavage = get_peptides_after_cleavage(
      sequences[[i]],
      cleavage_site_mapping[[i]]
    )
    cleavage_site_mapping[[i]]$PriorPeptide = peptides_before_cleavage
    cleavage_site_mapping[[i]]$PosteriorPeptide = peptides_after_cleavage
    out[[i]] = data.frame(peptide = c(peptides_before_cleavage[1], peptides_after_cleavage), start_site = c(1, cleavage_site_mapping[[i]]$Site), abundance = sequence_abundances)
  }
  return(list(abundances = out, cleavage_site_mapping = cleavage_site_mapping))
}

digest_proteins_PK = function(sequences, pk_mu){
  lapply(sequences, \(x) simulate_proteinase_k_cleavage(x, pk_mu))
}

digest_proteins_trypsin = function(sequences, tr_mu){
  lapply(sequences, \(x) simulate_trypsin_cleavage(x, tr_mu))
}

get_masked_sequences = function(masking_list){
  lapply(masking_list, \(x) x$masked_sequence)
}
get_sequences = function(masking_list){
  lapply(masking_list, \(x) x$sequence)
}

simulate_protein_for_sample = function(n_reps, mask_out, pk_missed_prop, tr_missed_prop, sequences, pk_prob = NULL, tr_prob = NULL, pk_mu, tr_mu){
  peptide_abundances = vector("list", n_reps)
  cleavage_mappings = vector("list", n_reps)
  for(i in seq_len(n_reps)){
    cleavage_mapping_example = get_cleavage_site_mapping(
      mask_out,
      pk_missed_prop,
      tr_missed_prop,
      pk_prob,
      tr_prob,
      pk_mu,
      tr_mu)
    tmp = get_peptide_abundances(
      sequences,
      cleavage_mapping_example)
    tmp$abundances = lapply(seq_along(tmp$abundances), function(x) {
      tmp$abundances[[x]]$protein <- mask_out[[x]]$sequence
      tmp$abundances[[x]]  # Return the modified data frame
    })

    tmp$cleavage_site_mapping = lapply(seq_along(tmp$abundances), function(x) {
      tmp$cleavage_site_mapping[[x]]$protein <- mask_out[[x]]$sequence
      tmp$cleavage_site_mapping[[x]]  # Return the modified data frame
    })

    peptide_abundances[[i]] = tmp$abundances %>% purrr::reduce(., rbind)
    cleavage_mappings[[i]] = tmp$cleavage_site_mapping %>% purrr::reduce(., rbind)

  }
  full_peptide_df = peptide_abundances %>%
    purrr::reduce(., rbind) %>%
    group_by(start_site, peptide, protein) %>%
    summarise(abundance = sum(abundance)) %>%
    ungroup()
  full_cleavage_mapping = cleavage_mappings %>%
    purrr::reduce(., rbind) %>%
    distinct() #TODO Maybe we dont want this, I forget though
  # arrange(start_site)
  return(list(peptide_df = full_peptide_df, cleavage_mapping = full_cleavage_mapping))
}


get_cleavage_site_mapping = function(mask_out, pk_missed_prop, tr_missed_prop, pk_prob = NULL, tr_prob = NULL, pk_mu = 1, tr_mu = 1){
  # digest_results_PK = adjust_digestions(digest_proteins_PK(get_masked_sequences(mask_out)), prop_to_miss = pk_missed_prop, prob = pk_prob)
  # digest_results_tryp = adjust_digestions(digest_proteins_trypsin(get_sequences(mask_out)), prop_to_miss = tr_missed_prop, prob = tr_prob)
  digest_results_PK = digest_proteins_PK(get_masked_sequences(mask_out), pk_mu)
  digest_results_tryp = digest_proteins_trypsin(get_sequences(mask_out), tr_mu)

  proteinaseK_cleavage_sites = lapply(digest_results_PK, \(x) x$start[-1])
  trypsin_cleavage_sites = lapply(digest_results_tryp, \(x) x$start[-1])

  PK_dat = lapply(proteinaseK_cleavage_sites, \(x) data.frame(Site = x, PK = 1))
  TR_dat = lapply(trypsin_cleavage_sites, \(x) data.frame(Site = x, TR = -1))
  full_dat = lapply(
    seq_along(PK_dat),
    \(x) full_join(PK_dat[[x]], TR_dat[[x]], by = "Site") %>%
      tidyr::replace_na(., list(PK = 0, TR = 0)) %>%
      mutate(
        sum_col = PK + TR,
        Enzyme = case_match(
          sum_col,
          -1 ~ "Trypsin",
          0 ~ "Proteinase K", # 0 is proteinase K because a trypsin cleavage cannot occur if a PK one already did
          1 ~ "Proteinase K"
        )) %>%
      select(Site, Enzyme) %>%
      arrange(Site)
  )

  return(full_dat)
}


#TODO: set this argument in the rest of the code
simulate_proteinase_k_cleavage = function(sequence, pk_mu = 8) {
  # Define amino acids considered for cleavage (aliphatic, aromatic, hydrophobic)
  cleavage_aa = c("A", "V", "L", "I", "F", "Y", "W", "M", "P")  # Aliphatic, aromatic, and hydrophobic residues

  # Convert the sequence to a character vector
  sequence_vector = unlist(strsplit(sequence, ""))

  # Initialize vectors to store peptide information
  peptides = c()
  start_indices = c()
  end_indices = c()

  # Initialize start index
  start = 1

  # Loop through the sequence to find potential cleavage sites
  buffer = 0
  min_size = rpois(1, pk_mu)
  for (i in 1:(length(sequence_vector) - 1)) {
    if (sequence_vector[i] %in% cleavage_aa & buffer > min_size) {
      peptides = c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices = c(start_indices, start)
      end_indices = c(end_indices, i)
      buffer = 0
      start = i + 1
      min_size = rpois(1, pk_mu)
    }
    buffer = buffer + 1
  }

  # Add the final peptide after the last cleavage site
  if (start <= length(sequence_vector)) {
    peptides = c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices = c(start_indices, start)
    end_indices = c(end_indices, length(sequence_vector))
  }

  # Create a dataframe
  cleavage_df = data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}


simulate_trypsin_cleavage = function(sequence, tr_mu = 1) {
  cleavage_aa = c("K", "R")
  skip_aa = "P"

  cleavage_rule = function(prev_aa, curr_aa) {
    return(curr_aa %in% cleavage_aa && prev_aa != skip_aa)
  }

  sequence_vector = unlist(strsplit(sequence, ""))

  peptides = c()
  start_indices = c()
  end_indices = c()

  start = 1
  buffer = 0
  min_size = rpois(1, tr_mu)
  for (i in 1:(length(sequence_vector) - 1)) {
    if (cleavage_rule(sequence_vector[i], sequence_vector[i + 1]) & buffer > min_size) {
      peptides = c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices = c(start_indices, start)
      end_indices = c(end_indices, i)
      buffer = 0
      start = i + 1
      min_size = rpois(1, tr_mu)
    }
    buffer = buffer + 1
  }

  if (start <= length(sequence_vector)) {
    peptides = c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices = c(start_indices, start)
    end_indices = c(end_indices, length(sequence_vector))
  }

  cleavage_df = data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}

#lets just assume this is for one protein
simulate_proteins_for_group = function(
    n_subjects,
    group_id,
    n_reps,
    masking_set,
    masking_proportions,
    full_masking,
    pk_missed_prop,
    tr_missed_prop,
    sequences,
    pk_prob = NULL,
    tr_prob = NULL,
    pk_mu = 1,
    tr_mu = 1) {
  full_pep_list = vector("list", length(masking_set))
  full_cleavage_mapping = vector("list", length(masking_set) * n_subjects)
  #loop over masking patterns
  #create mask_out from current masking pattern
  #adjust n_reps based on target proportion
  #continue as normal
  sequences = list(sequences)
  for(m in seq_along(masking_set)){
    sub_pep_list = vector("list", n_subjects)

    curr_masking_pattern = list(get_masking_pattern(full_masking, masking_pattern_set[[m]]))
    adjusted_n_reps = round(n_reps * masking_proportions[m])

    for(i in seq_len(n_subjects)){
      out = simulate_protein_for_sample(adjusted_n_reps, curr_masking_pattern, pk_missed_prop, tr_missed_prop, sequences, pk_prob, tr_prob, pk_mu, tr_mu)
      sub_pep_list[[i]] = out$peptide_df
      full_cleavage_mapping[[i + (m - 1) * n_subjects]] = out$cleavage_mapping
    }
    full_pep_list[[m]] = purrr::reduce(sub_pep_list, full_join, by = c("start_site", "peptide", "protein"))
  }
  #join data across masking patterns
  full_pep_list = purrr::reduce(full_pep_list, rbind) %>% group_by(start_site, peptide, protein) %>% summarise(across(everything(), sum, na.rm = T))
  full_cleavage_mapping = purrr::reduce(full_cleavage_mapping, rbind) %>% distinct()
  names(full_pep_list)[-c(1:3)] = paste("subject", seq_len(n_subjects), "group", group_id, sep = '_')

  return(list(peptide_df = full_pep_list, cleavage_mapping = full_cleavage_mapping))
}

get_masking_pattern = function(
    full_masking,
    init_masking_pattern
){
  full_masking$start_vec = init_masking_pattern$start_vec
  full_masking$stop_vec = init_masking_pattern$stop_vec
  full_masking$masked_sequence = remask_sequence(
    full_masking$sequence,
    full_masking$start_vec,
    full_masking$stop_vec
  )
  return(full_masking)
}

create_masking_set = function(full_masking, n_patterns) {
  masking_set = full_masking
  n_hash_regions <- count_hash_regions(full_masking$masked_sequence)

  masking_set = lapply(
    seq_len(n_patterns),
    \(x) {
      # mask_region_indices = sample(seq_len(n_hash_regions), sample(seq_len(n_hash_regions), 1))
      mask_region_indices = sample(seq_len(n_hash_regions), 2)
      list(start_vec = full_masking$start_vec[mask_region_indices], stop_vec = full_masking$stop_vec[mask_region_indices])
    }
  )
  return(masking_set)
}

count_hash_regions <- function(input_string) {
  modified_string <- gsub("#+", "#", input_string)

  hash_regions <- strsplit(modified_string, "[^#]")[[1]]
  hash_regions <- hash_regions[hash_regions != ""]

  return(length(hash_regions))
}

generate_f_data = function(full_abundance_df){
  subject_ids = names(full_abundance_df)[-which(names(full_abundance_df) %in% c('start_site', 'peptide', 'protein'))]
  f_data = data.frame(Subject_IDs = subject_ids) %>%
    mutate(Group_IDs = sapply(strsplit(Subject_IDs, '_'), function(x) x[4]))
  return(f_data)
}

### Kind of messy plotting and results code





get_full_pep_data = function(g1_dat, g2_dat){
  #wrap up into a get_site_results function
  full_pep_data = full_join(g1_dat$peptide_df, g2_dat$peptide_df, by = c("peptide", "start_site", "protein"))

  f_data = generate_f_data(full_pep_data)

  pep_data = pmartR::as.pepData(
    e_data = full_pep_data %>% ungroup() %>% mutate(Peptide_Site = paste(start_site, peptide, protein, sep = '_')) %>% select(-start_site, -peptide, -protein),
    f_data = f_data,
    edata_cname = "Peptide_Site",
    fdata_cname = "Subject_IDs")

  pep_data = pep_data %>%
    pmartR::edata_transform(., 'log2') %>%
    pmartR::normalize_global(., "all", "median", apply_norm = T)

  full_pep_data = pep_data$e_data %>%
    mutate(
      peptide = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][2]),
      start_site = as.numeric(sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][1])),
      protein = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][3])
    ) %>%
    select(-Peptide_Site)

  return(full_pep_data)
}

get_full_cleavage_mapping = function(g1_dat, g2_dat){
  full_cleavage_mapping = rbind(
    g1_dat$cleavage_mapping %>% distinct,
    g2_dat$cleavage_mapping %>% distinct
  ) %>%
    distinct() %>%
    as_tibble()
  return(full_cleavage_mapping)
}


get_site_res = function(full_pep_data, full_cleavage_mapping, aggregation_fn = median){
  long_cleavage_mapping = full_cleavage_mapping %>%  pivot_longer(c("PriorPeptide", "PosteriorPeptide"), values_to = "Peptide", names_to = "Type")

  data_for_rollup = full_pep_data %>%
    rename(Site = start_site, Peptide = peptide) %>%
    left_join(., long_cleavage_mapping, by = c("Site", "Peptide")) %>%
    drop_na(., Enzyme, Type)

  rolled_up_data = data_for_rollup %>% group_by(Site, Enzyme) %>% select(-Type) %>% summarize_if(is.numeric, aggregation_fn, na.rm = T)
  # need to figure out the best way to put all this data together

  rolled_up_data = data_for_rollup %>% group_by(Site, Enzyme) %>% select(-Type) %>% summarize_if(is.numeric, median, na.rm = T)

  site_data = pmartR::as.pepData(
    e_data = rolled_up_data %>% ungroup() %>% mutate(Site_Type = paste(Site, Enzyme, sep = '_')) %>% select(-Site, -Enzyme),
    f_data = f_data,
    edata_cname = "Site_Type",
    fdata_cname = "Subject_IDs",
    data_scale = "log2")

  site_data = site_data %>% group_designation(., main_effects = "Group_IDs")
  site_res = site_data %>% imd_anova(., test_method = "anova")

  return(site_res)
}

create_site_plot = function(site_res, masking_pattern_set){
  #wrap up into a get_site_plots function?
  ggdat = site_res %>%
    select(
      Site_Type,
      Fold_change_1_vs_2
    ) %>%
    rename(FC = Fold_change_1_vs_2) %>%
    mutate(
      Site = as.numeric(sapply(Site_Type, \(x) strsplit(x, split = '_')[[1]][1])),
      Enzyme = sapply(Site_Type, \(x) strsplit(x, split = '_')[[1]][2])
    ) %>%
    select(-Site_Type)

  masking_data = c()
  for(idx in seq_along(masking_pattern_set)){
    curr_masking_pattern = list(get_masking_pattern(full_masking, masking_pattern_set[[idx]]))
    tmp_masking_data = data.frame(
      Site = seq_len(nchar(curr_masking_pattern[[1]]$sequence)),
      AA = strsplit(curr_masking_pattern[[1]]$masked_sequence, '')[[1]],
      pattern_id = idx
    ) %>%
      mutate(Site_Type = ifelse(AA == "#", "Masked", "Unmasked"))
    masking_data = rbind(masking_data, tmp_masking_data)
  }

  all_sites <- data.frame(
    Site = rep(seq(min(masking_data$Site), max(masking_data$Site))), length(masking_pattern_set),
    pattern_id = rep(seq_along(masking_pattern_set), each = max(masking_data$Site)))

  plot_data_1 <- all_sites %>%
    left_join(masking_data, by = c("Site", "pattern_id")) %>%
    mutate(Enzyme = "Trypsin")
  plot_data_2 <- all_sites %>%
    left_join(masking_data, by = c("Site", "pattern_id")) %>%
    mutate(Enzyme = "Proteinase K")
  plot_data = rbind(plot_data_1, plot_data_2) %>%
    left_join(., ggdat, by = c("Site", "Enzyme"))

  plot_data = plot_data %>% filter(Enzyme != "Trypsin")
  site_line_data <- plot_data %>%
    filter(!is.na(FC))

  p = ggplot() +
    geom_tile(data = plot_data, aes(x = Site, y = (max(site_line_data$FC, na.rm = TRUE) + min(site_line_data$FC, na.rm = TRUE)) / 2, fill = Site_Type), height = max(site_line_data$FC, na.rm = TRUE) - min(site_line_data$FC, na.rm = TRUE), alpha = 0.2) +
    geom_line(data = site_line_data, aes(x = Site, y = FC)) +
    theme_bw() +
    scale_fill_manual(values = c("Unmasked" = "lightgray", "Masked" = "#FDE725FF")) +
    labs(fill = "Region Masked?") +
    ylab("log2 Fold Change") +
    facet_wrap(~ pattern_id, scales = "free", ncol = 2)
  return(list(p = p, plot_data = plot_data, line_data = site_line_data))
}




get_tryptic_res = function(full_pep_data, full_cleavage_mapping){
  d1 = full_cleavage_mapping %>%
    filter(Enzyme == "Trypsin") %>%
    arrange(Site) %>%
    select(-PriorPeptide) %>%
    rename(Start_Location = Site,
           Start_Location_Enzyme = Enzyme,
           Peptide = PosteriorPeptide) %>%
    mutate(End_Location = nchar(Peptide) + Start_Location) %>%
    distinct()


  d2 = full_cleavage_mapping %>%
    filter(Enzyme == "Trypsin") %>%
    arrange(Site) %>%
    select(-PosteriorPeptide) %>%
    rename(End_Location = Site,
           Start_Location_Enzyme = Enzyme,
           Peptide = PriorPeptide) %>%
    distinct()

  tryptic_peptides = inner_join(d1, d2, by = c("End_Location", "Peptide", "Start_Location_Enzyme", "protein"))

  tryptic_pep_edata = full_pep_data %>% filter(start_site %in% tryptic_peptides$Start_Location & peptide %in% tryptic_peptides$Peptide)

  tryptic_pep_data = pmartR::as.pepData(
    e_data = tryptic_pep_edata %>% ungroup() %>% mutate(Peptide_Site = paste(start_site, peptide, protein, sep = '_')) %>% select(-start_site, -peptide, -protein),
    f_data = f_data,
    edata_cname = "Peptide_Site",
    fdata_cname = "Subject_IDs",
    data_scale = "log2")


  tryptic_pep_data = tryptic_pep_data %>% group_designation(., main_effects = "Group_IDs")
  tryptic_res = tryptic_pep_data %>% imd_anova(., test_method = "anova")
  return(tryptic_res)
}

#tryptic peptide analysis



create_tryptic_plot = function(tryptic_res, masking_pattern_set){
  ggdat = tryptic_res %>%
    select(Peptide_Site,
           Fold_change_1_vs_2) %>%
    rename(FC = Fold_change_1_vs_2) %>%
    mutate(
      Peptide = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][2]),
      Site = as.numeric(sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][1])),
      # protein = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][3]),
      Enzyme = "Trypsin"
    ) %>%
    select(-Peptide_Site)

  masking_data = c()
  for(idx in seq_along(masking_pattern_set)){
    curr_masking_pattern = list(get_masking_pattern(full_masking, masking_pattern_set[[idx]]))
    tmp_masking_data = data.frame(
      Site = seq_len(nchar(curr_masking_pattern[[1]]$sequence)),
      AA = strsplit(curr_masking_pattern[[1]]$masked_sequence, '')[[1]],
      pattern_id = idx
    ) %>%
      mutate(Site_Type = ifelse(AA == "#", "Masked", "Unmasked"))
    masking_data = rbind(masking_data, tmp_masking_data)
  }

  all_sites <- data.frame(
    Site = rep(seq(min(masking_data$Site), max(masking_data$Site))), length(masking_pattern_set),
    pattern_id = rep(seq_along(masking_pattern_set), each = max(masking_data$Site)))

  plot_data_1 <- all_sites %>%
    left_join(masking_data, by = c("Site", "pattern_id")) %>%
    mutate(Enzyme = "Trypsin")
  plot_data_2 <- all_sites %>%
    left_join(masking_data, by = c("Site", "pattern_id")) %>%
    mutate(Enzyme = "Proteinase K")
  plot_data = rbind(plot_data_1, plot_data_2) %>%
    left_join(., ggdat, by = c("Site", "Enzyme"))

  tryptic_line_data <- plot_data %>%
    filter(!is.na(FC))

  plot_data = plot_data %>% filter(Enzyme == "Trypsin")


  # Creating the plot with faceting
  p = ggplot() +
    geom_tile(data = plot_data, aes(x = Site, y = (max(tryptic_line_data$FC, na.rm = TRUE) + min(tryptic_line_data$FC, na.rm = TRUE)) / 2, fill = Site_Type), height = max(tryptic_line_data$FC, na.rm = TRUE) - min(tryptic_line_data$FC, na.rm = TRUE), alpha = 0.2) +
    geom_line(data = tryptic_line_data, aes(x = Site, y = FC)) +
    theme_bw() +
    scale_fill_manual(values = c("Unmasked" = "lightgray", "Masked" = "#FDE725FF")) +
    labs(fill = "Region Masked?") +
    ylab("log2 Fold Change") +
    facet_wrap(~ pattern_id, scales = "free", ncol = 2)
  return(list(p = p, plot_data = plot_data, line_data = tryptic_line_data))
}


get_comparison_results = function(site_line_data, tryptic_line_data, plot_data, full_masking){
  region_id_data = data.frame(Site = seq_len(nchar(full_masking$sequence))) %>%
    left_join(
      data.frame(
        Site = c(
          full_masking$start_vec,
          full_masking$stop_vec + 1
        ),
        region_id = c(
          seq_along(full_masking$start_vec),
          rep(0, length(full_masking$stop_vec)))
      ),
      by = 'Site'
    ) %>%
    fill(
      region_id,
      .direction = 'down'
    ) %>%
    mutate(region_id = replace_na(region_id, 0))


  fc_data = data.frame(
    region_id = c(
      region_id_data %>% filter(Site %in% masking_pattern_set[[1]]$start_vec) %>% pull(region_id),
      region_id_data %>% filter(Site %in% masking_pattern_set[[2]]$start_vec) %>% pull(region_id)
    ),
    log2fc = rep(log2(expected_fc), each = length(masking_pattern_set[[1]]$start_vec))
  ) %>%
    group_by(region_id) %>%
    summarise(log2fc = mean(log2fc)) %>%
    ungroup()


  tryptic_region_dat = tryptic_line_data %>% left_join(
    region_id_data,
    by = "Site"
  ) %>% left_join(
    fc_data,
    by = "region_id"
  ) %>%
    mutate(
      log2fc = replace_na(log2fc, 0),
      masked_region_id = ifelse(log2fc != 0, 1, 0) * region_id,
      squared_error = (FC - log2fc)^2) %>%
    group_by(masked_region_id)

  tryptic_region_metrics = tryptic_region_dat %>%
    summarise(
      MSE = mean(squared_error),
      expected_FC = mean(log2fc),
      estimated_FC = mean(FC))

  site_region_dat = site_line_data %>% left_join(
    region_id_data,
    by = "Site"
  ) %>% left_join(
    fc_data,
    by = "region_id"
  ) %>%
    mutate(
      log2fc = replace_na(log2fc, 0) * -1,
      masked_region_id = ifelse(log2fc != 0, 1, 0) * region_id,
      squared_error = (FC - log2fc)^2) %>%
    group_by(masked_region_id)

  site_region_metrics = site_region_dat %>%
    summarise(
      MSE = mean(squared_error),
      expected_FC = mean(log2fc),
      estimated_FC = mean(FC))

  site_plot_data = region_id_data %>%
    left_join(
      site_region_dat %>% select(region_id, masked_region_id, log2fc) %>% distinct,
      by = "region_id"
    ) %>% left_join(
      site_region_metrics,
      by = "masked_region_id"
    )

  tryptic_plot_data = region_id_data %>%
    left_join(
      tryptic_region_dat %>% select(region_id, masked_region_id, log2fc) %>% distinct,
      by = "region_id"
    ) %>% left_join(
      tryptic_region_metrics,
      by = "masked_region_id"
    )

  ggdat = rbind(
    site_plot_data %>% mutate(Analysis_Method = "SiteRollup"),
    tryptic_plot_data %>% mutate(Analysis_Method = "TrypticPeptides")) %>%
    ungroup() %>%
    select(Site, MSE, log2fc, estimated_FC, Analysis_Method) %>%
    mutate(
      min = pmin(MSE, log2fc, estimated_FC),
      max = pmax(MSE, log2fc, estimated_FC)) %>%
    pivot_longer(c('log2fc', 'MSE', 'estimated_FC'), values_to = "Value", names_to = "Metric") %>%
    filter(Metric != "MSE")

  tile_data = plot_data %>%
    select(Site, Site_Type) %>%
    distinct() %>%
    mutate(tmp = ifelse(Site_Type == "Masked", 0, 1)) %>%
    select(-Site_Type) %>%
    group_by(Site) %>%
    summarize(tmp = prod(tmp)) %>%
    mutate(Site_Type = ifelse(tmp == 0, "Masked", "Unmasked")) %>%
    select(-tmp)

  compare_p = ggplot() +
    geom_tile(
      data = tile_data,
      aes(
        x = Site,
        y = (max(ggdat$max, na.rm = TRUE) + min(ggdat$min, na.rm = TRUE)) / 2,
        fill = Site_Type),
      height = max(ggdat$max, na.rm = TRUE) - min(ggdat$min, na.rm = TRUE),
      alpha = 0.4) +
    geom_line(data = ggdat, aes(x = Site, y = Value, color = Metric, linetype = Metric), linewidth = 1) +
    facet_wrap(~ Analysis_Method, scales = "free") +
    theme_bw() +
    scale_fill_manual(values = c("Unmasked" = "lightgray", "Masked" = "#FDE725FF")) +
    scale_color_manual(values = c("estimated_FC" = "#22A884FF", "log2fc" = "black")) +
    ggplot2::scale_linetype_manual(values = c("estimated_FC" = "solid", "log2fc" = "dotted")) +
    labs(fill = "Region Masked?") +
    ylab("Within Region log2 Fold Change") +
    theme(legend.position = 'bottom')

  #comparison across runs
  res = rbind(
    tryptic_region_metrics %>% mutate(masked = ifelse(masked_region_id == 0, F, T)) %>% group_by(masked) %>% summarise(MSE = mean(MSE)) %>% mutate(Analysis_Method = "TrypticPeptides"),
    site_region_metrics %>% mutate(masked = ifelse(masked_region_id == 0, F, T)) %>% group_by(masked) %>% summarise(MSE = mean(MSE)) %>% mutate(Analysis_Method = "SiteRollup")
  )

  return(list(p = compare_p, res = res))
}

