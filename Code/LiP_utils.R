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
