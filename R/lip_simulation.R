#' Title
#' #'
#'Description
#'
#'@param $$
#'
#'@export
#'
#'@example
#'
simulate_proteins_for_group <- function(
    n_subjects,
    group_id,
    n_reps,
    mask_out,
    pk_missed_prop,
    tr_missed_prop,
    sequences,
    pk_prob = NULL,
    tr_prob = NULL
    ) {

  #Initialize output objects
  full_pep_list <- vector("list", n_subjects)
  full_cleavage_mapping <- vector("list", n_subjects)

  #Loop through all subjects in group
  for (i in seq_len(n_subjects)) {
    out <- simulate_protein_for_sample(
        n_reps,
        mask_out,
        pk_missed_prop,
        tr_missed_prop,
        sequences,
        pk_prob,
        tr_prob
        )
    full_pep_list[[i]] <- out$peptide_df
    full_cleavage_mapping[[i]] <- out$cleavage_mapping
  }

  #Combine peptide and cleavage data across subjects
  full_pep_list <- purrr::reduce(full_pep_list, dplyr::full_join, by = c("start_site", "peptide", "protein"))
  full_cleavage_mapping <- purrr::reduce(full_cleavage_mapping, rbind)
  names(full_pep_list)[-c(1:3)] <- paste("subject", seq_len(n_subjects), "group", group_id, sep = "_")

  return(list(peptide_df = full_pep_list, cleavage_mapping = full_cleavage_mapping))
}

#TODO: This is limited to the same number of replicates for each protein in sequences, might want to relax this
simulate_protein_for_sample <- function(n_reps, mask_out, pk_missed_prop, tr_missed_prop, sequences, pk_prob = NULL, tr_prob = NULL) {
  #initialize output objects
  peptide_abundances <- vector("list", n_reps)
  cleavage_mappings <- vector("list", n_reps)

  #loop through and simulate each instance of protein(s) separately
  for (i in seq_len(n_reps)) {

    #get peptide abundance and cleavage objects
    cleavage_mapping_for_replicate <- get_cleavage_site_mapping(mask_out, pk_missed_prop, tr_missed_prop, pk_prob, tr_prob)
    peptide_abundances_for_replicate <- get_peptide_abundances(sequences, cleavage_mapping_for_replicate) #TODO: does this need to return the cleavage mapping too?

    #add protein label to abundance and cleavage objects
    peptide_abundances_for_replicate$abundances <- lapply(
        seq_along(peptide_abundances_for_replicate$abundances),
        function(x) {
            peptide_abundances_for_replicate$abundances[[x]]$protein <- mask_out[[x]]$sequence
            peptide_abundances_for_replicate$abundances[[x]]
            }
        )
    peptide_abundances_for_replicate$cleavage_site_mapping <- lapply(
        seq_along(peptide_abundances_for_replicate$abundances),
        function(x) {
            peptide_abundances_for_replicate$cleavage_site_mapping[[x]]$protein <- mask_out[[x]]$sequence
            peptide_abundances_for_replicate$cleavage_site_mapping[[x]]
            }
        )

    #join abundance and cleavage objects across all proteins
    peptide_abundances[[i]] <- peptide_abundances_for_replicate$abundances %>% purrr::reduce(., rbind)
    cleavage_mappings[[i]] <- peptide_abundances_for_replicate$cleavage_site_mapping %>% purrr::reduce(., rbind)
  }

  #format abundance and cleavage objects before returning
  full_peptide_df <- peptide_abundances %>%
    purrr::reduce(., rbind) %>%
    dplyr::group_by(start_site, peptide, protein) %>%
    dplyr::summarise(abundance = sum(abundance)) %>%
    dplyr::ungroup()
  full_cleavage_mapping <- cleavage_mappings %>%
    purrr::reduce(., rbind) %>%
    dplyr::distinct() #TODO Maybe we dont want this, I forget though

  return(list(peptide_df = full_peptide_df, cleavage_mapping = full_cleavage_mapping))
}

get_peptide_abundances <- function(sequences, cleavage_site_mapping, sequence_abundances = 1) {
  #initialize output objects
  abundances <- vector("list", length(sequences))

  #loop through proteins
  for (i in seq_along(sequences)) {
    #get each peptide before and after cleavage sites
    peptides_before_cleavage <- get_peptides_before_cleavage(sequences[[i]], cleavage_site_mapping[[i]])
    peptides_after_cleavage <- get_peptides_after_cleavage(sequences[[i]], cleavage_site_mapping[[i]])

    #append prior and posterior pepetides to cleavage map
    cleavage_site_mapping[[i]]$PriorPeptide <- peptides_before_cleavage
    cleavage_site_mapping[[i]]$PosteriorPeptide <- peptides_after_cleavage

    #create abundance object
    abundances[[i]] <- data.frame(
        peptide = c(peptides_before_cleavage[1], peptides_after_cleavage),
        start_site = c(1, cleavage_site_mapping[[i]]$Site),
        abundance = sequence_abundances
        )
  }

  return(list(abundances = abundances, cleavage_site_mapping = cleavage_site_mapping))
}

#given a sequence and cleavage mapping, returns the peptide before each clevage site
get_peptides_before_cleavage <- function(sequence, cleavage_site_mapping) {
  curr_start <- 1
  peptides <- c()
  for (i in seq_len(nrow(cleavage_site_mapping))) {
    curr_stop <- cleavage_site_mapping$Site[i] - 1
    peptides <- c(peptides, substr(sequence, start = curr_start, stop = curr_stop))
    curr_start <- curr_stop + 1
  }
  return(peptides)
}

#given a sequence and cleavage mapping, returns the peptide after each clevage site
get_peptides_after_cleavage <- function(sequence, cleavage_site_mapping) {
  peptides <- c()
  for (i in seq_len(nrow(cleavage_site_mapping) - 1)) {
    curr_start <- cleavage_site_mapping$Site[i]
    curr_stop <- cleavage_site_mapping$Site[i + 1] - 1
    peptides <- c(peptides, substr(sequence, start = curr_start, stop = curr_stop))
  }
  curr_start <- cleavage_site_mapping$Site[i + 1]
  curr_stop <- nchar(sequence)
  peptides <- c(peptides, substr(sequence, start = curr_start, stop = curr_stop))

  return(peptides)
}

get_cleavage_site_mapping <- function(mask_out, pk_missed_prop, tr_missed_prop, pk_prob = NULL, tr_prob = NULL) {
  #simulate imperfect digestion given missed cleavage probabilities
  digest_results_pk <- adjust_digestions(
    digest_results_pk(get_masked_sequences(mask_out)),
    prop_to_miss = pk_missed_prop,
    prob = pk_prob)
  digest_proteins_trypsin <- adjust_digestions(
    digest_proteins_trypsin(get_sequences(mask_out)),
    prop_to_miss = tr_missed_prop,
    prob = tr_prob)

  #get cleavage sites for proteinase k and trypsin
  pk_cleavage_sites <- lapply(digest_results_pk, \(x) x$start[-1])
  trypsin_cleavage_sites <- lapply(digest_proteins_trypsin, \(x) x$start[-1])

  #label pk sites as -1, and trypsin sites as +1, add them and label accordingly (for handling "shared" cleavage sites)
  pk_dat <- lapply(pk_cleavage_sites, \(x) data.frame(Site = x, PK = 1))
  trypsin_dat <- lapply(trypsin_cleavage_sites, \(x) data.frame(Site = x, TR = -1))
  full_dat <- lapply(
    seq_along(pk_dat),
    function(x) {
        dplyr::full_join(pk_dat[[x]], trypsin_dat[[x]], by = "Site") %>%
        tidyr::replace_na(., list(PK = 0, TR = 0)) %>%
        dplyr::mutate(
            sum_col = PK + TR,
            Enzyme = dplyr::case_match(
                dplyr::sum_col,
                -1 ~ "Trypsin",
                0 ~ "Proteinase K", # 0 if proteinase K because a trypsin cleavage cannot occur if a PK one already did
                1 ~ "Proteinase K"
                )
            ) %>%
        dplyr::select(Site, Enzyme) %>%
        dplyr::arrange(Site)
    }
  )

  return(full_dat)
}

adjust_digestions <- function(full_digest_results, ...) {
  lapply(full_digest_results, \(x) imperfect_digest(x, ...))
}

imperfect_digest <- function(ptm_site_mapping, prop_to_miss, prob = NULL) {
  n_peps <- nrow(ptm_site_mapping)

  # TODO does the probability of a successful split depend on the peptide length?
  ptm_site_mapping$merge_with_next <- seq_len(n_peps) %in% sample(seq_len(n_peps)[-n_peps], round(prop_to_miss * (n_peps - 1)), prob = prob)

  remove_vec <- rep(FALSE, nrow(ptm_site_mapping))

  for (i in seq_len(nrow(ptm_site_mapping) - 1)) {
    if (ptm_site_mapping$merge_with_next[i]) {
      peptide <- paste0(ptm_site_mapping$peptide[i], ptm_site_mapping$peptide[i + 1])
      start <- ptm_site_mapping$start[i]
      stop <- ptm_site_mapping$stop[i + 1]
      remove_vec[i] <- TRUE
      merge_with_next <- ptm_site_mapping$merge_with_next[i + 1]

      ptm_site_mapping$peptide[i + 1] <- peptide
      ptm_site_mapping$start[i + 1] <- start
      ptm_site_mapping$stop[i + 1] <- stop
      ptm_site_mapping$merge_with_next[i + 1] <- merge_with_next
    }
  }
  return(ptm_site_mapping[!remove_vec, ])
}

digest_proteins_pk <- function(sequences) {
  lapply(sequences, \(x) simulate_proteinase_k_cleavage(x))
}

digest_proteins_trypsin <- function(sequences) {
  lapply(sequences, \(x) simulate_trypsin_cleavage(x))
}

get_masked_sequences <- function(masking_list) {
  lapply(masking_list, \(x) x$masked_sequence)
}
get_sequences <- function(masking_list) {
  lapply(masking_list, \(x) x$sequence)
}

get_mask_indices <- function(sequence, start, delay_fn, max_masked_region_size, min_num_trypsin_sites) {
  if (start > length(amino_acids)) stop("Start position is larger than sequence length")
  adjusted_start <- start + delay_fn()
  curr_ind <- adjusted_start
  amino_acids <- unlist(lapply(sequence, \(x) strsplit(x, split = "")[[1]]))
  n_trypsin_sites <- 0
  region_size <- 0

  # Check if the adjusted start is valid
  if (adjusted_start > length(amino_acids)) {
    return(list(start = start, stop = Inf)) #stop at inf to flag finished procedure (this one will be tossed out)
  }

  #check if one of the conditions can be satisfied
  remaining_amino_acids <- amino_acids[adjusted_start:length(amino_acids)]
  remaining_spit_sites <- sum(grepl("R", remaining_amino_acids) | grepl("K", remaining_amino_acids))
  if (length(remaining_amino_acids) < max_masked_region_size & remaining_spit_sites < min_num_trypsin_sites) {
    return(list(start = adjusted_start, stop = Inf)) #stop at inf to flag finished procedure (this one will be tossed out)
  }

  while (n_trypsin_sites < min_num_trypsin_sites & region_size < max_masked_region_size) {
    curr_amino_acid <- amino_acids[curr_ind]
    if (curr_amino_acid == "K" | curr_amino_acid == "R") {
      n_trypsin_sites <- n_trypsin_sites + 1
    }
    region_size <- region_size + 1
    curr_ind <- curr_ind + 1
  }
  end_ind <- curr_ind + delay_fn()
  return(list(start = adjusted_start, stop = end_ind))
}

mask_sequence <- function(sequence, delay_fn, gap_fn, max_masked_region_size, min_num_trypsin_sites) {
  start_vec <- c()
  stop_vec <- c()
  curr_start <- 1 + gap_fn()
  last_stop <- 0
  masked_sequence <- sequence
  while (last_stop < nchar(sequence)) {
    mask_indices <- get_mask_indices(
      sequence = sequence,
      start = curr_start,
      delay_fn = delay_fn,
      max_masked_region_size = max_masked_region_size,
      min_num_trypsin_sites = min_num_trypsin_sites
    )
    last_stop <- mask_indices$stop
    curr_start <- mask_indices$stop + gap_fn()
    if (!(last_stop > nchar(sequence) | curr_start > nchar(sequence))) {
      start_vec <- c(start_vec, mask_indices$start)
      stop_vec <- c(stop_vec, mask_indices$stop)
      masked_sequence <- replace_between_indices(masked_sequence, mask_indices$start, mask_indices$stop, "#")
    }
  }
  return(list(sequence = sequence, masked_sequence = masked_sequence, start_vec = start_vec, stop_vec = stop_vec))
}

replace_between_indices <- function(x, start, stop, replacement) {
  before <- substr(x, 1, start - 1)
  after <- substr(x, stop + 1, nchar(x))
  replacement_string <- strrep(replacement, stop - start + 1)
  result <- paste0(before, replacement_string, after)
  return(result)
}

mask_all_proteins <- function(sequences, delay_fn, gap_fn, max_masked_region_size, min_num_trypsin_sites) {
  out_list <- vector("list", length = length(sequences))
  for (i in seq_along(sequences)) {
    out_list[[i]] <- mask_sequence(
      sequence = sequences[[i]],
      delay_fn = sample_delay,
      gap_fn = sample_gap,
      max_masked_region_size = max_masked_region_size,
      min_num_trypsin_sites = min_num_trypsin_sites
    )
  }
  return(out_list)
}

downsample_masking <- function(masking_output, downsample_fn, ...) {
  indices <- downsample_fn(seq_along(masking_output$start_vec), ...)

  masking_output$start_vec <- masking_output$start_vec[indices]
  masking_output$stop_vec <- masking_output$stop_vec[indices]
  masking_output$masked_sequence <- remask_sequence(masking_output$sequence, masking_output$start_vec, masking_output$stop_vec)
  return(masking_output)
}

remask_sequence <- function(sequence, start_vec, stop_vec) {
  masked_sequence <- sequence
  for (i in seq_along(start_vec)) {
    masked_sequence <- replace_between_indices(masked_sequence, start_vec[i], stop_vec[i], "#")
  }
  return(masked_sequence)
}

downsample_maskings <- function(masking_list, downsample_fn, ...) {
  out_list <- vector("list", length = length(masking_list))
  for (i in seq_along(masking_list)) {
    out_list[[i]] <- downsample_masking(masking_output = masking_list[[i]], downsample_fn = downsample_fn, ...)
  }
  return(out_list)
}

simulate_proteinase_k_cleavage <- function(sequence) {
  # Define amino acids considered for cleavage (aliphatic, aromatic, hydrophobic)
  cleavage_aa <- c("A", "V", "L", "I", "F", "Y", "W", "M", "P")  # Aliphatic, aromatic, and hydrophobic residues
  sequence_vector <- unlist(strsplit(sequence, ""))
  peptides <- c()
  start_indices <- c()
  end_indices <- c()

  start <- 1
  for (i in 1:(length(sequence_vector) - 1)) {
    if (sequence_vector[i] %in% cleavage_aa) {
      peptides <- c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices <- c(start_indices, start)
      end_indices <- c(end_indices, i)
      start <- i + 1
    }
  }

  # Add the final peptide after the last cleavage site
  if (start <= length(sequence_vector)) {
    peptides <- c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices <- c(start_indices, start)
    end_indices <- c(end_indices, length(sequence_vector))
  }

  cleavage_df <- data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}

insert_char_at_indices <- function(string, indices, char = "|") {
  indices <- sort(indices)
  indices <- indices[indices >= 1 & indices <= nchar(string)]

  chars <- unlist(strsplit(string, ""))

  offset <- 0
  for (index in indices) {
    chars <- append(chars, char, after = index + offset - 1)
    offset <- offset + 1
  }

  new_string <- paste(chars, collapse = "")

  return(new_string)
}

add_cleavage_characters <- function(masked_sequence_dat, digest_results, char = "|") {
  masked_sequence_dat$masked_sequence <- insert_char_at_indices(masked_sequence_dat$masked_sequence, digest_results$start, char = char)
  masked_sequence_dat$sequence <- insert_char_at_indices(masked_sequence_dat$sequence, digest_results$start, char = char)

  return(masked_sequence_dat)
}

simulate_trypsin_cleavage <- function(sequence) {
  cleavage_aa <- c("K", "R")
  skip_aa <- "P"

  cleavage_rule <- function(prev_aa, curr_aa) {
    return(curr_aa %in% cleavage_aa && prev_aa != skip_aa)
  }

  sequence_vector <- unlist(strsplit(sequence, ""))
  peptides <- c()
  start_indices <- c()
  end_indices <- c()

  start <- 1
  for (i in 1:(length(sequence_vector) - 1)) {
    if (cleavage_rule(sequence_vector[i], sequence_vector[i + 1])) {
      peptides <- c(peptides, paste(sequence_vector[start:i], collapse = ""))
      start_indices <- c(start_indices, start)
      end_indices <- c(end_indices, i)
      start <- i + 1
    }
  }

  if (start <= length(sequence_vector)) {
    peptides <- c(peptides, paste(sequence_vector[start:length(sequence_vector)], collapse = ""))
    start_indices <- c(start_indices, start)
    end_indices <- c(end_indices, length(sequence_vector))
  }

  cleavage_df <- data.frame(
    peptide = peptides,
    start = start_indices,
    stop = end_indices,
    stringsAsFactors = FALSE
  )

  return(cleavage_df)
}