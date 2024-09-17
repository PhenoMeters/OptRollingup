library(ggplot2)
library(dplyr)
library(OrgMassSpecR)
library(pmartR)
library(tidyverse)

######################################################
### Helper Functions
######################################################

here::i_am("Code/protein_sequence_simulation.R")

source(here::here("Code", "ptm_utils.R"))
source(here::here("Code", "LiP_utils.R"))


#Digestion
digest_proteins_PK = function(sequences){
  lapply(sequences, \(x) simulate_proteinase_k_cleavage(x))
}

digest_proteins_trypsin = function(sequences){
  lapply(sequences, \(x) simulate_trypsin_cleavage(x))
}

get_masked_sequences = function(masking_list){
  lapply(masking_list, \(x) x$masked_sequence)
}
get_sequences = function(masking_list){
  lapply(masking_list, \(x) x$sequence)
}

get_cleavage_site_mapping = function(mask_out, pk_missed_prop, tr_missed_prop, pk_prob = NULL, tr_prob = NULL){
  digest_results_PK = adjust_digestions(digest_proteins_PK(get_masked_sequences(mask_out)), prop_to_miss = pk_missed_prop, prob = pk_prob)
  digest_results_tryp = adjust_digestions(digest_proteins_trypsin(get_sequences(mask_out)), prop_to_miss = tr_missed_prop, prob = tr_prob)
  
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

imperfect_digest = function(
    ptm_site_mapping, 
    prop_to_miss, 
    prob = NULL
){
  n_peps = nrow(ptm_site_mapping)
  # coef_dat = data.frame(peptide = names(peptide_coef), Mult_Factor = exp(peptide_coef))
  # ptm_site_mapping = ptm_site_mapping %>% left_join(., coef_dat, by = "peptide")
  
  # TODO does the probability of a successful split depend on the peptide length?
  ptm_site_mapping$merge_with_next = seq_len(n_peps) %in% sample(seq_len(n_peps)[-n_peps], round(prop_to_miss * (n_peps - 1)), prob = prob)
  protein = ptm_site_mapping$protein %>% unique()
  
  remove_vec = rep(FALSE, nrow(ptm_site_mapping))
  # ptm_site_mapping = ptm_site_mapping %>% as_tibble()
  
  for(i in seq_len(nrow(ptm_site_mapping) - 1)){
    if(ptm_site_mapping$merge_with_next[i]){
      peptide = paste0(ptm_site_mapping$peptide[i], ptm_site_mapping$peptide[i + 1])
      start = ptm_site_mapping$start[i]
      stop = ptm_site_mapping$stop[i + 1]
      remove_vec[i] = TRUE
      merge_with_next = ptm_site_mapping$merge_with_next[i+1]
      
      
      ptm_site_mapping$peptide[i + 1] = peptide
      ptm_site_mapping$start[i + 1] = start
      ptm_site_mapping$stop[i + 1] = stop
      ptm_site_mapping$merge_with_next[i + 1] = merge_with_next
    }
  }
  return(ptm_site_mapping[!remove_vec,])
}

adjust_digestions = function(full_digest_results, ...){
  lapply(full_digest_results, \(x) imperfect_digest(x, ...))
}

get_cleavage_site_mapping = function(mask_out, pk_missed_prop, tr_missed_prop, pk_prob = NULL, tr_prob = NULL){
  digest_results_PK = adjust_digestions(digest_proteins_PK(get_masked_sequences(mask_out)), prop_to_miss = pk_missed_prop, prob = pk_prob)
  digest_results_tryp = adjust_digestions(digest_proteins_trypsin(get_sequences(mask_out)), prop_to_miss = tr_missed_prop, prob = tr_prob)
  
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

simulate_protein_for_sample = function(n_reps, mask_out, pk_missed_prop, tr_missed_prop, sequences, pk_prob = NULL, tr_prob = NULL){
  peptide_abundances = vector("list", n_reps)
  cleavage_mappings = vector("list", n_reps)
  for(i in seq_len(n_reps)){
    cleavage_mapping_example = get_cleavage_site_mapping(
      mask_out,
      pk_missed_prop,
      tr_missed_prop,
      pk_prob,
      tr_prob)
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

simulate_protein_for_group = function(n_subjects, group_id, n_reps, mask_out, pk_missed_prop, tr_missed_prop, sequences, pk_prob = NULL, tr_prob = NULL){
  full_pep_list = vector("list", n_subjects)
  full_cleavage_mapping = vector("list", n_subjects)
  for(i in seq_len(n_subjects)){
    out = simulate_protein_for_sample(n_reps, mask_out, pk_missed_prop, tr_missed_prop, sequences, pk_prob, tr_prob)
    full_pep_list[[i]] = out$peptide_df
    full_cleavage_mapping[[i]] = out$cleavage_mapping
  }
  full_pep_list = purrr::reduce(full_pep_list, full_join, by = c("start_site", "peptide", "protein"))
  full_cleavage_mapping = purrr::reduce(full_cleavage_mapping, rbind)
  names(full_pep_list)[-c(1:3)] = paste("subject", seq_len(n_subjects), "group", group_id, sep = '_')
  return(list(peptide_df = full_pep_list, cleavage_mapping = full_cleavage_mapping))
}


generate_f_data = function(full_abundance_df){
  subject_ids = names(full_abundance_df)[-which(names(full_abundance_df) %in% c('start_site', 'peptide', 'protein'))]
  f_data = data.frame(Subject_IDs = subject_ids) %>% 
    mutate(Group_IDs = sapply(strsplit(Subject_IDs, '_'), function(x) x[4]))
  return(f_data)
}


### multi protein rollup
# These functions could use some cleaning up, but they work for now and for 2 groups
site_rollup_analysis = function(g1_dat, g2_dat){
  
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
  
  full_cleavage_mapping = rbind(
    g1_dat$cleavage_mapping %>% distinct,
    g2_dat$cleavage_mapping %>% distinct
  ) %>% 
    distinct() %>%
    pivot_longer(c("PriorPeptide", "PosteriorPeptide"), values_to = "Peptide", names_to = "Type")
  
  
  data_for_rollup = full_pep_data %>% 
    rename(Site = start_site, Peptide = peptide, protein = protein) %>% 
    left_join(., full_cleavage_mapping, by = c("Site", "Peptide", "protein")) %>% 
    drop_na(., Enzyme, Type)
  
  det_nondet_data = data_for_rollup %>%
    mutate(across(-c(Peptide, protein, Site, Enzyme, Type), ~ ifelse(is.na(.), 0, 1)))
  
  det_nondet_rolled_up_data = det_nondet_data %>% group_by(Site, Enzyme, protein) %>% select(-Type) %>% summarize_if(is.numeric, mean, na.rm = T)
  
  long_data <- det_nondet_rolled_up_data %>%
    pivot_longer(cols = -c(Site, Enzyme, protein), names_to = "SampleID", values_to = "Y")
  
  long_data <- long_data %>%
    mutate(group = ifelse(grepl("group_1", SampleID), 1, 2))
  
  res = long_data %>% 
    mutate(Y_bin = ifelse(Y > 0.5, 1, 0)) %>%
    ungroup() %>% 
    nest_by(Site, Enzyme, protein) %>%
    mutate(res = list(glm(Y_bin ~ group, data = data, family = binomial))) %>%
    mutate(logodds = sum(coef(res) * c(0, 1)),
           est_prop = exp(logodds) / (1 + exp(logodds)))
  
  return(res)
}

tryptic_peptide_analysis = function(g1_dat, g2_dat){
  
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
  
  
  full_cleavage_mapping = rbind(
    g1_dat$cleavage_mapping %>% distinct,
    g2_dat$cleavage_mapping %>% distinct
  ) %>% 
    distinct() %>%
    as_tibble()
  
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
  res = tryptic_pep_data %>% imd_anova(., test_method = "anova")
  
  return(res)
}


create_site_plot = function(res, g1_masking, g2_masking){
  ggdat = res %>% 
    select(
      Site,
      Enzyme,
      logodds
    )
  
  masking_data = c()
  for(idx in seq_along(g1_masking)){
    tmp_masking_data = data.frame(
      Site = seq_len(nchar(g1_masking[[idx]]$sequence)),
      group_1 = strsplit(g1_masking[[idx]]$masked_sequence, '')[[1]],
      group_2 = strsplit(g2_masking[[idx]]$masked_sequence, '')[[1]]) %>%
      mutate(
        cat = paste0(group_1, group_2),
        category = case_when(
          cat == "##" ~ "Both",
          grepl("#[A-Za-z]", cat) ~ "Group 1 Only",
          grepl("[A-Za-z]#", cat) ~ "Group 2 Only",
          grepl("[A-Za-z][A-Za-z]", cat) ~ "Neither"
        ),
        protein = g1_masking[[idx]]$sequence
      )
    masking_data = rbind(masking_data, tmp_masking_data)
  }
  all_sites <- data.frame(Site = seq(min(masking_data$Site), max(masking_data$Site)))
  
  
  plot_data_1 <- all_sites %>%
    left_join(masking_data, by = "Site") %>%
    mutate(Enzyme = "Trypsin")
  plot_data_2 <- all_sites %>%
    left_join(masking_data, by = "Site") %>%
    mutate(Enzyme = "Proteinase K")
  plot_data = rbind(plot_data_1, plot_data_2) %>%
    select(-group_1, -group_2, -cat) %>%
    left_join(., ggdat, by = c("Site", "Enzyme", "protein"))
  
  
  #recode proteins
  unique_values <- unique(plot_data$protein)
  recode_map <- setNames(paste0("Sequence ", seq_along(unique_values)), unique_values)
  plot_data = plot_data %>%
    mutate(protein = recode(protein, !!!recode_map))
  ggdat = ggdat %>%
    mutate(protein = recode(protein, !!!recode_map))
  
  
  # Separate the data for the line plot (where FC is not NA)
  line_data <- plot_data %>%
    filter(!is.na(logodds))
  
  # Creating the plot with faceting
  p = ggplot() + 
    geom_tile(data = plot_data, aes(x = Site, y = (max(line_data$logodds, na.rm = TRUE) + min(line_data$logodds, na.rm = TRUE)) / 2, fill = category), height = max(line_data$logodds, na.rm = TRUE) - min(line_data$logodds, na.rm = TRUE), alpha = 0.2) + 
    geom_line(data = line_data, aes(x = Site, y = logodds)) + 
    theme_bw() +
    scale_fill_manual(values = c("Neither" = "gray", "Both" = "blue", "Group 1 Only" = "orange", "Group 2 Only" = "green")) +
    labs(fill = "Region Masked?") + 
    ylab("Difference in Proportion of Samples with Cleavage") +
    facet_grid(Enzyme ~ protein, scales = "free")
  
  return(p)
}


create_trypic_plot = function(res, g1_masking, g2_masking){
  ggdat = res %>% 
    select(Peptide_Site, 
           Fold_change_1_vs_2) %>% 
    rename(FC = Fold_change_1_vs_2) %>%
    mutate(
      Peptide = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][2]),
      Site = as.numeric(sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][1])),
      protein = sapply(Peptide_Site, \(x) strsplit(x, split = '_')[[1]][3]),
      Enzyme = "Trypsin"
    ) %>%
    select(-Peptide_Site)
  
  masking_data = c()
  for(idx in seq_along(g1_masking)){
    tmp_masking_data = data.frame(
      Site = seq_len(nchar(g1_masking[[idx]]$sequence)),
      group_1 = strsplit(g1_masking[[idx]]$masked_sequence, '')[[1]],
      group_2 = strsplit(g2_masking[[idx]]$masked_sequence, '')[[1]]) %>%
      mutate(
        cat = paste0(group_1, group_2),
        category = case_when(
          cat == "##" ~ "Both",
          grepl("#[A-Za-z]", cat) ~ "Group 1 Only",
          grepl("[A-Za-z]#", cat) ~ "Group 2 Only",
          grepl("[A-Za-z][A-Za-z]", cat) ~ "Neither"
        ),
        protein = g1_masking[[idx]]$sequence
      )
    masking_data = rbind(masking_data, tmp_masking_data)
  }
  all_sites <- data.frame(Site = seq(min(masking_data$Site), max(masking_data$Site)))
  
  
  plot_data_1 <- all_sites %>%
    left_join(masking_data, by = "Site") %>%
    mutate(Enzyme = "Trypsin")
  plot_data_2 <- all_sites %>%
    left_join(masking_data, by = "Site") %>%
    mutate(Enzyme = "Proteinase K")
  plot_data = rbind(plot_data_1, plot_data_2) %>%
    select(-group_1, -group_2, -cat) %>%
    left_join(., ggdat, by = c("Site", "Enzyme", "protein"))
  
  
  plot_data = plot_data %>% filter(Enzyme == "Trypsin")
  
  #recode proteins
  unique_values <- unique(plot_data$protein)
  recode_map <- setNames(paste0("Sequence ", seq_along(unique_values)), unique_values)
  plot_data = plot_data %>%
    mutate(protein = recode(protein, !!!recode_map))
  ggdat = ggdat %>%
    mutate(protein = recode(protein, !!!recode_map))
  
  
  # Separate the data for the line plot (where FC is not NA)
  line_data <- plot_data %>%
    filter(!is.na(FC))
  
  # Creating the plot with faceting
  p = ggplot() + 
    geom_tile(data = plot_data, aes(x = Site, y = (max(line_data$FC, na.rm = TRUE) + min(line_data$FC, na.rm = TRUE)) / 2, fill = category), height = max(line_data$FC, na.rm = TRUE) - min(line_data$FC, na.rm = TRUE), alpha = 0.2) + 
    geom_line(data = line_data, aes(x = Site, y = FC)) + 
    theme_bw() +
    scale_fill_manual(values = c("Neither" = "gray", "Both" = "blue", "Group 1 Only" = "orange", "Group 2 Only" = "green")) +
    labs(fill = "Region Masked?") + 
    ylab("Difference in Proportion of Samples with Cleavage") +
    facet_grid(Enzyme ~ protein, scales = "free")
  return(p)
}


######################################################
### Simulation
######################################################

# Create amino acid distribution
lip_data = readr::read_tsv(here::here("Data", "double_pept_lytic_sites.tsv"))
all_peptides = unlist(lapply(lip_data$Sequence, \(x) strsplit(x, split = '; ')[[1]]))
unique_peptides = unique(all_peptides)
amino_acids = unlist(lapply(unique_peptides, \(x) strsplit(x, split = "")[[1]]))
amino_acid_distribution = create_amino_acid_distribution(amino_acids)

#Generate one synthetic protein
synthetic_proteins = list(
  generate_protein(1000, amino_acid_distribution),
  generate_protein(1000, amino_acid_distribution),
  generate_protein(1000, amino_acid_distribution),
  generate_protein(1000, amino_acid_distribution),
  generate_protein(1000, amino_acid_distribution),
  generate_protein(1000, amino_acid_distribution)
)

# Tuning parameters for protein masking
sample_gap = function() rpois(1, 20)
sample_delay = function() rpois(1, 3)
downsample_fn = function(vec, prop) sample(vec, floor(prop * length(vec)))
max_masked_region_size = 100
min_num_trypsin_sites = 2

# Perform full masking
full_masking = mask_all_proteins(
  sequences = synthetic_proteins,
  delay_fn = sample_delay,
  gap_fn = sample_gap,
  max_masked_region_size = max_masked_region_size,
  min_num_trypsin_sites = min_num_trypsin_sites
)

#Create groups by downsampling
g1_masking = downsample_maskings(full_masking, downsample_fn, prop = 0.5)

g2_masking = downsample_maskings(full_masking, downsample_fn, prop = 0.5)

n_g1 = 20
n_g2 = 20
trypsin_missed_prob = 0.2
proteinase_k_missed_prob = 0.2
n_proteins = 200

g1_dat = simulate_protein_for_group(
  n_g1, 
  1, 
  n_proteins, 
  g1_masking, 
  proteinase_k_missed_prob,
  trypsin_missed_prob,  
  synthetic_proteins)

g2_dat = simulate_protein_for_group(
  n_g2, 
  2, 
  n_proteins, 
  g2_masking, 
  proteinase_k_missed_prob, 
  trypsin_missed_prob, 
  synthetic_proteins)

sw_res = site_rollup_analysis(g1_dat, g2_dat)
tryptic_pep_res = tryptic_peptide_analysis(g1_dat, g2_dat)

p1 = create_site_plot(sw_res, g1_masking, g2_masking)
ggsave(here::here("Out", "site_test_plot.png"), p1, width = 30, height = 4)

p2 = create_trypic_plot(tryptic_pep_res, g1_masking, g2_masking)
ggsave(here::here("Out", "tryptic_test_plot.png"), p2, width = 30, height = 4)



# TODO: synthetic TMT missingness?
# TODO: wrap up in for loop
# TODO: more sophisticated imperfect digestion


