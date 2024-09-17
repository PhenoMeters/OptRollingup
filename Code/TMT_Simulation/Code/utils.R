# Below are functions that help simulate and visualize missingness in pmartR omicsData objects. 

# A function to get the plex distribution for an omicsdata
    # obj: A pmartR omicsData
    # plex_cname: The plex column name in the omicsData fdata
    # samps_per_ples: Integer: The number of samples per plex in omicsData
get_plex_dist = function(obj, plex_cname, pool_id = NULL, samps_per_plex = NULL) {
  plexes = unique(obj$f_data[[plex_cname]])
  
  if (is.null(samps_per_plex)) {
    samps_per_plex = length(which(obj$f_data[[plex_cname]] == plexes[1]))
  }
  
  fdata_cname = get_fdata_cname(obj)
  fdata_no_pool = obj$f_data
  
  pool_samp = 0
  if (!is.null(pool_id)) {
    fdata_no_pool = fdata_no_pool[-which(is.na(fdata_no_pool[[pool_id]])),]
    pool_samp = 1
  }
  
  results = data.frame(miss_count = 0:(samps_per_plex - pool_samp),
                        count = rep(0, samps_per_plex + 1 - pool_samp))
  
  
  for (plex in plexes) {
    miss_counts = apply(
      obj$e_data[,
                 fdata_no_pool[fdata_no_pool[[plex_cname]] == plex,][[fdata_cname]]
      ], 1, \(x) length(which(is.na(x)))
    )
    
    for (count in 0:(samps_per_plex - pool_samp)) {
      results$count[count + 1] = results$count[count + 1] + length(which(miss_counts == count))
    }
  }
  
  class(results) = c("plexdist", class(results))
  
  return(results)
}

# A function to plot the plex distribution of an omicsData
    # x: The plex distribution from an omicsData. Obtained via get_plex_dist
plot.plexdist = function(x) {
  x = x %>% 
    mutate(
      miss_count = as.factor(miss_count),
      prop = count / sum(.$count) * 100
    ) %>%
    mutate(
      prop_lab = paste0(round(prop, 2), "%")
    )
  ggplot(x, aes(x = miss_count, y = prop)) + 
    geom_bar(stat = "identity") +
    xlab("Number of Observations Missing per Plex") +
    ylab("Percentage of Plexes") +
    geom_text(aes(label = prop_lab), vjust = -0.5) + 
    theme_bw() +
    theme(text = element_text(size = 18))
}

#Function to generate missingness bins based on the current missing data structure of omicsData
    # omicsData: A pmartR omicsData object which has missingness. 
get_missingness_bins = function(omicsData, num_bins = 10){
  edata_cname = get_edata_cname(omicsData)
  
  raw_vals = omicsData$e_data[-which(names(omicsData$e_data) == edata_cname)]
  
  # Generate means and missingness percents per each peptide
  means = apply(raw_vals, 1, \(x) mean(x, na.rm = TRUE))
  missingness_pct = apply(raw_vals, 1, \(x) 100 * length(which(is.na(x))) / length(names(raw_vals)))
  
  # Cut into 10 bins with approximately equal elements in each bin
  q_bins = cut_number(means, num_bins)
  
  # For plotting: change bins to be a numeric value
  n_bins = q_bins
  levels(n_bins) = seq_len(num_bins)
  n_bins = as.numeric(n_bins)
  
  #since our data is simulated and we dont have missing data, we approximate the missing data mechanism as a function of the mean abundances
  bin_means = sapply(seq_len(num_bins), \(x) mean(means[which(q_bins == levels(q_bins)[x])]))
  bin_probs = exp(-log2(bin_means)) / sum(exp(-log2(bin_means)))
  
  binned_peps = list()
  for (i in seq_len(num_bins)) {
    binned_peps[[i]] = omicsData$e_data[[edata_cname]][which(n_bins == i)]
  }
  
  return(list(binned_peps = binned_peps, bin_probs = bin_probs, bin_means = bin_means))
}


# Function to simulate missingness for either TMT or label free missing data mechanisms
    # pepObj: A pmartR peptide object which has been subsetted down to the complete case
    # prop_missingness: Float. Proportion of missing data
    # binned_peps: Located in results from get_missingness_bins
    # bin_probs: Located in results from get_missingness_bins
    # verbose: Boolean. If true verbose output
    # TMT: Boolean. If true simulated with TMT missingness, if false, label-free
simulate_missingness = function(pepObj, prop_missingness, binned_peps, bin_probs, verbose = T, TMT = T, p_adjust = F){
  # if(p_adjust & (min(bin_probs) == 0)){
  #   bin_probs[which.min(bin_probs)] = min(bin_probs[-which.min(bin_probs)]) / 10
  #   bin_probs = bin_probs / sum(bin_probs)
  # }
  fdata_cname = get_fdata_cname(pepObj)
  names(pepObj$f_data)[which(names(pepObj$f_data) == fdata_cname)] = "Subject_ID"
  
  #estimate the number of times a peptide must go completely missing within a plex
  n_obs = nrow(pepObj$e_data) * (ncol(pepObj$e_data) - 1)
  n_obs_go_missing = n_obs * prop_missingness
  
  binned_peps_tmp = binned_peps
  while(sum(is.na(pepObj$e_data[-which(names(pepObj$e_data) == get_edata_cname(pepObj))])) < n_obs_go_missing){
    
    if(verbose) cat("Current missingness: ", sum(is.na(pepObj$e_data[-which(names(pepObj$e_data) == get_edata_cname(pepObj))])), "| Target Missingness :", n_obs_go_missing, "\r")
    
    #sample bin
    #repeat sampling until we choose a bin that has peptides that can go missing
    n_peps_bin = 0
    while(n_peps_bin == 0){
      
      #adjust probability of selecting a bin
      which_bins_empty = (sapply(binned_peps_tmp, length) == 0)
      adj_probs = ifelse(which_bins_empty, 0, bin_probs)
      adj_probs = adj_probs / sum(adj_probs)
      
      #only sample from the bins that are not empty
      bin = sample(seq_len(10), 1, replace = F, prob = adj_probs)
      n_peps_bin = length(binned_peps_tmp[[bin]])
    }
    
    #randomly choose peptide in that bin
    pep_go_missing = sample(binned_peps_tmp[[bin]], 1, replace = F)
    
    if(TMT){
      #get the plexes that dont have missing data for chosen peptide
      curr_pep_vals = pepObj$e_data[which(pepObj$e_data[,1] == pep_go_missing),-which(names(pepObj$e_data) == get_edata_cname(pepObj))]
      valid_plexes = suppressMessages(left_join(pepObj$f_data %>% select(Plex, Subject_ID), tibble(vals = c(t(curr_pep_vals)), Subject_ID = names(curr_pep_vals)))) %>% drop_na %>% pull(Plex) %>% unique
      
      #remove all data from that peptide from randomly selected plex
      rm_plex = sample(valid_plexes, 1)
      
      #set all values within that plex for that peptide to NA
      rm_subjs = pepObj$f_data$Subject_ID[pepObj$f_data$Plex == rm_plex]
      pepObj$e_data[which(pepObj$e_data[,1] == pep_go_missing),which(names(pepObj$e_data) %in% rm_subjs)] = NA
      
      #remove from binned_peps_tmp if no more plexes to sample
      if(length(valid_plexes) == 1){
        binned_peps_tmp[[bin]] = binned_peps_tmp[[bin]][-which(binned_peps_tmp[[bin]] == pep_go_missing)]
      }
    } else {
      subj_id = names(pepObj$e_data)[names(pepObj$e_data) != get_edata_cname(pepObj)]
      curr_pep_vals = pepObj$e_data[which(pepObj$e_data[,1] == pep_go_missing),which(names(pepObj$e_data) != get_edata_cname(pepObj))]
      valid_subj = subj_id[which(!is.na(curr_pep_vals))]
      
      
      rm_subj = sample(valid_subj, 1)
      pepObj$e_data[which(pepObj$e_data[,1] == pep_go_missing),which(names(pepObj$e_data) == rm_subj)] = NA
      
      #remove from binned_peps_tmp if no more plexes to sample
      if(length(valid_subj) == 1){
        binned_peps_tmp[[bin]] = binned_peps_tmp[[bin]][-which(binned_peps_tmp[[bin]] == pep_go_missing)]
      }
    }
  }
  
  names(pepObj$f_data)[which(names(pepObj$f_data) == "Subject_ID")] = fdata_cname
  return(pepObj)
}