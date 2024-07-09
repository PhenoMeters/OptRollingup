get_binned_abundaces = function(omicsData, bins){
  browser()
  edata_cname = get_edata_cname(omicsData)
  
  raw_vals <- omicsData$e_data[-which(names(omicsData$e_data) == edata_cname)]
  
  # Generate means and missingness percents per each peptide
  means <- apply(raw_vals, 1, \(x) mean(x, na.rm = TRUE))
  missingness_pct <- apply(raw_vals, 1, \(x) 100 * length(which(is.na(x))) / length(names(raw_vals)))
  
  # Cut into 10 bins with approximately equal elements in each bin
  q_bins <- cut_number(means, bins)
  
  # For plotting: change bins to be a numeric value
  n_bins <- q_bins
  levels(n_bins) <- seq_len(bins)
  n_bins <- as.numeric(n_bins)
  
  # Calculate the bin medians and use this to calculate the probability of selection
  bin_medians <- sapply(seq_len(bins), \(x) median(raw_vals[which(q_bins == levels(q_bins)[x])]))
  
  return(list(binned_peps = binned_peps, bin_medians = bin_medians))
}




simulate_missingness = function(pepObj, prop_missingness, binned_peps, bin_probs, verbose = T, p_adjust = F, TMT = T){
  if(p_adjust & (min(bin_probs) == 0)){
    bin_probs[which.min(bin_probs)] = min(bin_probs[-which.min(bin_probs)]) / 10
    bin_probs = bin_probs / sum(bin_probs)
  }
  
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
  
  return(pepObj)
}