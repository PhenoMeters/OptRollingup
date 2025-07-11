


custom_rrollup = function (pepData, combine_fn, parallel = TRUE) 
{
  browser()
  if (!inherits(pepData, "pepData")) {
    stop("pepData is not an object of the appropriate class")
  }
  if (is.null(pepData$e_meta)) {
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname
  temp <- merge(x = pepData$e_meta[, c(pep_id, pro_id)], y = pepData$e_data, 
                by = pep_id, all.x = FALSE, all.y = TRUE) %>% dplyr::select(-dplyr::sym(pep_id)) %>% 
    data.frame(check.names = FALSE)
  unique_proteins <- unique(temp[[pro_id]])
  if (parallel == TRUE) {
    cores <- parallelly::availableCores()
    cl <- parallelly::makeClusterPSOCK(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  }
  else {
    foreach::registerDoSEQ()
  }
  r <- foreach::foreach(i = 1:length(unique_proteins)) %dopar% 
    {
      row_ind <- which(temp[, pro_id] == unique_proteins[i])
      current_subset <- temp[row_ind, ]
      current_subset <- current_subset[, -which(names(temp) == 
                                                  pro_id)]
      num_peps = nrow(current_subset)
      res = matrix(NA, nrow = 1, ncol = ncol(current_subset))
      if (num_peps == 1) {
        protein_val = unlist(current_subset)
      }
      else {
        na.cnt = apply(is.na(current_subset), 1, sum)
        least.na = which(na.cnt == min(na.cnt))
        if (length(least.na) > 1) {
          mds = apply(current_subset, 1, median, na.rm = TRUE)[least.na]
          least.na = least.na[which(mds == max(mds))]
        }
        prot_val = unlist(current_subset[least.na, ])
        scaling_factor = apply(rep(as.numeric(prot_val), 
                                   each = nrow(current_subset)) - current_subset, 
                               1, median, na.rm = TRUE)
        x_scaled = current_subset + rep(scaling_factor, 
                                        ncol(current_subset))
        protein_val = apply(x_scaled, 2, combine_fn)
      }
      res[1, ] = protein_val
      res <- data.frame(res)
      names(res) <- names(current_subset)
      res
    }
  final_result <- data.frame(unique_proteins, data.table::rbindlist(r), 
                             check.names = FALSE)
  names(final_result)[1] <- pro_id
  return(list(e_data = final_result, e_meta = NULL))
}


