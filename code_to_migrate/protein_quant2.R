protein_quant2 = function (pepData, method, isoformRes = NULL, qrollup_thresh = NULL, 
          single_pep = FALSE, single_observation = FALSE, combine_fn = "median", 
          parallel = TRUE, emeta_cols = NULL) 
{
  if (get_data_scale(pepData) == "abundance") {
    stop(paste("The data must be log transformed. Use edata_transform to", 
               "convert to the log scale.", sep = " "))
  }
  if (!inherits(pepData, "pepData")) 
    stop("pepData must be an object of class pepData")
  if (!(method %in% c("rollup", "rrollup", "qrollup", "zrollup"))) 
    stop("method must be one of, rollup, rrollup, qrollup, zrollup")
  if (!(combine_fn %in% c("median", "mean", "sum"))) 
    stop("combine_fn must be either 'mean', 'sum', or 'median'")
  if (method != "zrollup" && (single_pep == TRUE || single_observation == 
                              TRUE)) 
    message("single_pep and single_observation will be ignored, as they are only applicable if method is zrollup")
  if (method != "qrollup" && !is.null(qrollup_thresh)) 
    message("qrollup_thresh argument will be ignored, as it is only applicable if method is qrollup")
  if (!is.null(isoformRes) && !inherits(isoformRes, "isoformRes")) {
    stop("The input for isoformRes must be of class 'isoformRes'.")
  }
  if (!is.logical(single_pep)) 
    stop("sinlge_pep must be either TRUE or FALSE.")
  if (!is.logical(single_observation)) {
    stop("sinlge_observation must be either TRUE or FALSE.")
  }
  if (!is.null(emeta_cols) && !is.character(emeta_cols)) {
    stop("emeta_cols must be a character vector.")
  }
  if (combine_fn == "median") {
    chosen_combine_fn <- pmartR:::combine_fn_median
  }
  else if(combine_fn == "sum"){
    chosen_combine_fn <- combine_fn_sum
  } else {
    chosen_combine_fn <- pmartR:::combine_fn_mean
  }
  edata_cname <- attr(pepData, "cnames")$edata_cname
  fdata_cname <- attr(pepData, "cnames")$fdata_cname
  emeta_cname <- attr(pepData, "cnames")$emeta_cname
  edata_cname_id <- which(names(pepData$e_data) == edata_cname)
  data_scale <- get_data_scale(pepData)
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized
  if (!is.null(isoformRes)) {
    e_meta <- pepData$e_meta
    filtas <- attr(pepData, "filters")
    groupies <- attr(pepData, "group_DF")
    scales <- get_data_scale_orig(pepData)
    inovas <- attr(pepData, "imdanova")
    isoformRes2 <- attr(isoformRes, "isoformRes_subset")
    peptides <- which(pepData$e_data[, edata_cname] %in% 
                        isoformRes2[, edata_cname])
    pepData <- as.pepData(e_data = pepData$e_data[peptides, 
    ], f_data = pepData$f_data, e_meta = isoformRes2, 
    edata_cname = edata_cname, fdata_cname = fdata_cname, 
    emeta_cname = "Protein_Isoform", data_scale = data_scale, 
    is_normalized = is_normalized)
    attr(pepData, "filters") <- filtas
    attr(pepData, "group_DF") <- groupies
    attr(pepData, "data_info")$data_scale_orig <- scales
    attr(pepData, "imdanova") <- inovas
    e_meta_iso <- pepData$e_meta %>% dplyr::select(dplyr::any_of(c(emeta_cname, 
                                                                   "Protein_Isoform")))
    emeta_cname_iso <- "Protein_Isoform"
  }
  if (method == "rollup") {
    results <- pmartR:::pquant(pepData = pepData, combine_fn = chosen_combine_fn)
  }
  if (method == "rrollup") {
    results <- pmartR:::rrollup(pepData, combine_fn = chosen_combine_fn, 
                       parallel = parallel)
  }
  if (method == "qrollup") {
    if (is.null(qrollup_thresh)) {
      stop("qrollup_thresh parameter value must be specified")
    }
    results <- pmartR:::qrollup(pepData, qrollup_thresh = qrollup_thresh, 
                       combine_fn = chosen_combine_fn, parallel = parallel)
  }
  if (method == "zrollup") {
    proteomicsfilt <- pmartR::proteomics_filter(pepData)
    moleculefilt <- pmartR::molecule_filter(pepData)
    if (single_pep == FALSE && single_observation == FALSE) {
      if (any(moleculefilt$Num_Observations == 1) && any(proteomicsfilt$counts_by_pro == 
                                                         1)) 
        stop("There are peptides with less than 2 observations and proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides/proteins and then run zrollup, set both 'single_observation' and 'single_pep' arguments to TRUE")
      if (any(moleculefilt$Num_Observations == 1)) 
        stop("There are peptides with less than 2 observations. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
      if (any(proteomicsfilt$counts_by_pro == 1)) 
        stop("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE")
    }
    if (single_pep == TRUE && single_observation == FALSE) {
      pepData = applyFilt(proteomicsfilt, pepData, min_num_peps = 2)
      if (any(moleculefilt$Num_Observations == 1)) 
        stop("There are peptides with less than 2 observations. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
    }
    if (single_pep == FALSE && single_observation == TRUE) {
      pepData = applyFilt(moleculefilt, pepData)
      if (any(proteomicsfilt$counts_by_pro == 1)) 
        stop("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE.")
    }
    if (single_pep == TRUE && single_observation == TRUE) {
      pepData0 = applyFilt(moleculefilt, pepData, min_num = 2)
      tmp_filt = proteomics_filter(pepData0)
      pepData = applyFilt(tmp_filt, pepData0, min_num_peps = 2)
    }
    results <- pmartR:::zrollup(pepData, combine_fn = chosen_combine_fn, 
                       parallel = parallel)
  }
  if (!is.null(emeta_cols)) {
    n_row_edata <- nrow(results$e_data)
    n_row_emeta <- if (is.null(results$e_meta)) {
      if (is.null(isoformRes)) {
        nrow(unique(pepData$e_meta))
      }
      else {
        nrow(e_meta_iso)
      }
    }
    else {
      nrow(unique(results$e_meta))
    }
    if (n_row_edata != n_row_emeta) {
      emeta_cols <- NULL
    }
  }
  if (is.null(results$e_meta)) {
    if (!is.null(isoformRes)) {
      results$e_meta <- e_meta_iso
    }
    else {
      results$e_meta <- pepData$e_meta
    }
  }
  if (is.null(isoformRes)) {
    results$e_meta <- results$e_meta %>% dplyr::group_by(!!dplyr::sym(emeta_cname)) %>% 
      dplyr::mutate(peps_per_pro = dplyr::n()) %>% dplyr::mutate(n_peps_used = if ("n_peps_used" %in% 
                                                                                   colnames(results$e_meta)) {
        n_peps_used
      }
      else {
        peps_per_pro
      }) %>% dplyr::relocate(n_peps_used, .after = dplyr::last_col()) %>% 
      dplyr::select(dplyr::any_of(c(emeta_cname, "peps_per_pro", 
                                    "n_peps_used", emeta_cols))) %>% dplyr::distinct() %>% 
      data.frame(check.names = FALSE)
  }
  else {
    peps_used <- if ("n_peps_used" %in% colnames(results$e_meta)) {
      results$e_meta %>% dplyr::select(!!dplyr::sym(emeta_cname), 
                                       n_peps_used) %>% dplyr::distinct()
    }
    else {
      isoformRes2 %>% dplyr::group_by(Protein_Isoform) %>% 
        dplyr::mutate(n_peps_used = dplyr::n()) %>% dplyr::slice(1) %>% 
        dplyr::ungroup() %>% dplyr::select(!!dplyr::sym(emeta_cname), 
                                           n_peps_used)
    }
    peps_per_protein = e_meta %>% dplyr::group_by(!!dplyr::sym(emeta_cname)) %>% 
      dplyr::mutate(peps_per_pro = dplyr::n()) %>% dplyr::select(dplyr::any_of(c(emeta_cname, 
                                                                                 "peps_per_pro", emeta_cols)))
    if ("n_peps_used" %in% colnames(results$e_meta)) {
      results$e_meta <- results$e_meta %>% dplyr::select(-n_peps_used)
    }
    results$e_meta <- results$e_meta %>% dplyr::left_join(peps_per_protein, 
                                                          by = emeta_cname, multiple = "all", relationship = "many-to-many") %>% 
      dplyr::left_join(peps_used, by = emeta_cname, multiple = "all", 
                       relationship = "many-to-many") %>% dplyr::relocate(dplyr::any_of(emeta_cols), 
                                                                          .after = n_peps_used) %>% dplyr::distinct(.) %>% 
      data.frame(check.names = FALSE)
  }
  if (!is.null(isoformRes)) {
    emeta_cname <- emeta_cname_iso
  }
  prodata <- as.proData(e_data = results$e_data, f_data = pepData$f_data, 
                        e_meta = results$e_meta, edata_cname = emeta_cname, fdata_cname = fdata_cname, 
                        emeta_cname = emeta_cname, data_scale = get_data_scale(pepData), 
                        is_normalized = is_normalized)
  attr(prodata, "data_info")$data_scale_orig <- get_data_scale_orig(pepData)
  attr(prodata, "group_DF") <- attr(pepData, "group_DF")
  attr(prodata, "pro_quant_info")$method <- method
  return(prodata)
}

combine_fn_sum = function (x) 
{
  if (all(is.na(x))) 
    sum(x)
  else sum(x, na.rm = TRUE)
}






