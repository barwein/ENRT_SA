
source("src/bias_adjustment.r")
source("src/outcomes_models.R")
source("src/sensitivity_analysis.R")


# Ensure all necessary source files and libraries are loaded
library(data.table)

#' Run Sensitivity Analysis with Non-Parametric Bootstrap
#'
#' This function wraps the single-iteration sensitivity analysis (`SA_one_iter`)
#' in a non-parametric bootstrap procedure to estimate standard errors and
#' confidence intervals for bias-corrected indirect and direct effects.
#'
#' The bootstrap procedure resamples ego-networks (i.e., egos and all their
#' recruited alters) with replacement.
#'
#' @param Y_e A numeric vector of outcomes for the egos.
#' @param Y_a A numeric vector of outcomes for the alters.
#' @param X_e A numeric matrix of covariates for the egos.
#' @param X_a A numeric matrix of covariates for the alters.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of observed exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's index.
#'        For example, `ego_id_a[j] = i` means alter `j` belongs to ego `i`.
#'        Indices should correspond to the rows in the ego datasets (Y_e, X_e).
#' @param reg_model_egos A function for the egos' outcome regression model.
#' @param reg_model_alters A function for the alters' outcome regression model.
#' @param formula_egos A formula for the egos' regression model.
#' @param formula_alters A formula for the alters' regression model.
#' @param pi_list_ego_ego A list where each element is a vector of exposure
#'        probabilities (`pi_i^e`) for egos, corresponding to a sensitivity
#'        parameter value.
#' @param pi_list_alter_ego A list where each element is a vector of exposure
#'        probabilities (`pi_i^a`) for alters, corresponding to a sensitivity
#'        parameter value.
#' @param kappa_vec A numeric vector of values for the kappa sensitivity parameter.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param B An integer specifying the number of bootstrap iterations.
#' @param probs A numeric vector of probabilities for which to calculate quantiles,
#'        typically `c(0.025, 0.975)`.
#' @param ... Additional arguments passed to the regression model functions.
#'
#' @return A list containing two data.tables:
#' \describe{
#'   \item{IE_results}{A data.table with point estimates, standard errors, and
#'         quantiles for the indirect effects (RD and RR) for each sensitivity
#'         parameter. The first row (`pi_param = 0`) contains the naive estimate.}
#'   \item{DE_results}{A data.table with point estimates, standard errors, and
#'         quantiles for the direct effects (RD and RR) for each combination of
#'         sensitivity parameters. The first row (`pi_param = 0`, `kappa = 0`)
#'         contains the naive estimate.}
#' }
#' @export
run_sensitivity_bootstrap <- function(Y_e,
                                      Y_a,
                                      X_e = NULL,
                                      X_a = NULL,
                                      Z_e,
                                      F_a,
                                      ego_id_a,
                                      reg_model_egos,
                                      reg_model_alters,
                                      formula_egos,
                                      formula_alters,
                                      pi_list_ego_ego,
                                      pi_list_alter_ego,
                                      kappa_vec,
                                      pz = 0.5,
                                      B = 500,
                                      probs = c(0.025, 0.975),
                                      ...) {
  
  # --- Input Validation ---
  if (!is.vector(ego_id_a) || !is.numeric(ego_id_a) || length(ego_id_a) != length(Y_a)) {
    stop("'ego_id_a' must be a numeric vector with the same length as the number of alters.")
  }
  if (max(ego_id_a) > length(Y_e)) {
    stop("An ID in 'ego_id_a' references an ego index that does not exist.")
  }
  
  n_e <- length(Y_e)
  
  # --- 1. Point Estimates (Full Sample) ---
  
  # Run SA on the original data to get the main point estimates
  message("Calculating point estimates on the full dataset...")
  full_sample_results <- SA_one_iter(
    Y_e = Y_e, Y_a = Y_a, X_e = X_e, X_a = X_a, Z_e = Z_e, F_a = F_a,
    reg_model_egos = reg_model_egos, reg_model_alters = reg_model_alters,
    formula_egos = formula_egos, formula_alters = formula_alters,
    pi_list_ego_ego = pi_list_ego_ego, pi_list_alter_ego = pi_list_alter_ego,
    kappa_vec = kappa_vec, pz = pz, 
    # ...
    family = binomial(link = "logit")
  )
  
  IE_point_estimates <- full_sample_results$IE_corrected
  DE_point_estimates <- full_sample_results$DE_corrected
  
  kappa_vec = unique(DE_point_estimates$kappa)
  
  # --- 2. Non-Parametric Bootstrap ---
  
  boot_ie_results_list <- list()
  boot_de_results_list <- list()
  
  message(paste0("Starting ", B, " bootstrap iterations..."))
  
  for (b in 1:B) {
    
    # Resample ego indices with replacement
    boot_ego_indices <- sample(1:n_e, size = n_e, replace = TRUE)
    
    # Find corresponding alters
    # `match` returns the position of the first match
    boot_alter_indices <- which(ego_id_a %in% boot_ego_indices)
    # We need to remap the ego_id's to the new bootstrap indices (1...n_e)
    # This is not strictly needed for the SA function but is good practice
    
    # Find corresponding alters, respecting the "with replacement" sampling.
    # For each ego in the bootstrap sample (including duplicates), find its alters'
    # original indices and combine them.
    boot_alter_indices <- unlist(lapply(boot_ego_indices, function(i) which(ego_id_a == i)))
    
    # Create bootstrap datasets
    Y_e_boot <- Y_e[boot_ego_indices]
    Z_e_boot <- Z_e[boot_ego_indices]
    X_e_boot <- if (!is.null(X_e)) X_e[boot_ego_indices, , drop = FALSE] else NULL
    
    Y_a_boot <- Y_a[boot_alter_indices]
    F_a_boot <- F_a[boot_alter_indices]
    X_a_boot <- if (!is.null(X_a)) X_a[boot_alter_indices, , drop = FALSE] else NULL
    
    # Subset heterogeneous pi vectors if they are not scalar
    pi_list_ego_boot <- lapply(pi_list_ego_ego, function(pi_vec) {
      if (length(pi_vec) > 1) pi_vec[boot_ego_indices] else pi_vec
    })
    
    pi_list_alter_boot <- lapply(pi_list_alter_ego, function(pi_vec) {
      if (length(pi_vec) > 1) pi_vec[boot_alter_indices] else pi_vec
    })
    
    # Run SA for the bootstrap sample
    boot_results <- tryCatch({
      SA_one_iter(
        Y_e = Y_e_boot, Y_a = Y_a_boot, X_e = X_e_boot, X_a = X_a_boot,
        Z_e = Z_e_boot, F_a = F_a_boot,
        reg_model_egos = reg_model_egos, reg_model_alters = reg_model_alters,
        formula_egos = formula_egos, formula_alters = formula_alters,
        pi_list_ego_ego = pi_list_ego_boot, pi_list_alter_ego = pi_list_alter_boot,
        kappa_vec = kappa_vec, bound_kappa = FALSE, pz = pz,
        # ...
        family = binomial(link = "logit")
      )
    }, error = function(e) {
      # In case of convergence issues or other errors in a bootstrap sample
      warning(paste("Bootstrap iteration", b, "failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(boot_results)) {
      boot_ie_results_list[[b]] <- boot_results$IE_corrected
      boot_de_results_list[[b]] <- boot_results$DE_corrected
    }
  }
  
  # --- 3. Aggregate Bootstrap Results ---
  
  message("Aggregating bootstrap results...")
  
  # Combine all bootstrap results into single data.tables
  boot_ie_dt <- rbindlist(boot_ie_results_list)
  boot_de_dt <- rbindlist(boot_de_results_list)

  # Summarize IE results
  IE_summary <- boot_ie_dt[, .(
    ie_rd_se = sd(ie_rd, na.rm = TRUE),
    ie_rr_se = sd(ie_rr, na.rm = TRUE),
    ie_rd_q_low = quantile(ie_rd, probs = probs[1], na.rm = TRUE),
    ie_rd_q_high = quantile(ie_rd, probs = probs[2], na.rm = TRUE),
    ie_rr_q_low = quantile(ie_rr, probs = probs[1], na.rm = TRUE),
    ie_rr_q_high = quantile(ie_rr, probs = probs[2], na.rm = TRUE)
  ), by = pi_param]
  
  # Summarize DE results
  DE_summary <- boot_de_dt[, .(
    de_rd_se = sd(de_rd, na.rm = TRUE),
    de_rr_se = sd(de_rr, na.rm = TRUE),
    de_rd_q_low = quantile(de_rd, probs = probs[1], na.rm = TRUE),
    de_rd_q_high = quantile(de_rd, probs = probs[2], na.rm = TRUE),
    de_rr_q_low = quantile(de_rr, probs = probs[1], na.rm = TRUE),
    de_rr_q_high = quantile(de_rr, probs = probs[2], na.rm = TRUE)
  ), by = .(pi_param, kappa)]

  # --- 4. Combine results ---
  
  # Combine IE 
  IE_results <- merge(IE_point_estimates, IE_summary, by = "pi_param")
  
  # Combine DE 
  DE_results <- merge(DE_point_estimates, DE_summary, by = c("pi_param", "kappa"))

  
  # --- 5. Return Final Results ---
  
  return(list(
    IE_results = IE_results,
    DE_results = DE_results
  ))
}



