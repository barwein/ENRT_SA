source("src/bias_adjustment.r")
source("src/outcomes_models.R")
source("src/sa_single_iter.R")
source("src/sensitivity_params.r") 

library(data.table)
library(parallel)

#' Perform Probabilistic Bias Analysis (PBA) for contamination in ENRT
#'
#' This function performs a PBA by integrating over uncertainty in sensitivity
#' parameters and sampling uncertainty (via non-parametric bootstrap). It uses a
#' single Monte Carlo loop to efficiently generate distributions for two types
#' of uncertainty:
#' 1.  **Bias-Only Uncertainty:** Reflects uncertainty from the sensitivity
#'     parameters, holding the dataset fixed.
#' 2.  **Total Uncertainty:** Reflects both sampling and sensitivity
#'     parameter uncertainty.
#'
#' @param Y_e A numeric vector of outcomes for the egos.
#' @param Y_a A numeric vector of outcomes for the alters.
#' @param X_e A numeric matrix of covariates for the egos.
#' @param X_a A numeric matrix of covariates for the alters.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of `observed` exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's index.
#'        Indices should correspond to the rows in the ego datasets (Y_e, X_e).
#' @param reg_model_egos A function for the egos' outcome regression model, e.g., `glm`.
#' @param reg_model_alters A function for the alters' outcome regression model, e.g., `glm`.
#' @param formula_egos A formula for the egos' regression model.
#' @param formula_alters A formula for the alters' regression model.
#' @param B An integer specifying the number of Monte Carlo iterations.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param probs A numeric vector (length=`2`) of probabilities for quantiles.
#' @param n_cores Number of cores for parallel execution.
#' @param verbose A logical indicating whether to print progress messages.
#' @param prior_func_ie A function that returns a single sampled value for the IE
#'        sensitivity parameter (e.g., one 'rho_a' or 'm_a').
#'        Example: `function() runif(1, 0, 0.1)`
#' @param pi_func_ie The R function to calculate the IE pi vector (e.g., `pi_homo`).
#' @param pi_args_ie A list of all arguments for `pi_func_ie`, except for the
#'        one sampled by `prior_func_ie`.
#' @param pi_param_name_ie The string name of the argument in `pi_func_ie` being
#'        sampled. Example: `"rho_vec"` or `"m_vec"`
#' @param prior_func_de A function that returns a list with two named elements:
#'        'pi_param' (sampled DE pi parameter) and 'kappa'.
#'        Example: `function() list(pi_param = runif(1, 0, 0.2), kappa = rlnorm(1, 0, 0.5))`
#' @param pi_func_de The R function to calculate the DE pi vector (e.g., `pi_homo`).
#' @param pi_args_de A list of arguments for `pi_func_de` (excluding the sampled one).
#' @param pi_param_name_de The string name of the pi parameter in `pi_func_de`,
#'       e.g., `"rho_vec"` or `"m_vec"`.
#'
#' @param ... Additional arguments passed to the regression model functions.
#'
#' @return A list containing two data.tables: `IE_results` and `DE_results`.
#'         Each table provides the mean and quantiles for both the
#'         'bias_only' and 'total' uncertainty distributions.
#'
enrt_pba <- function(Y_e,
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
                     # PBA / Bootstrap params
                     B = 1e3,
                     pz = 0.5,
                     probs = c(0.025, 0.975),
                     n_cores = 1,
                     verbose = TRUE,
                     # IE Prior / Pi specifications
                     prior_func_ie,
                     pi_func_ie,
                     pi_args_ie,
                     pi_param_name_ie,
                     # DE Prior / Pi specifications
                     prior_func_de,
                     pi_func_de,
                     pi_args_de,
                     pi_param_name_de,
                     ...) {
  
  # --- Input Validation ---
  n_e <- length(Y_e)
  n_a <- length(Y_a)
  
  # --- Main Iteration Function (Single Loop) ---
  
  pba_single_loop_iteration <- function(iter) {
    # 1. Sample sensitivity parameters *once* for this iteration
    s_a <- prior_func_ie()
    s_de <- prior_func_de()
    
    # 2. Create one bootstrap dataset for this iteration
    boot_ego_indices <- sample(1:n_e, size = n_e, replace = TRUE)
    boot_alter_indices <- unlist(lapply(boot_ego_indices, function(i) which(ego_id_a == i)))
    
    # --- A. Calculate Bias-Only Estimates (Full Data D, Sampled s) ---
    tryCatch({
      # IE Pi List (full data)
      pi_args_alter_full <- pi_args_ie
      pi_args_alter_full[[pi_param_name_ie]] <- s_a
      pi_list_alter_full <- do.call(pi_func_ie, pi_args_alter_full)
      
      # DE Pi List (full data)
      pi_args_ego_full <- pi_args_de
      pi_args_ego_full[[pi_param_name_de]] <- s_de$pi_param
      pi_list_ego_full <- do.call(pi_func_de, pi_args_ego_full)
      
      res_bias_only <- SA_one_iter(
        Y_e = Y_e, Y_a = Y_a,
        X_e = X_e, X_a = X_a,
        Z_e = Z_e, F_a = F_a,
        reg_model_egos = reg_model_egos, 
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos, 
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_list_ego_full,
        pi_list_alter_ego = pi_list_alter_full,
        kappa_vec = s_de$kappa,
        bound_kappa = FALSE, 
        pz = pz, 
        ...
      )
    }, error = function(e) {
      warning(paste("Bias-only calculation failed in iteration", iter, ":", e$message))
      res_bias_only <- NULL
    })
    
    # --- B. Calculate Total Uncertainty Estimates (Bootstrap Data D*, Sampled s) ---
    tryCatch({
      # Bootstrap datasets
      Y_e_boot <- Y_e[boot_ego_indices]
      Z_e_boot <- Z_e[boot_ego_indices]
      X_e_boot <- if (!is.null(X_e)) X_e[boot_ego_indices, , drop = FALSE] else NULL
      Y_a_boot <- Y_a[boot_alter_indices]
      F_a_boot <- F_a[boot_alter_indices]
      X_a_boot <- if (!is.null(X_a)) X_a[boot_alter_indices, , drop = FALSE] else NULL
      
      # IE Pi List (bootstrap data)
      pi_args_alter_boot <- pi_args_ie
      pi_args_alter_boot[[pi_param_name_ie]] <- s_a
      if ("X_e" %in% names(pi_args_alter_boot)) pi_args_alter_boot$X_e <- X_e_boot
      if ("X_a" %in% names(pi_args_alter_boot)) pi_args_alter_boot$X_a <- X_a_boot
      if ("n_a" %in% names(pi_args_alter_boot)) pi_args_alter_boot$n_a <- length(Y_a_boot)
      if ("ego_index" %in% names(pi_args_alter_boot)) pi_args_alter_boot$ego_index <- match(ego_id_a[boot_alter_indices], boot_ego_indices)
      pi_list_alter_boot <- do.call(pi_func_ie, pi_args_alter_boot)
      
      # DE Pi List (bootstrap data)
      pi_args_ego_boot <- pi_args_de
      pi_args_ego_boot[[pi_param_name_de]] <- s_de$pi_param
      if ("X_e" %in% names(pi_args_ego_boot)) pi_args_ego_boot$X_e <- X_e_boot
      pi_list_ego_boot <- do.call(pi_func_de, pi_args_ego_boot)
      
      res_total <- SA_one_iter(
        Y_e = Y_e_boot, Y_a = Y_a_boot,
        X_e = X_e_boot, X_a = X_a_boot,
        Z_e = Z_e_boot, F_a = F_a_boot,
        reg_model_egos = reg_model_egos,
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos,
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_list_ego_boot,
        pi_list_alter_ego = pi_list_alter_boot,
        kappa_vec = s_de$kappa,
        bound_kappa = FALSE, 
        pz = pz,
        ...
      )
    }, error = function(e) {
      warning(paste("Total uncertainty calculation failed in iteration", iter, ":", e$message))
      res_total <- NULL
    })
    
    return(list(bias_only = res_bias_only, total = res_total))
  }
  
  # --- Run and Aggregate ---
  if (verbose) message(paste0("Starting ", B, " PBA iterations on ", n_cores, " core(s)..."))
  
  # Run the single loop
  results_list <- if (n_cores > 1 && .Platform$OS.type != "windows") {
    mclapply(1:B, pba_single_loop_iteration, mc.cores = n_cores)
  } else {
    lapply(1:B, pba_single_loop_iteration)
  }
  
  if (verbose) message("Aggregating results...")
  
  # Separate the results
  bias_only_res <- lapply(results_list, `[[`, "bias_only")
  total_res <- lapply(results_list, `[[`, "total")
  
  # Filter out NULLs from failed iterations
  bias_only_res <- bias_only_res[!sapply(bias_only_res, is.null)]
  total_res <- total_res[!sapply(total_res, is.null)]
  
  # Helper for summarization
  summarize_pba_results <- function(results, type) {
    ie_dt <- rbindlist(lapply(results, `[[`, "IE_corrected"))
    de_dt <- rbindlist(lapply(results, `[[`, "DE_corrected"))
    
    ie_summary <- ie_dt[, .(
      ie_rd_mean = mean(ie_rd, na.rm = TRUE),
      ie_rd_q_low = quantile(ie_rd, probs = probs[1], na.rm = TRUE),
      ie_rd_q_high = quantile(ie_rd, probs = probs[2], na.rm = TRUE),
      ie_rr_mean = mean(ie_rr, na.rm = TRUE),
      ie_rr_q_low = exp(quantile(log(ie_rr), probs = probs[1], na.rm = TRUE)),
      ie_rr_q_high = exp(quantile(log(ie_rr), probs = probs[2], na.rm = TRUE))
    )]
    
    de_summary <- de_dt[, .(
      de_rd_mean = mean(de_rd, na.rm = TRUE),
      de_rd_q_low = quantile(de_rd, probs = probs[1], na.rm = TRUE),
      de_rd_q_high = quantile(de_rd, probs = probs[2], na.rm = TRUE),
      de_rr_mean = mean(de_rr, na.rm = TRUE),
      de_rr_q_low = exp(quantile(log(de_rr), probs = probs[1], na.rm = TRUE)),
      de_rr_q_high = exp(quantile(log(de_rr), probs = probs[2], na.rm = TRUE))
    )]
    
    return(list(
      IE_results = data.table(uncertainty_type = type, ie_summary),
      DE_results = data.table(uncertainty_type = type, de_summary)
    ))
  }
  
  # Summarize both sets of results
  summary_bias_only <- summarize_pba_results(bias_only_res, "bias_only")
  summary_total <- summarize_pba_results(total_res, "total")
  
  # Combine into final data.tables
  final_ie_results <- rbind(summary_bias_only$IE_results, summary_total$IE_results)
  final_de_results <- rbind(summary_bias_only$DE_results, summary_total$DE_results)
  
  return(list(
    IE_results = final_ie_results,
    DE_results = final_de_results
  ))
}


