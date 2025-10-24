# Load necessary libraries
library(data.table)
library(ggplot2)
library(parallel)

# Source dependencies
# Make sure the paths to these files are correct
source("src/sa_bootstrap_wrap.R")
source("src/sa_single_iter.R")
source("src/outcomes_models.R")

#' End-to-End Sensitivity Analysis for Egocentric-Network Randomized Trials
#'
#' This is the main user-facing function. It calculates bias-corrected
#' point estimates for IE and DE, along with variance estimates, and constructs
#' confidence intervals using a normal approximation.
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
#' @param pi_lists_ego_ego A named list of lists. Each inner list is a valid
#'        `pi_list_ego_ego` argument.
#' @param pi_lists_alter_ego A named list of lists, corresponding to `pi_lists_ego_ego`.
#' @param kappa_vec A numeric vector of values for the kappa sensitivity parameter.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param bootstrap A logical value. If `TRUE`, variance is estimated via
#'        non-parametric bootstrap. If `FALSE` (default), the empirical
#'        (cluster-robust for IE) variance is used.
#' @param B An integer specifying the number of bootstrap iterations (if `bootstrap = TRUE`).
#' @param n_folds An integer specifying the number of folds for cross-fitting.
#' @param n_cores An integer specifying the number of CPU cores to use for
#'        parallel processing.
#' @param alpha The significance level for confidence intervals (default 0.05).
#' @param plot A logical value. If `TRUE` (default), graphical displays are
#'        generated and returned.
#' @param verbose A logical value. If `TRUE` (default), progress messages are
#'        printed to the console.
#' @param ... Additional arguments passed to the regression model functions (e.g., `family`).
#'
#' @return A list with four elements:
#' \describe{
#'   \item{null_results}{A list containing two data.tables (`IE` and `DE`)
#'         with the naive estimates and CIs.}
#'   \item{sa_results}{A list containing two data.tables (`IE` and `DE`) with
#'         the aggregated SA results and CIs.}
#'   \item{ie_rd_plot}{A ggplot object for the Indirect Effect (RD), or `NULL`.}
#'   \item{de_rd_plot}{A ggplot object for the Direct Effect (RD), or `NULL`.}
#' }
#' @export
enrt_sa <- function(Y_e,
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
                    pi_lists_ego_ego,
                    pi_lists_alter_ego,
                    kappa_vec,
                    pz = 0.5,
                    bootstrap = FALSE,
                    B = 1e3,
                    n_folds = 2,
                    n_cores = 1,
                    alpha = 0.05,
                    plot = TRUE,
                    verbose = TRUE,
                    ...) {
  
  # --- 1. Argument checks ---
  if (!is.list(pi_lists_ego_ego) || is.null(names(pi_lists_ego_ego))) {
    stop("'pi_lists_ego_ego' must be a named list of lists.")
  }
  if (length(pi_lists_ego_ego) != length(pi_lists_alter_ego)) {
    stop("Ego-ego and alter-ego pi lists must have the same number of specifications.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.")
  }
  
  z_val <- qnorm(1 - alpha / 2)
  
  # --- 2. Define Naive (Null) Case ---
  pi_list_null_ee <- list('0' = 0)
  pi_list_null_ae <- list('0' = pz) # Naive case for alters is pi = pz
  
  # Combine naive case with other specifications
  spec_names_to_run <- c("Naive", names(pi_lists_ego_ego))
  pi_ee_to_run <- c(list(Naive = pi_list_null_ee), pi_lists_ego_ego)
  pi_ae_to_run <- c(list(Naive = pi_list_null_ae), pi_lists_alter_ego)
  names(pi_ee_to_run) <- spec_names_to_run
  names(pi_ae_to_run) <- spec_names_to_run
  
  
  # --- 3. Run Sensitivity Analyses ---
  
  # Define the core function to be called by mapply
  run_spec_func <- function(pi_ee, pi_ae, spec_name) {
    if(verbose) message(paste("Running specification:", spec_name))
    
    # Use full kappa_vec for SA specs, but *only c(0)* for Naive DE spec
    k_vec_to_use <- if (spec_name == "Naive") c(0) else kappa_vec
    
    if (bootstrap) {
      if(verbose) message(paste0("...using ", B, " Bootstrap Iterations."))
      run_sensitivity_bootstrap(
        Y_e = Y_e,
        Y_a = Y_a, 
        X_e = X_e, 
        X_a = X_a,
        Z_e = Z_e, 
        F_a = F_a, 
        ego_id_a = ego_id_a,
        reg_model_egos = reg_model_egos, 
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos, 
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_ee,
        pi_list_alter_ego = pi_ae,
        kappa_vec = k_vec_to_use, 
        pz = pz, 
        B = B,
        n_folds = n_folds, 
        n_cores = n_cores,
        verbose = FALSE,
        ...
      )
    } else {
      if(verbose) message("...using Empirical Variance.")
      # Call SA_one_iter directly and get empirical variance
      SA_one_iter(
        Y_e = Y_e, 
        Y_a = Y_a, 
        X_e = X_e, 
        X_a = X_a,
        Z_e = Z_e,
        F_a = F_a,
        ego_id_a = ego_id_a,
        reg_model_egos = reg_model_egos,
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos, 
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_ee,
        pi_list_alter_ego = pi_ae,
        kappa_vec = k_vec_to_use, # Use correct k_vec
        pz = pz, n_folds = n_folds,
        estimate_var = TRUE,
        ...
      )
    }
  }
  
  # Run mapply over all specifications
  all_results_list <- mapply(
    FUN = run_spec_func,
    pi_ee_to_run,
    pi_ae_to_run,
    names(pi_ee_to_run),
    SIMPLIFY = FALSE
  )
  
  # --- 4. Aggregate Results and Calculate CIs ---
  # Aggregate IE
  if (bootstrap) {
    ie_results_dt <- rbindlist(lapply(all_results_list, `[[`, "IE_results"), idcol = "spec")
    ie_results_dt[, var_to_use := ie_rd_boot_var]
  } else {
    ie_results_dt <- rbindlist(lapply(all_results_list, `[[`, "IE_corrected"), idcol = "spec")
    ie_results_dt[, var_to_use := ie_rd_var]
  }
  
  # Aggregate DE
  if (bootstrap) {
    de_results_dt <- rbindlist(lapply(all_results_list, `[[`, "DE_results"), idcol = "spec")
    de_results_dt[, var_to_use := de_rd_boot_var]
  } else {
    de_results_dt <- rbindlist(lapply(all_results_list, `[[`, "DE_corrected"), idcol = "spec")
    de_results_dt[, var_to_use := de_rd_var]
  }
  
  # Calculate CIs for IE
  ie_results_dt[, std_err := sqrt(var_to_use)]
  ie_results_dt[, ci_low := ie_rd - z_val * std_err]
  ie_results_dt[, ci_high := ie_rd + z_val * std_err]
  
  # Calculate CIs for DE
  de_results_dt[, std_err := sqrt(var_to_use)]
  de_results_dt[, ci_low := de_rd - z_val * std_err]
  de_results_dt[, ci_high := de_rd + z_val * std_err]
  
  # Separate Naive from SA results
  null_results_ie <- ie_results_dt[spec == "Naive"]
  null_results_de <- de_results_dt[spec == "Naive"] # Will be 1 row
  sa_results_ie <- ie_results_dt[spec != "Naive"]
  sa_results_de <- de_results_dt[spec != "Naive"]
  
  null_results <- list(IE = null_results_ie, DE = null_results_de)
  sa_results <- list(IE = sa_results_ie, DE = sa_results_de)
  
  # --- Initialize plot variables to NULL ---
  ie_rd_plot <- de_rd_plot <- NULL
  
  # --- 5. Generate Plots (if requested) ---
  if (plot) {
    # --- IE Plot (RD only) ---
    if (nrow(sa_results_ie) > 0) {
      ie_plot_data_rd <- sa_results_ie
      naive_ie_rd <- null_results_ie
      
      subtitle_ie_rd <- paste0(
        sprintf("Naive Estimate (%.0f%% CI): %.2f [%.2f, %.2f]",
                (1 - alpha) * 100, naive_ie_rd$ie_rd, naive_ie_rd$ci_low, naive_ie_rd$ci_high),
        ". Dashed line at 0 indicates no indirect effect."
      )
      
      ie_rd_plot <- ggplot(ie_plot_data_rd, aes(x = as.numeric(pi_param), y = ie_rd, color = spec)) +
        # Naive point
        geom_point(data = naive_ie_rd, aes(x = 0, y = ie_rd), color = "black", size = 3, inherit.aes = FALSE) +
        geom_errorbar(data = naive_ie_rd, aes(x = 0, ymin = ci_low, ymax = ci_high),
                      color = "black", width = 0.01, linewidth = 0.8, inherit.aes = FALSE) +
        # SA lines and CIs
        geom_line(aes(group = spec), linewidth = 1) +
        geom_point(size = 2.5) +
        geom_errorbar(aes(ymin = ci_low, ymax = ci_high, group = spec), width = 0.01, linewidth = 0.8) +
        # Null line
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey48", alpha=0.6, linewidth = 1) +
        labs(title = "Sensitivity of Indirect Effect (Risk Difference)", subtitle = subtitle_ie_rd,
             x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = "Estimated IE (RD)", color = "Specification") +
        theme_bw(base_size = 14) + theme(legend.position = "bottom")
    }
    
    # --- DE Plot (RD only) ---
    if (nrow(sa_results_de) > 0) {
      naive_de_rd <- null_results_de # This is now the single (0,0) point
      subtitle_de_rd <- paste0(
        sprintf("Naive Estimate (pi=0, kappa=0, %.0f%% CI): %.2f [%.2f, %.2f]",
                (1-alpha)*100, naive_de_rd$de_rd, naive_de_rd$ci_low, naive_de_rd$ci_high),
        ". Dashed red line at 0 indicates no direct effect."
      )
      
      de_plot_data_rd <- sa_results_de
      de_rd_plot <- ggplot(de_plot_data_rd, aes(x = as.numeric(pi_param), y = kappa, z = de_rd)) +
        geom_contour_filled(alpha = 0.8) +
        geom_contour(color = "white", linewidth = 0.2) +
        geom_contour(breaks = 0, color = "red", linetype = "dashed", linewidth = 1.2) + # Null for RD is 0
        facet_wrap(~ spec) +
        scale_fill_viridis_d(direction = -1, option = "plasma") +
        labs(title = "Sensitivity of Direct Effect (Risk Difference)", subtitle = subtitle_de_rd,
             x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = expression(paste(kappa, " Parameter")), fill = "Estimated DE (RD)") +
        theme_bw(base_size = 14) + theme(legend.position = "bottom", strip.background = element_rect(fill="white"))
    }
  }
  
  # --- 6. Return list of results ---
  return(invisible(list(
    null_results = null_results,
    sa_results = sa_results,
    ie_rd_plot = ie_rd_plot,
    de_rd_plot = de_rd_plot
  )))
}