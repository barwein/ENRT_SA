# Load necessary libraries
library(data.table)
library(ggplot2)
library(parallel)

# source("src/sa_bootstrap_wrap.R") # Assumes the bootstrap function is in this file

#' End-to-End Sensitivity Analysis for Egocentric-Network Randomized Trials
#'
#' This is the main user-facing function that orchestrates the sensitivity
#' analysis. It runs the bootstrap procedure for a naive (null) scenario and for
#' one or more specifications of the sensitivity parameters, returning the
#' results in both tabular and graphical formats.
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
#'        `pi_list_ego_ego` argument for `run_sensitivity_bootstrap`,
#'        representing a different specification (e.g., "Homogeneous", "Heterogeneous").
#' @param pi_lists_alter_ego A named list of lists, corresponding to `pi_lists_ego_ego`.
#' @param kappa_vec A numeric vector of values for the kappa sensitivity parameter.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param B An integer specifying the number of bootstrap iterations.
#' @param n_cores An integer specifying the number of CPU cores to use for
#'        parallel processing via `parallel::mclapply`. Defaults to 1 (sequential `lapply`).
#'        Note: `mclapply` is not supported on Windows.
#' @param plot A logical value. If `TRUE` (default), graphical displays are
#'        generated and returned.
#' @param verbose A logical value. If `TRUE` (default), progress messages are
#'        printed to the console.
#' @param ... Additional arguments passed to the regression model functions (e.g., `family`).
#'
#' @return A list with six elements:
#' \describe{
#'   \item{null_results}{A list containing two data.tables (`IE` and `DE`)
#'         with the naive estimates.}
#'   \item{sa_results}{A list containing two data.tables (`IE` and `DE`) with
#'         the aggregated results from all sensitivity analysis specifications.}
#'   \item{ie_rd_plot}{A ggplot object for the Indirect Effect (RD), or `NULL`.}
#'   \item{ie_rr_plot}{A ggplot object for the Indirect Effect (RR), or `NULL`.}
#'   \item{de_rd_plot}{A ggplot object for the Direct Effect (RD), or `NULL`.}
#'   \item{de_rr_plot}{A ggplot object for the Direct Effect (RR), or `NULL`.}
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
                    B = 1e3,
                    n_cores = parallel::detectCores() %/% 2,
                    plot = TRUE,
                    verbose = TRUE,
                    ...) {
  
  # --- Argument checks ---
  if (!is.list(pi_lists_ego_ego) || is.null(names(pi_lists_ego_ego))) {
    stop("'pi_lists_ego_ego' must be a named list of lists.")
  }
  if (length(pi_lists_ego_ego) != length(pi_lists_alter_ego)) {
    stop("Ego-ego and alter-ego pi lists must have the same number of specifications.")
  }
  
  # --- 1. Calculate Naive Estimates ---
  if(verbose) message("--- Step 1: Calculating Naive Estimates ---")
  
  pi_list_null_ee <- list('0' = 0)
  pi_list_null_ae <- list('0' = pz) # Naive case for alters is pi = pz
  
  null_bootstrap_res <- run_sensitivity_bootstrap(
    Y_e = Y_e, Y_a = Y_a, 
    X_e = X_e, X_a = X_a, 
    Z_e = Z_e, F_a = F_a,
    ego_id_a = ego_id_a, 
    reg_model_egos = reg_model_egos,
    reg_model_alters = reg_model_alters,
    formula_egos = formula_egos, 
    formula_alters = formula_alters,
    pi_list_ego_ego = pi_list_null_ee,
    pi_list_alter_ego = pi_list_null_ae,
    kappa_vec = c(0), pz = pz, 
    B = B, n_cores = n_cores,
    verbose = FALSE,
    ...
  )
  
  null_results <- list(
    IE = null_bootstrap_res$IE_results,
    DE = null_bootstrap_res$DE_results
  )
  
  # --- 2. Run Sensitivity Analyses for all specifications ---
  if(verbose) message("\n--- Step 2: Running Sensitivity Analyses for all Specifications ---")
  
  sa_results_list <- mapply(
    FUN = function(pi_ee, pi_ae, spec_name) {
      if(verbose) message(paste("Running specification:", spec_name))
      run_sensitivity_bootstrap(
        Y_e = Y_e, Y_a = Y_a,
        X_e = X_e, X_a = X_a,
        Z_e = Z_e, F_a = F_a,
        ego_id_a = ego_id_a, 
        reg_model_egos = reg_model_egos,
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos, 
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_ee,
        pi_list_alter_ego = pi_ae,
        kappa_vec = kappa_vec,
        pz = pz, 
        B = B,
        n_cores = n_cores,
        verbose = FALSE,
        ...
      )
    },
    pi_lists_ego_ego,
    pi_lists_alter_ego,
    names(pi_lists_ego_ego),
    SIMPLIFY = FALSE
  )
  
  # --- 3. Aggregate SA Results ---
  if(verbose) message("\n--- Step 3: Aggregating All Results ---")
  
  ie_results_dt <- rbindlist(lapply(sa_results_list, `[[`, "IE_results"), idcol = "spec")
  de_results_dt <- rbindlist(lapply(sa_results_list, `[[`, "DE_results"), idcol = "spec")
  
  sa_results <- list(IE = ie_results_dt, DE = de_results_dt)
  
  # --- Initialize plot variables to NULL ---
  ie_rd_plot <- ie_rr_plot <- de_rd_plot <- de_rr_plot <- NULL
  
  # --- 4. Generate Plots (if requested) ---
  if (plot) {
    if(verbose) message("--- Step 4: Generating Plots ---")
    # Define position dodge for IE plots
    pd <- position_dodge(1.5)
    
    # --- IE Plots ---
    ie_plot_data_rd <- ie_results_dt[, .(spec, pi_param = as.numeric(pi_param), estimate = ie_rd, q_low = ie_rd_q_low, q_high = ie_rd_q_high)]
    naive_ie_rd <- null_results$IE[, .(est = ie_rd, q_low = ie_rd_q_low, q_high = ie_rd_q_high)]
    subtitle_ie_rd <- paste0(sprintf("Naive Estimate (95%% CI): %.2f [%.2f, %.2f]", naive_ie_rd$est, naive_ie_rd$q_low, naive_ie_rd$q_high),
                            ". Dashed line at 0 indicates no indirect effect.")
    
    ie_rd_plot <- ggplot(ie_plot_data_rd, aes(x = pi_param, y = estimate, color = spec)) +
      geom_point(data = naive_ie_rd, aes(x = 0, y = est), color = "black", size = 3, inherit.aes = FALSE) +
      geom_errorbar(data = naive_ie_rd, aes(x = 0, ymin = q_low, ymax = q_high), color = "black", width = 0.01, linewidth = 0.8, inherit.aes = FALSE) +
      geom_line(aes(group = spec), linewidth = 1, position = pd) +
      geom_point(size = 2.5, position = pd) +
      geom_errorbar(aes(ymin = q_low, ymax = q_high, group = spec), width = 0.01, linewidth = 0.8, position = pd) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey48", alpha=0.6, linewidth = 1) + # Null for RD is 0
      labs(title = "Sensitivity of Indirect Effect (Risk Difference)", subtitle = subtitle_ie_rd, 
           x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = "Estimated IE (RD)", color = "Specification") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
    
    ie_plot_data_rr <- ie_results_dt[, .(spec, pi_param = as.numeric(pi_param), estimate = ie_rr, q_low = ie_rr_q_low, q_high = ie_rr_q_high)]
    naive_ie_rr <- null_results$IE[, .(est = ie_rr, q_low = ie_rr_q_low, q_high = ie_rr_q_high)]
    subtitle_ie_rr <- paste0(sprintf("Naive Estimate (95%% CI): %.2f [%.2f, %.2f]", naive_ie_rr$est, naive_ie_rr$q_low, naive_ie_rr$q_high),
                             ". Dashed line at 1 indicates no indirect effect.")
    
    ie_rr_plot <- ggplot(ie_plot_data_rr, aes(x = pi_param, y = estimate, color = spec)) +
      geom_point(data = naive_ie_rr, aes(x = 0, y = est), color = "black", size = 3, inherit.aes = FALSE) +
      geom_errorbar(data = naive_ie_rr, aes(x = 0, ymin = q_low, ymax = q_high), color = "black", width = 0.01, linewidth = 0.8, inherit.aes = FALSE) +
      geom_line(aes(group = spec), linewidth = 1, position = pd) +
      geom_point(size = 2.5, position = pd) +
      geom_errorbar(aes(ymin = q_low, ymax = q_high, group = spec), width = 0.01, linewidth = 0.8, position = pd) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey48", alpha=0.6, linewidth = 1) + # Null for RR is 1
      labs(title = "Sensitivity of Indirect Effect (Risk Ratio)", subtitle = subtitle_ie_rr, 
           x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = "Estimated IE (RR)", color = "Specification") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
    
    # --- DE Plots ---
    naive_de_rd <- null_results$DE[, .(est = de_rd, q_low = de_rd_q_low, q_high = de_rd_q_high)]
    subtitle_de_rd <- paste0(sprintf("Naive Estimate (95%% CI): %.2f [%.2f, %.2f]", naive_de_rd$est, naive_de_rd$q_low, naive_de_rd$q_high),
                             ". Dashed red line at 0 indicates no direct effect.")
    
    de_plot_data_rd <- de_results_dt[, .(spec, pi_param = as.numeric(pi_param), kappa, estimate = de_rd)]
    de_rd_plot <- ggplot(de_plot_data_rd, aes(x = pi_param, y = kappa, z = estimate)) +
      geom_contour_filled(alpha = 0.8) +
      geom_contour(color = "white", linewidth = 0.2) +
      geom_contour(breaks = 0, color = "red", linetype = "dashed", linewidth = 1.2) +
      facet_wrap(~ spec) +
      scale_fill_viridis_d(direction = -1, option = "plasma") +
      labs(title = "Sensitivity of Direct Effect (Risk Difference)", subtitle = subtitle_de_rd, 
           x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = expression(paste(kappa, " Parameter")), fill = "Estimated DE (RD)") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom", strip.background = element_rect(fill="white"))
    
    naive_de_rr <- null_results$DE[, .(est = de_rr, q_low = de_rr_q_low, q_high = de_rr_q_high)]
    subtitle_de_rr <- paste0(sprintf("Naive Estimate (95%% CI): %.2f [%.2f, %.2f]", naive_de_rr$est, naive_de_rr$q_low, naive_de_rr$q_high),
                             ". Dashed red line at 1 indicates no direct effect.")
    
    de_plot_data_rr <- de_results_dt[, .(spec, pi_param = as.numeric(pi_param), kappa, estimate = de_rr)]
    de_rr_plot <- ggplot(de_plot_data_rr, aes(x = pi_param, y = kappa, z = estimate)) +
      geom_contour_filled(alpha = 0.8) +
      geom_contour(color = "white", linewidth = 0.2) +
      geom_contour(breaks = 1, color = "red", linetype = "dashed", linewidth = 1.2) + # Null for RR is 1
      facet_wrap(~ spec) +
      scale_fill_viridis_d(direction = -1, option = "plasma") +
      labs(title = "Sensitivity of Direct Effect (Risk Ratio)", subtitle = subtitle_de_rr, 
           x = expression(paste(pi, " Parameter (m or ", rho, ")")), y = expression(paste(kappa, " Parameter")), fill = "Estimated DE (RR)") +
      theme_bw(base_size = 14) + theme(legend.position = "bottom", strip.background = element_rect(fill="white"))
  }
  
  # --- 5. Return list of results ---
  if(verbose) message("\n--- Analysis Complete ---")
  
  return(invisible(list(
    null_results = null_results,
    sa_results = sa_results,
    ie_rd_plot = ie_rd_plot,
    ie_rr_plot = ie_rr_plot,
    de_rd_plot = de_rd_plot,
    de_rr_plot = de_rr_plot
  )))
}