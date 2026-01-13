
library(data.table)

#' @title Create Population with Independent Latent Edges and True Estimands
#' @description
#' Generates the true population network A, the observed network A_tilde,
#' the potential outcomes (using specified naming), and calculates the
#' true sample average IE and DE (Risk Difference).
#'
#' @param X_e Matrix of ego covariates (n_e x p).
#' @param X_a Matrix of alter covariates (n_a x p).
#' @param alters_per_ego Integer, the number of alters *recruited* by each ego.
#' @param m_e Integer, the *expected* number of latent ego-ego edges.
#'        Used for Scenario 1 (homogeneous).
#' @param m_a Integer, the *expected* number of latent alter-ego edges.
#'        Used for Scenario 1 (homogeneous).
#' @param rho_egos Optional. An (n_e x n_e) matrix of heterogeneous
#'        ego-ego edge probabilities. Used for Scenario 2.
#' @param rho_alters Optional. An (n_e x n_a) matrix of heterogeneous
#'        alter-ego edge probabilities. Used for Scenario 2.
#' @param params_po List of parameters for potential outcomes, e.g.,
#'        list(b0_e, b1_e, b2_e, b3_e, b0_a, b2_a).
#' @param params_covar List of covariate effect vectors, e.g.,
#'        list(g_e, g_a).
#' @param seed Integer for reproducibility.
#'
#' @return A list containing A, A_tilde, po_egos, po_alters, X_e, X_a,
#'         ego_id_a_map, n_e, n_a, IE_RD, and DE_RD.
#'         
#' @keywords simulation
#'          
create_population <- function(X_e,
                              X_a,
                              alters_per_ego,
                              m_e = NULL,
                              m_a = NULL,
                              rho_egos = NULL,
                              rho_alters = NULL,
                              params_po,
                              params_covar,
                              seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_e <- nrow(X_e)
  n_a <- nrow(X_a)
  n <- n_e + n_a
  
  if (n_a != n_e * alters_per_ego) {
    stop("Mismatch in alter counts. Assumes n_a = n_e * alters_per_ego.")
  }
  
  # --- 1. Create Ego-Alter Mapping and Observed Network (A_tilde) ---
  
  ego_id_a_map <- rep(1:n_e, each = alters_per_ego)
  
  A_tilde <- matrix(0, nrow = n, ncol = n)
  alters_global_idx <- (n_e + 1):n
  
  for (j in 1:n_a) {
    i_ego <- ego_id_a_map[j]
    j_alter_global <- alters_global_idx[j]
    A_tilde[i_ego, j_alter_global] <- 1
    A_tilde[j_alter_global, i_ego] <- 1
  }
  
  A <- A_tilde
  
  # --- 2. Add Latent Ego-Ego Edges to A ---
  
  m_e_sampled <- 0
  
  if (is.null(rho_egos)) {
    if(is.null(m_e)) stop("Must provide m_e for homogeneous case.")
    rho_e_homo <- m_e / choose(n_e, 2)
    if(rho_e_homo > 1) rho_e_homo <- 1
    rho_ee_mat <- matrix(rho_e_homo, n_e, n_e)
  } else {
    rho_ee_mat <- rho_egos
  }
  
  for (i in 1:(n_e - 1)) {
    for (j in (i + 1):n_e) {
      prob <- rho_ee_mat[i, j]
      edge_sampled <- rbinom(1, 1, prob)
      
      A[i, j] <- edge_sampled
      A[j, i] <- edge_sampled
      
      # Increment counter if edge was added
      m_e_sampled <- m_e_sampled + edge_sampled
    }
  }
  
  # --- 3. Add Latent Alter-Ego Edges to A ---
  
  m_a_sampled <- 0
  
  if (is.null(rho_alters)) {
    if(is.null(m_a)) stop("Must provide m_a for homogeneous case.")
    rho_a_homo <- m_a / (n_a * (n_e - 1))
    if(rho_a_homo > 1) rho_a_homo <- 1
    rho_ae_mat <- matrix(rho_a_homo, n_e, n_a)
  } else {
    rho_ae_mat <- rho_alters
  }
  
  for (k in 1:n_e) { 
    for (j in 1:n_a) {
      if (k == ego_id_a_map[j]) {
        next
      }
      j_alter_global <- alters_global_idx[j]
      prob <- rho_ae_mat[k, j]
      edge_sampled <- rbinom(1, 1, prob)
      
      A[k, j_alter_global] <- edge_sampled
      A[j_alter_global, k] <- edge_sampled
      
      m_a_sampled <- m_a_sampled + edge_sampled
    }
  }
  
  # --- 4. Generate Potential Outcomes ---
  
  # Calculate covariate linear predictors
  Xg_e <- if(is.null(params_covar$g_e)) 0 else as.numeric(X_e %*% params_covar$g_e)
  Xg_a <- if(is.null(params_covar$g_a)) 0 else as.numeric(X_a %*% params_covar$g_a)
  
  # Generate POs using baseline parameters and covariate effects
  po_egos <- with(params_po, {
    lp_e_00 <- b0_e + Xg_e
    lp_e_01 <- b0_e + Xg_e + b2_e
    lp_e_10 <- b0_e + Xg_e + b1_e
    lp_e_11 <- b0_e + Xg_e + b1_e + b2_e + b3_e
    ego_epsilon <- rnorm(n_e, mean = 0, sd = 1)
    
    cbind(
      Y_e_00 = lp_e_00 + ego_epsilon,
      Y_e_01 = lp_e_01 + ego_epsilon,
      Y_e_10 = lp_e_10 + ego_epsilon,
      Y_e_11 = lp_e_11 + ego_epsilon
      # Y_e_00 = rbinom(n_e, 1, plogis(lp_e_00)),
      # Y_e_01 = rbinom(n_e, 1, plogis(lp_e_01)),
      # Y_e_10 = rbinom(n_e, 1, plogis(lp_e_10)),
      # Y_e_11 = rbinom(n_e, 1, plogis(lp_e_11))
    )
  })
  
  po_alters <- with(params_po, {
    lp_a_00 <- b0_a + Xg_a
    lp_a_01 <- b0_a + Xg_a + b2_a
    alter_epsilon <- rnorm(n_a, mean = 0, sd = 1)
    
    cbind(
      Y_a_00 = lp_a_00 + alter_epsilon,
      Y_a_01 = lp_a_01 + alter_epsilon
      # Y_a_00 = rbinom(n_a, 1, plogis(lp_a_00)),
      # Y_a_01 = rbinom(n_a, 1, plogis(lp_a_01))
    )
  })
  
  # --- 5. Calculate True Causal Estimands (RD) ---
  
  # IE = E[Y(0, 1) - Y(0, 0)] for alters 
  IE_RD <- mean(po_alters[, "Y_a_01"] - po_alters[, "Y_a_00"])
  
  # DE = E[Y(1, 0) - Y(0, 0)] for egos 
  DE_RD <- mean(po_egos[, "Y_e_10"] - po_egos[, "Y_e_00"])
  
  # --- 6. Return all data ---
  return(list(
    A = A,
    A_tilde = A_tilde,
    po_egos = po_egos,
    po_alters = po_alters,
    X_e = X_e,
    X_a = X_a,
    ego_id_a_map = ego_id_a_map,
    n_e = n_e,
    n_a = n_a,
    IE_RD = IE_RD,
    DE_RD = DE_RD,
    m_e_sampled = m_e_sampled,
    m_a_sampled = m_a_sampled
  ))
}

#' Run One Randomized Trial on a Fixed Population
#'
#' Takes a fixed population (from create_population()) and simulates a
#' single Bernoulli trial (treatment assignment Z) to produce one
#' instance of observed data (Z_e, F_a_obs, Y_e, Y_a).
#'
#' @param population A list object generated by create_population().
#' @param pz The probability of ego treatment, Pr(Z=1).
#' @param seed An optional random seed for this specific trial.
#'
#' @return A list formatted for the enrt_sa() function:
#'   (Y_e, Y_a, X_e, X_a, Z_e, F_a, ego_id_a)
#'
run_trial <- function(population, pz = 0.5, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- 1. Extract Fixed Quantities from Population ---
  A <- population$A
  po_egos <- population$po_egos
  po_alters <- population$po_alters
  ego_id_a_map <- population$ego_id_a_map
  n_e <- population$n_e
  n_a <- population$n_a
  n <- n_e + n_a
  
  
  # --- 2. Run Trial: Assign Treatment (Z) ---
  Z_e <- rbinom(n_e, 1, pz)
  Z <- c(Z_e, rep(0, n_a)) # Full n x 1 treatment vector
  
  
  # --- 3. Determine True (F_i) and Observed (\tilde{F}_i) Exposures ---
  
  # 3a. True Exposure (F_i) - based on true network A
  neighbors_treated <- A %*% Z
  F_true <- as.numeric(neighbors_treated > 0)
  F_e_true <- F_true[1:n_e]
  F_a_true <- F_true[(n_e + 1):n]
  
  # 3b. Observed Exposure (\tilde{F}_a) - based on observed network
  # F_a_obs = 1 iff alter's *own* ego is treated
  F_a_obs <- Z_e[ego_id_a_map]
  
  
  # --- 4. Determine Observed Outcomes (Y_i) ---
  # Based on Consistency: Y_i = Y_i(Z_i, F_i_true)
  
  # 4a. Egos
  Y_e <- po_egos[,"Y_e_00"] # Start with Z=0, F=0
  
  idx_e_10 <- (Z_e == 1 & F_e_true == 0)
  idx_e_01 <- (Z_e == 0 & F_e_true == 1)
  idx_e_11 <- (Z_e == 1 & F_e_true == 1)
  
  Y_e[idx_e_10] <- po_egos[,"Y_e_10"][idx_e_10]
  Y_e[idx_e_01] <- po_egos[,"Y_e_01"][idx_e_01]
  Y_e[idx_e_11] <- po_egos[,"Y_e_11"][idx_e_11]
  
  # 4b. Alters (Z_i = 0 for all alters)
  Y_a <- po_alters[,"Y_a_00"] # Start with Z=0, F=0
  
  idx_a_01 <- (F_a_true == 1)
  Y_a[idx_a_01] <- po_alters[,"Y_a_01"][idx_a_01]
  
  
  # --- 5. Return Data Formatted for enrt_sa ---
  
  return(list(
    Y_e = Y_e,
    Y_a = Y_a,
    X_e = population$X_e,
    X_a = population$X_a,
    Z_e = Z_e,
    F_a = F_a_obs,  # Use *OBSERVED* exposure
    ego_id_a = ego_id_a_map
  ))
}


#' Generate Covariate Matrices for Egos and Alters
#'
create_covariates <- function(n_e,
                              n_a,
                              seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  x_e_ber1 <- rbinom(n_e, 1, 0.6)
  x_e_ber2 <- rbinom(n_e, 1, 0.2)
  x_e_norm <- rnorm(n_e, mean = 0, sd = 2)
  X_e <- cbind(x_e_ber1, x_e_ber2, x_e_norm)
  
  x_a_ber1 <- rbinom(n_a, 1, 0.5)
  x_a_ber2 <- rbinom(n_a, 1, 0.3)
  x_a_norm <- rnorm(n_a, mean = 0, sd = 2)
  X_a <- cbind(x_a_ber1, x_a_ber2, x_a_norm)
  
  return(list(X_e = X_e, X_a = X_a))
}


