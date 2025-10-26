
#' Create a Fixed Population for ENRT Simulation
#'
#' Generates all theoretical, fixed quantities for a simulation, including
#' the true network (A), all potential outcomes (POs), and the
#' true causal estimands (IE, DE). This output does not involve
#' any randomization from treatment assignment.
#'
#' @param X_e A numeric matrix of covariates for the egos (n_e rows).
#' @param X_a A numeric matrix of covariates for the alters (n_a rows).
#' @param alters_per_ego Number of alters recruited by each ego.
#' @param m_e Number of *contaminating* ego-ego edges to add.
#' @param m_a Number of *contaminating* alter-ego edges to add (to egos
#'            other than the alter's own).
#' @param params_po A list of intercept/effect coefficients for the
#'                  potential outcome models.
#' @param params_covar A list of covariate coefficients (gamma) for the
#'                     potential outcome models.
#' @param seed An optional random seed for reproducibility of the
#'             population and POs.
#'
#' @return A list containing:
#'   \item{n_e, n_a, n}{Unit counts.}
#'   \item{X_e, X_a}{The input covariate matrices.}
#'   \item{A}{The true n x n adjacency matrix (including contamination).}
#'   \item{ego_id_a_map}{Vector mapping alters (1...n_a) to ego indices (1...n_e).
#'                       This represents the observed network structure.}
#'   \item{po_egos}{A list of 4 vectors: Y_e_00, Y_e_10, Y_e_01, Y_e_11.}
#'   \item{po_alters}{A list of 2 vectors: Y_a_00, Y_a_01.}
#'   \item{true_estimands}{A list with the true $IE and $DE.}
#'
create_population <- function(X_e,
                              X_a,
                              alters_per_ego = 2,
                              m_e = 50,
                              m_a = 100,
                              params_po = list(
                                b0_e = -2.0, b1_e = 0.7, b2_e = 0.5, b3_e = 0.2,
                                b0_a = -2.5, b2_a = 0.6
                              ),
                              params_covar = list(
                                g_e = c(0.5),
                                g_a = c(0.4)
                              ),
                              seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- 1. Setup Units & Validate Inputs ---
  n_e <- nrow(X_e)
  n_a <- nrow(X_a)
  n_covar <- ncol(X_e)
  
  if (n_a != n_e * alters_per_ego) {
    stop("nrow(X_a) does not match nrow(X_e) * alters_per_ego.")
  }
  if (ncol(X_a) != n_covar) {
    stop("X_e and X_a must have the same number of columns.")
  }
  if (length(params_covar$g_e) != n_covar) {
    stop("Length of params_covar$g_e must equal ncol(X_e).")
  }
  if (length(params_covar$g_a) != n_covar) {
    stop("Length of params_covar$g_a must equal ncol(X_a).")
  }
  
  n <- n_e + n_a # Total units
  ego_indices <- 1:n_e
  alter_indices <- (n_e + 1):n
  ego_id_a_map <- rep(ego_indices, times = alters_per_ego)
  
  # --- 1b. Calculate linear covariate term ---
  lin_pred_X_e <- X_e %*% params_covar$g_e
  lin_pred_X_a <- X_a %*% params_covar$g_a
  
  
  # --- 2. Create True Adjacency Matrix (A) ---
  A <- matrix(0, nrow = n, ncol = n)
  
  # 2a. Observed ego-alter "star" edges (represents \tilde{A})
  for (j in 1:n_a) {
    ego_idx <- ego_id_a_map[j]
    alter_idx <- alter_indices[j]
    A[ego_idx, alter_idx] <- 1
    A[alter_idx, ego_idx] <- 1
  }
  
  # 2b. Contaminating ego-ego edges (m_e)
  if (m_e > 0) {
    possible_ee_pairs <- which(upper.tri(matrix(NA, n_e, n_e)), arr.ind = TRUE)
    m_e <- min(m_e, nrow(possible_ee_pairs))
    sampled_rows <- sample(1:nrow(possible_ee_pairs), m_e)
    ee_edges <- possible_ee_pairs[sampled_rows, , drop = FALSE]
    A[ee_edges] <- 1
    A[ee_edges[, c(2,1), drop = FALSE]] <- 1 # Symmetric
  }
  
  # 2c. Contaminating alter-ego edges (m_a)
  if (m_a > 0) {
    all_ae_pairs <- expand.grid(alter_idx_global = alter_indices,
                                ego_idx_global = ego_indices)
    alter_pos <- all_ae_pairs$alter_idx_global - n_e
    all_ae_pairs$own_ego_idx <- ego_id_a_map[alter_pos]
    valid_ae_pairs <- all_ae_pairs[all_ae_pairs$ego_idx_global != all_ae_pairs$own_ego_idx, ]
    m_a <- min(m_a, nrow(valid_ae_pairs))
    sampled_rows <- sample(1:nrow(valid_ae_pairs), m_a)
    ae_edges <- valid_ae_pairs[sampled_rows, c("alter_idx_global", "ego_idx_global")]
    A[as.matrix(ae_edges)] <- 1
    A[as.matrix(ae_edges[, c(2, 1)])] <- 1 # Symmetric
  }
  
  diag(A) <- 0
  
  
  # --- 3. Generate Fixed Potential Outcomes ---
  
  # 3a. Egos (n_e)
  Y_e_00 <- rbinom(n_e, 1, plogis(params_po$b0_e + lin_pred_X_e))
  Y_e_10 <- rbinom(n_e, 1, plogis(params_po$b0_e + params_po$b1_e + lin_pred_X_e))
  Y_e_01 <- rbinom(n_e, 1, plogis(params_po$b0_e + params_po$b2_e + lin_pred_X_e))
  Y_e_11 <- rbinom(n_e, 1, plogis(params_po$b0_e + params_po$b1_e + 
                                    params_po$b2_e + params_po$b3_e + 
                                    lin_pred_X_e))
  
  # 3b. Alters (n_a)
  Y_a_00 <- rbinom(n_a, 1, plogis(params_po$b0_a + lin_pred_X_a))
  Y_a_01 <- rbinom(n_a, 1, plogis(params_po$b0_a + params_po$b2_a + lin_pred_X_a))
  
  
  # --- 4. Define True Causal Estimands ---
  IE_true <- mean(Y_a_01 - Y_a_00)
  DE_true <- mean(Y_e_10 - Y_e_00)
  
  # --- 5. Return All Fixed Population Quantities ---
  
  return(list(
    n_e = n_e,
    n_a = n_a,
    n = n,
    X_e = X_e,
    X_a = X_a,
    A = A,
    ego_id_a_map = ego_id_a_map,
    po_egos = list(
      Y_e_00 = Y_e_00, Y_e_10 = Y_e_10, Y_e_01 = Y_e_01, Y_e_11 = Y_e_11
    ),
    po_alters = list(
      Y_a_00 = Y_a_00, Y_a_01 = Y_a_01
    ),
    true_estimands = list(
      IE = IE_true,
      DE = DE_true
    )
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
  n <- population$n
  
  
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
  Y_e <- po_egos$Y_e_00 # Start with Z=0, F=0
  
  idx_e_10 <- (Z_e == 1 & F_e_true == 0)
  idx_e_01 <- (Z_e == 0 & F_e_true == 1)
  idx_e_11 <- (Z_e == 1 & F_e_true == 1)
  
  Y_e[idx_e_10] <- po_egos$Y_e_10[idx_e_10]
  Y_e[idx_e_01] <- po_egos$Y_e_01[idx_e_01]
  Y_e[idx_e_11] <- po_egos$Y_e_11[idx_e_11]
  
  # 4b. Alters (Z_i = 0 for all alters)
  Y_a <- po_alters$Y_a_00 # Start with Z=0, F=0
  
  idx_a_01 <- (F_a_true == 1)
  Y_a[idx_a_01] <- po_alters$Y_a_01[idx_a_01]
  
  
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


