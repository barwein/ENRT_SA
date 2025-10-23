

# Source + Packages --------------------------------------------------------------

source("src/sensitivity_params.r")

# Indirect effects --------------------------------------------------------

#' Calculate augmented IE_RD for a single pi_vec
ie_aug_point_one_param_ <- function(Y_a,
                                    F_a,
                                    m_a_1,
                                    m_a_0,
                                    pi_vec,
                                    pz) {
  
  # Ensure pi_vec is valid
  if (any(pi_vec >= 1)) {
    warning("pi_vec contains values >= 1. IE estimate will be NA.")
    pi_vec[pi_vec >= 1] <- NA
  }
  
  # Weights for bias correction
  weights <- (1 - pz) / (1 - pi_vec)
  
  I_1 <- (F_a == 1)
  I_0 <- (F_a == 0)
  
  # Augmented estimator components
  term_1 <- (I_1 * (Y_a - m_a_1)) / pz
  term_2 <- (I_0 * (Y_a - m_a_0)) / (1 - pz)
  term_3 <- m_a_1 - m_a_0
  
  inside_sum <- term_1 - term_2 + term_3
  
  # Final estimator for RD
  ie_rd_ <- mean(weights * inside_sum, na.rm = TRUE)
  
  return(c(ie_rd_))
}

#' Calculate augmented IE estimates for a grid of pi_params
#'
#' @param Y_a Vector of alter outcomes.
#' @param F_a Vector of alter observed exposures.
#' @param m_a_1 Vector of cross-fit predictions E[Y|F=1,X].
#' @param m_a_0 Vector of cross-fit predictions E[Y|F=0,X].
#' @param pi_list List of pi_vec sensitivity parameters.
#' @param pz Scalar, Pr(Z=1) among egos.
#'
#' @return A data.table with pi_param, ie_rd, ie_rr
ie_aug_point_grid_ <- function(Y_a,
                               F_a,
                               m_a_1,
                               m_a_0,
                               pi_list,
                               pz) {
  
  if (length(Y_a) != length(F_a) ||
      length(Y_a) != length(m_a_1) ||
      length(Y_a) != length(m_a_0)) {
    stop("Y_a, F_a, m_a_1, and m_a_0 must have the same length (n_a).")
  }
  
  if (pz <= 0 | pz >= 1) {
    stop("Invalid treatment probability 'pz' input.")
  }
  
  ie_grid_ <- vapply(pi_list, function(pi_v) {
    # Ensure pi_v has correct length if it's not scalar
    if (length(pi_v) > 1 && length(pi_v) != length(Y_a)) {
      stop(paste("Heterogeneous pi_vec length", length(pi_v), 
                 "does not match alter sample size", length(Y_a)))
    }
    
    # 'out' is now a scalar
    out <- ie_aug_point_one_param_(
      Y_a = Y_a,
      F_a = F_a,
      m_a_1 = m_a_1,
      m_a_0 = m_a_0,
      pi_vec = pi_v,
      pz = pz
    )
    out
  }, FUN.VALUE = numeric(1))
  
  res_dt <- data.table(
    pi_param = as.numeric(names(pi_list)),
    ie_rd = ie_grid_
  )
  
  return(res_dt)
}


# Direct effects ----------------------------------------------------------

egos_adjust_pi <- function(pi_list,
                           mu_00,
                           mu_01,
                           epsilon = 1e-6){ # Increased precision slightly
  
  # Calculate ratio and handle numerical edge cases robustly
  ratio <- mu_00 / mu_01
  
  # If mu_01 is 0, ratio is Inf or NaN. We want pi_v < Inf, so this is okay.
  # If both are 0, ratio is NaN. We treat this case as unconstrained (Inf).
  ratio[is.nan(ratio)] <- Inf 
  
  # Adjust pi_vec to ensure the denominator of RR_i^e(0) is positive
  pi_list_adjusted <- lapply(pi_list, function(pi_v){
    # Ensure pi_v is broadcastable to the length of ratio for the comparison
    if (length(pi_v) == 1) {
      pi_v_broadcast <- rep(pi_v, length(ratio))
    } else {
      pi_v_broadcast <- pi_v
    }
    
    # Adjust pi where the condition is violated
    pi_v_adjusted <- ifelse(pi_v_broadcast >= ratio, ratio - epsilon, pi_v_broadcast)
    
    # Ensure pi values are non-negative
    pi_v_adjusted[pi_v_adjusted < 0] <- 0
    
    return(pi_v_adjusted)
  })
  
  return(pi_list_adjusted)
}

de_ego_rr_zero <- function(mu_00,
                           mu_01,
                           pi_vec){
  # Aux functions that computes RR_i^e(0) given estimates and pi vec
  rr_i_0 <- ((1 - pi_vec)*mu_01) / (mu_00 - pi_vec*mu_01)
  if(any(rr_i_0 < 0)){
    warning("Some values of RR_i^e(0) are negative. Check inputs.")
  }
  return(rr_i_0)
}

de_point_one_pi_kappa <- function(mu_10,
                                  mu_00,
                                  mu_01,
                                  pi_vec,
                                  kappa_, 
                                  rr_i_0){
  # Estimate bias-corrected DE estimate for a single (pi_param, kappa) combination
  # Weights vectors
  w_vec_10 <- 1 / (1 + pi_vec*(kappa_*rr_i_0 - 1))
  w_vec_00 <- 1 / (1 + pi_vec*(rr_i_0 - 1))
  
  if(any(w_vec_10 < 0) | any(w_vec_00 < 0)){
    warning("Some DE weights are negative. Check inputs.")
  }
  
  # Direct effects on differences scale
  de_rd_ <- mean((mu_10*w_vec_10) - (mu_00*w_vec_00))
  
  # Direct effects on ratio scale
  de_rr_ <- sum(mu_10*w_vec_10) / sum(mu_00*w_vec_00)
  
  return(c(de_rd_, de_rr_))
}

de_grid_one_pi_multi_kappa <- function(mu_10,
                                       mu_00,
                                       mu_01,
                                       pi_vec,
                                       kappa_vec_,
                                       bound_kappa = TRUE){
  # Estimate bias-corrected DE estimates for a single pi_param and a grid of kappa values
  pi_vec_numeric <- unlist(pi_vec)
  rr_i_0_ <- de_ego_rr_zero(mu_00 = mu_00,
                            mu_01 = mu_01,
                            pi_vec = pi_vec_numeric)
  
  # data-driven upper-bound for kappa
  if (bound_kappa){
    upp_bound <- (1 - pi_vec_numeric) / (rr_i_0_*(mu_10 - pi_vec_numeric))
    u_i <- ifelse(mu_10 > pi_vec_numeric,
                  upp_bound,
                  Inf)
    k_bound <- min(u_i[u_i > 0])
    # Bound kappa_vec accordingly
    kappa_vec_[kappa_vec_ >= k_bound] <- k_bound
  }
  
  kappa_list <- as.list(kappa_vec_)
  names(kappa_list) <- round(kappa_vec_, 3)
  
  de_kappa_mat <- t(vapply(kappa_list, function(kappa_){
    out <- de_point_one_pi_kappa(mu_10 = mu_10,
                                 mu_00 = mu_00,
                                 mu_01 = mu_01,
                                 pi_vec = pi_vec_numeric,
                                 kappa_ = kappa_,
                                 rr_i_0 = rr_i_0_)
    out
  }, FUN.VALUE = numeric(2)))
  
  res_dt <- data.table(
    kappa = rownames(de_kappa_mat),
    de_rd = de_kappa_mat[,1],
    de_rr = de_kappa_mat[,2]
  )
  
  return(res_dt)
}


de_grid_multi_pi_kappa <- function(mu_10,
                                   mu_00,
                                   mu_01,
                                   pi_list,
                                   kappa_vec,
                                   bound_kappa = FALSE){
  
  # Estimate bias-corrected DE estimates for a grid of (pi_param, kappa) combinations
  
  if(length(mu_01) != length(mu_00) | length(mu_01) != length(mu_10)){
    stop("Estimates vectors are with different lengths.")
  }
  if (all(lapply(pi_list, function(p){all(p>=0)}) == FALSE)){
    stop("Some pi_vec contains negative values.")
  }
  
  de_list <- lapply(pi_list, function(pi_v){
    de_grid_one_pi_multi_kappa(mu_10 = mu_10,
                               mu_00 = mu_00,
                               mu_01 = mu_01,
                               pi_vec = pi_v,
                               kappa_vec_ = kappa_vec,
                               bound_kappa = bound_kappa)
  })
  res_dt <- rbindlist(de_list, idcol = "pi_param")
  res_dt$pi_param <- as.numeric(res_dt$pi_param)
  res_dt$kappa <- as.numeric(res_dt$kappa)
  
  return(res_dt)
}




