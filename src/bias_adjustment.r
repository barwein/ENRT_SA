

# Source + Packages --------------------------------------------------------------

source("src/sensitivity_params.r")

# TODO: add wrapper for both IE and DE functions
# For both, add function that uses non-parameteric bootstrap
# Such that we'll have bias-corrected estimated + 95% CI for each parameters values
# For IE, it is CI for each param
# For DE, it is CI for each (pi_param, kappa) combination.

# TODO: after completing the wrappers, add function that plot the results
# Simple grid for IE, and contour plot for DE
# The plot should include line/point for the naive estimate

# Indirect effects --------------------------------------------------------

ie_point_one_param_ <- function(esti_mat,
                                pi_vec,
                                n_a,
                                pz){
  # Estimate bias-corrected IE estimate for a single pi_param
  
  # Indirect effects on differences scale
  w_vec <- (1 - pz) / (1 - pi_vec)
  ie_rd_ <- mean((esti_mat %*% c(1,-1)) * w_vec)
  
  # IE on ratio scale 
  sum_01 <- sum(esti_mat %*% c(1,0))
  ie_rr_ <- sum_01 / (sum_01 - n_a*ie_rd_)
  
  return(c(ie_rd_, ie_rr_))
}

ie_pi_homo_point_grid_ <- function(mu_01,
                                   mu_00,
                                   pi_list,
                                   pz){
  # Estimate bias-corrected IE estimates for a grid of pi_params
  if(length(mu_01) != length(mu_00)){
    stop("Estimates vectors are with different lengths.")
  }
  esti_mat <- as.matrix(cbind(mu_01,mu_00))
  
  if (pz <= 0 | pz >= 1){
    stop("Invalid treatment probability 'pz' input.")
  }
  
  n_a_ <- nrow(esti_mat)
  
  ie_grid_ <- t(vapply(pi_list, function(pi_v){
    out <- ie_point_one_param_(esti_mat = esti_mat,
                               pi_vec = pi_v,
                               n_a = n_a_,
                               pz = pz)
    out
  }, FUN.VALUE = numeric(2)))
  
  res_dt <- data.table(
    pi_param = as.numeric(rownames(ie_grid_)),
    ie_rd = ie_grid_[,1],
    ie_rr = ie_grid_[,2]
  )
  
  return(res_dt)
}

# Direct effects ----------------------------------------------------------

egos_adjust_pi <- function(pi_list,
                           mu_00,
                           mu_01,
                           epsilon = 1e-3){
  ratio <- mu_00 / mu_01
  # Adjust pi_vec to ensure RR_i^e(0) >= 0
  pi_list_adjusted <- lapply(pi_list, function(pi_v){
    pi_v <- ifelse(pi_v >= ratio, ratio - epsilon, pi_v)
    pi_v[pi_v < 0] <- 0
    return(pi_v)
  })
  return(pi_list_adjusted)
}


de_ego_rr_zero <- function(mu_00,
                           mu_01,
                           pi_vec){
  # Aux functions that computes RR_i^e(0) given estimates and pi vec
  rr_i_0 <- ((1 - pi_vec)*mu_01) / (mu_00 - pi_vec*mu_01)
  # print(paste("range of RR_i^e(0):", paste(round(min(rr_i_0),3), round(max(rr_i_0),3), sep = " to ")))
  if(any(rr_i_0 < 0)){
    print(paste("number of i with pi >= ratio:", sum(pi_vec >= (mu_00 / mu_01))))
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
  
  # print(paste("range of w_vec_10:", paste(round(min(w_vec_10),3), round(max(w_vec_10),3), sep = " to ")))
  
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
                                   bound_kappa = TRUE){
  
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




