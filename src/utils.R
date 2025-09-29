# Write note on this script -- utility function for the sensitivty analysis
# Should be a note documentation

#' @title Sensitivity Analysis Utility Function
#' @description This function performs a sensitivity analysis on a given model.
#' @param model A fitted model object.
#' @param data A data frame containing the data used for the model.
#' @param params A list of parameters to vary in the sensitivity analysis.
#' @return A data frame containing the results of the sensitivity analysis.


# TODO: Think about the general structure of functions:
# 1. Functions that compute estimates of IE and DE given pi (and kappa) values
# 2. Functions that compute pi values given sensitivity parameters
# 2a. Homogeneous probs
# 2b. Homogeneous number of edges
# 2c. Heterogeneous probs of edges
# 2d. Heterogeneous number of edges
# 3. Function in 1+2 should be efficent enough to run with bootstrap iterations
# 4. Wrapper function that either do 
# 4a. Sensitivty analysis for a grid of rho and kappa values (bias-corrected estimates and their associated CI)
# 4b. PBA that takes a prior distribution over rho and kappa and returns a single bias-corrected estimate
# 5. Wrapper function for IE and DE that take as input all the relevant requests and return the relevant result

### We focus only on identification of IE and DE and not the DE bounds!!! ###

# Libraries ---------------------------------------------------------------

library(data.table)
library(proxy)

# General functions -------------------------------------------------------

dist_lp_norm <- function(X_e, X_a = NULL, p){
  # L_p Norm distance
  
  # --- ONE MATRIX (output dim: n_e x n_e) ---
  if (is.null(X_a)) {
    method <- if (is.infinite(p)) "maximum" else
      if (p == 2) "euclidean" else
        if (p == 1) "manhattan" else "minkowski"
    return(as.matrix(stats::dist(X_e,
                                 method = method,
                                 p = if (method == "minkowski") p else 2)))
  }
  
  # --- TWO MATRICES (output dim: n_e x n_a) ---
  lp_fun <- function(a, b, p) {
    if (is.infinite(p)) max(abs(a - b)) else (sum(abs(a - b)^p))^(1 / p)
  }
  return(as.matrix(proxy::dist(X_e, X_a, method = lp_fun, p = p)))
  
  
  
  # if (is.null(X_a)) {X_a <- X_e}
  # EE <- rowSums(X_e^2) # length n_e
  # AA <- rowSums(X_a^2) # length n_a
  # DD_ <- -2 * X_e %*% t(X_a) + outer(EE, AA, `+`)
  # DD_[DD_ < 0] <- 0 # numerical stability
  # return(sqrt(DD_))
  # # return(as.matrix(dist(X)))
}

dist_inner_cosine <- function(X_e, X_a = NULL, dist){
  # Inner product: X_i^TX_j
  # And possible Cosine distance: 1 - (X_i^TX_j) / (||X_i|| ||X_j||)
  center <- FALSE
  if (is.null(X_a)) {
    X_a <- X_e
    center <- TRUE
  }
  inner_prod <- X_e %*% t(X_a)
  # Center the inner product matrix by removing the diagonal
  if (center){
    inner_prod_centered <- inner_prod - diag(diag(inner_prod))
  } else {
    inner_prod_centered <- inner_prod
  }
  # Return the appropriate distance metric
  if (dist == "inner"){
    return(inner_prod_centered)
  }
  if (dist == "cosine"){
    norms_e <- sqrt(rowSums(X_e^2)) # length n_e
    norms_a <- sqrt(rowSums(X_a^2)) # length n_a
    denom <- outer(norms_e, norms_a)
    # To avoid division by zero, set zero norms to a small value
    denom[denom == 0] <- .Machine$double.eps
    return(1 - inner_prod / denom)
    # return(1 - inner_prod_centered / denom)
    # return(inner_prod_centered / outer(norms, norms))
  }
  stop("Invalid distance metric specified.")
}

dist_func <- function(X_e, X_a = NULL, dist = "norm", p = 1){
  stopifnot(is.matrix(X_e))
  if (!is.null(X_a)) stopifnot(is.matrix(X_a), ncol(X_a) == ncol(X_e))
  if (!(is.infinite(p) || (is.numeric(p) && p >= 1)))
    stop("p must be >= 1 or Inf (unweighted).")
  
  if(dist == "norm"){
    return(dist_lp_norm(X_e, X_a, p))
  }
  if(dist %in% c("cosine","inner")){
    return(dist_inner_cosine(X_e, X_a, dist))
  }
  stop("Invalid distance metric specified.")
}


get_dist_matrix <- function(X_e,
                            X_a = NULL,
                            dist = "norm",
                            p = 1){
  
  if (!dist %in% c("norm","cosine","inner")){
    stop(paste0("Type=", dist, " of covariates distance is not supported."))
  }
  n_e <- nrow(X_e)
  if(is.null(X_a)){n_a <- n_e}
  else{n_a <- nrow(X_a)}
  # Build the n_e x n_a distance matrix
  # Note that when X_a=NULL, this is the n_e x n_e distance matrix (with 0 diagonal)
  # and when X_a is given, this is the n_e x n_a distance matrix
  # with no restriction on the diagonal
  D <- dist_func(X_e, X_a, dist, p)
  return(D)
}

# Sensitivity parameter ----------------------------------------------------

pi_a_homo <- function(rho,
                      n_a,
                      n_e,
                      pz){
  
  # edge prob with homogeneous prob or number of missing edges
  rho_ij <- 0
  if(rho <= 1 && rho >= 0 && !is.integer(rho)){
    rho_ij <- rho
  }
  else if (rho >= 1 && is.integer(rho)){
    rho_ij <- rho / (n_a*(n_e-1))
  }
  else{
    stop("Invalid rho value. rho should be either a probability or an integer.")
  }
  # Prob an alter is not exposed to additional treated ego
  p_not_expos <- (1 - pz*rho_ij)^(n_e-1)
  # Prob an alter is exposed
  pi_a <- pz + (1 - pz)*(1 - p_not_expos)
  return(pi_a)
}


.col_logsumexp <- function(L){
  # Helper function to compute the column-wise LogSumExp of matrix L
  # L is either n_e x n_e or n_e x n_a
  col_max <- apply(L, 2, max, na.rm=TRUE)
  Lc <- sweep(L, 2, m, "-")
  Lc[!is.finite(Lc)] <- -Inf                   # handles -Inf - (-Inf) = NaN
  S  <- colSums(exp(Lc))
  logsumexp <- col_max + log(S)
  logsumexp[!is.finite(col_max)] <- -Inf
  return(logsumexp)
}

hetero_pi_weight_from_dist_ <- function(D, 
                                        gamma = 1,
                                        self_zero = FALSE){
  # D: n x m (or n x n) distance matrix
  # gamma: scalar "temperature" (often negative for distances)
  # If self_zero=TRUE and D is square, force P[ii] = 0 (excluded in the softmax)

  L <- gamma * D
  weights <- exp(L)
  
  L[is.na(L)] <- -Inf # exclude missing distances
  # For ego-ego distance, make the diagonal -inf -> prob = 0
  if (self_zero && nrow(D) == ncol(D)) diag(L) <- -Inf  # exclude self
  
  # Compute probs matrix via softmax (by column) and LogSumExp trick 
  lse <- .col_logsumexp(L)
  log_prob <- sweep(L, 2, lse, "-") # log-softmax by column
  log_prob[, is.infinite(lse)] <- -Inf # if all -Inf in col, set all probs to 0
  
  prob <- exp(log_prob)
  prob[!is.finite(prob)] <- 0
  
  return(list(weights = weights,
              prob = prob))
}

# TODO: continue here. I have function that computes the weights W and probs matrix P
# This for one gamma value.
# Can compute for multiple using apply function
# Given a list of these W and P matrices, and point estimate µ..,
# I can compute the IE and DE estimates for each gamma value very fast
# Repeat this process using non-parameteric bootstrap 
# For PBA, we can do similar process, but in the end take an average over gamma values



hetero_pi <- function(X_1,
                      X_2,
                      egonet_index,
                      rho,
                      dist = "euclid",
                      type = "number",
                      gamma = 1,
                      pz = 0.5){
  #' @title Heterogeneous pi values
  #' 
  
  X_1 <- as.matrix(X_1)
  X_2 <- as.matrix(X_2)
  
  # input checks
  hetero_pi_input_check_(X_1, 
                         X_2,
                         egonet_index,
                         rho, 
                         dist, 
                         type, 
                         gamma)
  
  # get weights matrix 
  W <- hetero_pi_weight_matrix_(X_1,
                                X_2,
                                egonet_index,
                                dist,
                                gamma)
  
  # Get normalized weights
  if (type == "number"){
    # number of missing edges
    rho_ve <- rowSums(W) / (ncol(X_2) - 1)
  } else if (type == "prob"){
    # probability of missing edges
    pi_vec <- rowSums(W) / ncol(X_2)
  }
  
  # TODO: think on to design this function.
  # In both 'prob' and 'number' scenarios, I need to compute the weights matrix W once
  # However, in each scenario, the 'm' or 'gamma' grids change the way I compute the 'pi' probs
  # Take this into consideration when modifiying it.
  
}

hetero_pi_input_check_ <- function(X_1,
                                   X_2,
                                   egonet_index,
                                   rho,
                                   dist = "euclid",
                                   type = "number",
                                   gamma = 1){
  # input checks
  
  if (ncol(X_1) != ncol(X_2)){
    stop("X matrices must have the same number of columns.")
  }
  n_1 <- nrow(X_1); n_2 <- nrow(X_2)
  
  if (length(egonet_index) != n_1){
    stop("egonet_index must have length nrow(X_1).")
  }
  if (any(is.na(egonet_index)) || any(egonet_index < 1) || any(egonet_index > n_2)){
    stop("egonet_index entries must be integers in 1,...,nrow(X_2).")
  }
  
  if (!type %in% c("number","prob")){
    stop("type must be either 'number' or 'prob'.")
  }
  
  if(!is.integer(rho) && type == "number"){
    stop("Invalid rho and 'type' combination.")
  }
  if(rho <= 1 && rho >= 0 && type == "number"){
    stop("Invalid rho and 'type' combination.")
  }
  if(is.integer(rho) && type == "prob"){
    stop("Invalid rho and 'type' combination.")
  }
  
  if (gamma <= 0){
    stop("gamma value must be larger than zero.")
  }
}

# Indirect effects --------------------------------------------------------

ie_pi_point_ <- function(esti_mat,
                         pi_vec,
                         pz){
  #' @title Point estimates of indirect effects given pi_vec
  
  # Get weights
  w_vec <- (1 - pz) / (1 - pi_vec)
  
  # Indirect effects on differences scale
  ie_rd_ <- mean((esti_mat %*% c(1,-1)) * w_vec)
  
  # Indirect effects on ratio scale
  n_a <- nrow(esti_mat)
  sum_01 <- sum(esti_mat %*% c(1,0))
  ie_rr_ <- sum_01 / (sum_01 - n_a*ie_rd_)
  
  # Return results
  return(data.table(ie_rd = ie_rd_,
                    ie_rr = ie_rr_,
                    pi_vec = pi_vec))
}

ie_pi_homo_point_grid_ <- function(esti_mat,
                                   pi_list,
                                   pz){
  
  esti_mat <- as.matrix(esti_mat)
  # check if all estimates vec have the same length
  if (ncol(esti_mat) != 2){
    stop("Invalid dimension of estimates matrix.")
  }
  if (pz <= 0 | pz >= 1){
    stop("Invalid treatment probability 'pz' input.")
  }
  
  # TODO: fill this function that return data.table of IE_RD and IE_RR for each pi value
  # should probably convert the 'pi_list' to 'rho_list' or something like that for better
  # clarity and further assistance in the displaying the results.
  
}

# Direct effects ----------------------------------------------------------


