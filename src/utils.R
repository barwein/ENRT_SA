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

# General functions -------------------------------------------------------

dist_euclid <- function(X_1, X_2){
  # Euclidean fast-path: sqrt( ||X_1||^2 + ||X_2||^2 - 2 X_1·X_2 )
  AA <- rowSums(X_1^2) # length n_1
  EE <- rowSums(X_2^2) # length n_2
  D2 <- outer(AA, EE, `+`) - 2 * (X_1 %*% t(X_2))
  DD  <- sqrt(pmax(D2, 0)) # numeric safety
  return(DD)
}

dist_cosine <- function(X_1, X_2){
  # Cosine fast-path: (X_1·X_1) / (||X_2|| ||X_2||)
  AA <- sqrt(rowSums(X_1^2)) # length n_1
  EE <- sqrt(rowSums(X_2^2)) # length n_2
  DD_ <- X_1 %*% t(X_2) / outer(AA, EE, `*`)
  DD <- 1 - DD_
  return(DD)
}

dist_inner <- function(X_1, X_2){
  # Inner product: X_1·X_2
  DD <- X_1 %*% t(X_2)
  return(DD)
}

dist_func <- function(X_1, X_2, dist){
  if(dist == "euclid"){
    return(dist_euclid(X_1, X_2))
  }
  if(dist == "cosine"){
    return(dist_cosine(X_1, X_2))
  }
  if(dist == "inner"){
    return(dist_inner(X_1, X_2))
  }
  stop("Invalid distance metric specified.")
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
  if (rho >= 1 && is.integer(rho)){
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

get_dist_matrix <- function(X_1,
                            X_2,
                            egonet_index,
                            dist = "euclid"){
  
  if (!dist %in% c("euclid","cosine","inner")){
    stop(paste0("Type=", dist, " of covariates distance is not supported."))
  }
  
  X_1 <- as.matrix(X_1)
  X_2 <- as.matrix(X_2)
  
  n_1 <- nrow(X_1)
  n_2 <- nrow(X_2)
  
  # Build the n_1 x n_2 distance matrix
  D <- dist_func(X_1, X_2, dist)
  # Enforce the convention: own-ego distance is 0
  D[cbind(seq_len(n_1), egonet_index)] <- 0
  dimnames(D) <- list(alter = seq_len(n_a), ego = seq_len(n_e))
  
  return(D)
}

hetero_pi_weight_matrix_ <- function(X_1,
                                     X_2,
                                     egonet_index,
                                     dist = "euclid",
                                     gamma = 1){
  # input checks
  X_1 <- as.matrix(X_1)
  X_2 <- as.matrix(X_2)
  
  if (ncol(X_1) != ncol(X_2)){
    stop("X matrices must have the same number of columns.")
  }
  n_1 <- nrow(x_alters); n_2 <- nrow(x_egos)
  
  if (length(egonet_index) != n_1){
    stop("egonet_index must have length nrow(X_1).")
  }
  if (any(is.na(egonet_index)) || any(egonet_index < 1) || any(egonet_index > n_2)){
    stop("egonet_index entries must be integers in 1,...,nrow(X_2).")
  }
  
  if (gamma <= 0){
    stop("gamma value must be larger than zero.")
  }
  
  # Get distances
  D <- get_dist_matrix(X_1,
                       X_2,
                       egonet_index,
                       dist)
  
  # ensure valid interpretation of 'gamma' parameter
  gamma_ <- ifelse(dist %in% c("euclid", "cosine"),
                   -gamma,
                   gamma )
  
  # Convert distances to weights ensuring numerical stability
  L <- gamma_ * D
  # own ego should NOT receive weight -> set to -Inf before exp
  L[cbind(seq_len(nrow(D)), egonet_index)] <- -Inf
  # Normalize weights by subtracting the row max (e.g., each alter max weight)
  row_max <- apply(L, 1, max)
  Lc <- sweep(L, 1, row_max, "-")
  W <- exp(Lc)
  W[!is.finite(W)] <- 0 # turns -Inf to 0
  
  if (dim(W)[1] != n_1 || dim(W)[2] != n_2){
    stop("Error: invalid dimensions of the weights matrix.")
  }
  return(W)
}


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


