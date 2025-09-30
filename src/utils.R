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
  # Build the n_e x n_a distance matrix
  # Note that when X_a=NULL, this is the n_e x n_e distance matrix (with 0 diagonal)
  # and when X_a is given, this is the n_e x n_a distance matrix
  # with no restriction on the diagonal
  return(dist_func(X_e, X_a, dist, p))
}

# Sensitivity parameters (edges probs) ----------------------------------------------------

pi_homo <- function(rho_vec = NULL,
                    m_vec = NULL,
                    n_e,
                    n_a = NULL,
                    type,
                    pz = 0.5
                    ){
  
  # Assert that 'type' in c("ego","alter")
  if (!type %in% c("ego","alter")){
    stop("Invalid prob type; 'type' should be one of 'ego' or 'alter'")
  }
  # Edges prob with homogeneous prob or number of missing edges
  # If rho is a probability (Homogeneous probs case)
  # if(rho <= 1 && rho >= 0 && !is.integer(rho)){
  if(!is.null(rho_vec)){
    rho_ij_vec <- rho_vec
  } else if (!is.null(m_vec)){
  # If rho is an integer (Homogeneous number of missing edges case)
    if(any(m_vec < 0)){
      stop("All values in 'm_vec' must be non-negative.")
    }
    rho_ij_vec <- lapply(m_vec,
                         function(m){
                           ifelse(type == "ego",
                                   m / choose(n_e, 2),
                                   m / (n_a*(n_e-1)))
                       })
  } else{
    stop("Invalid inputs. One of 'rho_vec' or 'm_vec' should be supplied.")
  }
  # Prob unit is not exposed to additional treated ego
  # p_not_expos <- (1 - pz*rho_ij)^(n_e-1)
  # pi_ <- ifelse(type == "ego",
  #               1 - p_not_expos,
  #               pz + (1 - pz)*(1 - p_not_expos))
  # return(pi_)
  pi_list <- lapply(rho_ij_vec, function(rho_ij){
    p_not_expos <- (1 - pz*rho_ij)^(n_e-1)
    if (type == "ego"){
      1 - p_not_expos
    } else {
      pz + (1 - pz)*(1 - p_not_expos)
    }
  })
  if (!is.null(rho_vec)){
    names(pi_list) <- paste0("rho=", rho_vec)
  } else {
    names(pi_list) <- paste0("m=", m_vec)
  }
  return(pi_list)
}


.col_logsumexp <- function(L){
  # Helper function to compute the column-wise LogSumExp of matrix L
  # L is either n_e x n_e or n_e x n_a
  # Were doing Softmax by column in L
  # For one column X: LSE(X) = log(sum(exp(X))) = max(X) + log(sum(exp(X - max(X))))
  # This is more numerically stable
  col_max <- apply(L, 2, max, na.rm=TRUE) # max by column
  Lc <- sweep(L, 2, col_max, "-") # X - max(X)
  Lc[!is.finite(Lc)] <- -Inf 
  S  <- colSums(exp(Lc)) # sum(exp(X - max(X)))
  logsumexp <- col_max + log(S) # max(X) + log(sum(exp(X - max(X))))
  logsumexp[!is.finite(col_max)] <- -Inf
  return(logsumexp)
}

hetero_pi_weight_from_dist_ <- function(D, 
                                        gamma = -1,
                                        ego_index = NULL
                                        ){
  # D: n_e x n_a (or n_e x n_e) distances matrix
  # gamma: scalar "temperature" (negative for distances; positive for inner prod)
  # If self_zero=TRUE and D is square, force P[ii] = 0 (excluded in the softmax)

  L <- gamma * D
  L[is.na(L)] <- -Inf # exclude missing distances
  # For ego-ego distance, make the diagonal -inf -> prob = 0
  if (!is.null(ego_index)){
    L[cbind(ego_index, seq_along(ego_index))] <- -Inf
  } 
  if (is.null(ego_index) && nrow(D) == ncol(D)){
    diag(L) <- -Inf  # exclude self distances in ego-ego case
  } 
  
  # Weights matrix: W_ij = exp(gamma * d_ij)
  # Required for Heterogeneous number of missing edges case
  weights <- exp(L)
  
  # Probs will be computed via Softmax by columns of W
  # We use LogSumExp trick for numerical stability
  # P_ij = exp(gamma*D_ij - LSE(D_j))
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

pi_hetero <- function(X_e,
                      X_a = NULL,
                      m_vec = NULL,
                      gamma = -1,
                      dist = "norm",
                      ego_index = NULL,
                      pz = 0.5,
                      p = 1){
  #' @title Heterogeneous pi values
  #' 
  
  
  type_ <- ifelse(is.null(X_a), "ego", "alter")
  
  if (type_ == "alter" && ncol(X_e) != ncol(X_a)){
    stop("X_e and X_a must have the same number of columns.")
  }
  
  if (!is.null(m_vec) && length(gamma) > 1){
    stop("When 'm_vec' is given, 'gamma' must be a scalar.")
  }
  
  if (!is.null(m_vec) && any(m_vec < 0)){
    stop("All values in 'm_vec' must be non-negative.")
  }
  
  # Distance matrix
  D <- get_dist_matrix(X_e, X_a, dist, p)
  # TODO: in the "alter" type case, need to convert D_ij = - Inf 
  # for each ego i that its alter j is in its ego-network (e(j)=i).
  # Should take an input vector of length n_a such that e_vec[j] = i
  # where i will be value in (1,...,n_e) corresponding to the row of ego e(j).
  
  # get weights and probs matrices (for each gamma value)
  W_P_list <- lapply(gamma, function(g){
    hetero_pi_weight_from_dist_(D, 
                                g, 
                                ego_index)
  })
    
  if(is.null(m_vec)){
    # Heterogeneous probs case. 
    pi_by_gamma <- lapply(W_P_list, function(wp){
      P <- wp$prob
      p_exposed <- 1 - apply(pz*P, 2, function(x) prod(1 - x))
      if (type_ == "ego"){
        p_exposed
      } else {
        pz + (1 - pz)*p_exposed
      }
    })
    names(pi_by_gamma) <- paste0("gamma=", gamma)
    if (length(pi_by_gamma) == 1){
      return(pi_by_gamma[[1]])
    } else {
      return(pi_by_gamma)
    }
  } else{
    # Heterogeneous number of missing edges case
    # Compute all 'pi' probs by all 'm_vec' values at once
    W <- W_P_list[[1]]$weights
    W_sum <- ifelse(type_ == "ego",
                    sum(W[lower.tri(W)]),
                    sum(W))
    
    rho_by_m <- lapply(m_vec, function(m){m*W / W_sum})
    names(rho_by_m) <- paste0("m=", m_vec)
    
    # rho_ij <- missing_num*W / W_sum
    if (type_ == "ego"){
      # rho_ij <- rho_ij - diag(diag(rho_ij)) # force diagonal to 0
      rho_by_m <- lapply(rho_by_m, function(rho_mat){
        rho_mat - diag(diag(rho_mat))
      })
    }
    # Compute pi probs for each m
    pi_by_m <- lapply(rho_by_m, function(rho_ij){
      p_exposed <- 1 - apply(pz*rho_ij, 2, function(x) prod(1 - x))
      if (type_ == "ego"){
        p_exposed
      } else {
        pz + (1 - pz)*p_exposed
      }
    })
    return(pi_by_m)
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

