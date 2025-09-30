###
# Functions that compute the sensitivity parameters 
# In our case, it is about the probability of exposure to at least one treated ego
# \pi_i using ego-ego or alter-ego edges probabilities
# Consider all option of homogeneous or heterogeneous probs/number of missing edges
###


# Source ------------------------------------------------------------------

source("src/utils.R")



# Homogeneous case --------------------------------------------------------


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


# Heterogeneous case ------------------------------------------------------

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


