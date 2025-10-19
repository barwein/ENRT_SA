



source("src/bias_adjustment.r")
source("src/outcomes_models.R")





SA_one_iter <- function(Y_e,
                    Y_a,
                    X_e = NULL,
                    X_a = NULL,
                    Z_e,
                    F_a,
                    reg_model_egos,
                    reg_model_alters,
                    formula_egos,
                    formula_alters,
                    pi_list_ego_ego,
                    pi_list_alter_ego,
                    kappa_vec,
                    bound_kappa = TRUE,
                    pz,
                    ...){
  
  # Check inputs
  if (!is.vector(Y_e) || !is.numeric(Y_e)) stop("Y_e must be a numeric vector.")
  if (!is.vector(Y_a) || !is.numeric(Y_a)) stop("Y_a must be a numeric vector.")
  if (!is.null(X_e) && (!is.matrix(X_e) || !is.numeric(X_e))) stop("X_e must be a numeric matrix.")
  if (!is.null(X_e) && (!is.matrix(X_a) || !is.numeric(X_a))) stop("X_a must be a numeric matrix.")
  if (!is.vector(Z_e) || !is.numeric(Z_e)) stop("Z_e must be a numeric vector.")
  if (!all(Z_e %in% c(0, 1))) stop("Z_e must be a binary vector (0 or 1).")
  if (!is.vector(F_a) || !is.numeric(F_a)) stop("F_a must be a numeric vector.")
  if (!all(F_a %in% c(0, 1))) stop("F_a must be a binary vector (0 or 1).")
  if (length(Y_e) != length(Z_e)) stop("Y_e and Z_e must have the same length (n_e).")
  if (length(Y_a) != length(F_a)) stop("Y_a and F_a must have the same length (n_a).")
  if (!is.function(reg_model_egos)) stop("reg_model_egos must be a function.")
  if (!is.function(reg_model_alters)) stop("reg_model_alters must be a function.")
  if (!inherits(formula_egos, "formula")) stop("formula_egos must be a formula object.")
  if (!inherits(formula_alters, "formula")) stop("formula_alters must be a formula object.")
  if (!is.numeric(kappa_vec)) stop("kappa_vec must be a numeric vector.")
  if (pz <= 0 | pz >= 1) stop("Invalid treatment probability 'pz' input.")
  

  # IE of alters
  # Alters regression model with observed data (naive estimates)
  alters_reg <- reg_func_alters(
                  Y_a = Y_a,
                  X_a = X_a,
                  F_a = F_a,
                  X_e = X_e,
                  reg_model = reg_model_alters,
                  formula = formula_alters,
                  ...
                )
  
  # Save point estimates for alters
  alters_mu01 <- alters_reg$mu_1_alters
  alters_mu00 <- alters_reg$mu_0_alters
  
  # Point estimates for Ego of E[Y(0,1)] (needed for DE adjustment)
  egos_mu01 <- alters_reg$mu_01_egos
  
  # Bias-corrected IE estimates
  # 'IE_corrected' is data.table with columns: pi_param, ie_rd, ie_rr
  # number of rows: length(pi_list_alter_ego)
  IE_corrected <- ie_pi_homo_point_grid_(
                    mu_01 = alters_mu01,
                    mu_00 = alters_mu00,
                    pi_list = pi_list_alter_ego,
                    pz = pz
                    )
  
  # Egos regression model with observed data (naive estimates)
  egos_reg <- reg_func_egos(
                Y_e = Y_e,
                X_e = X_e,
                Z_e = Z_e,
                reg_model = reg_model_egos,
                formula = formula_egos,
               ...
              )
  
  # Save point estimates for egos
  egos_mu10 <- egos_reg$mu_10_e
  egos_mu00 <- egos_reg$mu_00_e
  
  # TODO: add condition: pi_i^e = min(pi_i^e, egos_mu00/egos_mu01 - epsilon)
  
  # Bias-corrected DE estimates (2D grid of (\pi, kappa))
  # 'DE_corrected' is data.table with columns: pi_param, kappa, de_rd, de_rr
  # number of rows: length(pi_list_ego_ego)*length(kappa_vec)
  DE_corrected <- de_grid_multi_pi_kappa(
                      mu_10 = egos_mu10,
                      mu_00 = egos_mu00,
                      mu_01 = egos_mu01,
                      pi_list = pi_list_ego_ego, 
                      kappa_vec = kappa_vec,
                      bound_kappa = bound_kappa
                      )
  
  return(list(
    IE_corrected = IE_corrected,
    DE_corrected = DE_corrected
  ))
}



