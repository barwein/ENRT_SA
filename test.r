
source("src/sensitivity_params.r")
source("src/bias_adjustment.r")
source("src/outcomes_models.R")
source("src/sa_single_iter.R")
source("src/sa_bootstrap_wrap.R")
source("src/sensitivity_analysis.R")
source("src/pba.R")

# n_e <- 100
n_e <- 200
# n_a <- 200
n_a <- 400
pz <- 0.5

X_e <- matrix(rbinom(n=n_e*3, size = 1, prob = 0.4), nrow=n_e, ncol = 3)
X_a <- matrix(rbinom(n=n_a*3, size = 1, prob = 0.3), nrow = n_a, ncol = 3)

ego_index <- rep(seq(n_e), each = 2)

Z_e <- rbinom(n_e, 1, pz)
F_a_tilde <- Z_e[ego_index] # observed exposure
F_a <- F_a_tilde 
F_a[F_a==0] <- rbinom(sum(F_a==0), 1, 0.3)

ego_expos <- rbinom(n_e, 1, 0.2)

probs_ego <- 1 / (1 + exp(-(-3 + Z_e*1 + 1.5*ego_expos - 0.5*Z_e*ego_expos + X_e %*% c(0.25, -0.25, -0.3))))
Y_e <- rbinom(n_e, 1, probs_ego)

true_de_rd <- mean(1 / (1 + exp(-(-3 + 1 + 1.5*ego_expos - 0.5*ego_expos + X_e %*% c(0.25, -0.25, -0.3))))) - 
  mean(1 / (1 + exp(-(-3 + 0 + 1.5*ego_expos + X_e %*% c(0.25, -0.25, -0.3)))))
true_de_rr <- sum(1 / (1 + exp(-(-3 + 1 + 1.5*ego_expos - 0.5*ego_expos + X_e %*% c(0.25, -0.25, -0.3))))) / 
  sum(1 / (1 + exp(-(-3 + 0 + 1.5*ego_expos + X_e %*% c(0.25, -0.25, -0.3)))))

probs_alters <- 1 / (1 + exp(-(-3 + F_a*1.5 + X_a %*% c(0.25, -0.25, -0.3))))
Y_a <- rbinom(n_a, 1, probs_alters)

true_ie_rd <- mean(1 / (1 + exp(-(-3 + 1.5 + X_a %*% c(0.25, -0.25, -0.3))))) - 
  mean(1 / (1 + exp(-(-3 + 0 + X_a %*% c(0.25, -0.25, -0.3)))))
true_ie_rr <- sum(1 / (1 + exp(-(-3 + 1.5 + X_a %*% c(0.25, -0.25, -0.3))))) / 
  sum(1 / (1 + exp(-(-3 + 0 + X_a %*% c(0.25, -0.25, -0.3)))))

rho_vec <- seq(1e-3,1e-2, 1e-3)
m_vec_ee <- c(1, seq(5,200,10))
m_vec_ae <- c(1, seq(5,200,10))

# EGO-EGO
pi_homo_ee <- pi_homo(rho_vec = rho_vec,
                      n_e = n_e,
                      pz = pz,
                      type = "ego")

pi_num_homo_ee <- pi_homo(m_vec = m_vec_ee,
                          n_e = n_e,
                          pz = pz,
                          type = "ego")

pi_hetero_ee <- pi_hetero(X_e = X_e,
                          gamma = seq(0.1,3,0.25), 
                          dist = "norm", 
                          p = 1,
                          pz = pz)

pi_num_hetero_ee <- pi_hetero(X_e = X_e,
                          gamma = -1,
                          m_vec = m_vec_ee,
                          dist = "norm", 
                          p = 1,
                          pz = pz)


# ALTER-EGO
pi_homo_ae <- pi_homo(rho_vec = rho_vec,
                      n_e = n_e,
                      n_a = n_a,
                      type = "alter",
                      pz = pz)

pi_num_homo_ae <- pi_homo(m_vec = m_vec_ae,
                          n_e = n_e,
                          n_a =n_a,
                          pz = pz,
                          type = "alter")

pi_hetero_ae <- pi_hetero(X_e = X_e,
                          X_a = X_a,
                          gamma = seq(0.1,3,0.25), 
                          dist = "norm",
                          p=1,
                          ego_index = ego_index,
                          pz = pz)

pi_num_hetero_ae <- pi_hetero(X_e=X_e,
                              X_a=X_a,
                              m_vec = m_vec_ae, 
                              gamma=-1,
                              dist = "norm",
                              p=2,
                              ego_index = ego_index,
                              pz=pz)

# TEST BIAS ADJUSTMENT FOR IE
ui_alters <- runif(n_a)
mu_01 <- 1 / (1 + exp(-(ui_alters + 0.5)))
mu_00 <- 1 / (1 + exp(-(ui_alters)))

esti_mat <- as.matrix(cbind(mu_01, mu_00))

ie_rd_naive <- mean(esti_mat %*% c(1,-1))
ie_rr_naive <- sum(esti_mat %*% c(1,0)) / sum(esti_mat %*% c(0,1))

ie_homo_num <- ie_pi_point_grid_(mu_01 =  mu_01,
                                      mu_00 = mu_00,
                                      pi_list = pi_num_homo_ae,
                                      pz = pz)

ie_hetero_num <- ie_pi_point_grid_(mu_01 =  mu_01,
                                        mu_00 = mu_00,
                                        pi_list = pi_num_hetero_ae,
                                        pz = pz)



# TEST BIAS ADJUSTMENT FOR DE
ui_egos <- runif(n_e)
mu_10 <- 1 / (1 + exp(-(ui_egos + 0.5)))
mu_00 <- 1 / (1 + exp(-(ui_egos + 0.1)))
mu_01 <- 1 / (1 + exp(-(ui_egos + 0.8)))

# kappa_vec <- seq(1, 1.25, 0.025)
kappa_vec <- seq(0.5, 1, 0.1)
# kappa_vec <- seq(1, 1.5, 0.05)

de_hetero_num <- de_grid_multi_pi_kappa(mu_10 = mu_10,
                                        mu_00 = mu_00,
                                        mu_01 = mu_01,
                                        pi_list = pi_num_hetero_ee, 
                                        kappa_vec = kappa_vec)

de_homo_num <- de_grid_multi_pi_kappa(mu_10 = mu_10,
                                        mu_00 = mu_00,
                                        mu_01 = mu_01,
                                        pi_list = pi_num_homo_ee, 
                                        kappa_vec = kappa_vec)



# TEST OUTCOME MODEL

frmal_alter <- as.formula(Y ~ F + X1 + X2 + X3)
preds_alters <- reg_func_alters(
  Y_a = Y_a,
  X_a = X_a,
  F_a = F_a,
  X_e = X_e,
  reg_model = glm,
  formula = frmal_alter,
  family = binomial(link = "logit") # Additional arg for glm
)



preds_egos <- reg_func_egos(
  Y_e = Y_e,
  X_e = X_e,
  Z_e = Z_e,
  reg_model = glm,
  formula = as.formula(Y ~ Z + X1 + X2 + X3),
  family = binomial(link = "logit") # Additional arg for glm
)


# one iter of SA
one_iter_num_hetero <- SA_one_iter(
  Y_e = Y_e,
  Y_a = Y_a,
  X_e = X_e,
  X_a = X_a,
  Z_e = Z_e,
  F_a = F_a,
  reg_model_egos = glm,
  reg_model_alters = glm,
  formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
  formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
  pi_list_ego_ego = pi_num_hetero_ee,
  pi_list_alter_ego = pi_num_hetero_ae,
  kappa_vec = kappa_vec,
  pz = pz, 
  family = binomial(link = "logit") # Additional arg for glm
)


# TEST the bootstrap wrapper

sa_bootstrap_res <- run_sensitivity_bootstrap(
  Y_e = Y_e,
  Y_a = Y_a,
  X_e = X_e,
  X_a = X_a,
  Z_e = Z_e,
  F_a = F_a_tilde,
  ego_id_a = ego_index,
  reg_model_egos = glm,
  reg_model_alters = glm,
  formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
  formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
  pi_list_ego_ego = pi_num_hetero_ee,
  pi_list_alter_ego = pi_num_hetero_ae,
  kappa_vec = kappa_vec,
  B = 1e3,
  n_cores = 8,
  family = binomial(link = "logit") # Additional arg for glm
)


null_sa_bootstrap_res <- run_sensitivity_bootstrap(
  Y_e = Y_e,
  Y_a = Y_a,
  X_e = X_e,
  X_a = X_a,
  Z_e = Z_e,
  F_a = F_a_tilde,
  ego_id_a = ego_index,
  reg_model_egos = glm,
  reg_model_alters = glm,
  formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
  formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
  pi_list_ego_ego = list('0'=0),
  pi_list_alter_ego = list('0' = 0.5),
  kappa_vec = c(0),
  B = 1e3,
  n_cores = 8,
  family = binomial(link = "logit") # Additional arg for glm
)




# # TEST THE FULL SENSITIVITY ANALYSIS FUNCTION
full_sa_res <- enrt_sa(Y_e = Y_e,
                        Y_a = Y_a,
                        X_e = X_e,
                        X_a = X_a,
                        Z_e = Z_e,
                        F_a = F_a_tilde,
                        ego_id_a = ego_index,
                        reg_model_egos = glm,
                        reg_model_alters = glm,
                        formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                        formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                        pi_lists_ego_ego = list("hetero" = pi_num_hetero_ee, "homo" = pi_num_homo_ee),
                        pi_lists_alter_ego = list("hetero" = pi_num_hetero_ae, "homo" = pi_num_homo_ae),
                        kappa_vec = kappa_vec,
                        pz = pz,
                        B = 5e3,
                        n_cores = 8,
                        family = binomial(link = "logit") # Additional arg for glm
  )
full_sa_res$null_results



# TEST the PBA

prior_ie_poiss <- function() {
  rpois(1, lambda = 50)
}

# pi_ie_args <- list(n_e = n_e, n_a = n_a, type = "alter", pz = 0.5)
pi_ie_args <- list(X_e=X_e,
                   X_a=X_a,
                   gamma=-1,
                   dist = "norm",
                   p=2,
                   ego_index = ego_index,
                   pz=pz)

prior_de_poiss_lognorm <- function() {
  list(
    pi_param = rpois(1, lambda = 50), # This will be m^e
    kappa = rlnorm(1, meanlog = 0, sdlog = 0.05)
  )
}

# pi_de_args <- list(n_e = n_e, type = "ego", pz = 0.5)
pi_de_args <- list(X_e = X_e,
                   gamma = -1,
                   dist = "norm", 
                   p = 2,
                   pz = pz)

pba_res <- enrt_pba(Y_e = Y_e, 
                    Y_a = Y_a, 
                    X_e = X_e,
                    X_a = X_a,
                    Z_e = Z_e,
                    F_a = F_a_tilde,
                    ego_id_a = ego_index,
                    reg_model_egos = glm,
                    reg_model_alters = glm,
                    formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                    formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                    B = 1e4,
                    n_cores = 8,
                    prior_func_ie = prior_ie_poiss,
                    # pi_func_ie = pi_homo,
                    pi_func_ie = pi_hetero,
                    pi_args_ie = pi_ie_args,
                    pi_param_name_ie = "m_vec",
                    prior_func_de = prior_de_poiss_lognorm,
                    # pi_func_de = pi_homo,
                    pi_func_de = pi_hetero,
                    pi_args_de = pi_de_args,
                    pi_param_name_de = "m_vec",
                    family = binomial(link = "logit") # Additional arg for glm
    )



# TODO: PBA with pi_hetero computional extensive (compared to homo), but works fine!
# TODO: 1. Test on more realistic DGP that mimics the network sampling design.
# TODO: the DGP can use data from HPTN (covariates, treatment, etc.)
# TODO: just simulate contamination, and outcomes.




