
source("src/sensitivity_params.r")
source("src/bias_adjustment.r")
source("src/outcomes_models.R")
source("src/sa_single_iter.R")
source("src/sa_bootstrap_wrap.R")
source("src/sensitivity_analysis.R")
source("src/pba.R")
source("simulations/dgp.R")
source("simulations/sim_one_iter.R")


# Sim data
# TODO: Think what kind of numerical results to show. Possibilities are:
#       1. Impact of ego-centric sampling on bias in a stylized example
#       2. Bias-correction (and coverage) of adjusted estimators with *true* sensitivity params
# TODO: Run the data analysis!!!


set.seed(515115)
n_e <- 200
n_a <- 400
pz <- 0.5
covar_ <- create_covariates(n_e=n_e, n_a=n_a, seed = 938312)
X_e <- covar_$X_e
X_a <- covar_$X_a
m_e <- 150
m_a <- 200
kappa_ <- 2

rho_hetero_ee <- pi_hetero(X_e = X_e,
                              gamma = -1,
                              m_vec = m_e,
                              dist = "norm", 
                              p = 2,
                              pz = pz, 
                              return_rho_ij = TRUE)$rho$`150`

rho_hetero_ae <- pi_hetero(X_e = X_e,
                           X_a = X_a,
                           m_vec = m_a, 
                           gamma=-1,
                           dist = "norm",
                           p=2,
                           ego_index = rep(1:n_e, each = 2),
                           pz=pz,
                           return_rho_ij = TRUE)$rho$`200`

pop_fixed_data <- create_population(
  X_e = X_e,
  X_a = X_a,
  alters_per_ego = 2,
  # m_e = m_e,
  # m_a = m_a,
  rho_egos = rho_hetero_ee,
  rho_alters = rho_hetero_ae,
  params_po = list(
    b0_e = -2.0,
    b1_e = 2.0, 
    b2_e = 0.5, 
    b3_e = (kappa_-1)*2.0,
    b0_a = -2.0,
    b2_a = 2.0
  ),
  params_covar = list(
    g_e = c(-0.5, -0.3, 0.2),
    g_a = c(-0.4, -0.2, 0.1)
  ),
  seed = 258411
)

# SA params
# m_vec_ee <- seq(10, 200, 10)
# m_vec_ae <- seq(10, 200, 10)
# kappa_vec <- seq(0.9, 2, 0.1)
m_vec_ee <- m_e
m_vec_ae <- m_a
kappa_vec <- kappa_

# EGO-EGO
pi_num_homo_ee <- pi_homo(m_vec = m_vec_ee,
                          n_e = n_e,
                          pz = pz,
                          type = "ego")

pi_num_hetero_ee <- pi_hetero(X_e = pop_fixed_data$X_e,
                          gamma = -1,
                          m_vec = m_vec_ee,
                          dist = "norm", 
                          p = 2,
                          pz = pz)


# ALTER-EGO
pi_num_homo_ae <- pi_homo(m_vec = m_vec_ae,
                          n_e = n_e,
                          n_a =n_a,
                          pz = pz,
                          type = "alter")


pi_num_hetero_ae <- pi_hetero(X_e = pop_fixed_data$X_e,
                              X_a = pop_fixed_data$X_a,
                              m_vec = m_vec_ae, 
                              gamma=-1,
                              dist = "norm",
                              p=2,
                              ego_index = pop_fixed_data$ego_id_a_map,
                              pz=pz)

# PBA priors

prior_ie_poiss <- function() {
  # rpois(1, lambda = m_a)
  rnorm(1, m_a, sqrt(5))
}

pi_ie_args_homo <- list(n_e = n_e, n_a = n_a, type = "alter", pz = 0.5)
pi_ie_args_hetero <- list(X_e=pop_fixed_data$X_e,
                   X_a=pop_fixed_data$X_a,
                   gamma=-1,
                   dist = "norm",
                   p=2,
                   ego_index = pop_fixed_data$ego_id_a_map,
                   pz=pz)

prior_de_poiss_lognorm <- function() {
  list(
    # pi_param = rpois(1, lambda = m_e), # This will be m^e
    pi_param = rnorm(1, m_e, sqrt(5)), # This will be m^e
    # kappa = rlnorm(1, meanlog = 0, sdlog = 0.05)
    kappa = rnorm(1, kappa_, .1)
  )
}

pi_de_args_homo <- list(n_e = n_e, type = "ego", pz = 0.5)
pi_de_args_hetero <- list(X_e = pop_fixed_data$X_e,
                   gamma = -1,
                   dist = "norm",
                   p = 2,
                   pz = pz)

# Check
all_results <- list()
for (i in seq(2)){
  
  iter_res <- sim_estimate_trial(pop_fixed_data = pop_fixed_data,
                     pz = pz,
                     setup_name = "homo",
                     m_e = m_e,
                     m_a = m_a,
                     kappa_ = kappa_, 
                     iter = i,
                     seed = 2548 + i,
                     true_ie = pop_fixed_data$IE_RD,
                     true_de = pop_fixed_data$DE_RD,
                     reg_model_egos = glm,
                     reg_model_alters = glm,
                     formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                     formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                     pi_lists_ego_ego = list("hetero" = pi_num_hetero_ee, "homo" = pi_num_homo_ee),
                     pi_lists_alter_ego = list("hetero" = pi_num_hetero_ae, "homo" = pi_num_homo_ae),
                     kappa_vec = kappa_vec, 
                     n_cores = 8,
                     B_sa = 1e2,
                     B_pba = 1e2, 
                     prior_func_ie = prior_ie_poiss,
                     pi_args_ie_homo =  pi_ie_args_homo,
                     pi_args_ie_hetero = pi_ie_args_hetero,
                     pi_param_name_ie = "m_vec", 
                     prior_func_de = prior_de_poiss_lognorm,
                     pi_args_de_homo = pi_de_args_homo,
                     pi_args_de_hetero = pi_de_args_hetero,
                     pi_param_name_de = "m_vec",
                     family = binomial(link = "logit") # Additional arg for glm
  )
  if (length(all_results) == 0){
    all_results <- iter_res
  } else{
    all_results <- Map(
      f = function(old_dt, new_dt) {
        # rbind the old accumulated data.table with the new one
        rbindlist(list(old_dt, new_dt), use.names = TRUE, fill = TRUE)
      },
      all_results, # The list of accumulated results
      iter_res  # The list from the current run
    )
  }
}
                 

true_ie <- pop_fixed_data$IE_RD
true_de <- pop_fixed_data$DE_RD
true_ie
true_de


