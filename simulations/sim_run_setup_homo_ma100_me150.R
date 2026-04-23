
# Source files
source("simulations/dgp.R")
source("simulations/sim_one_iter.R")

# load libraries
library(ENRTsensitivity)
library(data.table)

# Params setup
n_e <- 200
n_a <- 400
pz <- 0.5
covar_ <- create_covariates(n_e=n_e, n_a=n_a, seed = 938312)
X_e <- covar_$X_e
X_a <- covar_$X_a
m_e <- 150
m_a <- 100
kappa_ <- 1.5
n_iter <- 5e3
B_pba <- 5e3
B_sa <- 5e3
N_CORES <- 8

# Results path
results_path <- "simulations/sim_results/"

# rho_ij in heterogeneous case
rho_hetero_ee <- pi_hetero(X_e = X_e,
                           gamma = -1,
                           m_vec = m_e,
                           dist = "norm", 
                           p = 2,
                           pz = pz)$`150`$rho

rho_hetero_ae <- pi_hetero(X_e = X_e,
                           X_a = X_a,
                           m_vec = m_a, 
                           gamma=-1,
                           dist = "norm",
                           p=2,
                           ego_index = rep(1:n_e, each = 2),
                           pz=pz)$`100`$rho

pop_fixed_data <- create_population(
  X_e = X_e,
  X_a = X_a,
  alters_per_ego = 2,
  m_e = m_e,
  m_a = m_a,
  # rho_egos = rho_hetero_ee,
  # rho_alters = rho_hetero_ae,
  params_po = list(
    b0_e = -.5,
    b1_e = 2.0, 
    b2_e = 0.5, 
    b3_e = (kappa_-1)*2.0,
    b0_a = -.5,
    b2_a = 2.0
  ),
  params_covar = list(
    g_e = c(-0.5, -0.3, 0.2),
    g_a = c(-0.4, -0.2, 0.1)
  ),
  param_covar_expos_inter = 0,
  binary_po = FALSE,
  seed = 258411
)

# Verify Ego PO
y_11_e <- pop_fixed_data$po_egos[,"Y_e_11"]
y_10_e <- pop_fixed_data$po_egos[,"Y_e_10"]
y_00_e <- pop_fixed_data$po_egos[,"Y_e_00"]
y_01_e <- pop_fixed_data$po_egos[,"Y_e_01"]

delta_1_e <- y_11_e - y_01_e
delta_0_e <- y_10_e - y_00_e

kappa_val <- mean(delta_1_e) / mean(delta_0_e)


# SA params
m_vec_ee <- m_e
m_vec_ae <- m_a
# kappa_vec <- kappa_
kappa_vec <- kappa_val

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

# Check cov(pi, delta)
pi_used_ee <- pi_num_homo_ee$`150`$pi
print(paste("Cov(pi, delta(1) =", cov(rep(pi_used_ee,times = length(delta_1_e)), delta_1_e),
            "; Cov(pi, delta(0) =", cov(rep(pi_used_ee,times = length(delta_1_e)), delta_0_e)))



# PBA priors

prior_ie_norm <- function(){
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

prior_de_norm <- function(){
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
for (i in seq(n_iter)){
  
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
                                 reg_model_egos = lm,
                                 reg_model_alters = lm,
                                 # reg_model_egos = glm,
                                 # reg_model_alters = glm,
                                 formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                                 formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                                 pi_lists_ego_ego = list("hetero" = pi_num_hetero_ee, "homo" = pi_num_homo_ee),
                                 pi_lists_alter_ego = list("hetero" = pi_num_hetero_ae, "homo" = pi_num_homo_ae),
                                 kappa_vec = kappa_vec, 
                                 n_cores = N_CORES,
                                 B_sa = B_sa,
                                 B_pba = B_pba, 
                                 prior_func_ie = prior_ie_norm,
                                 pi_args_ie_homo =  pi_ie_args_homo,
                                 pi_args_ie_hetero = pi_ie_args_hetero,
                                 pi_param_name_ie = "m_vec", 
                                 prior_func_de = prior_de_norm,
                                 pi_args_de_homo = pi_de_args_homo,
                                 pi_args_de_hetero = pi_de_args_hetero,
                                 pi_param_name_de = "m_vec",
                                 # family = binomial(link = "logit") # Additional arg for glm
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


# Save results
setup_desc <- paste0("pi_homo_","ma", m_a,"_m_e",m_e, "_niter", n_iter,".csv")

for (dt in names(all_results)) {
  res_path <- paste0(results_path,dt,"_",setup_desc)
  write.csv(x = all_results[[dt]],
            res_path)
}




