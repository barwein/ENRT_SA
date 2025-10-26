
source("src/sensitivity_params.r")
source("src/bias_adjustment.r")
source("src/outcomes_models.R")
source("src/sa_single_iter.R")
source("src/sa_bootstrap_wrap.R")
source("src/sensitivity_analysis.R")
source("src/pba.R")
source("simulations/dgp.R")


# Sim data
# TODO: Think what kind of numerical results to show. Possibilities are:
#       1. Impact of ego-centric sampling on bias in a stylized example
#       2. Bias-correction (and coverage) of adjusted estimators with *true* sensitivity params
# TODO: Run the data analysis!!!


set.seed(515115)
n_e <- 200
n_a <- 400
pz <- 0.5
X_e <- matrix(rbinom(n=n_e*3, size = 1, prob = 0.2), nrow=n_e, ncol = 3)
X_a <- matrix(rbinom(n=n_a*3, size = 1, prob = 0.2), nrow = n_a, ncol = 3)
m_e <- 50
m_a <- 120
kappa <- 1.5

pop_fixed_data <- create_population(
  X_e = X_e,
  X_a = X_a,
  alters_per_ego = 2,
  m_e = m_e,
  m_a = m_a,
  params_po = list(
    b0_e = -2.0,
    b1_e = 1.5, 
    b2_e = 0.5, 
    b3_e = (kappa-1)*1.5,
    b0_a = -2.0,
    b2_a = 2.0
  ),
  params_covar = list(
    g_e = c(0.5, -0.3, 0.2),
    g_a = c(0.4, -0.2, 0.1)
  ),
  seed = 56462
)

# SA params
m_vec_ee <- seq(10, 200, 10)
m_vec_ae <- seq(10, 200, 10)
kappa_vec <- seq(0.9, 2, 0.1)

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
                              ego_index = cur_data$ego_id_a,
                              pz=pz)

# PBA priors

prior_ie_poiss <- function() {
  # rpois(1, lambda = m_a)
  rnorm(1, m_a, sqrt(5))
}

pi_ie_args <- list(n_e = n_e, n_a = n_a, type = "alter", pz = 0.5)
# pi_ie_args <- list(X_e=cur_data$X_e,
#                    X_a=cur_data$X_a,
#                    gamma=-1,
#                    dist = "norm",
#                    p=2,
#                    ego_index = cur_data$ego_id_a,
#                    pz=pz)

prior_de_poiss_lognorm <- function() {
  list(
    # pi_param = rpois(1, lambda = m_e), # This will be m^e
    pi_param = rnorm(1, m_e, sqrt(5)), # This will be m^e
    # kappa = rlnorm(1, meanlog = 0, sdlog = 0.05)
    kappa = rnorm(1, kappa, .2)
  )
}

pi_de_args <- list(n_e = n_e, type = "ego", pz = 0.5)
# pi_de_args <- list(X_e = cur_data$X_e,
#                    gamma = -1,
#                    dist = "norm",
#                    p = 2,
#                    pz = pz)


null_ie_res <- data.table()
null_de_res <- data.table()
ie_res <- data.table()
de_res <- data.table()
ie_res_pba <- data.table()
de_res_pba <- data.table()

for (i in seq(1e3)){
  cur_data <- run_trial(population = pop_fixed_data, 
                        pz = pz,
                        seed = i)
  cur_sa_res <- enrt_sa(Y_e = cur_data$Y_e,
                        Y_a = cur_data$Y_a,
                        X_e = cur_data$X_e,
                        X_a = cur_data$X_a,
                        Z_e = cur_data$Z_e,
                        F_a = cur_data$F_a,
                        ego_id_a = cur_data$ego_id_a,
                        reg_model_egos = glm,
                        # reg_model_egos = lm,
                        # reg_model_egos = NULL,
                        reg_model_alters = glm,
                        # reg_model_alters = lm,
                        formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                        formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                        pi_lists_ego_ego = list("hetero" = pi_num_hetero_ee, "homo" = pi_num_homo_ee),
                        pi_lists_alter_ego = list("hetero" = pi_num_hetero_ae, "homo" = pi_num_homo_ae),
                        kappa_vec = kappa_vec, 
                        verbose = FALSE,
                        pz = pz,
                        n_cores = 8,
                        family = binomial(link = "logit") # Additional arg for glm
  )
  
  null_ie_res <- rbindlist(list(
    null_ie_res,
    cur_sa_res$null_results$IE
  ))
  null_de_res <- rbindlist(list(
    null_de_res,
    cur_sa_res$null_results$DE
  ))
  ie_res <- rbindlist(list(
    ie_res,
    # cur_sa_res$sa_results$IE[spec=="homo" & pi_param == m_a, ]
    cur_sa_res$sa_results$IE[spec=="hetero" & pi_param == m_a, ]
  ))
  de_res <- rbindlist(list(
    de_res,
    # cur_sa_res$sa_results$DE[spec=="homo" & pi_param == m_e & kappa == "1.5", ]
    cur_sa_res$sa_results$DE[spec=="hetero" & pi_param == m_e & kappa == "1.5", ]
  ))
  
  cur_pba_res <- enrt_pba(Y_e = cur_data$Y_e,
                      Y_a = cur_data$Y_a,
                      X_e = cur_data$X_e,
                      X_a = cur_data$X_a,
                      Z_e = cur_data$Z_e,
                      F_a = cur_data$F_a,
                      ego_id_a = cur_data$ego_id_a,
                      reg_model_egos = glm,
                      # reg_model_egos = NULL,
                      reg_model_alters = glm,
                      formula_egos = as.formula(Y ~ Z + X1 + X2 + X3),
                      formula_alters = as.formula(Y ~ F + X1 + X2 + X3),
                      bootstrap = FALSE,
                      # bootstrap = TRUE,
                      verbose = FALSE,
                      B = 1e3,
                      # B = 10,
                      n_cores = 8,
                      prior_func_ie = prior_ie_poiss,
                      pi_func_ie = pi_homo,
                      # pi_func_ie = pi_hetero,
                      pi_args_ie = pi_ie_args,
                      pi_param_name_ie = "m_vec",
                      prior_func_de = prior_de_poiss_lognorm,
                      pi_func_de = pi_homo,
                      # pi_func_de = pi_hetero,
                      pi_args_de = pi_de_args,
                      pi_param_name_de = "m_vec",
                      family = binomial(link = "logit") # Additional arg for glm
  )
  
  ie_res_pba <- rbindlist(list(
    ie_res_pba,
    cur_pba_res$IE_results
  ))
  
  de_res_pba <- rbindlist(list(
    de_res_pba,
    cur_pba_res$DE_results
  ))
  
  if(i %% 10 == 0){
    print(paste("finished iter: ", i))
  }
}

true_ie <- pop_fixed_data$true_estimands$IE
true_de <- pop_fixed_data$true_estimands$DE
true_ie
true_de

ie_res[,.(mean = mean(ie_rd), cover = mean(ci_low <= true_ie & ci_high >= true_ie))]
de_res[,.(mean = mean(de_rd), cover = mean(ci_low <= true_de & ci_high >= true_de))]
null_ie_res[,.(mean = mean(ie_rd), cover = mean(ci_low <= true_ie & ci_high >= true_ie))]
null_de_res[,.(mean = mean(de_rd), cover = mean(ci_low <= true_de & ci_high >= true_de))]

ie_res_pba[,.(mean = mean(ie_rd_mean), 
              cover = mean(ie_rd_q_low <= true_ie & ie_rd_q_high >= true_ie)),
           by = uncertainty_type]

de_res_pba[,.(mean = mean(de_rd_mean), 
              cover = mean(de_rd_q_low <= true_de & de_rd_q_high >= true_de)),
           by = uncertainty_type]


