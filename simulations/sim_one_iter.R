
# Source
source("simulations/dgp.R")
library(ENRTsensitivity)

#' Function that run one iter of the simulation
#' First, it simulate the trial (randomize treaments Z assignment to egos)
#' Then, it estimate IE and DE using bias-corrected and naive estimators
#' We use all setups (homo/iter, bootstrap yes/no, augmented yes/no) for SA and PBA
#' For SA we use the *true* sensitivity parameter values
#' 
sim_estimate_trial <- function(pop_fixed_data,
                               pz,
                               setup_name,
                               m_e,
                               m_a,
                               kappa_,
                               iter,
                               seed,
                               true_ie,
                               true_de,
                               reg_model_egos,
                               reg_model_alters,
                               formula_egos,
                               formula_alters,
                               pi_lists_ego_ego,
                               pi_lists_alter_ego,
                               kappa_vec = kappa_vec,
                               n_cores,
                               B_sa,
                               B_pba,
                               prior_func_ie,
                               pi_args_ie_homo,
                               pi_args_ie_hetero,
                               pi_param_name_ie = "m_vec",
                               prior_func_de,
                               pi_args_de_homo,
                               pi_args_de_hetero,
                               pi_param_name_de = "m_vec",
                               ...){
  # Simulate the trial
  cur_data <- run_trial(population = pop_fixed_data, 
                        pz = pz,
                        seed = seed)
  
  # Run the SA in all setups
  # 1. SA with empirical var; not augmented
  sa_res <- enrt_sa(Y_e = cur_data$Y_e,
                    Y_a = cur_data$Y_a,
                    X_e = cur_data$X_e,
                    X_a = cur_data$X_a,
                    Z_e = cur_data$Z_e,
                    F_a = cur_data$F_a,
                    ego_id_a = cur_data$ego_id_a,
                    reg_model_egos = NULL,
                    reg_model_alters = NULL,
                    formula_egos = NULL,
                    formula_alters = NULL,
                    pi_lists_ego_ego = pi_lists_ego_ego,
                    pi_lists_alter_ego = pi_lists_alter_ego,
                    kappa_vec = kappa_vec, 
                    verbose = FALSE,
                    pz = pz,
                    n_cores = n_cores,
                        ...)
  sa_res_IE <- sa_res$sa_results$IE
  sa_res_DE <- sa_res$sa_results$DE
  sa_res_IE[, `:=`(setup = setup_name,
                   bootstrap = FALSE,
                   augmented = FALSE)]
  sa_res_DE[, `:=`(setup = setup_name,
                   bootstrap = FALSE,
                   augmented = FALSE)]
  
  # null results
  sa_null_IE <- sa_res$null_results$IE
  sa_null_DE <- sa_res$null_results$DE
  sa_null_IE[, `:=`(setup = setup_name,
                    bootstrap = FALSE,
                    augmented = FALSE)]
  sa_null_DE[, `:=`(setup = setup_name,
                    bootstrap = FALSE,
                    augmented = FALSE)]
  
  # 2. SA with empirical var; augmented
  sa_res_aug <- enrt_sa(Y_e = cur_data$Y_e,
                        Y_a = cur_data$Y_a,
                        X_e = cur_data$X_e,
                        X_a = cur_data$X_a,
                        Z_e = cur_data$Z_e,
                        F_a = cur_data$F_a,
                        ego_id_a = cur_data$ego_id_a,
                        reg_model_egos = reg_model_egos,
                        reg_model_alters = reg_model_alters,
                        formula_egos = formula_egos,
                        formula_alters = formula_alters,
                        pi_lists_ego_ego = pi_lists_ego_ego,
                        pi_lists_alter_ego = pi_lists_alter_ego,
                        kappa_vec = kappa_vec, 
                        verbose = FALSE,
                        pz = pz,
                        n_cores = n_cores,
                        ...)
  sa_res_IE_aug <- sa_res_aug$sa_results$IE
  sa_res_DE_aug <- sa_res_aug$sa_results$DE
  sa_res_IE_aug[, `:=`(setup = setup_name,
                       bootstrap = FALSE,
                       augmented = TRUE)]
  sa_res_DE_aug[, `:=`(setup = setup_name,
                       bootstrap = FALSE,
                       augmented = TRUE)]
  # null results
  sa_null_IE_aug <- sa_res_aug$null_results$IE
  sa_null_DE_aug <- sa_res_aug$null_results$DE
  sa_null_IE_aug[, `:=`(setup = setup_name,
                        bootstrap = FALSE,
                        augmented = TRUE)]
  sa_null_DE_aug[, `:=`(setup = setup_name,
                        bootstrap = FALSE,
                        augmented = TRUE)]
  
  # 3. SA with bootstrap var; not augmented
  sa_res_boot <- enrt_sa(Y_e = cur_data$Y_e,
                         Y_a = cur_data$Y_a,
                         X_e = cur_data$X_e,
                         X_a = cur_data$X_a,
                         Z_e = cur_data$Z_e,
                         F_a = cur_data$F_a,
                         ego_id_a = cur_data$ego_id_a,
                         reg_model_egos = NULL,
                         reg_model_alters = NULL,
                         formula_egos = NULL,
                         formula_alters = NULL,
                         pi_lists_ego_ego = pi_lists_ego_ego,
                         pi_lists_alter_ego = pi_lists_alter_ego,
                         kappa_vec = kappa_vec, 
                         verbose = FALSE,
                         pz = pz,
                         n_cores = n_cores,
                         bootstrap = TRUE,
                         B = B_sa,
                         ...)
  sa_res_IE_boot <- sa_res_boot$sa_results$IE
  sa_res_DE_boot <- sa_res_boot$sa_results$DE
  sa_res_IE_boot[, `:=`(setup = setup_name,
                        bootstrap = TRUE,
                        augmented = FALSE)]
  sa_res_IE_boot <- sa_res_IE_boot[,.SD, .SDcols = names(sa_res_IE)]
  sa_res_DE_boot[, `:=`(setup = setup_name,
                        bootstrap = TRUE,
                        augmented = FALSE)]
  sa_res_DE_boot <- sa_res_DE_boot[,.SD, .SDcols = names(sa_res_DE)]
  
  # null results
  sa_null_IE_boot <- sa_res_boot$null_results$IE
  sa_null_DE_boot <- sa_res_boot$null_results$DE
  sa_null_IE_boot[, `:=`(setup = setup_name,
                         bootstrap = TRUE,
                         augmented = FALSE)]
  sa_null_DE_boot[, `:=`(setup = setup_name,
                         bootstrap = TRUE,
                         augmented = FALSE)]
  sa_null_IE_boot <- sa_null_IE_boot[,.SD, .SDcols = names(sa_null_IE)]
  sa_null_DE_boot <- sa_null_DE_boot[,.SD, .SDcols = names(sa_null_DE)]
  
  # 4. SA with bootstrap var; augmented
  sa_res_boot_aug <- enrt_sa(Y_e = cur_data$Y_e,
                             Y_a = cur_data$Y_a,
                             X_e = cur_data$X_e,
                             X_a = cur_data$X_a,
                             Z_e = cur_data$Z_e,
                             F_a = cur_data$F_a,
                             ego_id_a = cur_data$ego_id_a,
                             reg_model_egos = reg_model_egos,
                             reg_model_alters = reg_model_alters,
                             formula_egos = formula_egos,
                             formula_alters = formula_alters,
                             pi_lists_ego_ego = pi_lists_ego_ego,
                             pi_lists_alter_ego = pi_lists_alter_ego,
                             kappa_vec = kappa_vec,
                             verbose = FALSE,
                             pz = pz,
                             n_cores = n_cores,
                             bootstrap = TRUE,
                             B = B_sa,
                             ...)
  sa_res_IE_boot_aug <- sa_res_boot_aug$sa_results$IE
  sa_res_DE_boot_aug <- sa_res_boot_aug$sa_results$DE
  sa_res_IE_boot_aug[, `:=`(setup = setup_name,
                             bootstrap = TRUE,
                             augmented = TRUE)]
  sa_res_DE_boot_aug[, `:=`(setup = setup_name,
                             bootstrap = TRUE,
                             augmented = TRUE)]
  sa_res_IE_boot_aug <- sa_res_IE_boot_aug[,.SD, .SDcols = names(sa_res_IE)]
  sa_res_DE_boot_aug <- sa_res_DE_boot_aug[,.SD, .SDcols = names(sa_res_DE)]
  
  # null results
  sa_null_IE_boot_aug <- sa_res_boot_aug$null_results$IE
  sa_null_DE_boot_aug <- sa_res_boot_aug$null_results$DE
  sa_null_IE_boot_aug[, `:=`(setup = setup_name,
                            bootstrap = TRUE,
                            augmented = TRUE)]
  sa_null_DE_boot_aug[, `:=`(setup = setup_name,
                            bootstrap = TRUE,
                            augmented = TRUE)]
  sa_null_IE_boot_aug <- sa_null_IE_boot_aug[,.SD, .SDcols = names(sa_null_IE)]
  sa_null_DE_boot_aug <- sa_null_DE_boot_aug[,.SD, .SDcols = names(sa_null_DE)]
  
  # Combine SA res
  sa_ie <- rbindlist(list(sa_res_IE,
                             sa_res_IE_aug,
                             sa_res_IE_boot,
                             sa_res_IE_boot_aug))
  sa_de <- rbindlist(list(sa_res_DE,
                             sa_res_DE_aug,
                             sa_res_DE_boot,
                             sa_res_DE_boot_aug))
  sa_null_ie <- rbindlist(list(sa_null_IE,
                                  sa_null_IE_aug,
                                  sa_null_IE_boot,
                                  sa_null_IE_boot_aug))
  sa_null_de <- rbindlist(list(sa_null_DE,
                                  sa_null_DE_aug,
                                  sa_null_DE_boot,
                                  sa_null_DE_boot_aug))
                           
  # return results
  sa_ie <- rbindlist(list(
    sa_ie,
    sa_null_ie
  ))
  
  sa_ie[, `:=`(iter = iter,
               m_a = m_a,
               true_ie = true_ie)]
  
  sa_de <- rbindlist(list(
    sa_de,
    sa_null_de
  ))
  sa_de[, `:=`(iter = iter,
               m_e = m_e,
               kappa_ = kappa_,
               true_de = true_de)]

  return(list(
    sa_ie = sa_ie,
    sa_de = sa_de
  ))
  
}








