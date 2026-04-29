###
# Example code for analysis ENRTs with the proposed senstivity analysis (GSA and PBA)
# The data mimics HPTN 037 settings
###

# Load relevant scripts (similar to those used in the simulation study)
source("simulations/dgp.R")

# Load relevant packages
library(ggplot2)
library(latex2exp)
library(kableExtra)
library(data.table)

library(ENRTsensitivity) # The paper package!


# 1. Define parameters and simulate data
n_e <- 150 # same as HTPN 037
n_a <- 300 # similar to HTPN 037 (HPTN 037 had n_a=263 but we round it here for ease of simulation)
pz <- 0.5 # same as HTPN 037

# Create covariates
x_e_ber1 <- rbinom(n_e, 1, 0.3)
x_e_ber2 <- rbinom(n_e, 1, 0.1)
x_e_norm <- rnorm(n_e, mean = 0, sd = 1)
X_e <- cbind(x_e_ber1, x_e_ber2, x_e_norm)

x_a_ber1 <- rbinom(n_a, 1, 0.4)
x_a_ber2 <- rbinom(n_a, 1, 0.1)
x_a_norm <- rnorm(n_a, mean = 0, sd = 1)
X_a <- cbind(x_a_ber1, x_a_ber2, x_a_norm)


# Contamination parameters (TRUE VALUES)
m_e <- 35 # Number of missing ego-ego edges
m_a <- 105 # Number of missing alter-ego edges
kappa_ <- 1.5 # Proportionality parameter for DE

# Hyperparams for the run
N_CORES <- 8
B <- 1e4 # number of PBA iterations

# Simulate contamination (heterogeneous)

# rho_ij in heterogeneous case
rho_hetero_ee <- pi_hetero(X_e = X_e,
                           gamma = -1,
                           m_vec = m_e,
                           dist = "norm", 
                           p = 2,
                           pz = pz)$`35`$rho

rho_hetero_ae <- pi_hetero(X_e = X_e,
                           X_a = X_a,
                           m_vec = m_a, 
                           gamma=-1,
                           dist = "norm",
                           p=2,
                           ego_index = rep(1:n_e, each = 2),
                           pz=pz)$`105`$rho

# Sample potential outcomes (binary) -- similar DGP to simulations
pop_fixed_data <- create_population(
  X_e = X_e,
  X_a = X_a,
  alters_per_ego = 2, # two alters per ego-network
  rho_egos = rho_hetero_ee,
  rho_alters = rho_hetero_ae,
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
  binary_po = TRUE,
  seed = 42
)

# Get observed outcomes
obs_data <- run_trial(population = pop_fixed_data,
                      pz = pz,
                      seed = 421)


# Verify EGOS potential outcomes
y_11_e <- pop_fixed_data$po_egos[,"Y_e_11"]
y_10_e <- pop_fixed_data$po_egos[,"Y_e_10"]
y_00_e <- pop_fixed_data$po_egos[,"Y_e_00"]
y_01_e <- pop_fixed_data$po_egos[,"Y_e_01"]

delta_1_e <- y_11_e - y_01_e
delta_0_e <- y_10_e - y_00_e

true_kappa <- mean(delta_1_e) / mean(delta_0_e)


# 2. Define sensitivity parameter range for GSA 

m_e_vec <- seq(10, 150, 10)
m_a_vec <- seq(10, 300, 10)
kappa_vec <- seq(1,2, by = 0.05)

# 3. Specify exposure probs (homo + hetero)

# Alters
pi_list_ae_homo <- pi_homo(m_vec = m_a_vec, 
                           n_e = n_e,
                           n_a = n_a,
                           type = "alter",
                           pz = pz)

pi_list_ae_hetero <- pi_hetero(X_e = obs_data$X_e, 
                               X_a = obs_data$X_a,
                               m_vec = m_a_vec,
                               gamma = 1,
                               dist = "norm",
                               p = 2,
                               ego_index = pop_fixed_data$ego_id_a_map,
                               pz = pz)

pi_list_ae <- list(
  Homo = pi_list_ae_homo,
  Hetero = pi_list_ae_hetero
)

# Egos
pi_list_ee_homo <- pi_homo(m_vec = m_e_vec,
                           n_e = n_e,
                           type = "ego",
                           pz = pz)

pi_list_ee_hetero <- pi_hetero(X_e = obs_data$X_e, 
                               m_vec = m_e_vec,
                               gamma = 1,
                               dist = "norm",
                               p = 2, 
                               pz = pz)


pi_list_ee <- list(
  Homo = pi_list_ee_homo,
  Hetero = pi_list_ee_hetero
)

# 4. Run GSA (with augmented outcome model -- logistic regression)

set.seed(42)
gsa_results <- enrt_sa(
    Y_e = obs_data$Y_e,
    Y_a = obs_data$Y_a,
    X_e = obs_data$X_e,
    X_a = obs_data$X_a,
    Z_e = obs_data$Z_e,
    F_a = obs_data$F_a,
    ego_id_a = obs_data$ego_id_a,
    augmented = TRUE,
    reg_model_egos = glm,
    reg_model_alters = glm,
    formula_egos = as.formula(Y ~ Z + .),
    formula_alters = as.formula(Y ~ F + .),
    pi_lists_ego_ego = pi_list_ee,
    pi_lists_alter_ego = pi_list_ae,
    kappa_vec = kappa_vec,
    pz = pz,
    n_cores = N_CORES,
    n_folds = 2,
    family = binomial(link = "logit") # Passed to glm()
  )

gsa_results$ie_rd_plot
ggsave("ENRT_example_code/gsa_ie_results.png",
  gsa_results$ie_rd_plot,
  width = 10, height = 6, dpi = 150)

gsa_results$de_rd_plot
ggsave("ENRT_example_code/gsa_de_results.png",
  gsa_results$de_rd_plot,
  width = 10, height = 6, dpi = 150)

# 5. Define PBA prior functions

# Alters
prior_ie_uniform <- function() {
  sample(seq(1, max(m_a_vec), by=1), 1)
}

prior_ie_poisson <- function() {
  rpois(1, lambda = max(m_a_vec)/2)
}

# Egos
min_kappa <- min(kappa_vec)
max_kappa <- max(kappa_vec)
mean_kappa <- mean(kappa_vec)

prior_de_uniform <- function() {
  list(
    pi_param = sample(seq(1, max(m_e_vec), by=1), 1),
    kappa = runif(1, min = min_kappa, max = max_kappa)
  )
}

prior_de_poisson_uniform <- function() {
  list(
    pi_param = rpois(1, lambda = max(m_e_vec)/2),
    kappa = runif(1, min = min_kappa, max = max_kappa)
  )
}


# Args
pi_ie_homo_args <- list(n_e = n_e,
                        n_a = n_a, 
                        type = "alter", 
                        pz = pz)

pi_ie_hetero_args <- list(X_e=X_e,
                          X_a=X_a,
                          gamma=-1,
                          dist = "norm",
                          p=2,
                          ego_index = obs_data$ego_id_a,
                          pz=pz)

pi_de_homo_args <- list(n_e = n_e, type = "ego", pz = 0.5)

pi_de_hetero_args <- list(X_e = X_e,
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          pz = pz)


# 6. Run PBA


set.seed(142)
pba_homo_unif <-  enrt_pba(Y_e = obs_data$Y_e, 
                               Y_a = obs_data$Y_a, 
                               X_e = obs_data$X_e,
                               X_a = obs_data$X_a,
                               Z_e = obs_data$Z_e,
                               F_a = obs_data$F_a,
                               ego_id_a = obs_data$ego_id_a,
                               augmented = TRUE,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               B = B,
                               n_cores = N_CORES,
                               n_folds = 2,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_uniform,
                               pi_func_ie = pi_homo,
                               pi_args_ie = pi_ie_homo_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_uniform,
                               pi_func_de = pi_homo,
                               pi_args_de = pi_de_homo_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_homo_unif$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "homo"
)]

pba_homo_unif$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "homo"
)]



set.seed(342)
pba_homo_pois <-  enrt_pba(Y_e = obs_data$Y_e, 
                               Y_a = obs_data$Y_a, 
                               X_e = obs_data$X_e,
                               X_a = obs_data$X_a,
                               Z_e = obs_data$Z_e,
                               F_a = obs_data$F_a,
                               ego_id_a = obs_data$ego_id_a,
                               augmented = TRUE,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               B = B,
                               n_cores = N_CORES,
                               n_folds = 2,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_poisson,
                               pi_func_ie = pi_homo,
                               pi_args_ie = pi_ie_homo_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_poisson_uniform,
                               pi_func_de = pi_homo,
                               pi_args_de = pi_de_homo_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_homo_pois$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson",
  model = "homo"
)]

pba_homo_pois$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson + Uniform",
  model = "homo"
)]


set.seed(1425)
pba_hetero_unif <-  enrt_pba(Y_e = obs_data$Y_e, 
                                 Y_a = obs_data$Y_a, 
                                 X_e = obs_data$X_e,
                                 X_a = obs_data$X_a,
                                 Z_e = obs_data$Z_e,
                                 F_a = obs_data$F_a,
                                 ego_id_a = obs_data$ego_id_a,
                                 augmented = TRUE,
                                 reg_model_egos = glm,
                                 reg_model_alters = glm, 
                                 formula_egos = as.formula(Y ~ Z + .),
                                 formula_alters = as.formula(Y ~ F + .),
                                 B = B,
                                 n_cores = N_CORES,
                                 n_folds = 2,
                                 pz = pz,
                                 verbose = TRUE,
                                 prior_func_ie = prior_ie_uniform,
                                 pi_func_ie = pi_hetero,
                                 pi_args_ie = pi_ie_hetero_args,
                                 pi_param_name_ie = "m_vec",
                                 prior_func_de = prior_de_uniform,
                                 pi_func_de = pi_hetero,
                                 pi_args_de = pi_de_hetero_args,
                                 pi_param_name_de = "m_vec", 
                                 family = binomial(link = "logit") # Additional arg for glm
)

pba_hetero_unif$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "hetero"
)]

pba_hetero_unif$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "hetero"
)]



set.seed(3425)
pba_hetero_pois <-  enrt_pba(Y_e = obs_data$Y_e, 
                                 Y_a = obs_data$Y_a, 
                                 X_e = obs_data$X_e,
                                 X_a = obs_data$X_a,
                                 Z_e = obs_data$Z_e,
                                 F_a = obs_data$F_a,
                                 ego_id_a = obs_data$ego_id_a,
                                 augmented = TRUE,
                                 reg_model_egos = glm,
                                 reg_model_alters = glm, 
                                 formula_egos = as.formula(Y ~ Z + .),
                                 formula_alters = as.formula(Y ~ F + .),
                                 B = B,
                                 n_cores = N_CORES,
                                 n_folds = 2,
                                 pz = pz,
                                 verbose = TRUE,
                                 prior_func_ie = prior_ie_poisson,
                                 pi_func_ie = pi_hetero,
                                 pi_args_ie = pi_ie_hetero_args,
                                 pi_param_name_ie = "m_vec",
                                 prior_func_de = prior_de_poisson_uniform,
                                 pi_func_de = pi_hetero,
                                 pi_args_de = pi_de_hetero_args,
                                 pi_param_name_de = "m_vec", 
                                 family = binomial(link = "logit") # Additional arg for glm
)

pba_hetero_pois$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson",
  model = "hetero"
)]

pba_hetero_pois$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson + Uniform",
  model = "hetero"
)]


# combine results
pba_ie_results <- rbindlist(
  list(
    pba_homo_unif$IE_results,
    pba_homo_pois$IE_results,
    pba_hetero_unif$IE_results,
    pba_hetero_pois$IE_results
  )
)

pba_ie_results <- pba_ie_results[augmented == TRUE & uncertainty_type == "total",]

pba_de_results <- rbindlist(
  list(
    pba_homo_unif$DE_results,
    pba_homo_pois$DE_results,
    pba_hetero_unif$DE_results,
    pba_hetero_pois$DE_results
  )
)

pba_de_results <- pba_de_results[augmented == TRUE & uncertainty_type == "total",]

# results of naive estimator (from GSA)
null_results_ie <- gsa_results$null_results$IE
null_results_de <- gsa_results$null_results$DE

null_results_ie[,`:=`(
  ie_rd_mean = ie_rd,
  ie_rd_q_low = ci_low,
  ie_rd_q_high = ci_high,
  augmented = TRUE,
  model = spec,
  uncertainty_type = "total",
  prior = "Naive"
  
)]

null_results_ie <- null_results_ie[,.(
  uncertainty_type, ie_rd_mean, ie_rd_q_low, ie_rd_q_high, augmented, prior, model)
]

null_results_de[,`:=`(
  de_rd_mean = de_rd,
  de_rd_q_low = ci_low,
  de_rd_q_high = ci_high,
  augmented = TRUE,
  model = spec,
  uncertainty_type = "total",
  prior = "Naive"
  
)]

null_results_de <- null_results_de[,.(
  uncertainty_type, de_rd_mean, de_rd_q_low, de_rd_q_high, augmented, prior, model)
]

# combine with pba results
prior_levels_ie <- c("Naive", "Poisson", "Uniform")
ie_pba_plot_data <- rbindlist(list(pba_ie_results, null_results_ie))
ie_pba_plot_data$prior <- factor(ie_pba_plot_data$prior, levels = prior_levels_ie)

prior_levels_de <- c("Naive", "Uniform","Poisson + Uniform")
de_pba_plot_data <- rbindlist(list(pba_de_results, null_results_de))
de_pba_plot_data$prior <- factor(de_pba_plot_data$prior, levels = prior_levels_de)

# Plot results

prior_colors <- c(
  "Uniform" = "#E41A1C", # Red
  "Poisson" = "#377EB8" # Blue
)

model_labels <- c(
  "homo"   = "Homogeneous",
  "hetero" = "Heterogeneous"
)

model_shapes <- c(
  "homo"   = 19,
  "hetero" = 17
)


pba_ie_plot <- ggplot(
    data = ie_pba_plot_data[model != "Naive",],
    aes(x = prior,
        y = ie_rd_mean,
        ymin = ie_rd_q_low,
        ymax = ie_rd_q_high,
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # Naive data
    geom_errorbar(
      data = ie_pba_plot_data[model == "Naive",],
      color = "black",
      width = 0.15,     
      linewidth = 0.8
    ) +
    geom_point(
      data = ie_pba_plot_data[model == "Naive",],
      color = "black",
      size = 4,
      shape = 15
    ) +
    geom_errorbar(
      aes(color = prior, group = model),
      position = position_dodge(.3),
      width = 0.3, 
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = prior, shape = model, group = model),
      position = position_dodge(.3),
      size = 4
    ) +
    scale_color_manual(values = prior_colors) +
    scale_shape_manual(
      name = "Model", # Set legend title
      values = model_shapes,
      labels = model_labels
    ) +
    # scale_y_continuous(breaks = seq(-0.5,0.1,0.1),
    #                    limits = c(-0.59,0.1)) +
    guides(
      color = "none" # Hide the color (prior) legend
    ) +
    labs(
      x = "",
      y = "Estimated IE",
      title = "Indirect Effect"
    ) +
  theme_bw(base_size = 14)  +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 14)) 

print(pba_ie_plot)
ggsave("ENRT_example_code/pba_ie_results.png",
       pba_ie_plot,
       width =8, height = 6, dpi=150)


# DE plot
prior_colors_de <- c(
  "Poisson + Uniform" = "#E41A1C", # Red
  "Uniform"      = "#4DAF4A"  # Green
)

pba_de_plot <- ggplot(
    data = de_pba_plot_data[model != "Naive",],
    aes(x = prior,
        y = de_rd_mean,
        ymin = de_rd_q_low,
        ymax = de_rd_q_high,
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # Naive data
    geom_errorbar(
      data = de_pba_plot_data[model == "Naive",],
      color = "black",
      width = 0.15,     
      linewidth = 0.8
    ) +
    geom_point(
      data = de_pba_plot_data[model == "Naive",],
      color = "black",
      size = 4,
      shape = 15
    ) +
    geom_errorbar(
      aes(color = prior, group = model),
      position = position_dodge(.3),
      width = 0.3, 
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = prior, shape = model, group = model),
      position = position_dodge(.3),
      size = 4
    ) +
    scale_color_manual(values = prior_colors_de) +
    scale_shape_manual(
      name = "Model", # Set legend title
      values = model_shapes,
      labels = model_labels
    ) +
    # scale_y_continuous(breaks = seq(-0.5,0.1,0.1),
    #                    limits = c(-0.59,0.1)) +
    guides(
      color = "none" # Hide the color (prior) legend
    ) +
    labs(
      x = "",
      y = "Estimated DE",
      title = "Direct Effect"
    ) +
  theme_bw(base_size = 14)  +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 14)) 

print(pba_de_plot)
ggsave("ENRT_example_code/pba_de_results.png",
       pba_de_plot,
       width =8, height = 6, dpi=150)

