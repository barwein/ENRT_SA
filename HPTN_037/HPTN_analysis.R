###
# Analyze HPTN-037 data using the SA and PBA for ENRT
###

# --- Libraries ---
library(data.table)
library(ggplot2)
library(latex2exp)
library(kableExtra)

source("src/sensitivity_analysis.R")
source("src/pba.R")
source("src/sensitivity_params.R")

# Global parameters

N_CORES <- 8
# BOOTSTRAP_B <- 1e4
BOOTSTRAP_B <- 1e2  # For testing; change to 1e4 for final analysis

# --- Read data ---
hptn_df <- fread("HPTN_037/hptn_clean_df.csv")


# --- Summary stats and params ---
n_e <- sum(hptn_df$is_ego)
n_a <- nrow(hptn_df) - n_e
print(paste("Number of egos:", n_e))
print(paste("Number of alters:", n_a))

pz <- 0.5  # Randomized 1:1

alter_per_egonet_table <- table(hptn_df[is_ego == FALSE,]$NKID)
print(paste("Mean number of alters per ego network:", mean(alter_per_egonet_table)))
print(paste("Median number of alters per ego network:", median(alter_per_egonet_table)))
print(paste("Max number of alters per ego network:", max(alter_per_egonet_table)))
print(paste("Min number of alters per ego network:", min(alter_per_egonet_table)))      
# Create histogram with ggplot2
netsize_hist <- ggplot(data.frame(alter_per_egonet = as.numeric(alter_per_egonet_table)),
         aes(x = alter_per_egonet)) +
    geom_histogram(binwidth = 1, fill = "gray50", color = "black") +
    scale_x_continuous(breaks = seq(0, max(alter_per_egonet_table), by = 1)) +
    labs(title = "Number of Alters per Ego-network",
         x = "Number of Alters",
         y = "Frequency") +
    theme_minimal()
ggsave("HPTN_037/figures/alter_per_egonet_histogram.png",
       plot = netsize_hist, width = 6, height = 4,
       dpi = 300,
       bg = "white")


# --- Prepare data for SA + PBA ---

ego_dt <- hptn_df[is_ego == TRUE, ]
alter_dt <- hptn_df[is_ego == FALSE, ]

# Create ego_id_a mapping
ego_id_a_vec <- match(alter_dt$NKID, ego_dt$NKID)

# Outcomes
Y_e <- ego_dt$injrisk_any_6m
Y_a <- alter_dt$injrisk_any_6m

# Treatment and observed exposures
Z_e <- ego_dt$ego_treat
F_a_tilde <- alter_dt$ego_treat

# Covariates
# We have unit-level and network level summary
# We're using only the unit-level for the sensitivity parameters \pi 
unit_level_covar <- c("injrisk_any_base",
                      "drunk_base",
                      "heroin_and_cocaine_base",
                      "is_male",
                      "age",
                      "nonwhite",
                      "hispanic")

net_level_covar <- c("net_avg_age",
                     "net_prev_nonwhite",
                     "net_prop_male",
                     "net_prev_injrisk_any_base",
                     "net_prev_cocaine_base",
                     "net_prev_heroin_and_cocaine")

X_e_ind <- as.matrix(ego_dt[, ..unit_level_covar])
X_e_net <- as.matrix(ego_dt[, ..net_level_covar])
X_e <- cbind(X_e_ind, X_e_net)
X_a_ind <- as.matrix(alter_dt[, ..unit_level_covar])
X_a_net <- as.matrix(alter_dt[, ..net_level_covar])
X_a <- cbind(X_a_ind, X_a_net)


# --- Sensitivity parameters ---
m_vec_egos <- seq(10, 300, by=10)
m_vec_alters <- seq(10, 400, by=10)

kappa_vec <- seq(1.1, 3.0, by=0.1)

pi_ee_homo <- pi_homo(m_vec = m_vec_egos,
                      n_e = n_e,
                      pz = pz,
                      type = "ego")

pi_ee_hetero <- pi_hetero(X_e = X_e_ind,
                          m_vec = m_vec_egos,
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          pz = pz)

pi_ae_homo <- pi_homo(m_vec = m_vec_alters,
                      n_e = n_e,
                      n_a = n_a,
                      type = "alter",
                      pz =  pz)

pi_ae_hetero <- pi_hetero(X_e = X_e_ind, 
                          X_a = X_a_ind, 
                          m_vec = m_vec_alters, 
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          ego_index = ego_id_a_vec, 
                          pz = pz)


# --- Run Sensitivity Analysis ---

set.seed(142)
hptn_sa_bootstrap_aug <- enrt_sa(Y_e = Y_e,
                                 Y_a = Y_a,
                                 X_e = X_e,
                                 X_a = X_a,
                                 Z_e = Z_e,
                                 F_a = F_a_tilde, 
                                 ego_id_a = ego_id_a_vec,
                                 reg_model_egos = glm,
                                 reg_model_alters = glm, 
                                 formula_egos = as.formula(Y ~ Z + .),
                                 formula_alters = as.formula(Y ~ F + .),
                                 pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                                         "homo" = pi_ee_homo),
                                 pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                           "homo" = pi_ae_homo),
                                 kappa_vec = kappa_vec, 
                                 n_cores = N_CORES,
                                 pz = pz,
                                 bootstrap = TRUE,
                                 B = BOOTSTRAP_B, 
                                 plot = TRUE,
                                 family = binomial(link = "logit") # Additional arg for glm
)

set.seed(242)
hptn_sa_bootstrap_not_aug <- enrt_sa(Y_e = Y_e,
                                 Y_a = Y_a,
                                 X_e = X_e,
                                 X_a = X_a,
                                 Z_e = Z_e,
                                 F_a = F_a_tilde, 
                                 ego_id_a = ego_id_a_vec,
                                 # reg_model_egos = glm,
                                 # reg_model_alters = glm, 
                                 # formula_egos = as.formula(Y ~ Z + .),
                                 # formula_alters = as.formula(Y ~ F + .),
                                 pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                                         "homo" = pi_ee_homo),
                                 pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                           "homo" = pi_ae_homo),
                                 kappa_vec = kappa_vec, 
                                 n_cores = N_CORES,
                                 pz = pz,
                                 bootstrap = TRUE,
                                 B = BOOTSTRAP_B, 
                                 plot = TRUE)



set.seed(342)
hptn_sa_aug <- enrt_sa(Y_e = Y_e,
                                 Y_a = Y_a,
                                 X_e = X_e,
                                 X_a = X_a,
                                 Z_e = Z_e,
                                 F_a = F_a_tilde, 
                                 ego_id_a = ego_id_a_vec,
                                 reg_model_egos = glm,
                                 reg_model_alters = glm, 
                                 formula_egos = as.formula(Y ~ Z + .),
                                 formula_alters = as.formula(Y ~ F + .),
                                 pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                                         "homo" = pi_ee_homo),
                                 pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                           "homo" = pi_ae_homo),
                                 kappa_vec = kappa_vec, 
                                 n_cores = N_CORES,
                                 pz = pz,
                                 bootstrap = FALSE,
                                 # B = BOOTSTRAP_B, 
                                 plot = TRUE,
                                 family = binomial(link = "logit") # Additional arg for glm
)

set.seed(442)
hptn_sa_not_aug <- enrt_sa(Y_e = Y_e,
                                 Y_a = Y_a,
                                 X_e = X_e,
                                 X_a = X_a,
                                 Z_e = Z_e,
                                 F_a = F_a_tilde, 
                                 ego_id_a = ego_id_a_vec,
                                 # reg_model_egos = glm,
                                 # reg_model_alters = glm, 
                                 # formula_egos = as.formula(Y ~ Z + .),
                                 # formula_alters = as.formula(Y ~ F + .),
                                 pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                                         "homo" = pi_ee_homo),
                                 pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                           "homo" = pi_ae_homo),
                                 kappa_vec = kappa_vec, 
                                 n_cores = N_CORES,
                                 pz = pz,
                                 bootstrap = FALSE,
                                 # B = BOOTSTRAP_B, 
                                 plot = TRUE)



# Combine IE results (bootstrap var; with/w.o. outcome model augmentation)
sa_ie_aug <- hptn_sa_bootstrap_aug$sa_results$IE
sa_ie_aug[, aug := TRUE]

sa_ie_no_aug <- hptn_sa_bootstrap_not_aug$sa_results$IE
sa_ie_no_aug[, aug := FALSE]

sa_ie_aug_null <- hptn_sa_bootstrap_aug$null_results$IE
sa_ie_aug_null[, aug := TRUE]

sa_ie_no_aug_null <- hptn_sa_bootstrap_not_aug$null_results$IE
sa_ie_no_aug_null[, aug := FALSE]

sa_ie_res <- rbindlist(list(
  sa_ie_aug,
  sa_ie_no_aug,
  sa_ie_aug_null,
  sa_ie_no_aug_null
))

sa_ie_res[, spec_name := ifelse(spec == "homo", "Homo",
                                ifelse(spec == "hetero",
                                       "Hetero", spec))]
# Plot

sa_ie_naive_data <- sa_ie_res[spec == "Naive"]
sa_ie_data <- sa_ie_res[spec != "Naive"]

# Create a new interaction variable for the 4 SA models
sa_ie_data[, model_type := paste0(spec_name," & ", ifelse(aug, "Aug", "Not-Aug"))]

# Define the order for the legend and colors
model_levels <- c("Homo & Not-Aug", 
                  "Homo & Aug", 
                  "Hetero & Not-Aug", 
                  "Hetero & Aug")
sa_ie_data[, model_type := factor(model_type, levels = model_levels)]

# --- 2. Define Colors and Labels ---

# A logical color palette: 
# Blues for "homo", Oranges for "hetero"
# Light shade for Unaugmented, Dark shade for Augmented
palette <- c(
  "Homo & Not-Aug"   = "#1C8BD9", # Light Blue
  # "Homo & Not-Aug"   = "#7DAFD1", # Light Blue
  # "Homo & Aug"     = "#1F78B4", # Dark Blue
  "Homo & Aug"     = "#085380", # Dark Blue
  # "Hetero & Not-Aug" = "#FDBF6F", # Light Orange
  "Hetero & Not-Aug" = "#DBA663", # Light Orange
  # "Hetero & Aug"   = "#FF7F00"  # Dark Orange
  "Hetero & Aug"   = "#874D08"  # Dark Orange
)

# Labels for the legend
plot_labels <- c(
  "Homo & Not-Aug" = "Homo",
  "Homo & Aug"     = "Homo (Aug)",
  "Hetero & Not-Aug" = "Hetero",
  "Hetero & Aug"   = "Hetero (Aug)"
)

# --- 3. Create the Plot ---

ie_plot <- ggplot() + 
  geom_ribbon(data = sa_ie_data,
  aes(
    x = pi_param,
    ymin = ci_low,
    ymax = ci_high,
    fill = model_type
  ),
  alpha = 0.3 
) +
  geom_line(
    data = sa_ie_data,
    aes(
      x = pi_param,
      y = ie_rd,
      color = model_type
    ),
    linewidth = 1.5
  ) +
  # Add the Naive point estimates
  geom_point(
    data = sa_ie_naive_data,
    aes(
      x = pi_param,
      y = ie_rd,
      shape = aug,
      group = aug
    ),
    color = "black",
    size = 4,
    position = position_dodge(12)
  ) +
  # Add the Naive error bars
  geom_errorbar(
    data = sa_ie_naive_data,
    aes(
      x = pi_param,
      ymin = ci_low,
      ymax = ci_high,
      group = aug
    ),
    color = "black",
    width = 10,
    linewidth = 0.7,
    position = position_dodge(12)
  ) +
  # Horizontal line at y=0 (no effect)
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey48"
  ) +
  # Use the custom color palette for lines and fills
  scale_color_manual(
    name = "SA Model:",
    values = palette,
    labels = plot_labels
  ) +
  scale_fill_manual(
    name = "SA Model:",
    values = palette,
    labels = plot_labels
  ) +
  # Customize the shape legend
  scale_shape_manual(
    name = "Naive Model:",
    values = c("TRUE" = 17, "FALSE" = 16), # Triangle and Circle
    labels = c("TRUE" = "Augmented", "FALSE" = "Unaugmented")
  ) +
  # Add labels (using latex2exp for the x-axis)
  labs(
    title = "Sensitivity Analysis for Indirect Effect",
    y = "Estimated IE",
    x = TeX("$m_a$ (Expected number of missing alter-ego edges)")
  ) +
  
  # Adjust x-axis scale
  scale_x_continuous(breaks = seq(0, max(m_vec_alters), 50),
                     labels = seq(0, max(m_vec_alters), 50)) +
  # scale_y_continuous(breaks = c(-0.45, -0.35, -0.25, -0.15, -0.05,0,0.05,0.15)) +
  scale_y_continuous(breaks = seq(-0.5,0.2,0.1)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(1.5, "cm"),
    axis.text = element_text(size = 12)
  )

# Display the plot (IE)
print(ie_plot)
ggsave("HPTN_037/figures/sa_ie_plot.png",
       plot = ie_plot, width = 8, height = 5,
       dpi = 300,
       bg = "white")

# TODO: make two plots: augmented & not augmented
#       Probably will show only the augmented in the main text.


# DE plot (bootstrap & augmented scenario)
spec_labels <- c(
  "homo" = "Homogeneous model",
  "hetero" = "Heterogeneous model"
  # "homo" = TeX("Homogeneous $\\pi_i^e$"),
  # "hetero" = TeX("Heterogeneous $\\pi_i^e$")
)
de_plot <- hptn_sa_bootstrap_aug$de_rd_plot
de_plot <- de_plot + 
        scale_x_continuous(breaks = seq(0, max(m_vec_egos), 50),
                           labels = seq(0, max(m_vec_egos), 50)) +
        labs(x = TeX("$m_e$ (Expected number of missing ego-ego edges)"),
             fill = "Estimated DE",
             title = "Sensitivity Analysis for Direct Effect") +
  facet_wrap(~spec,
             labeller = labeller(spec = spec_labels)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face="bold")
  )
print(de_plot)

ggsave("HPTN_037/figures/sa_de_plot.png",
       plot = de_plot, width = 10, height = 6,
       dpi = 300,
       bg = "white")



# --- PBA results ---

# Set priors.

# Priors for IE (for m_e)

prior_ie_uniform <- function() {
  sample(seq(1, 400, by=1), 1)
}

prior_ie_poisson <- function() {
  rpois(1, lambda = 200)
}

prior_ie_neg_binom <- function() {
  rnbinom(1, size = 10, mu = 200)
}

# Priors for DE (for m_a and kappa)

prior_de_uniform <- function() {
  list(
    pi_param = sample(seq(1, 300, by=1), 1),
    kappa = runif(1, min = 1, max = 3.0)
  )
}

prior_de_nb_lognormal <- function() {
  list(
    pi_param = rnbinom(1, size = 10, mu = 150),
    kappa = rlnorm(1, meanlog = log(2), sdlog = 0.2)
  )
}

prior_de_poisson_uniform <- function() {
  list(
    pi_param = rpois(1, lambda = 150),
    kappa = runif(1, min = 1, max = 3.0)
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
                          ego_index = ego_id_a_vec,
                          pz=pz)

pi_de_homo_args <- list(n_e = n_e, type = "ego", pz = 0.5)

pi_de_hetero_args <- list(X_e = X_e,
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          pz = pz)

# --- Run PBA ---
# We run with/without augmented model + bootstrap for each prior spec

set.seed(142)
pba_homo_unif_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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

pba_homo_unif_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "homo"
)]

pba_homo_unif_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "homo"
)]


set.seed(242)
pba_homo_unif <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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
  augmented = FALSE,
  prior = "Uniform",
  model = "homo"
)]

pba_homo_unif$DE_results[, `:=`(
  augmented = FALSE,
  prior = "Uniform",
  model = "homo"
)]


set.seed(342)
pba_homo_pois_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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

pba_homo_pois_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson",
  model = "homo"
)]

pba_homo_pois_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson + Uniform",
  model = "homo"
)]


set.seed(442)
pba_homo_pois <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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
  augmented = FALSE,
  prior = "Poisson",
  model = "homo"
)]

pba_homo_pois$DE_results[, `:=`(
  augmented = FALSE,
  prior = "Poisson + Uniform",
  model = "homo"
)]


set.seed(542)
pba_homo_nb_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_neg_binom,
                               pi_func_ie = pi_homo,
                               pi_args_ie = pi_ie_homo_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_nb_lognormal,
                               pi_func_de = pi_homo,
                               pi_args_de = pi_de_homo_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_homo_nb_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "NB",
  model = "homo"
)]

pba_homo_nb_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "NB + Lognormal",
  model = "homo"
)]

set.seed(642)
pba_homo_nb <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_neg_binom,
                               pi_func_ie = pi_homo,
                               pi_args_ie = pi_ie_homo_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_nb_lognormal,
                               pi_func_de = pi_homo,
                               pi_args_de = pi_de_homo_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_homo_nb$IE_results[, `:=`(
  augmented = FALSE,
  prior = "NB",
  model = "homo"
)]

pba_homo_nb$DE_results[, `:=`(
  augmented = FALSE,
  prior = "NB + Lognormal",
  model = "homo"
)]


set.seed(1425)
pba_hetero_unif_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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

pba_hetero_unif_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "hetero"
)]

pba_hetero_unif_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Uniform",
  model = "hetero"
)]


set.seed(2425)
pba_hetero_unif <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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
  augmented = FALSE,
  prior = "Uniform",
  model = "hetero"
)]

pba_hetero_unif$DE_results[, `:=`(
  augmented = FALSE,
  prior = "Uniform",
  model = "hetero"
)]



set.seed(3425)
pba_hetero_pois_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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

pba_hetero_pois_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson",
  model = "hetero"
)]

pba_hetero_pois_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "Poisson + Uniform",
  model = "hetero"
)]



set.seed(4425)
pba_hetero_pois <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
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
  augmented = FALSE,
  prior = "Poisson",
  model = "hetero"
)]

pba_hetero_pois$DE_results[, `:=`(
  augmented = FALSE,
  prior = "Poisson + Uniform",
  model = "hetero"
)]

set.seed(5425)
pba_hetero_nb_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               reg_model_egos = glm,
                               reg_model_alters = glm, 
                               formula_egos = as.formula(Y ~ Z + .),
                               formula_alters = as.formula(Y ~ F + .),
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_neg_binom,
                               pi_func_ie = pi_hetero,
                               pi_args_ie = pi_ie_hetero_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_nb_lognormal,
                               pi_func_de = pi_hetero,
                               pi_args_de = pi_de_hetero_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_hetero_nb_aug$IE_results[, `:=`(
  augmented = TRUE,
  prior = "NB",
  model = "hetero"
)]

pba_hetero_nb_aug$DE_results[, `:=`(
  augmented = TRUE,
  prior = "NB + Lognormal",
  model = "hetero"
)]


set.seed(6425)
pba_hetero_nb <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
                               bootstrap = TRUE,
                               B = BOOTSTRAP_B,
                               n_cores = N_CORES,
                               pz = pz,
                               verbose = TRUE,
                               prior_func_ie = prior_ie_neg_binom,
                               pi_func_ie = pi_hetero,
                               pi_args_ie = pi_ie_hetero_args,
                               pi_param_name_ie = "m_vec",
                               prior_func_de = prior_de_nb_lognormal,
                               pi_func_de = pi_hetero,
                               pi_args_de = pi_de_hetero_args,
                               pi_param_name_de = "m_vec", 
                               family = binomial(link = "logit") # Additional arg for glm
)

pba_hetero_nb$IE_results[, `:=`(
  augmented = FALSE,
  prior = "NB",
  model = "hetero"
)]

pba_hetero_nb$DE_results[, `:=`(
  augmented = FALSE,
  prior = "NB + Lognormal",
  model = "hetero"
)]


# combine IE results
pba_ie_results <- rbindlist(
  list(
    pba_homo_unif_aug$IE_results,
    pba_homo_unif$IE_results,
    pba_homo_pois_aug$IE_results,
    pba_homo_pois$IE_results,
    pba_homo_nb_aug$IE_results,
    pba_homo_nb$IE_results,
    pba_hetero_unif_aug$IE_results,
    pba_hetero_unif$IE_results,
    pba_hetero_pois_aug$IE_results,
    pba_hetero_pois$IE_results,
    pba_hetero_nb_aug$IE_results,
    pba_hetero_nb$IE_results
  )
)

pba_de_results <- rbindlist(
  list(
    pba_homo_unif_aug$DE_results,
    pba_homo_unif$DE_results,
    pba_homo_pois_aug$DE_results,
    pba_homo_pois$DE_results,
    pba_homo_nb_aug$DE_results,
    pba_homo_nb$DE_results,
    pba_hetero_unif_aug$DE_results,
    pba_hetero_unif$DE_results,
    pba_hetero_pois_aug$DE_results,
    pba_hetero_pois$DE_results,
    pba_hetero_nb_aug$DE_results,
    pba_hetero_nb$DE_results
  )
)


# Plot results
naive_ie_aug_boot <- hptn_sa_bootstrap_aug$null_results$IE
# naive_ie_aug_boot[,`;=`(
#   model = spec,
#   augmented = aug,
#   ie_rd_mean = ie_rd,
#   ie_rd_q_low = ci_low,
#   ie_rd_q_high = ci_high
#   
# )]

# TODO: make the plot with the naive estimator to the left as well.
#       Repeat for Augmented & non augmented and for DE (same)

pba_ie_plot_aug <- ggplot(pba_ie_results[augmented == TRUE & uncertainty_type == "total",],
                      aes(x = prior,
                          y = ie_rd_mean,
                          ymin = ie_rd_q_low,
                          ymax = ie_rd_q_high,
                          color = prior,
                          shape = model)) +
  geom_point(position = position_dodge(.5),
             size = 3) +
  geom_errorbar(position = position_dodge(.5),
                width = 0.3, 
                linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "",
    y = "Estimated IE",
    title = "Probabilistic Bias Analysis for Indirect Effect"
  ) +
  theme_bw(base_size = 14)  +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 15, face = "bold")) 

print(pba_ie_plot_aug)

