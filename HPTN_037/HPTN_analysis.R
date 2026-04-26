###
# Analyze HPTN-037 data using the SA and PBA for ENRT
###

# --- Libraries ---
library(data.table)
library(ggplot2)
library(latex2exp)
library(kableExtra)

library(ENRTsensitivity)

# source("src/sensitivity_analysis.R")
# source("src/pba.R")
# source("src/sensitivity_params.R")
source("HPTN_037/hptn_plot_functions.R")

# Global parameters

N_CORES <- 8
B <- 1e4
# B <- 1e2  # For testing; change to 1e4 for final analysis

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

# Egos covariates
X_e_ind <- ego_dt[, ..unit_level_covar]
X_e_net <- ego_dt[, ..net_level_covar]
# Standarize Age variable 
X_e_ind$age <- scale(X_e_ind$age)
X_e_net$net_avg_age <- scale(X_e_net$net_avg_age)
# As matrix
X_e_ind <- as.matrix(X_e_ind)
X_e_net <- as.matrix(X_e_net)
X_e <- cbind(X_e_ind, X_e_net)

# Alters covariates
X_a_ind <- alter_dt[, ..unit_level_covar]
X_a_net <- alter_dt[, ..net_level_covar]
# Standarize Age variable
X_a_ind$age <- scale(X_a_ind$age)
X_a_net$net_avg_age <- scale(X_a_net$net_avg_age)
# As matrix
X_a_ind <- as.matrix(X_a_ind)
X_a_net <- as.matrix(X_a_net)
X_a <- cbind(X_a_ind, X_a_net)


# --- Sensitivity parameters ---
m_vec_egos <- seq(10, 150, by=5)
m_vec_alters <- seq(10, 500, by=10)
# m_vec_egos <- seq(10, 150, by=10)
# m_vec_alters <- seq(10, 500, by=10)

# kappa_vec <- seq(1.1, 3.0, by=0.1)
# kappa_vec <- seq(0.5, 2.0, by=0.05)
kappa_vec <- seq(1.0, 2.0, by=0.05)

pi_ee_homo <- pi_homo(m_vec = m_vec_egos,
                      n_e = n_e,
                      pz = pz,
                      type = "ego")

pi_ee_hetero <- pi_hetero(X_e = X_e_ind,
                          m_vec = m_vec_egos,
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          # p = 1,
                          pz = pz)

pi_ee_hetero_gamma2 <- pi_hetero(X_e = X_e_ind,
                          m_vec = m_vec_egos,
                          gamma = -2,
                          dist = "norm",
                          p = 2,
                          # p = 1,
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
                          # p = 2,
                          p = 1,
                          ego_index = ego_id_a_vec, 
                          pz = pz)

pi_ae_hetero_gamma2 <- pi_hetero(X_e = X_e_ind, 
                          X_a = X_a_ind, 
                          m_vec = m_vec_alters, 
                          gamma = -2,
                          dist = "norm",
                          # p = 2,
                          p = 1,
                          ego_index = ego_id_a_vec, 
                          pz = pz)


# --- Run Sensitivity Analysis ---

set.seed(342)
hptn_sa_aug <- enrt_sa(Y_e = Y_e,
                       Y_a = Y_a,
                       X_e = X_e,
                       X_a = X_a,
                       Z_e = Z_e,
                       F_a = F_a_tilde, 
                       ego_id_a = ego_id_a_vec,
                       augmented = TRUE,
                       reg_model_egos = glm,
                       reg_model_alters = glm, 
                       formula_egos = as.formula(Y ~ Z + .),
                       formula_alters = as.formula(Y ~ F + .),
                       pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                               # "hetero_gamma2" = pi_ee_hetero_gamma2,
                                               "homo" = pi_ee_homo),
                       pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                 # "hetero_gamma2" = pi_ae_hetero_gamma2,
                                                 "homo" = pi_ae_homo),
                       kappa_vec = kappa_vec, 
                       n_cores = N_CORES,
                       n_folds = 2,
                       pz = pz,
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
                           augmented = FALSE,
                           # reg_model_egos = glm,
                           # reg_model_alters = glm, 
                           # formula_egos = as.formula(Y ~ Z + .),
                           # formula_alters = as.formula(Y ~ F + .),
                           pi_lists_ego_ego = list("hetero" = pi_ee_hetero,
                                                     # "hetero_gamma2" = pi_ee_hetero_gamma2,
                                                   "homo" = pi_ee_homo),
                           pi_lists_alter_ego = list("hetero" = pi_ae_hetero,
                                                     # "hetero_gamma2" = pi_ae_hetero_gamma2,
                                                     "homo" = pi_ae_homo),
                           kappa_vec = kappa_vec, 
                           n_cores = N_CORES,
                           n_folds = 2,
                           pz = pz,
                           plot = TRUE)



# Combine IE results (bootstrap var; with/w.o. outcome model augmentation)
sa_ie_aug <- hptn_sa_aug$sa_results$IE
sa_ie_aug[, aug := TRUE]

sa_ie_no_aug <- hptn_sa_not_aug$sa_results$IE
sa_ie_no_aug[, aug := FALSE]

sa_ie_aug_null <- hptn_sa_aug$null_results$IE
sa_ie_aug_null[, aug := TRUE]

sa_ie_no_aug_null <- hptn_sa_not_aug$null_results$IE
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

sa_ie_naive_data[,pi_param := as.numeric(pi_param)]

# Create a new interaction variable for the 4 SA models
sa_ie_data[, model_type := paste0(spec_name," & ", ifelse(aug, "Aug", "Not-Aug"))]

# Define the order for the legend and colors
model_levels <- c("Homo & Not-Aug", 
                  "Homo & Aug", 
                  "Hetero & Not-Aug", 
                  "Hetero & Aug")
sa_ie_data[, model_type := factor(model_type, levels = model_levels)]
sa_ie_data[,pi_param := as.numeric(pi_param)]
# --- 2. Define Colors and Labels ---

# A logical color palette: 
# Blues for "homo", Oranges for "hetero"
# Light shade for Unaugmented, Dark shade for Augmented
palette <- c(
  "Homo & Not-Aug"   = "#085380", # Light Blue
  # "Homo & Not-Aug"   = "#7DAFD1", # Light Blue
  # "Homo & Aug"     = "#1F78B4", # Dark Blue
  "Homo & Aug"     = "#085380", # Dark Blue
  # "Hetero & Not-Aug" = "#FDBF6F", # Light Orange
  "Hetero & Not-Aug" = "#BA7013", # Light Orange
  # "Hetero & Aug"   = "#FF7F00"  # Dark Orange
  "Hetero & Aug"   = "#BA7013"  # Dark Orange
)

# Labels for the legend
plot_labels <- c(
  "Homo & Not-Aug" = "Homogeneous", 
  "Homo & Aug"     = "Homogeneous",
  "Hetero & Not-Aug" = "Heterogeneous",
  "Hetero & Aug"   = "Heterogeneous"
)

# --- 3. Create the Plot ---

ie_plot_aug <- sa_ie_plot(sa_ie_data = sa_ie_data[aug == TRUE,],
                          sa_ie_naive_data = sa_ie_naive_data[aug == TRUE,],
                          palette = palette,
                          plot_labels = plot_labels,
                          m_vec_alters = m_vec_alters)
ie_plot_aug <- ie_plot_aug + geom_vline(xintercept = 105.5, linetype = 3, 
                         color = "grey48")
print(ie_plot_aug)
ggsave("HPTN_037/figures/sa_ie_plot_augmented.png",
       plot = ie_plot_aug, width = 8, height = 5,
       dpi = 300,
       bg = "white")

ie_plot_not_aug <- sa_ie_plot(sa_ie_data = sa_ie_data[aug == FALSE,],
                          sa_ie_naive_data = sa_ie_naive_data[aug == FALSE,],
                          palette = palette,
                          plot_labels = plot_labels,
                          m_vec_alters = m_vec_alters)
ie_plot_not_aug <- ie_plot_not_aug + geom_vline(xintercept = 105.5, linetype = 3, 
                                        color = "grey48")
print(ie_plot_not_aug)
ggsave("HPTN_037/figures/sa_ie_plot_not_augmented.png",
       plot = ie_plot_not_aug, width = 8, height = 5,
       dpi = 300,
       bg = "white")

# DE plot 
spec_labels <- c(
  "homo" = "Homogeneous model",
  "hetero" = "Heterogeneous model"
  # "homo" = TeX("Homogeneous $\\pi_i^e$"),
  # "hetero" = TeX("Heterogeneous $\\pi_i^e$")
)
# de_plot <- hptn_sa_aug$de_rd_plot
de_plot_not_aug <- hptn_sa_not_aug$de_rd_plot
de_plot_not_aug <- de_plot_not_aug + 
        scale_x_continuous(breaks = seq(0, max(m_vec_egos), 25),
                           labels = seq(0, max(m_vec_egos), 25)) +
        labs(x = TeX("$m^e$ (Expected number of missing ego-ego edges)"),
             fill = "Estimated DE",
             title = "Direct Effect") +
    geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48") +
  facet_wrap(~spec,
             labeller = labeller(spec = spec_labels)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face="bold")
  )
print(de_plot_not_aug)

ggsave("HPTN_037/figures/sa_de_not_aug_plot.png",
       plot = de_plot_not_aug, width = 10, height = 6,
       dpi = 300,
       bg = "white")

de_plot_aug <- hptn_sa_aug$de_rd_plot
de_plot_aug <- de_plot_aug + 
        scale_x_continuous(breaks = seq(0, max(m_vec_egos), 25),
                           labels = seq(0, max(m_vec_egos), 25)) +
        labs(x = TeX("$m^e$ (Expected number of missing ego-ego edges)"),
             fill = "Estimated DE",
             title = "Direct Effect") +
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48") +
  facet_wrap(~spec,
             labeller = labeller(spec = spec_labels)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face="bold")
  )
print(de_plot_aug)

ggsave("HPTN_037/figures/sa_de_aug_plot.png",
       plot = de_plot_aug, width = 10, height = 6,
       dpi = 300,
       bg = "white")




# DE plot given kappa=1.5
palette <- c(
  "homo"   = "#085380", # Light Blue
  "hetero"   = "#BA7013"  # Dark Orange
)

# Labels for the legend
plot_labels <- c(
  "homo" = "Homogeneous", 
  "hetero"= "Heterogeneous"
)


hptn_sa_not_aug$sa_results$DE[, pi_param := as.numeric(pi_param)]
hptn_sa_not_aug$null_results$DE[, pi_param := as.numeric(pi_param)]
de_plot_kappa1.5_not_aug <- sa_de_plot_given_kappa(
                                         sa_de_data = hptn_sa_not_aug$sa_results$DE[kappa == 1.5,],
                                         sa_de_naive_data = hptn_sa_not_aug$null_results$DE,
                                         # sa_de_data = hptn_sa_aug$sa_results$DE[kappa == 1.5,],
                                         # sa_de_naive_data = hptn_sa_aug$null_results$DE,
                                         kappa = 1.5,
                                         palette = palette,
                                         plot_labels = plot_labels, 
                                         m_vec_egos = m_vec_egos)
de_plot_kappa1.5_not_aug <- de_plot_kappa1.5_not_aug + 
        geom_vline(xintercept = 33.4, linetype = 3, 
               color = "grey48") 
print(de_plot_kappa1.5_not_aug)
ggsave("HPTN_037/figures/sa_de_plot_not_aug_kappa1.5.png",
       plot = de_plot_kappa1.5_not_aug, width = 8, height = 5,
       dpi = 300,
       bg = "white")

hptn_sa_aug$sa_results$DE[, pi_param := as.numeric(pi_param)]
hptn_sa_aug$null_results$DE[, pi_param := as.numeric(pi_param)]
de_plot_kappa1.5_aug <- sa_de_plot_given_kappa(
                                         sa_de_data = hptn_sa_aug$sa_results$DE[kappa == 1.5,],
                                         sa_de_naive_data = hptn_sa_aug$null_results$DE,
                                         # sa_de_data = hptn_sa_aug$sa_results$DE[kappa == 1.5,],
                                         # sa_de_naive_data = hptn_sa_aug$null_results$DE,
                                         kappa = 1.5,
                                         palette = palette,
                                         plot_labels = plot_labels, 
                                         m_vec_egos = m_vec_egos)
de_plot_kappa1.5_aug <- de_plot_kappa1.5_aug + 
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48") 
print(de_plot_kappa1.5_aug) 
ggsave("HPTN_037/figures/sa_de_plot_aug_kappa1.5.png",
       plot = de_plot_kappa1.5_aug, width = 8, height = 5,
       dpi = 300,
       bg = "white")

# Run analysis with gamma=2

set.seed(342)
hptn_sa_aug_gamma2 <- enrt_sa(Y_e = Y_e,
                       Y_a = Y_a,
                       X_e = X_e,
                       X_a = X_a,
                       Z_e = Z_e,
                       F_a = F_a_tilde, 
                       ego_id_a = ego_id_a_vec,
                       augmented = TRUE,
                       reg_model_egos = glm,
                       reg_model_alters = glm, 
                       formula_egos = as.formula(Y ~ Z + .),
                       formula_alters = as.formula(Y ~ F + .),
                       pi_lists_ego_ego = list(
                         "hetero_gamma1" = pi_ee_hetero,
                         "hetero_gamma2" = pi_ee_hetero_gamma2,
                         "homo" = pi_ee_homo),
                       pi_lists_alter_ego = list(
                         "hetero_gamma1" = pi_ae_hetero,
                         "hetero_gamma2" = pi_ae_hetero_gamma2,
                         "homo" = pi_ae_homo),
                       kappa_vec = kappa_vec, 
                       n_cores = N_CORES,
                       n_folds = 2,
                       pz = pz,
                       plot = TRUE,
                       family = binomial(link = "logit") # Additional arg for glm
)

set.seed(442)
hptn_sa_not_aug_gamma2<- enrt_sa(Y_e = Y_e,
                           Y_a = Y_a,
                           X_e = X_e,
                           X_a = X_a,
                           Z_e = Z_e,
                           F_a = F_a_tilde, 
                           ego_id_a = ego_id_a_vec,
                           augmented = FALSE,
                           # reg_model_egos = glm,
                           # reg_model_alters = glm, 
                           # formula_egos = as.formula(Y ~ Z + .),
                           # formula_alters = as.formula(Y ~ F + .),
                           pi_lists_ego_ego = list("hetero_gamma1" = pi_ee_hetero,
                                                   "hetero_gamma2" = pi_ee_hetero_gamma2,
                                                   "homo" = pi_ee_homo),
                           pi_lists_alter_ego = list("hetero_gamma1" = pi_ae_hetero,
                                                     "hetero_gamma2" = pi_ae_hetero_gamma2,
                                                     "homo" = pi_ae_homo),
                           kappa_vec = kappa_vec, 
                           n_cores = N_CORES,
                           n_folds = 2,
                           pz = pz,
                           plot = TRUE)



# Combine IE results (bootstrap var; with/w.o. outcome model augmentation)
sa_ie_aug_g2 <- hptn_sa_aug_gamma2$sa_results$IE
sa_ie_aug_g2[, aug := TRUE]

sa_ie_no_aug_g2 <- hptn_sa_not_aug_gamma2$sa_results$IE
sa_ie_no_aug_g2[, aug := FALSE]

sa_ie_aug_null_g2 <- hptn_sa_aug_gamma2$null_results$IE
sa_ie_aug_null_g2[, aug := TRUE]

sa_ie_no_aug_null_g2 <- hptn_sa_not_aug_gamma2$null_results$IE
sa_ie_no_aug_null_g2[, aug := FALSE]

sa_ie_res_g2 <- rbindlist(list(
  sa_ie_aug_g2,
  sa_ie_no_aug_g2,
  sa_ie_aug_null_g2,
  sa_ie_no_aug_null_g2
))

sa_ie_res_g2[, spec_name := ifelse(spec == "homo", "Homo",
                                ifelse(spec == "hetero_gamma1",
                                       "Hetero_g1", 
                                       ifelse(spec == "hetero_gamma2",
                                              "Hetero_g2", spec)))]

sa_ie_naive_data_g2 <- sa_ie_res_g2[spec == "Naive"]
sa_ie_data_g2 <- sa_ie_res_g2[spec != "Naive"]

sa_ie_naive_data_g2[,pi_param := as.numeric(pi_param)]



# Create a new interaction variable for the 6 SA models
sa_ie_data_g2[, model_type := paste0(spec_name," & ", ifelse(aug, "Aug", "Not-Aug"))]

# Define the order for the legend and colors
model_levels_g2 <- c("Homo & Not-Aug", 
                  "Homo & Aug", 
                  "Hetero_g1 & Not-Aug", 
                  "Hetero_g1 & Aug",
                  "Hetero_g2 & Not-Aug", 
                  "Hetero_g2 & Aug")
sa_ie_data_g2[, model_type := factor(model_type, levels = model_levels_g2)]
sa_ie_data_g2[,pi_param := as.numeric(pi_param)]

# A logical color palette: 
# Blues for "homo", Oranges for "hetero"
# Light shade for Unaugmented, Dark shade for Augmented
palette_g2 <- c(
  "Homo & Not-Aug"   = "#085380", # Light Blue
  # "Homo & Not-Aug"   = "#7DAFD1", # Light Blue
  # "Homo & Aug"     = "#1F78B4", # Dark Blue
  "Homo & Aug"     = "#085380", # Dark Blue
  # "Hetero & Not-Aug" = "#FDBF6F", # Light Orange
  "Hetero_g1 & Not-Aug" = "#BA7013", # Light Orange
  # "Hetero & Aug"   = "#FF7F00"  # Dark Orange
  "Hetero_g1 & Aug"   = "#BA7013",  # Dark Orange,
  "Hetero_g2 & Not-Aug" = "#941D07", 
  "Hetero_g2 & Aug"   = "#941D07"  
)

# Labels for the legend
plot_labels_g2 <- c(
  "Homo & Not-Aug" = "Homogeneous", 
  "Homo & Aug"     = "Homogeneous",
  "Hetero_g1 & Not-Aug" = "Heterogeneous (gamma=1)",
  "Hetero_g1 & Aug"   = "Heterogeneous (gamma=1)",
  "Hetero_g2 & Not-Aug" = "Heterogeneous (gamma=3)",
  "Hetero_g2 & Aug"   = "Heterogeneous (gamma=3)"
)

# --- 3. Create the Plot ---

ie_plot_aug_g2 <- sa_ie_plot(sa_ie_data = sa_ie_data_g2[aug == TRUE,],
                          sa_ie_naive_data = sa_ie_naive_data_g2[aug == TRUE,],
                          palette = palette_g2,
                          plot_labels = plot_labels_g2,
                          m_vec_alters = m_vec_alters)
ie_plot_aug_g2 <- ie_plot_aug_g2 + geom_vline(xintercept = 105.5, linetype = 3, 
                                                color = "grey48")
print(ie_plot_aug_g2)
ggsave("HPTN_037/figures/sa_ie_plot_augmented_with_gamma2.png",
       plot = ie_plot_aug_g2, width = 10, height = 5,
       dpi = 300,
       bg = "white")

ie_plot_not_aug_g2 <- sa_ie_plot(sa_ie_data = sa_ie_data_g2[aug == FALSE,],
                              sa_ie_naive_data = sa_ie_naive_data_g2[aug == FALSE,],
                              palette = palette_g2,
                              plot_labels = plot_labels_g2,
                              m_vec_alters = m_vec_alters)
ie_plot_not_aug_g2 <- ie_plot_not_aug_g2 + geom_vline(xintercept = 105.5, linetype = 3, 
                                              color = "grey48")
print(ie_plot_not_aug_g2)
ggsave("HPTN_037/figures/sa_ie_plot_not_augmented_with_gamma2.png",
       plot = ie_plot_not_aug_g2, width = 10, height = 5,
       dpi = 300,
       bg = "white")


# DE plot 
spec_labels_g2 <- c(
  "homo" = "Homogeneous",
  "hetero_gamma1" = "Heterogeneous (gamma=1)",
  "hetero_gamma2" = "Heterogeneous (gamma=2)"
  # "homo" = TeX("Homogeneous $\\pi_i^e$"),
  # "hetero" = TeX("Heterogeneous $\\pi_i^e$")
)
# de_plot <- hptn_sa_aug$de_rd_plot
de_plot_not_aug_g2 <- hptn_sa_not_aug_gamma2$de_rd_plot
de_plot_not_aug_g2 <- de_plot_not_aug_g2 + 
  scale_x_continuous(breaks = seq(0, max(m_vec_egos), 25),
                     labels = seq(0, max(m_vec_egos), 25)) +
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48") +
  labs(x = TeX("$m^e$ (Expected number of missing ego-ego edges)"),
       fill = "Estimated DE",
       title = "Direct Effect") +
  
  facet_wrap(~spec,
             labeller = labeller(spec = spec_labels_g2)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face="bold")
  )
print(de_plot_not_aug_g2)

ggsave("HPTN_037/figures/sa_de_not_aug_plot_with_gamma2.png",
       plot = de_plot_not_aug_g2, width = 10, height = 6,
       dpi = 300,
       bg = "white")

de_plot_aug_g2 <- hptn_sa_aug_gamma2$de_rd_plot
de_plot_aug_g2 <- de_plot_aug_g2 + 
  scale_x_continuous(breaks = seq(0, max(m_vec_egos), 25),
                     labels = seq(0, max(m_vec_egos), 25)) +
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48") +
  labs(x = TeX("$m^e$ (Expected number of missing ego-ego edges)"),
       fill = "Estimated DE",
       title = "Direct Effect") +
  facet_wrap(~spec,
             labeller = labeller(spec = spec_labels_g2)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face="bold")
  )
print(de_plot_aug_g2)

ggsave("HPTN_037/figures/sa_de_aug_plot_with_gamma2.png",
       plot = de_plot_aug_g2, width = 10, height = 6,
       dpi = 300,
       bg = "white")



# DE plot given kappa=1.5
palette_g2 <- c(
  "homo"   = "#085380", # Light Blue
  "hetero_gamma1"   = "#BA7013",  # Dark Orange
  "hetero_gamma2"   = "#941D07"
)

# Labels for the legend
plot_labels_g2 <- c(
  "homo" = "Homogeneous", 
  "hetero_gamma1"= "Heterogeneous (gamma=1)",
  "hetero_gamma2"= "Heterogeneous (gamma=2)"
)


hptn_sa_not_aug_gamma2$sa_results$DE[, pi_param := as.numeric(pi_param)]
hptn_sa_not_aug_gamma2$null_results$DE[, pi_param := as.numeric(pi_param)]
de_plot_kappa1.5_not_aug_g2 <- sa_de_plot_given_kappa(
  sa_de_data = hptn_sa_not_aug_gamma2$sa_results$DE[kappa == 1.5,],
  sa_de_naive_data = hptn_sa_not_aug_gamma2$null_results$DE,
  # sa_de_data = hptn_sa_aug$sa_results$DE[kappa == 1.5,],
  # sa_de_naive_data = hptn_sa_aug$null_results$DE,
  kappa = 1.5,
  palette = palette_g2,
  plot_labels = plot_labels_g2, 
  m_vec_egos = m_vec_egos)
de_plot_kappa1.5_not_aug_g2 <- de_plot_kappa1.5_not_aug_g2 + 
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48")
print(de_plot_kappa1.5_not_aug_g2)
ggsave("HPTN_037/figures/sa_de_plot_not_aug_kappa1.5_with_gamma2.png",
       plot = de_plot_kappa1.5_not_aug_g2, width = 10, height = 5,
       dpi = 300,
       bg = "white")

hptn_sa_aug_gamma2$sa_results$DE[, pi_param := as.numeric(pi_param)]
hptn_sa_aug_gamma2$null_results$DE[, pi_param := as.numeric(pi_param)]
de_plot_kappa1.5_aug_g2 <- sa_de_plot_given_kappa(
  sa_de_data = hptn_sa_aug_gamma2$sa_results$DE[kappa == 1.5,],
  sa_de_naive_data = hptn_sa_aug_gamma2$null_results$DE,
  # sa_de_data = hptn_sa_aug$sa_results$DE[kappa == 1.5,],
  # sa_de_naive_data = hptn_sa_aug$null_results$DE,
  kappa = 1.5,
  palette = palette_g2,
  plot_labels = plot_labels_g2, 
  m_vec_egos = m_vec_egos)
de_plot_kappa1.5_aug_g2 <- de_plot_kappa1.5_aug_g2 + 
  geom_vline(xintercept = 33.4, linetype = 3, 
             color = "grey48")
print(de_plot_kappa1.5_aug_g2)
ggsave("HPTN_037/figures/sa_de_plot_aug_kappa1.5_with_gamma2.png",
       plot = de_plot_kappa1.5_aug_g2, width = 10, height = 5,
       dpi = 300,
       bg = "white")

# Plot given m^e = 35 (approximately the estimate from the validation data)
de_plot_m_e_valid_aug <- sa_de_plot_given_m_e(
  sa_de_data = hptn_sa_aug_gamma2$sa_results$DE[pi_param == 35,],
  sa_de_naive_data = hptn_sa_aug_gamma2$null_results$DE,
  m_e = 35,
  palette = palette_g2,
  plot_labels = plot_labels_g2)
print(de_plot_m_e_valid_aug)
ggsave("HPTN_037/figures/sa_de_plot_aug_m_e_valid_with_gamma2.png",
       plot = de_plot_m_e_valid_aug, width = 10, height = 5,
       dpi = 300,
       bg = "white")

# --- PBA results ---

# Set priors.

# Priors for IE (for m_e)

prior_ie_uniform <- function() {
  sample(seq(1, max(m_vec_alters), by=1), 1)
}

prior_ie_poisson <- function() {
  rpois(1, lambda = max(m_vec_alters)/2)
}

prior_ie_neg_binom <- function() {
  rnbinom(1, size = 10, mu = max(m_vec_alters)/2)
}

# Priors for DE (for m_a and kappa)
min_kappa <- min(kappa_vec)
max_kappa <- max(kappa_vec)
mean_kappa <- mean(kappa_vec)

prior_de_uniform <- function() {
  list(
    pi_param = sample(seq(1, max(m_vec_egos), by=1), 1),
    # kappa = runif(1, min = 1, max = 3)
    kappa = runif(1, min = min_kappa, max = max_kappa)
  )
}

prior_de_nb_lognormal <- function() {
  list(
    pi_param = rnbinom(1, size = 10, mu = max(m_vec_egos)/2),
    # kappa = rlnorm(1, meanlog = log(2), sdlog = 0.2)
    kappa = rlnorm(1, meanlog = log(mean_kappa), sdlog = 0.2)
  )
}

prior_de_poisson_uniform <- function() {
  list(
    pi_param = rpois(1, lambda = max(m_vec_egos)/2),
    # kappa = runif(1, min = 1, max = 3)
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
                          ego_index = ego_id_a_vec,
                          pz=pz)

pi_de_homo_args <- list(n_e = n_e, type = "ego", pz = 0.5)

pi_de_hetero_args <- list(X_e = X_e,
                          gamma = -1,
                          dist = "norm",
                          p = 2,
                          pz = pz)

# --- Run PBA ---
# We run with/without augmented model 

set.seed(142)
pba_homo_unif_aug <-  enrt_pba(Y_e = Y_e, 
                               Y_a = Y_a, 
                               X_e = X_e,
                               X_a = X_a,
                               Z_e = Z_e,
                               F_a = F_a_tilde,
                               ego_id_a = ego_id_a_vec,
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
                           augmented = FALSE,
                           B = B,
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
                           augmented = FALSE,
                           B = B,
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
                               augmented = FALSE,
                               B = B,
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
                               augmented = FALSE,
                               B = B,
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
                               augmented = FALSE,
                               B = B,
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
                               augmented = FALSE,
                               B = B,
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
prior_levels <- c("Naive", "Poisson", "NB", "Uniform")

prior_colors <- c(
  "Uniform" = "#E41A1C", # Red
  "Poisson" = "#377EB8", # Blue
  "NB"      = "#4DAF4A"  # Green
)

model_labels <- c(
  "homo"   = "Homogeneous",
  "hetero" = "Heterogeneous"
)

model_shapes <- c(
  "homo"   = 19,
  "hetero" = 17
)


pba_ie_aug_figure <- pba_ie_plot(
  pba_ie_data = pba_ie_results[augmented == TRUE & uncertainty_type == "total",],
  null_results = hptn_sa_aug$null_results$IE,
  prior_levels = prior_levels,
  prior_colors = prior_colors,
  model_labels = model_labels,
  model_shapes = model_shapes,
  aug = TRUE
)

print(pba_ie_aug_figure)

ggsave("HPTN_037/figures/pba_ie_augmented.png",
       plot = pba_ie_aug_figure, width = 8, height = 5,
       dpi = 300,
       bg = "white")

pba_ie_not_aug_figure <- pba_ie_plot(
  pba_ie_data = pba_ie_results[augmented == FALSE & uncertainty_type == "total",],
  null_results = hptn_sa_not_aug$null_results$IE,
  prior_levels = prior_levels,
  prior_colors = prior_colors,
  model_labels = model_labels,
  model_shapes = model_shapes,
  aug = FALSE
)

print(pba_ie_not_aug_figure)
ggsave("HPTN_037/figures/pba_ie_not_augmented.png",
       plot = pba_ie_not_aug_figure, width = 8, height = 5,
       dpi = 300,
       bg = "white")


# DE plot
prior_levels_de <- c("Naive", "Uniform","Poisson + Uniform", "NB + Lognormal")
prior_colors_de <- c(
  "Poisson + Uniform" = "#E41A1C", # Red
  "NB + Lognormal" = "#377EB8", # Blue
  "Uniform"      = "#4DAF4A"  # Green
)
pba_de_aug_figure <- pba_de_plot(
  pba_de_data = pba_de_results[augmented == TRUE & uncertainty_type == "total",],
  null_results = hptn_sa_aug$null_results$DE,
  prior_levels = prior_levels_de,
  prior_colors = prior_colors_de,
  model_labels = model_labels,
  model_shapes = model_shapes,
  aug = TRUE
)
print(pba_de_aug_figure)
ggsave("HPTN_037/figures/pba_de_augmented.png",
       plot = pba_de_aug_figure, width = 8, height = 5,
       dpi = 300,
       bg = "white")

pba_de_not_aug_figure <- pba_de_plot(
  pba_de_data = pba_de_results[augmented == FALSE & uncertainty_type == "total",],
  null_results = hptn_sa_not_aug$null_results$DE,
  prior_levels = prior_levels_de,
  prior_colors = prior_colors_de,
  model_labels = model_labels,
  model_shapes = model_shapes,
  aug = FALSE
)
print(pba_de_not_aug_figure)
ggsave("HPTN_037/figures/pba_de_not_augmented.png",
       plot = pba_de_not_aug_figure, width = 8, height = 5,
       dpi = 300,
       bg = "white")



