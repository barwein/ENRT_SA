###
# Read and analyze the simulation results
###


# Packages ---------------------------------------------------------------

library(data.table)
library(kableExtra)


# Correct specified settings (continuous) ---------------------------------

# --- Setup: read data and combine files ---
folder_path <- "simulations/sim_results/"

prefixes_map <- c(
  "sa_ie_homo"   = "sa_ie_pi_homo_ma",
  "sa_ie_hetero" = "sa_ie_pi_hetero_ma",
  "sa_de_homo"   = "sa_de_pi_homo_ma",
  "sa_de_hetero" = "sa_de_pi_hetero_ma"
)

all_data_tables_list <- lapply(prefixes_map, function(file_prefix) {
  
  pattern_to_search <- paste0("^", file_prefix, ".*\\.csv$")
  file_list <- list.files(
    path = folder_path,
    pattern = pattern_to_search,
    full.names = TRUE
  )

  if (length(file_list) == 0) {
    # If no files are found, print a warning and return an empty data.table
    # This prevents errors from trying to read an empty list
    warning("No files found with prefix: '", file_prefix, "' in folder: '", folder_path, "'")
    return(data.table())
  }
  
  # Read all matching files into a list of data.tables
  list_of_dts <- lapply(file_list, fread)
  
  # Combine the list of data.tables into one
  combined_dt <- rbindlist(list_of_dts, use.names = TRUE, fill = TRUE)
  return(combined_dt)
})

list2env(all_data_tables_list, envir = .GlobalEnv)


# --- Analysis: summarize results ---

sa_ie_homo_summary <- sa_ie_homo[,.(
          bias = mean(ie_rd - true_ie),
          # rel_bias = mean((ie_rd - true_ie)/true_ie),
          coverage = mean((ci_low <= true_ie) & (true_ie <= ci_high)),
          bias_sd = sd(ie_rd),
          mean_se = mean(sqrt(var_to_use)),
          ese_ase = sd(ie_rd) / mean(sqrt(var_to_use)),
          true_ie = mean(true_ie)
        ),
        by = c("m_a","spec", "augmented")]
        

sa_ie_hetero_summary <- sa_ie_hetero[,.(
                bias = mean(ie_rd - true_ie),
                # rel_bias = mean((ie_rd - true_ie)/true_ie),
                coverage = mean((ci_low <= true_ie) & (true_ie <= ci_high)),
                bias_sd = sd(ie_rd),
                mean_se = mean(sqrt(var_to_use)),
                ese_ase = sd(ie_rd) / mean(sqrt(var_to_use)),
                true_ie = mean(true_ie)
              ),
                by = c("m_a","spec", "augmented")]


sa_de_homo_summary <- sa_de_homo[,.(
          bias = mean(de_rd - true_de),
          # rel_bias = mean((de_rd - true_de)/true_de),
          coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
          bias_sd = sd(de_rd),
          mean_se = mean(sqrt(var_to_use)),
          ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
          true_de = mean(true_de)
        ),
        by = c("m_e","spec", "augmented", "kappa_")]
        

sa_de_hetero_summary <- sa_de_hetero[,.(
            bias = mean(de_rd - true_de),
            # rel_bias = mean((de_rd - true_de)/true_de),
            coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
            bias_sd = sd(de_rd),
            mean_se = mean(sqrt(var_to_use)),
            ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
            true_de = mean(true_de)
          ),
          by = c("m_e","spec", "augmented", "kappa_")]



# --- Main Results: Heterogeneous true contamination + Augmented estimators ---

# remove 'augmented' column
# ie_hetero_summ_augmented[, augmented := NULL]
# ie_hetero_summ[, Scenario := paste0("$m_a=", m_a, ", IE=", round(true_ie,3),"$")]
sa_ie_hetero_summary[, Scenario := paste0("$m_a=", m_a,"$")]
# ie_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
sa_ie_hetero_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(sa_ie_hetero_summary, c("m_a","spec","augmented"))
# change variable order
sa_ie_hetero_summary <- sa_ie_hetero_summary[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(sa_ie_hetero_summary, 
      caption = "IE Results - Heterogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_ie_hetero_summary, 
    caption = "IE Results - Heterogeneous Exposure", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lcccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )


sa_de_hetero_summary[, Scenario := paste0("$m_e=", m_e,"$")]
# de_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
sa_de_hetero_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(sa_de_hetero_summary, c("m_e","spec","augmented"))
# change variable order
sa_de_hetero_summary <- sa_de_hetero_summary[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(sa_de_hetero_summary, 
      caption = "DE Results - Heterogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_de_hetero_summary, 
    caption = "DE Results - Heterogeneous Exposure - Augmented", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )

# --- Appendix tables ---

#  Homogeneous true contamination 

ie_homo_summ <- sa_ie_homo_summary[m_a %in% c(100, 200),]
# remove 'augmented' column
# ie_homo_summ_augmented[, augmented := NULL]
ie_homo_summ[, Scenario := paste0("$m_a=", m_a,"$")]
# ie_homo_summ[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
ie_homo_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(ie_homo_summ, c("m_a","spec","augmented"))
# change variable order
ie_homo_summ <- ie_homo_summ[,.(Scenario, spec, augmented, bias, coverage, ese_ase)]
kable(ie_homo_summ, 
      caption = "IE Results - homogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(ie_homo_summ, 
    caption = "IE Results - homogeneous Exposure", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )


de_homo_summ <- sa_de_homo_summary[m_e %in% c(150, 250),]
# de_homo_summ[, augmented := NULL]
de_homo_summ[, Scenario := paste0("$m_e=", m_e,"$")]
# de_homo_summ[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
de_homo_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(de_homo_summ, c("m_e","spec","augmented"))
# change variable order
de_homo_summ <- de_homo_summ[,.(Scenario, spec, augmented, bias, coverage, ese_ase)]
kable(de_homo_summ, 
      caption = "DE Results - homogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(de_homo_summ, 
    caption = "DE Results - homogeneous Exposure", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )





# Misspecified kappa (continuous) ---------------------------------

# --- Setup: read data and combine files ---
folder_path <- "simulations/sim_results/"

sa_de_hetero_misspec_kappa05 <- fread("simulations/sim_results/sa_de_pi_hetero_truekappa1.5_usedkappa0.5_ma100_m_e150_niter5000.csv")
sa_de_hetero_misspec_kappa05[,kappa := 0.5]

sa_de_hetero_misspec_kappa11 <- fread("simulations/sim_results/sa_de_pi_hetero_truekappa1.5_usedkappa1.1_ma100_m_e150_niter5000.csv")
sa_de_hetero_misspec_kappa11[,kappa := 1.1]

sa_de_hetero_misspec_kappa <- rbindlist(list(
  sa_de_hetero_misspec_kappa11,
  sa_de_hetero_misspec_kappa05
))

# --- Analysis: summarize results ---

sa_de_homo_misspec_kappa_summary <- sa_de_hetero_misspec_kappa[,.(
          bias = mean(de_rd - true_de),
          # rel_bias = mean((de_rd - true_de)/true_de),
          coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
          bias_sd = sd(de_rd),
          mean_se = mean(sqrt(var_to_use)),
          ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
          true_de = mean(true_de)
        ),
        by = c("kappa","spec", "augmented")]
        

sa_de_homo_misspec_kappa_summary[, Scenario := paste0("Using $kappa=", kappa,"$")]
sa_de_homo_misspec_kappa_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(sa_de_homo_misspec_kappa_summary, c("kappa","spec","augmented"))
# change variable order
sa_de_homo_misspec_kappa_summary <- sa_de_homo_misspec_kappa_summary[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(sa_de_homo_misspec_kappa_summary, 
      caption = "DE Results - Heterogeneous Exposure with misspecified kappa. True DE = 2. True kappa=1.5", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_de_homo_misspec_kappa_summary, 
    caption = "DE Results - Heterogeneous Exposure with misspecified kappa. True DE = 2. True kappa=1.5", , 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lcccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )






# Violation of ass.4 (continuous) ---------------------------------

# --- Setup: read data and combine files ---
folder_path <- "simulations/sim_results/"

sa_de_hetero_viol_15 <- fread("simulations/sim_results/sa_de_pi_hetero_violass4_inter15_ma100_m_e150_niter5000.csv")
sa_de_hetero_viol_15[,inter_param := 1.5]
sa_de_hetero_viol_05 <- fread("simulations/sim_results/sa_de_pi_hetero_violass4_inter05_ma100_m_e150_niter5000.csv")
sa_de_hetero_viol_05[,inter_param := 0.5]

sa_de_hetero_viol <- rbindlist(list(
  sa_de_hetero_viol_15,
  sa_de_hetero_viol_05
))

# --- Analysis: summarize results ---

sa_de_hetero_viol_summary <- sa_de_hetero_viol[,.(
          bias = mean(de_rd - true_de),
          # rel_bias = mean((de_rd - true_de)/true_de),
          coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
          bias_sd = sd(de_rd),
          mean_se = mean(sqrt(var_to_use)),
          ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
          true_de = mean(true_de)
        ),
        by = c("kappa_","spec", "augmented","inter_param")]
        

sa_de_hetero_viol_summary[, Scenario := paste0("$beta_4=",inter_param,"$; ","$kappa=", kappa_,"$")]
sa_de_hetero_viol_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(sa_de_hetero_viol_summary, c("kappa_","spec","augmented"))
# change variable order
sa_de_hetero_viol_summary <- sa_de_hetero_viol_summary[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(sa_de_hetero_viol_summary, 
      caption = "DE Results - Heterogeneous Exposure with violation of ass.4. True DE = 2 and m_e=150", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_de_hetero_viol_summary, 
    caption = "DE Results - Heterogeneous Exposure with violation of ass.4. True DE = 2 and m_e=150", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lcccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )






# Correct specified settings (binary) ---------------------------------

# --- Setup: read data and combine files ---
folder_path <- "simulations/sim_results/"

prefixes_map <- c(
  "sa_ie_homo"   = "sa_ie_pi_homo_binary",
  "sa_ie_hetero" = "sa_ie_pi_hetero_binary",
  "sa_de_homo"   = "sa_de_pi_homo_binary",
  "sa_de_hetero" = "sa_de_pi_hetero_binary"
)

all_data_tables_list <- lapply(prefixes_map, function(file_prefix) {
  
  pattern_to_search <- paste0("^", file_prefix, ".*\\.csv$")
  file_list <- list.files(
    path = folder_path,
    pattern = pattern_to_search,
    full.names = TRUE
  )
  
  if (length(file_list) == 0) {
    # If no files are found, print a warning and return an empty data.table
    # This prevents errors from trying to read an empty list
    warning("No files found with prefix: '", file_prefix, "' in folder: '", folder_path, "'")
    return(data.table())
  }
  
  # Read all matching files into a list of data.tables
  list_of_dts <- lapply(file_list, fread)
  
  # Combine the list of data.tables into one
  combined_dt <- rbindlist(list_of_dts, use.names = TRUE, fill = TRUE)
  return(combined_dt)
})

list2env(all_data_tables_list, envir = .GlobalEnv)


# --- Analysis: summarize results ---

sa_ie_homo_summary <- sa_ie_homo[,.(
  bias = mean(ie_rd - true_ie),
  # rel_bias = mean((ie_rd - true_ie)/true_ie),
  coverage = mean((ci_low <= true_ie) & (true_ie <= ci_high)),
  bias_sd = sd(ie_rd),
  mean_se = mean(sqrt(var_to_use)),
  ese_ase = sd(ie_rd) / mean(sqrt(var_to_use)),
  true_ie = mean(true_ie)
),
by = c("m_a","spec", "augmented")]


sa_ie_hetero_summary <- sa_ie_hetero[,.(
  bias = mean(ie_rd - true_ie),
  # rel_bias = mean((ie_rd - true_ie)/true_ie),
  coverage = mean((ci_low <= true_ie) & (true_ie <= ci_high)),
  bias_sd = sd(ie_rd),
  mean_se = mean(sqrt(var_to_use)),
  ese_ase = sd(ie_rd) / mean(sqrt(var_to_use)),
  true_ie = mean(true_ie)
),
by = c("m_a","spec", "augmented")]


sa_de_homo_summary <- sa_de_homo[,.(
  bias = mean(de_rd - true_de),
  # rel_bias = mean((de_rd - true_de)/true_de),
  coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
  bias_sd = sd(de_rd),
  mean_se = mean(sqrt(var_to_use)),
  ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
  true_de = mean(true_de)
),
by = c("m_e","spec", "augmented", "kappa_")]


sa_de_hetero_summary <- sa_de_hetero[,.(
  bias = mean(de_rd - true_de),
  # rel_bias = mean((de_rd - true_de)/true_de),
  coverage = mean((ci_low <= true_de) & (true_de <= ci_high)),
  bias_sd = sd(de_rd),
  mean_se = mean(sqrt(var_to_use)),
  ese_ase = sd(de_rd) / mean(sqrt(var_to_use)),
  true_de = mean(true_de)
),
by = c("m_e","spec", "augmented", "kappa_")]



# --- Main Results: IE ---

sa_ie_hetero_summary[, Scenario := paste0("$m_a=", m_a,"$")]
# ie_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
sa_ie_hetero_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                      ifelse(spec == "hetero",
                                             "Heterogeneous", spec))]
setorderv(sa_ie_hetero_summary, c("m_a","spec","augmented"))
# change variable order
sa_ie_hetero_summary <- sa_ie_hetero_summary[,
                                             .(Scenario, spec, augmented, bias, coverage, ese_ase)
]
kable(sa_ie_hetero_summary, 
      caption = "IE Results - Heterogeneous Exposure with binary PO", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_ie_hetero_summary, 
    caption = "IE Results - Heterogeneous Exposure with binary PO", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lcccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )


ie_homo_summ <- sa_ie_homo_summary
# remove 'augmented' column
# ie_homo_summ_augmented[, augmented := NULL]
ie_homo_summ[, Scenario := paste0("$m_a=", m_a,"$")]
# ie_homo_summ[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
ie_homo_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                              ifelse(spec == "hetero",
                                     "Heterogeneous", spec))]
setorderv(ie_homo_summ, c("m_a","spec","augmented"))
# change variable order
ie_homo_summ <- ie_homo_summ[,.(Scenario, spec, augmented, bias, coverage, ese_ase)]
kable(ie_homo_summ, 
      caption = "IE Results - homogeneous Exposure with binary PO", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(ie_homo_summ, 
    caption = "IE Results - homogeneous Exposure with binary PO", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )


# --- Main Results: DE ---

sa_de_hetero_summary[, Scenario := paste0("$m_e=", m_e,"$")]
# de_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
sa_de_hetero_summary[, spec := ifelse(spec == "homo", "Homogeneous",
                                      ifelse(spec == "hetero",
                                             "Heterogeneous", spec))]
setorderv(sa_de_hetero_summary, c("m_e","spec","augmented"))
# change variable order
sa_de_hetero_summary <- sa_de_hetero_summary[,
                                             .(Scenario, spec, augmented, bias, coverage, ese_ase)
]
kable(sa_de_hetero_summary, 
      caption = "DE Results - Heterogeneous Exposure with binary PO", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(sa_de_hetero_summary, 
    caption = "DE Results - Heterogeneous Exposure with binary PO", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )

de_homo_summ <- sa_de_homo_summary[m_e %in% c(150, 250),]
# de_homo_summ[, augmented := NULL]
de_homo_summ[, Scenario := paste0("$m_e=", m_e,"$")]
# de_homo_summ[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
de_homo_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                              ifelse(spec == "hetero",
                                     "Heterogeneous", spec))]
setorderv(de_homo_summ, c("m_e","spec","augmented"))
# change variable order
de_homo_summ <- de_homo_summ[,.(Scenario, spec, augmented, bias, coverage, ese_ase)]
kable(de_homo_summ, 
      caption = "DE Results - homogeneous Exposure with binary PO", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(de_homo_summ, 
    caption = "DE Results - homogeneous Exposure with binary PO", 
    digits = 3,
    format = "latex",
    booktabs = TRUE,
    align = "lccccc",
    linesep = "")  %>%
  collapse_rows(
    columns = 1:2,      
    valign = "middle"   
  )




