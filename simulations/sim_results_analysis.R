###
# Read and analyze the simulation results
###

library(data.table)
library(kableExtra)

# --- Setup: read data and combine files ---
folder_path <- "simulations/sim_results/"

prefixes_map <- c(
  "sa_ie_homo"   = "sa_ie_pi_homo",
  "sa_ie_hetero" = "sa_ie_pi_hetero",
  "sa_de_homo"   = "sa_de_pi_homo",
  "sa_de_hetero" = "sa_de_pi_hetero"
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

ie_hetero_summ <- sa_ie_hetero_summary[m_a %in% c(100, 200),]
# remove 'augmented' column
# ie_hetero_summ_augmented[, augmented := NULL]
# ie_hetero_summ[, Scenario := paste0("$m_a=", m_a, ", IE=", round(true_ie,3),"$")]
ie_hetero_summ[, Scenario := paste0("$m_a=", m_a,"$")]
# ie_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
ie_hetero_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(ie_hetero_summ, c("m_a","spec","augmented"))
# change variable order
ie_hetero_summ <- ie_hetero_summ[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(ie_hetero_summ, 
      caption = "IE Results - Heterogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(ie_hetero_summ, 
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


de_hetero_summ <- sa_de_hetero_summary[m_e %in% c(150, 250),]
# de_hetero_summ[, augmented := NULL]
de_hetero_summ[, Scenario := paste0("$m_e=", m_e,"$")]
# de_hetero_summ_augmented[, Variance := ifelse(bootstrap, "Bootstrap", "Analytic")]
de_hetero_summ[, spec := ifelse(spec == "homo", "Homogeneous",
                                          ifelse(spec == "hetero",
                                                 "Heterogeneous", spec))]
setorderv(de_hetero_summ, c("m_e","spec","augmented"))
# change variable order
de_hetero_summ <- de_hetero_summ[,
                                     .(Scenario, spec, augmented, bias, coverage, ese_ase)
                                     ]
kable(de_hetero_summ, 
      caption = "DE Results - Heterogeneous Exposure", 
      digits = 3) %>%
  kable_styling(full_width = FALSE)


kbl(de_hetero_summ, 
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


