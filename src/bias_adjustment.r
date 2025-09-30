

# Source + Packages --------------------------------------------------------------

source("src/sensitivity_params.r")



# Indirect effects --------------------------------------------------------

IE_adjustment <- function(...){
  # TODO:
  # Take as input the relevant data and specification
  # return list of lists (or something similar) of bias corrected estimates
  # Should include the point estimates and the bootstrap CIs
  # Should work for all setups of SA params
  # Should work for both SA and PBA
  # Consist of lower level funcs:
  # 1. Function that just estimate bias correction given one pi vector
  # 2. For multiple pi vectors
  # 3. Bootstrap across units to repeat point 2
  # 4. Get point estimate + bootstrap CIs/
  # 5. similar logic for PBA or SA (average or grid over SA params)
  # Should work for both RD and RR estimands. 
}

# Direct effects ----------------------------------------------------------


