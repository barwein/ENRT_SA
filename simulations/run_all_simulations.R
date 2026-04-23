### 
# Run all simulations setup #
# All seeds are within the sourced scripts
###

# 1. Homogeneous exposure probabilities (no violations, continuous outcomes)
# m_a = 100, m_e = 150
source("simulations/sim_run_setup_homo_ma100_me150.R")

# m_a = 200, m_e = 250
source("simulations/sim_run_setup_homo_ma200_me250.R")

# 2. Heterogeneous exposure probabilities (no violations, continuous outcomes)
# m_a = 100, m_e = 150
source("simulations/sim_run_setup_hetero_ma100_me150.R")

# m_a = 200, m_e = 250
source("simulations/sim_run_setup_hetero_ma200_me250.R")


# 3. Heterogeneous exposure probabilities + violation of Assumption 4 (continuous outcomes)

# Interaction param = 0.5
source("simulations/sim_run_setup_hetero_violass4_inter05_ma100_me150.R")

# Interaction param = 1.5
source("simulations/sim_run_setup_hetero_violass4_inter15_ma100_me150.R")


# 4. Heterogeneous exposure probabilities + misspecification of kappa (continuous outcomes)

# Used kappa = 1.1, true kappa = 1.5
source("simulations/sim_run_setup_hetero_misspec_kappa11_ma100_me150.R")

# Used kappa = 0.5, true kappa = 1.5
source("simulations/sim_run_setup_hetero_misspec_kappa05_ma100_me150.R")


# 5. Binary outcomes

# Heterogeneous probs
source("simulations/sim_run_setup_hetero_binarypo_ma100_me150.R")

# Homogeneous probs
source("simulations/sim_run_setup_homo_binarypo_ma100_me150.R")









