### 
# Read and wrangle HTPN 037 data
###

library(haven)
library(data.table)
library(foreign)
library(stringi)

# Read data ---------------------------------------------------------------

file.dir <- dir("HPTN_037/SAS/")
file.name <- stri_extract(file.dir, regex = "[a-z]+")
path.name <- "HPTN_037/SAS/"

for (i in seq(length(file.dir))) {
  assign(file.name[i], 
         as.data.table(read_sas(paste0(path.name,file.dir[i]))))
}


# Keys for variables:
# site_id = "222" --> Pennsylvanian site
# `uid` --> unique id per subject
# NKID --> Network ID (one per ego-network)
# MRID --> member type in the network. MRID=="0" is ego, and alter otherwise
# `visnum` --> follow-up visit number. `0` for baseline, `6` for 6 months, etc.

# Datasets:
# dm: demographic data (sex, age, race, ethnicity, marriage, schooling,)
# ra: risk assessment (sexual behavior, drug use, etc.) -- the secondary outcomes! 
# il: Intervention Session Participation Log... Probably about the treatment but have to review it


# TODO:
# 1. Define baseline covariates for analysis
# 2. Define treatment variable for ego-networks
# 3. Create tabular data.frame where each row is unit, and each column is a variable
# 4. Create ego, alter, and ego-network indicating variables
# 5. Define the composite outcome of `any risky behavior`. Either 6 or 12 months to follow-up.
# 6. Define the sensitivity probabilities for alter-ego and ego-ego edges.

