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

PEN_SITE <- "222" # Pennsylvanian site

# Keys for variables:
# site_id = "222" --> Pennsylvanian site
# `uid` --> unique id per subject
# NKID --> Network ID (one per ego-network)
# MRID --> member type in the network. MRID=="0" is ego, and alter otherwise
# `visnum` --> follow-up visit number. `0` for baseline, `6` for 6 months, etc.


# Datasets:
# dm: demographic data (sex, age, race, ethnicity, marriage, schooling,)
# ra: risk assessment (sexual behavior, drug use, etc.) -- the secondary outcomes!
# ns: network summary
# il: Intervention Session Participation Log (session attendance; helps infer treatment)
# keys: key variables per participant (often contains arm assignment if available)

# Keep only Pennsylvania site + prepare visit windows -------------------------------------------------------------------------

dm_pen <- dm[site_id == PEN_SITE]
dm_pen[, `:=`(
  is_ego = MRID == "0",
  is_male = DM1sex == "Male",
  uid = as.character(uid),
  NKID = as.character(NKID)
)]

# Risk assessment — coerce core fields
ra_pen <- ra[site_id == PEN_SITE]
ra_pen[, `:=`(
  visnum = as.numeric(visnum),
  studyday = as.numeric(studyday),
  visit = as.character(visit),
  uid = as.character(uid)
)]

# 6‑month follow‑up: take the earliest completed assessment per uid
ra_pen_6m <- ra_pen[
  visnum == 6 & !is.na(studyday),
  .SD[which.min(studyday)],
  by = uid
]


# Restrict other tables to those with 6m outcomes
uid_6m <- unique(ra_pen_6m$uid)
ra_pen_sub <- ra_pen[uid %in% uid_6m]
dm_sub <- dm_pen[uid %in% uid_6m]


 
# Outcome: injection risk behaviors at 6 months (binary indicators + union) -------------------------------------------------------------------------

# Helper: treat blanks as NA
blank_to_na <- function(x) {x[x %in% c("")] <- NA; x}

ra_pen_6m[, c("RA3nwatr","RA3ncook","RA3ncotn","RA3nload","RA3nndl",
              "RA4othrs","RA4dkwel","RA3injlm","RA3inj6m",
              "RA3ndays") := lapply(.SD, blank_to_na),
          .SDcols = c("RA3nwatr","RA3ncook","RA3ncotn","RA3nload","RA3nndl",
                      "RA4othrs","RA4dkwel","RA3injlm","RA3inj6m","RA3ndays")]


# Numeric sharing counts (>0 = risky)
num_risk <- c("RA3nwatr","RA3ncook","RA3ncotn","RA3nload","RA3nndl")
ra_pen_6m[, (paste0(num_risk,"_bin")) := lapply(.SD,
                                                function(v) as.integer(as.numeric(v) > 0)),
          .SDcols = num_risk]


# Public/shooting gallery (Yes/No)
ra_pen_6m[, RA4othrs_bin := as.integer(RA4othrs == "Yes")]


# Injecting with people you don't know well (anything except "rarely or never")
ra_pen_6m[, RA4dkwel_bin := as.integer(!RA4dkwel %in% c(NA, "rarely or never"))]


# If no injection in last month (RA3injlm == 2 or "No") or last 6 months, set relevant risk to 0
no_inj_month <- ra_pen_6m$RA3injlm %in% c("2","No")
no_inj_6mo <- ra_pen_6m$RA3inj6m %in% c("0","No")


risk_cols <- c(paste0(num_risk,"_bin"), "RA4othrs_bin","RA4dkwel_bin")
# treat NA as zeros
for (cc in risk_cols) {
  ra_pen_6m[no_inj_month | no_inj_6mo, (cc) := 0L]
  ra_pen_6m[is.na(get(cc)), (cc) := 0L]
}


# Composite: any risky injection behavior at 6m
ra_pen_6m[, injrisk_any_6m := as.integer(
  do.call(pmax, c(.SD, list(na.rm = TRUE)))
), .SDcols = risk_cols]


# Also mirror SAS "Injected >14 days last month" (inject14)
ra_pen_6m[, inject14_6m := as.integer(as.numeric(RA3ndays) > 14)]


# Baseline covariates (visnum = 0), one record per uid -------------------------------------------------------------------------

ra_base <- ra_pen[visnum == 0 & !is.na(studyday), .SD[which.min(studyday)], by = uid]
ra_base[, c("RA3nwatr","RA3ncook","RA3ncotn","RA3nload","RA3nndl",
            "RA4othrs","RA4dkwel","RA3injlm","RA3inj6m",
            "RA1drunk","RA3ndays","RA3hero","RA3cocai") :=
          lapply(.SD, blank_to_na),
        .SDcols = c("RA3nwatr","RA3ncook","RA3ncotn","RA3nload","RA3nndl",
                    "RA4othrs","RA4dkwel","RA3injlm","RA3inj6m",
                    "RA1drunk","RA3ndays","RA3hero","RA3cocai")]


# Make RA3 substance flags 0 when no recent injection, per SAS logic
ra_base[RA3inj6m %in% c("0","No") | RA3injlm %in% c("0","2","No"), `:=`(
  RA3hero = "0", RA3cocai = "0"
)]


# Baseline risk union (same rules as 6m)
ra_base[, (paste0(num_risk,"_bin")) := lapply(.SD,
                                              function(v) as.integer(as.numeric(v) > 0)),
        .SDcols = num_risk]

ra_base[, RA4othrs_bin := as.integer(RA4othrs == "Yes")]
ra_base[, RA4dkwel_bin := as.integer(!RA4dkwel %in% c(NA, "rarely or never"))]
no_inj_month_b <- ra_base$RA3injlm %in% c("2","No","0")
no_inj_6mo_b <- ra_base$RA3inj6m %in% c("0","No")
# NA as zeros
for (cc in c(paste0(num_risk,"_bin"),"RA4othrs_bin","RA4dkwel_bin")) {
  ra_base[no_inj_month_b | no_inj_6mo_b, (cc) := 0L]
  ra_base[is.na(get(cc)), (cc) := 0L]
}
# Aggregate any baseline risk
ra_base[, injrisk_any_base := as.integer(do.call(pmax, c(.SD, list(na.rm=TRUE)))),
        .SDcols = c(paste0(num_risk,"_bin"),"RA4othrs_bin","RA4dkwel_bin")]


# Daily/frequent injection at baseline (>14 days last month)
ra_base[, inject14_base := as.integer(as.numeric(RA3ndays) > 14)]
ra_base[is.na(inject14_base), inject14_base := 0L] # treat NA as 0


# Alcohol ("got drunk") at baseline
ra_base[, drunk_base := as.integer(RA1drunk != "rarely or never")]
ra_base[is.na(drunk_base), drunk_base := 0L] # treat NA as 0

# Baseline dual use (injected heroin AND cocaine)
ra_base[, heroin_and_cocaine_base := 
          as.integer(RA3hero == "Yes" & RA3cocai == "Yes")]

# cocaine base
ra_base[, cocaine_base := as.integer(RA3cocai == "Yes")]

# Demographics + per‑network indicators -------------------------------------------------------------------------

dm_base <- unique(dm_sub[, .(uid = as.character(uid),
                             NKID = as.character(NKID),
                             is_ego,
                             is_male = is_male,
                             age = suppressWarnings(as.numeric(age)),
                             nonwhite = as.integer(DM1white != "X"),
                             hispanic = as.integer(DM1latin == "Yes"))])


# Network‑level aggregates at baseline
net_agg <- merge(dm_base[, .(uid, NKID, age, nonwhite, is_male)],
                 ra_base[, .(uid, injrisk_any_base, cocaine_base,
                             heroin_and_cocaine_base)], by = "uid", all.x = TRUE)


net_agg <- net_agg[, .(
  net_avg_age = mean(age, na.rm = TRUE),
  net_prev_nonwhite = mean(nonwhite, na.rm = TRUE),
  net_prop_male = mean(is_male, na.rm = TRUE),
  net_prev_injrisk_any_base = mean(injrisk_any_base, na.rm = TRUE),
  net_prev_cocaine_base = mean(cocaine_base, na.rm = TRUE),
  net_prev_heroin_and_cocaine = mean(heroin_and_cocaine_base, na.rm = TRUE)
), by = NKID]



# Treatment assignment for ego‑networks -------------------------------------------------------------------------

# Preferred: use explicit treatment from `keys` if available (1=treatment, 2=control)
arm_by_uid <- NULL
if (exists("keys")) {
  keys_pen <- as.data.table(keys)[site_id == PEN_SITE]
  if ("treat" %in% names(keys_pen)) {
    arm_by_uid <- keys_pen[, .(uid = as.character(uid), treat = as.integer(treat))]
  }
}


# Fallback: infer from Intervention Session Participation Log (IL-1)
if (is.null(arm_by_uid)) {
  if (exists("il")) {
    il_pen <- as.data.table(il)[site_id == PEN_SITE]
    # gather participant ID columns (ILP1ID..ILP12ID)
    pid_cols <- grep("^ILP[0-9]+ID$", names(il_pen), value = TRUE)
    if (length(pid_cols)) {
      il_long <- melt(il_pen, measure.vars = pid_cols, value.name = "uid", variable.name = "slot")
      il_long <- il_long[!is.na(uid) & uid != "", .(uid = as.character(uid))]
      # Map to NKID using dm_base and flag any network with attendance as treated
      il_long <- unique(il_long[dm_base, on = .(uid), nomatch = 0L])
      treat_by_nkid <- il_long[, .(treat = 1L), by = NKID]
      arm_by_uid <- treat_by_nkid[dm_base, on = .(NKID),
                                   .(uid, treat = fifelse(is.na(treat), 0L, treat))]
    }
  }
}


if (is.null(arm_by_uid)) {
  # If still missing, set to NA for transparency
  arm_by_uid <- dm_base[, .(uid, treat = NA_integer_)]
}


# `ego_treat`: network‑level assignment applied to all members
# If keys$treat used: treat==1 is intervention, otherwise control (0)
arm_by_uid[, ego_treat := fifelse(is.na(treat), NA_integer_, as.integer(treat == 1))]
arm_by_uid[, treat := NULL]



# Build final analysis table (one row per uid with 6m outcome) -------------------------------------------------------------------------

 
# Build final analysis table — merge all uid-based tables first, then add net_agg by NKID
analysis_uid <- Reduce(function(x, y) merge(x, y, by = "uid", all.x = TRUE), list(
  ra_pen_6m[, .(uid,
                injrisk_any_6m,
                inject14_6m)],
  dm_base,
  ra_base[, .(uid,
              injrisk_any_base,
              inject14_base,
              drunk_base,
              heroin_and_cocaine_base)],
  arm_by_uid
))


# Now bring in network-level features by NKID
analysis_dt <- merge(analysis_uid, net_agg, by = "NKID", all.x = TRUE)


# Remove alters without ego in the data
ego_nkids <- unique(analysis_dt[is_ego == TRUE]$NKID)
analysis_dt <- analysis_dt[NKID %in% ego_nkids]


# Order columns nicely
setcolorder(analysis_dt, c("uid","NKID","is_ego","ego_treat",
                           "injrisk_any_6m","inject14_6m",
                           "injrisk_any_base","inject14_base","drunk_base","heroin_and_cocaine_base",
                           "is_male","age","nonwhite","hispanic",
                           "net_avg_age","net_prev_nonwhite", "net_prop_male",
                           "net_prev_injrisk_any_base",
                           "net_prev_cocaine_base","net_prev_heroin_and_cocaine"))


# Print number of NAs per column
data.frame(apply(analysis_dt, 2, function(x){sum(is.na(x))}))

# Save as csv
write.csv(analysis_dt, "HPTN_037/hptn_clean_df.csv", row.names = FALSE)

