


#  DATA PREPROCESSING
#    First, parse the date fields
#    Then, construct the three clocks
#    Construct observed time and event indicator
#    Get covariate matrix and screening
#    Construct analysis dataframes

vax_col       <- "VAC1"  

ice_vars      <- c("de_2080","de_2050","de_1040","de_1070","dinf1618") # List of ICE date columns

A_col         <- "Treat"     # Treatment indicator (1 = Vax, 0 = Control)
event_ind_col <- "Y"
tpos_days_col <- "EVDVC1"    # Outcome date column

tfu_days_col  <- "FUTVC1"    # End of follow-up date column

base_covars   <- c("AGE", "CTN", "RACE", "HORM_FLG", "CYTOD0", "CHLAD0", "NB_PAR", "STD_FLG","SMK","NB_PAR_I","CHLAD0_I", "STD_FLG_I","SMK_I") # Baseline covariates for matching/adjustment



# Dates 
# Vaccination date arrives as DDMONYYYY (e.g. "01MAY2004").
dat <- dat %>%
  mutate(

    VAC1_chr = toupper(trimws(as.character(.data[[vax_col]]))),   # normalise case/whitespace
    vax_day  = as.integer(str_extract(VAC1_chr, "^\\d{1,2}")),    # 1-2 digits at the start
    vax_monT = str_match(VAC1_chr, "^\\d{1,2}([A-Z]{3,})\\d{2,4}$")[, 2],  # the letters between
    vax_year = as.integer(str_extract(VAC1_chr, "\\d{2,4}$")),    # digits at the end
    vax_mon  = match(vax_monT, c("JAN","FEB","MAR","APR","MAY","JUN",
                                 "JUL","AUG","SEP","OCT","NOV","DEC")),  # NA if unrecognised
    vax_date = as.Date(sprintf("%04d-%02d-%02d", vax_year, vax_mon, vax_day)),
    
    # The ICE and positive-test dates are already ISO, so they only need parsing.
    across(all_of(ice_vars), ~ as.Date(.x)),
    # Earliest ICE of any type: only the first one enters the analysis.
    # pmin can drop the Date class, hence the as.Date around it.
    ice_date = as.Date(do.call(pmin, c(across(all_of(ice_vars)), na.rm = TRUE))),
    pos_date = as.Date(.data[[tpos_days_col]]),
    
    # Days from vaccination, with Inf where the event never happened, so the
    # pmin() that builds the observed time treats "never" as something that
    # can never be the minimum rather than letting NA propagate.
    # A Date difference is a difftime, so as.numeric() makes it a plain count.
    t_ice = ifelse(is.na(ice_date), Inf, as.numeric(ice_date - vax_date)),
    t_pos = ifelse(is.na(pos_date), Inf, as.numeric(pos_date - vax_date)),
    # tfu_days_col is already a day count, not a date
    t_fu  = as.numeric(.data[[tfu_days_col]])
  )

# The three clocks 
t_pos_cc <- dat$t_pos                     # time to CIN2+   (Inf if never)
t_ice_cc <- dat$t_ice                     # time to first ICE (Inf if none)
t_fu_cc  <- dat$t_fu                      # administrative end of follow-up
# Censoring time: whichever of the ICE and the end of follow-up comes first.
# This is also the capped ICE clock the PLOT estimators use.
t_ice_capped <- pmin(t_ice_cc, t_fu_cc)

# Observed time and event
# Observed time is the earlier of the event and the censoring time. The event
# counts only if it was actually recorded and happened at or before censoring;
# a tie counts as an event.
raw_event_ind <- as.numeric(as.character(dat[[event_ind_col]]))
Y_obs     <- pmin(t_pos_cc, t_ice_capped)
Delta_obs <- as.numeric(raw_event_ind == 1 & t_pos_cc <= t_ice_capped)
A_treat   <- as.numeric(as.character(dat[[A_col]]))   # 0 = placebo, 1 = vaccine

#  Covariates 
W_numeric <- model.matrix(~ . - 1, data = as.data.frame(dat[, base_covars, drop = FALSE]))
# Near-zero-variance screening at caret's defaults: a column goes when its
# most frequent value occurs more than 19 times as often as its second.
bad_cols_idx <- caret::nearZeroVar(W_numeric)
W_clean <- if (length(bad_cols_idx)) W_numeric[, -bad_cols_idx, drop = FALSE] else W_numeric
# Collinearity check: largest over smallest singular value after centring and
# scaling. Below 30 indicates no problem.
condition_number <- kappa(scale(W_clean), exact = TRUE)

W_covars_clean <- as.data.frame(W_clean)
colnames(W_covars_clean) <- make.names(colnames(W_covars_clean), unique = TRUE)
# The chlamydia main effect was screened out, so its missingness indicator
# goes with it.
W_covars_clean$CHLAD0_I1 <- NULL

# Analysis frames 
valid_time_idx <- Y_obs > 0        # drop anyone with no positive observation time
raw_keep       <- dat[valid_time_idx, ]   # necessary at a later point

analysis_df <- cbind(data.frame(Event        = Delta_obs,
                                Treat        = A_treat,
                                PersonTime   = Y_obs,
                                t_pos        = t_pos_cc,
                                t_ice_capped = t_ice_capped),
                     W_covars_clean)[valid_time_idx, ]
rownames(analysis_df) <- NULL

# The same rows as vectors: CFsurvival takes time, event, treat and
# confounders separately rather than as a data frame.
Y_surv     <- analysis_df$PersonTime
Delta_surv <- analysis_df$Event
A_surv     <- analysis_df$Treat
W_surv     <- analysis_df[, names(W_covars_clean), drop = FALSE]

# alglm does take a frame, under its own column names.
df <- data.frame(Y_time = Y_surv, Delta_event = Delta_surv, A_treat = A_surv, W_surv)
