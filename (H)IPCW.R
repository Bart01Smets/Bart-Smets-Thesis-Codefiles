#  IPCW AND HIPCW analysis (Part I: IPCW, Part II: HIPCW, Part III: Appendix)
RAW <- dat

# PART I: IPCW POISSON
##  Who was censored specifically by an intercurrent event
# Four conditions: an ICE occurred, at or before the
# administrative end of follow-up, it determined the observed time, and
# no event was recorded.
Censored_ICE <- as.integer(
  is.finite(t_ice_cc) &
    t_ice_cc <= t_fu_cc &
    t_ice_cc == Y_obs   &
    Delta_obs == 0
)
cat("\nICE censoring indicator:\n"); print(table(Censored_ICE))
cat("\nICE censoring by treatment:\n")
print(table(Treat = A_treat, Censored_ICE = Censored_ICE))

## Analysis dataset
# Same analysis sample as analysis_df, differing only in column selection. It
# adds Censored_ICE (the IPCW censoring indicator) and omits t_pos/t_ice_capped,
# which the weighting never uses. 
analysis_ipcw <- data.frame(
  Event        = Delta_obs[valid_time_idx],
  Treat        = A_treat[valid_time_idx],
  PersonTime   = Y_obs[valid_time_idx],
  Censored_ICE = Censored_ICE[valid_time_idx]
)
analysis_ipcw   <- cbind(analysis_ipcw,
                         W_covars_clean[valid_time_idx, , drop = FALSE])
covariate_names <- paste(colnames(W_covars_clean), collapse = " + ")
cat(sprintf("\nanalysis sample: %d | events: %d (vaccine %d, control %d)\n",
            nrow(analysis_ipcw), sum(analysis_ipcw$Event),
            sum(analysis_ipcw$Event == 1 & analysis_ipcw$Treat == 1),
            sum(analysis_ipcw$Event == 1 & analysis_ipcw$Treat == 0)))
cat(sprintf("censored by an intercurrent event: %d (%.1f%%)\n",
            sum(analysis_ipcw$Censored_ICE),
            100 * mean(analysis_ipcw$Censored_ICE)))

## Censoring model and weights
# Cox model for the hazard of being censored by an intercurrent event, on
# treatment and the baseline covariates.
cens_formula <- as.formula(paste("Surv(PersonTime, Censored_ICE) ~ Treat +",
                                 covariate_names))
cens_mod     <- coxph(cens_formula, data = analysis_ipcw)
cumhaz_den   <- predict(cens_mod, type = "expected")
analysis_ipcw$G_den       <- exp(-cumhaz_den)   # P(uncensored | Treat, W)
analysis_ipcw$ipcw_weight <- exp(cumhaz_den)    # = 1 / G_den
stopifnot(all(is.finite(analysis_ipcw$ipcw_weight)),
          all(analysis_ipcw$ipcw_weight >= 1))
cat("\n--- IPCW WEIGHTS ---\n")
print(round(quantile(analysis_ipcw$ipcw_weight,
                     c(0, .5, .9, .95, .99, .999, 1)), 4))
ess <- sum(analysis_ipcw$ipcw_weight)^2 / sum(analysis_ipcw$ipcw_weight^2)
cat(sprintf("mean %.4f | max %.3f | max/mean %.2f | effective n %.0f of %d (%.1f%%)\n",
            mean(analysis_ipcw$ipcw_weight), max(analysis_ipcw$ipcw_weight),
            max(analysis_ipcw$ipcw_weight) / mean(analysis_ipcw$ipcw_weight),
            ess, nrow(analysis_ipcw), 100 * ess / nrow(analysis_ipcw)))
cat(sprintf("smallest estimated probability of remaining ICE-free: %.4f\n",
            min(analysis_ipcw$G_den)))

## Treatment-only censoring model
# Not used for the primary weights. It is the numerator of the stabilised
# weights and the marginal term in HIPCW, so it is fitted here once.
cens_mod_num <- coxph(Surv(PersonTime, Censored_ICE) ~ Treat,
                      data = analysis_ipcw)
cumhaz_num   <- predict(cens_mod_num, type = "expected")

## Weighted Poisson model
formula_ipcw <- as.formula(paste("Event ~ Treat +", covariate_names,
                                 "+ offset(log(PersonTime))"))
mod_ipcw <- glm(formula_ipcw, family = poisson(link = "log"),
                data = analysis_ipcw, weights = ipcw_weight)
irr_ipcw   <- exp(coef(mod_ipcw)["Treat"])
robust_cov <- vcovHC(mod_ipcw, type = "HC0")
ci_ipcw    <- exp(coefci(mod_ipcw, vcov = robust_cov)["Treat", ])
ve_ipcw       <- as.numeric((1 - irr_ipcw) * 100)
ve_lower_ipcw <- as.numeric((1 - ci_ipcw[2]) * 100)   # endpoints reverse
ve_upper_ipcw <- as.numeric((1 - ci_ipcw[1]) * 100)
cat(sprintf("\nIPCW IRR: %.4f (95%% CI: %.4f to %.4f)\n",
            irr_ipcw, ci_ipcw[1], ci_ipcw[2]))
cat(sprintf("IPCW Vaccine Efficacy: %.2f%% (95%% CI: %.2f%% to %.2f%%)\n",
            ve_ipcw, ve_lower_ipcw, ve_upper_ipcw))

# PART II:  HIPCW POISSON
#  Stabilised HIPCW weight, equation (10) of Waterschoot et al.:
#     w = [G_num(t) / G_den(t)] x exp(lp_num - lp_den)
survival_ratio    <- exp(cumhaz_den - cumhaz_num)        # G_num / G_den
lp_den <- predict(cens_mod,     type = "lp", reference = "zero")  # conditional
lp_num <- predict(cens_mod_num, type = "lp", reference = "zero")  # marginal
hazard_correction <- exp(lp_num - lp_den)
analysis_ipcw$hipcw_weight <- survival_ratio * hazard_correction
wh <- analysis_ipcw$hipcw_weight
fit_ve <- function(w, d = analysis_ipcw) {
  d$wgt <- w
  m  <- glm(formula_ipcw, family = poisson(link = "log"), data = d, weights = wgt)
  ci <- exp(coefci(m, vcov = vcovHC(m, type = "HC0"))["Treat", ])
  irr <- as.numeric(exp(coef(m)["Treat"]))
  c(IRR = irr, VE = 100 * (1 - irr),
    VE_lo = as.numeric(100 * (1 - ci[2])),
    VE_hi = as.numeric(100 * (1 - ci[1])))
}
wsum <- function(w) c(mean = mean(w), max = max(w),
                      eff_n = sum(w)^2 / sum(w^2))
mod_hipcw   <- glm(formula_ipcw, family = poisson(link = "log"),
                   data = analysis_ipcw, weights = hipcw_weight)
irr_hipcw   <- exp(coef(mod_hipcw)["Treat"])
ci_hipcw    <- exp(coefci(mod_hipcw,
                          vcov = vcovHC(mod_hipcw, type = "HC0"))["Treat", ])
ve_hipcw    <- as.numeric((1 - irr_hipcw) * 100)
ve_lo_hipcw <- as.numeric((1 - ci_hipcw[2]) * 100)
ve_hi_hipcw <- as.numeric((1 - ci_hipcw[1]) * 100)
cat("\n--- IPCW AGAINST HIPCW ---\n")
# max/mean, not max: the HIPCW weight omits the baseline-hazard ratio, so
# its absolute magnitude is not on the same scale as the IPCW weight's.
# The ratio is, and it is what the comparison rests on.
print(data.frame(
  method       = c("IPCW", "HIPCW"),
  IRR          = round(c(irr_ipcw, irr_hipcw), 4),
  VE           = round(c(ve_ipcw, ve_hipcw), 2),
  VE_lower     = round(c(ve_lower_ipcw, ve_lo_hipcw), 2),
  VE_upper     = round(c(ve_upper_ipcw, ve_hi_hipcw), 2),
  max_over_mean = round(c(max(analysis_ipcw$ipcw_weight) /
                            mean(analysis_ipcw$ipcw_weight),
                          max(wh) / mean(wh)), 2)),
  row.names = FALSE)
cat(sprintf("HIPCW moves the estimate by %.2f pp relative to IPCW\n",
            abs(ve_hipcw - ve_ipcw)))
##############################
# PART III: Appendix (D.3)
w0 <- analysis_ipcw$ipcw_weight
n  <- nrow(analysis_ipcw)

# Distribution of the weights (D.3.1)
cat(sprintf("censored by an intercurrent event: %d of %d (%.1f%%)\n",
            sum(analysis_ipcw$Censored_ICE), n,
            100 * mean(analysis_ipcw$Censored_ICE)))
eff <- sum(w0)^2 / sum(w0^2)
cat(sprintf("weights: mean %.3f | max %.3f | effective n %.0f (%.1f%% of n)\n",
            mean(w0), max(w0), eff, 100 * eff / n))
cat(sprintf("smallest estimated probability of remaining ICE-free: %.4f\n",
            min(analysis_ipcw$G_den)))

# Stabilisation and truncation, both ultimately not included
trunc1 <- function(w, p = 0.01) {
  b <- quantile(w, c(p, 1 - p), na.rm = TRUE); pmin(pmax(w, b[1]), b[2])
}
variants_w <- list(
  "unstabilised, untruncated" = exp(cumhaz_den),
  "unstabilised, truncated"   = trunc1(exp(cumhaz_den)),
  "stabilised, untruncated"   = exp(cumhaz_den - cumhaz_num),
  "stabilised, truncated"     = trunc1(exp(cumhaz_den - cumhaz_num)))
tab_w <- do.call(rbind, lapply(names(variants_w), function(nm)
  data.frame(variant = nm,
             t(round(wsum(variants_w[[nm]]), 3)),
             t(round(fit_ve(variants_w[[nm]]), 3)))))
cat("\n"); print(tab_w, row.names = FALSE)
cat(sprintf("\nstabilisation moves VE by %.3f pp; truncation by %.3f pp\n",
            abs(tab_w$VE[3] - tab_w$VE[1]), abs(tab_w$VE[2] - tab_w$VE[1])))
write.csv(tab_w, "tables/ipcw_weight_construction.csv", row.names = FALSE)

# Proportional-hazards assumption (D.3.2, Table 21, Figure 8)
ph_den <- cox.zph(cens_mod)
print(ph_den)
g <- ph_den$table["GLOBAL", ]
cat(sprintf("\n>> GLOBAL: chi2 = %.1f on %.0f df, p = %.3g\n",
            g["chisq"], g["df"], g["p"]))
tab_ph   <- as.data.frame(ph_den$table)
tv_terms <- rownames(tab_ph)[rownames(tab_ph) != "GLOBAL" & tab_ph$p < 0.01]
tv_terms_pool <- tv_terms
tv_terms      <- setdiff(tv_terms, "AGE")
cat(sprintf(">> terms failing at the 1%% level: %d of %d\n",
            length(tv_terms_pool), nrow(tab_ph) - 1))
# Left panel: the treatment coefficient over time.
i_treat <- which(rownames(ph_den$table) == "Treat")
pdf("figures/ipcw_schoenfeld_denominator.pdf", width = 6, height = 5)
plot(ph_den[i_treat]); abline(h = 0, lty = 2); dev.off()
# Right panel: the term with the strongest evidence against proportionality.
i_worst <- which.max(tab_ph$chisq[rownames(tab_ph) != "GLOBAL"])
cat(sprintf(">> worst-violating term is column index %d (chi2 = %.1f)\n",
            i_worst, tab_ph$chisq[i_worst]))
pdf("figures/ipcw_schoenfeld_worst.pdf", width = 6, height = 5)
plot(ph_den[i_worst]); abline(h = 0, lty = 2); dev.off()

# Composition and timing of the censoring process (D.3.3)
# Which category produced each participant's earliest intercurrent event.
M <- sapply(ice_vars, function(v) as.numeric(raw_keep[[v]] - raw_keep$vax_date))
M[is.na(M)] <- Inf                    # no event of this category
first_t   <- apply(M, 1, min)
first_cat <- ice_vars[apply(M, 1, which.min)]
sel       <- is.finite(first_t) & analysis_ipcw$Censored_ICE == 1
cat_f   <- factor(first_cat[sel], levels = ice_vars)
tab_ice <- data.frame(category = ice_vars,
                      n      = as.vector(table(cat_f)),
                      median = as.vector(tapply(first_t[sel], cat_f, median)))
print(tab_ice, row.names = FALSE)
cat(sprintf("\ntotal ICE-censored: %d\n", sum(sel)))
cat(sprintf("vaccination-schedule: %d | infection: %d | remainder: %d\n",
            sum(tab_ice$n[tab_ice$category %in% c("de_1070", "de_2080")]),
            tab_ice$n[tab_ice$category == "dinf1618"],
            sum(tab_ice$n[tab_ice$category %in% c("de_1040", "de_2050")])))
write.csv(tab_ice, "tables/ice_timing.csv", row.names = FALSE)

# Influence diagnostics (D.3.4, Figure 9)
dfb <- residuals(cens_mod, type = "dfbetas")
mx  <- apply(abs(dfb), 2, max)
cat(sprintf(">> overall maximum |DFBETAS|: %.3f\n", max(mx)))
cat(sprintf(">> coefficients with max |DFBETAS| > 0.10: %d of %d\n",
            sum(mx > 0.10), length(mx)))
cat(sprintf(">> share of all values within +/- 0.2: %.1f%%\n",
            100 * mean(abs(dfb) <= 0.2)))
pdf("figures/ipcw_dfbetas.pdf", width = 7, height = 5)
matplot(dfb, type = "p", pch = 20, cex = 0.4, col = grey(0.3),
        xlab = "Participant", ylab = "DFBETAS")
abline(h = 0, lty = 2); dev.off()

# Sensitivity to non-proportional hazards (D.3.5, Table 22)
## 5a. No covariates
cens_km <- coxph(Surv(PersonTime, Censored_ICE) ~ 1, data = analysis_ipcw)
w_km    <- exp(predict(cens_km, type = "expected"))
## 5b. Stratified on the three worst-violating terms
ph_terms    <- tab_ph[rownames(tab_ph) != "GLOBAL", ]      # drop the global row
strat_terms <- rownames(ph_terms)[order(ph_terms$chisq, decreasing = TRUE)][1:3]
cat("stratifying on:", paste(strat_terms, collapse = ", "), "\n")
# those three move from the linear predictor into strata(), so each gets its
# own baseline hazard instead of a proportional one
rhs_strat <- c(setdiff(c("Treat", colnames(W_covars_clean)), strat_terms),
               paste0("strata(", strat_terms, ")"))
cens_strat <- coxph(as.formula(paste("Surv(PersonTime, Censored_ICE) ~",
                                     paste(rhs_strat, collapse = " + "))),
                    data = analysis_ipcw)
w_strat <- exp(predict(cens_strat, type = "expected"))     # 1/G under this model
## 5c. Time-varying coefficients
# Let the hazard vary linearly with log(time) for the eleven terms
# failing the test (AGE excluded). (corresponds to one extra term per coefficient)
tv_form <- as.formula(paste("Surv(PersonTime, Censored_ICE) ~ Treat +",
                            covariate_names, "+",
                            paste0("tt(", tv_terms, ")", collapse = " + ")))
cens_mod_tv <- coxph(tv_form, data = analysis_ipcw,
                     tt = function(x, t, ...) x * log(t))
cat(sprintf("time-varying terms: %d (%d coefficients against %d)\n",
            length(tv_terms), length(coef(cens_mod_tv)), length(coef(cens_mod))))
# predict(type = "expected") is not defined for tt models, so build the
# cumulative hazard by hand from Breslow increments.
cf   <- coef(cens_mod_tv)
Xfix <- as.matrix(analysis_ipcw[, c("Treat", colnames(W_covars_clean))])
a_i  <- as.numeric(Xfix %*% cf[colnames(Xfix)])                  # intercept
b_i  <- as.numeric(as.matrix(analysis_ipcw[, tv_terms, drop = FALSE]) %*%
                     cf[paste0("tt(", tv_terms, ")")])           # slope in log(t)
t_i  <- analysis_ipcw$PersonTime                      # each participant's time
ev_t <- t_i[analysis_ipcw$Censored_ICE == 1]          # times an ICE occurred
u    <- sort(unique(ev_t)); logu <- log(u)            # distinct event times
d_k  <- tabulate(match(ev_t, u), nbins = length(u))   # events tied at each one
# Baseline hazard, one increment per event time. The events observed at u_k
# divided by the summed hazard of everyone still at risk at u_k. Only those
# with t_i >= u_k are at risk, and each contributes exp(a_i + b_i*log(u_k)).
dL <- vapply(seq_along(u), function(k) {
  r <- t_i >= u[k]
  d_k[k] / sum(exp(a_i[r] + b_i[r] * logu[k]))
}, numeric(1))
# Each participant's cumulative hazard is their own hazard summed over the
# event times they lived through, weighted by those increments. findInterval
# gives how many event times fall at or before t_i. The weight is exp(H),
# the reciprocal of the estimated probability of staying uncensored.
n_at <- findInterval(t_i, u)
w_tv <- exp(vapply(seq_along(t_i), function(i)
  sum(exp(a_i[i] + b_i[i] * logu[seq_len(n_at[i])]) * dL[seq_len(n_at[i])]),
  numeric(1)))

## The comparison table
spec <- list("No covariates"             = w_km,
             "Proportional hazards"      = w0,
             "Stratified on three terms" = w_strat,
             "Time-varying coefficients" = w_tv)
tab_sens <- do.call(rbind, lapply(names(spec), function(nm)
  data.frame(model = nm,
             t(round(wsum(spec[[nm]]), 3)),
             t(round(fit_ve(spec[[nm]]), 3)))))
cat("\n"); print(tab_sens, row.names = FALSE)
cat(sprintf("\nrange across the four specifications: %.3f pp\n",
            max(tab_sens$VE) - min(tab_sens$VE)))
cat(sprintf("covariate adjustment is worth: %.3f pp\n",
            abs(tab_sens$VE[2] - tab_sens$VE[1])))
write.csv(tab_sens, "tables/ipcw_ph_sensitivity.csv", row.names = FALSE)

## Does the flexible model actually fit better?
cat("\n")
print(c(PH = AIC(cens_mod), time_varying = AIC(cens_mod_tv)))
dev <- as.numeric(2 * (logLik(cens_mod_tv) - logLik(cens_mod)))
ddf <- length(coef(cens_mod_tv)) - length(coef(cens_mod))
cat(sprintf(">> deviance difference %.1f on %d df (%d against %d parameters), p = %.3g\n",
            dev, ddf, length(coef(cens_mod_tv)), length(coef(cens_mod)),
            pchisq(dev, ddf, lower.tail = FALSE)))

## Sensitivity to the definition of an ICE (D.3.6, Table 23)
t_pos0 <- RAW$t_pos                                    # time to CIN2+
t_fu0  <- RAW$t_fu                                     # end of follow-up
ev_ind <- as.numeric(as.character(RAW[[event_ind_col]]))   # was CIN2+ recorded
# Rebuilds the whole analysis from the raw dates for a given ICE definition,
# so the only thing changing between the two runs is which categories count.
run_variant <- function(vars, label) {
  # days from vaccination to each ICE category, inf where it never happened
  Mv <- sapply(vars, function(v) as.numeric(RAW[[v]] - RAW$vax_date))
  Mv[is.na(Mv)] <- Inf
  t_ice <- apply(Mv, 1, min)                   # first ICE under this definition
  
  C <- pmin(t_ice, t_fu0)                      # censoring: ICE or end of study
  Y <- pmin(t_pos0, C)                         # observed follow-up time
  D <- as.numeric(ev_ind == 1 & t_pos0 <= C)   # CIN2+ seen before censoring
  # censored BY an ICE: it happened, within follow-up, it set Y, and no event
  I <- as.integer(is.finite(t_ice) & t_ice <= t_fu0 & t_ice == Y & D == 0)
  
  # same participants as the primary analysis, new clocks and indicator
  d <- data.frame(Event = D[valid_time_idx], Treat = A_treat[valid_time_idx],
                  PersonTime = Y[valid_time_idx], Censored_ICE = I[valid_time_idx])
  d <- cbind(d, W_covars_clean[valid_time_idx, , drop = FALSE])
  
  cm  <- coxph(as.formula(paste("Surv(PersonTime, Censored_ICE) ~ Treat +",
                                covariate_names)), data = d)   # censoring model
  d$w <- exp(predict(cm, type = "expected"))   # 1/G, as in the primary analysis
  v   <- fit_ve(d$w, d = d)                    # weighted Poisson -> VE and CI
  zt  <- as.data.frame(cox.zph(cm)$table)      # PH check for this variant
  
  data.frame(variant = label, n = nrow(d), events = sum(d$Event),
             ev_vac  = sum(d$Event == 1 & d$Treat == 1),
             ev_ctrl = sum(d$Event == 1 & d$Treat == 0),
             ice_cens = sum(d$Censored_ICE),
             mean_w = round(mean(d$w), 4), max_w = round(max(d$w), 3),
             VE = round(v["VE"], 3),
             VE_lo = round(v["VE_lo"], 2), VE_hi = round(v["VE_hi"], 2),
             ph_fail = sum(zt$p[rownames(zt) != "GLOBAL"] < 0.01),  # terms at 1%
             ph_terms = nrow(zt) - 1,                               # minus GLOBAL
             ph_p = signif(zt["GLOBAL", "p"], 3))
}
tab_def <- rbind(run_variant(ice_vars, "A: infection as ICE"),
                 run_variant(setdiff(ice_vars, "dinf1618"),
                             "B: infection as ordinary follow-up"))
print(tab_def, row.names = FALSE)
# A should rebuild the primary analysis exactly, which checks the whole
# reconstruction above against the main script.
cat(sprintf("\nvariant A reproduces the primary estimate: %s (%.3f vs %.3f)\n",
            abs(tab_def$VE[1] - ve_ipcw) < 1e-3, tab_def$VE[1], ve_ipcw))
# the events B gains are the ones A was censoring at infection
cat(sprintf("extra events under B: %d (vaccine %+d, control %+d)\n",
            tab_def$events[2] - tab_def$events[1],
            tab_def$ev_vac[2] - tab_def$ev_vac[1],
            tab_def$ev_ctrl[2] - tab_def$ev_ctrl[1]))
# those gained events over the participants B stopped censoring
cat(sprintf("event rate among the %d released from censoring: %.1f%%\n",
            tab_def$ice_cens[1] - tab_def$ice_cens[2],
            100 * (tab_def$events[2] - tab_def$events[1]) /
              (tab_def$ice_cens[1] - tab_def$ice_cens[2])))
# B rests on more events, so its interval should be narrower
cat(sprintf("CI width: %.1f pp (A) against %.1f pp (B)\n",
            tab_def$VE_hi[1] - tab_def$VE_lo[1],
            tab_def$VE_hi[2] - tab_def$VE_lo[2]))
# how much of the PH violation was infection-driven
cat(sprintf("Schoenfeld failures at 1%%: %d of %d (A) against %d of %d (B); global p %.3g against %.3g\n",
            tab_def$ph_fail[1], tab_def$ph_terms[1],
            tab_def$ph_fail[2], tab_def$ph_terms[2],
            tab_def$ph_p[1], tab_def$ph_p[2]))
write.csv(tab_def, "tables/ipcw_ice_definition.csv", row.names = FALSE)

# Functional form of age (D.3.7, Table 24, Figure 10)
cat(sprintf("range of age: %d to %d years\n",
            min(analysis_ipcw$AGE), max(analysis_ipcw$AGE)))
# Figure 10: martingale residuals against age. A LOWESS curve staying flat
# near zero means the linear term is not missing curvature.
mr <- residuals(cens_mod, type = "martingale")
pdf("figures/ipcw_age_martingale.pdf", width = 6, height = 5)
plot(analysis_ipcw$AGE, mr, pch = 20, cex = 0.3, col = grey(0.4),
     xlab = "Age (years)", ylab = "Martingale residual")
lines(lowess(analysis_ipcw$AGE, mr, iter = 0), lwd = 2)
abline(h = 0, lty = 2); dev.off()
# Two (other flexible) versions of the same censoring model, for Table 24
m_qua <- coxph(update(cens_formula, . ~ . + I(AGE^2)), data = analysis_ipcw)
m_spl <- coxph(update(cens_formula, . ~ . - AGE + ns(AGE, df = 3)),
               data = analysis_ipcw)
tab_age <- data.frame(
  specification = c("Linear", "Linear + quadratic", "Natural spline (3 df)"),
  logLik = round(c(logLik(cens_mod), logLik(m_qua), logLik(m_spl)), 2),
  AIC    = round(c(AIC(cens_mod),    AIC(m_qua),    AIC(m_spl)), 2))
cat("\n"); print(tab_age, row.names = FALSE)
write.csv(tab_age, "tables/ipcw_age_functional_form.csv", row.names = FALSE)
# Do the weights move the estimate?
w_spl <- exp(predict(m_spl, type = "expected"))
cat("\ntreatment effect under each age specification:\n")
print(rbind(data.frame(age = "Linear",         t(round(fit_ve(w0), 5))),
            data.frame(age = "Natural spline", t(round(fit_ve(w_spl), 5)))),
      row.names = FALSE)
