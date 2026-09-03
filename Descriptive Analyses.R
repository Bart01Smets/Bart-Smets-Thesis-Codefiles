#  DESCRIPTIVE ANALYSIS (thesis Section 6 and Appendix C)
#  Run immediately after the data have been processed.
#
#    Cohort derivation, follow-up and events (Tables 4, 5)
#    Covariate pre-processing (Table 6)
#    Intercurrent events (Tables 7, 8)
#    Censoring mechanism (Table 9)
#    Cumulative incidence (Figure 2)
#    Type-specific ICE models (Table 14)
#    Timing of intercurrent events (Figure 4)
#    Kaplan-Meier survival (Figure 5)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)
ice_vars <- c("de_1040",   # administration of vaccine(s) forbidden by protocol
              "de_1070",   # dose not administered according to protocol
              "de_2050",   # underlying medical condition forbidden by protocol
              "de_2080",   # noncompliance with vaccination schedule
              "dinf1618")  # incident HPV-16/18 before schedule completion

#  6.1.1  Cohort derivation, follow-up and events
# Participants whose ICE fell at or before the day of first vaccination
# contribute no positive observation time and are excluded. 
excl <- dat$t_ice <= 0

                              # raw rows kept
Wdf     <- W_covars_clean[valid_time_idx, , drop = FALSE] # covariates, analysis sample
t_ice_a <- t_ice_cc[valid_time_idx]                       # ICE clock, analysis sample

cat(sprintf("N = %d; excluded %d (%d vaccine, %d control); analysis sample %d\n",
            nrow(dat), sum(excl), sum(excl & dat$Treat == 1), sum(excl & dat$Treat == 0),
            nrow(raw_keep)))
cat(sprintf("vaccine %d (%.1f%%) | control %d (%.1f%%)\n",
            sum(analysis_df$Treat == 1), 100 * mean(analysis_df$Treat == 1),
            sum(analysis_df$Treat == 0), 100 * mean(analysis_df$Treat == 0)))

# Table 4: follow-up time and person-years by arm.
arm_col <- function(k) c(sum(k), median(analysis_df$PersonTime[k]),
                         max(analysis_df$PersonTime[k]),
                         sum(analysis_df$PersonTime[k]) / 365.25)
table4 <- data.frame(
  Quantity = c("Participants, n", "Median follow-up, days",
               "Maximum follow-up, days", "Person-time, person-years"),
  Vaccine  = arm_col(analysis_df$Treat == 1),
  Control  = arm_col(analysis_df$Treat == 0),
  Total    = arm_col(rep(TRUE, nrow(analysis_df))))
cat("\nTable 4\n"); print(table4, row.names = FALSE, digits = 6)

# Table 5: events and crude incidence rates, per 1000 person-years.
py <- c(sum(analysis_df$PersonTime[analysis_df$Treat == 1]),
        sum(analysis_df$PersonTime[analysis_df$Treat == 0])) / 365.25
ev <- c(sum(analysis_df$Event == 1 & analysis_df$Treat == 1),
        sum(analysis_df$Event == 1 & analysis_df$Treat == 0))
table5 <- data.frame(
  Quantity = c("Events (Delta = 1)", "Event-free (Delta = 0)", "Incidence rate (/1000 py)"),
  Vaccine  = c(ev[1], sum(analysis_df$Treat == 1) - ev[1], 1000 * ev[1] / py[1]),
  Control  = c(ev[2], sum(analysis_df$Treat == 0) - ev[2], 1000 * ev[2] / py[2]),
  Total    = c(sum(ev), nrow(analysis_df) - sum(ev), NA))
cat("\nTable 5\n"); print(table5, row.names = FALSE, digits = 4)
cat(sprintf("crude IRR %.3f -> unadjusted VE %.1f%%\n",
            (ev[1]/py[1]) / (ev[2]/py[2]), 100 * (1 - (ev[1]/py[1]) / (ev[2]/py[2]))))

#  6.1.2  Covariate pre-processing
# Near-zero-variance screening, at caret's defaults: a column goes when its
# most frequent value occurs more than 19 times as often as its second.
cat("\nremoved by near-zero-variance screening:\n")
print(colnames(W_numeric)[bad_cols_idx])
cat(sprintf("condition number: %.1f (below 30 indicates no problem)\n", condition_number))

# Table 6: baseline characteristics by arm, after screening, so it describes
# the columns the estimators actually use. Age is the only continuous
# covariate, so it gets a mean and everything else n (%).
fmt <- function(x) {
  if (length(unique(x)) > 2) sprintf("%.2f", mean(x))          # continuous: mean
  else sprintf("%d (%.1f%%)", sum(x), 100 * mean(x))           # binary: n (%)
}
vax <- analysis_df$Treat == 1
cat("\nTable 6\n")
print(data.frame(Characteristic = names(Wdf),
                 Vaccine = sapply(Wdf[vax, ],  fmt),
                 Control = sapply(Wdf[!vax, ], fmt)), row.names = FALSE)

#  6.1.3  Intercurrent events (Tables 7 and 8)
# Table 7: ICEs at or before time zero, split into strictly before and
# exactly at. Rows count events, so they exceed the 1204 excluded
# participants where someone recorded two types on the same date.
table7 <- data.frame(
  ICE    = ice_vars,
  Before = sapply(ice_vars, function(v) sum(!is.na(dat[[v]]) & dat[[v]] <  dat$vax_date)),
  At     = sapply(ice_vars, function(v) sum(!is.na(dat[[v]]) & dat[[v]] == dat$vax_date)))
table7$Total <- table7$Before + table7$At
cat("\nTable 7\n"); print(table7, row.names = FALSE)
cat(sprintf("any ICE: %d strictly before, %d exactly at, %d participants (types sum to %d)\n",
            sum(dat$t_ice < 0), sum(dat$t_ice == 0), sum(excl), sum(table7$Total)))
# the HPV-16/18 exclusions are balanced across arms, so this reads as an
# eligibility restriction rather than a treatment-induced one
hpv0 <- !is.na(dat[["dinf1618"]]) & dat[["dinf1618"]] <= dat$vax_date
cat(sprintf("HPV-16/18 at time zero: %d vaccine, %d control, p = %.2f\n",
            sum(hpv0 & dat$Treat == 1), sum(hpv0 & dat$Treat == 0),
            prop.test(c(sum(hpv0 & dat$Treat == 1), sum(hpv0 & dat$Treat == 0)),
                      c(sum(dat$Treat == 1), sum(dat$Treat == 0)),
                      correct = FALSE)$p.value))

# Table 8: ICEs strictly after time zero, in the analysis sample.
table8 <- data.frame(
  ICE     = ice_vars,
  Vaccine = sapply(ice_vars, function(v)
    sum(!is.na(raw_keep[[v]]) & raw_keep[[v]] > raw_keep$vax_date & raw_keep$Treat == 1)),
  Control = sapply(ice_vars, function(v)
    sum(!is.na(raw_keep[[v]]) & raw_keep[[v]] > raw_keep$vax_date & raw_keep$Treat == 0)))
table8$Total <- table8$Vaccine + table8$Control
cat("\nTable 8\n"); print(table8, row.names = FALSE)
any_ice <- is.finite(t_ice_a) & t_ice_a > 0     # t_ice is Inf when no ICE occurred
cat(sprintf("any ICE: %d participants (%.1f%%), %d vaccine, %d control (types sum to %d)\n",
            sum(any_ice), 100 * mean(any_ice), sum(any_ice & analysis_df$Treat == 1),
            sum(any_ice & analysis_df$Treat == 0), sum(table8$Total)))

#  6.2.1  Censoring mechanism (Table 9)
# Censored by an intercurrent event: all four conditions are needed. An ICE
# occurred, at or before the administrative end of follow-up, it determined
# the observed time, and no CIN2+ event was recorded.
Censored_ICE <- as.integer(is.finite(t_ice_cc) & t_ice_cc <= t_fu_cc &
                             t_ice_cc == Y_obs & Delta_obs == 0)[valid_time_idx]
# Cox model for the time to the first ICE of any type. This asks the
# question; is censoring selective with respect to treatment and
# measured covariates? not whether it is informative, which these data
# cannot settle.
cens_df  <- cbind(data.frame(PersonTime = analysis_df$PersonTime, Censored_ICE,
                             Treat = analysis_df$Treat), Wdf)
cens_mod <- coxph(as.formula(paste("Surv(PersonTime, Censored_ICE) ~ Treat +",
                                   paste(names(Wdf), collapse = " + "))),
                  data = cens_df)
s9 <- summary(cens_mod)
table9 <- data.frame(
  Predictor = rownames(s9$coefficients),
  HR        = s9$coefficients[, "exp(coef)"],
  CI        = sprintf("%.2f-%.2f", s9$conf.int[, "lower .95"], s9$conf.int[, "upper .95"]),
  p         = s9$coefficients[, "Pr(>|z|)"],
  sig       = ifelse(s9$coefficients[, "Pr(>|z|)"] < 0.01, "**",
                     ifelse(s9$coefficients[, "Pr(>|z|)"] < 0.05, "*", "")))
cat("\nTable 9\n"); print(table9, row.names = FALSE, digits = 3)

#  6.2.2  Cumulative incidence (Figure 2)
km_fit <- survfit(Surv(PersonTime, Event) ~ Treat, data = analysis_df)
lr     <- survdiff(Surv(PersonTime, Event) ~ Treat, data = analysis_df)
cat(sprintf("\nlog-rank: chisq = %.1f on 1 df, p = %.3g\n",
            lr$chisq, pchisq(lr$chisq, 1, lower.tail = FALSE)))
sm <- summary(km_fit, times = c(1461, 1522))
print(data.frame(arm = sm$strata, day = sm$time, at_risk = sm$n.risk,
                 cum_inc_pct = 100 * (1 - sm$surv)), row.names = FALSE, digits = 3)
cat(sprintf("%d of %d control-arm events fall after day 1522\n",
            sum(analysis_df$Event == 1 & analysis_df$Treat == 0 &
                  analysis_df$PersonTime > 1522), ev[2]))
cat(sprintf("%.0f%% of the %d events occur after day 1000\n",
            100 * mean(analysis_df$PersonTime[analysis_df$Event == 1] > 1000), sum(ev)))
fig2 <- ggsurvplot(
  km_fit, data = analysis_df, fun = "event",
  conf.int = TRUE, censor = FALSE, size = 1.2,
  palette = c("#d95f02", "#1b9e77"), break.time.by = 365,
  xlab = "Days since vaccination", ylab = "Cumulative probability of CIN2+",
  legend.title = "Treatment group", legend.labs = c("Placebo (0)", "Vaccine (1)"),
  risk.table = "nrisk_cumcensor",          # "number at risk (number censored)"
  ggtheme = theme_minimal())
pdf("figures/cumulative_incidence.pdf", width = 8, height = 6)
print(fig2); dev.off()

#  C.1  Type-specific ICE models (Table 14)
# One Cox model per ICE type, on the follow-up clock, with "did this ICE
# occur" as the event. Table 14 reports a hazard ratio only where the
# association is significant.
covs   <- names(Wdf)
base14 <- data.frame(Time = analysis_df$PersonTime, Treat = analysis_df$Treat, Wdf)
form14 <- as.formula(paste("Surv(Time, ice_k) ~ Treat +", paste(covs, collapse = " + ")))
ice_models <- lapply(ice_vars, function(v) {
  d <- base14; d$ice_k <- as.integer(!is.na(raw_keep[[v]]))
  coxph(form14, data = d)
})
names(ice_models) <- paste0("ICE", seq_along(ice_vars))
hr_cell <- function(m, term) {
  s <- summary(m)$coefficients
  if (!term %in% rownames(s)) return("")
  p <- s[term, "Pr(>|z|)"]
  if (p >= 0.05) return("")
  sprintf("%.2f (%s)", s[term, "exp(coef)"], if (p < 0.01) "**" else "*")
}
terms14 <- c("Treat", covs)
table14 <- sapply(ice_models, function(m) vapply(terms14, hr_cell, character(1), m = m))
rownames(table14) <- terms14
table14 <- table14[apply(table14 != "", 1, any), , drop = FALSE]   # drop never-significant rows
cat("\nTable 14 (blank = not significant)\n"); print(table14, quote = FALSE)

#  C.2  Timing of intercurrent events (Figure 4)
ice_long <- raw_keep %>%
  select(Treat, vax_date, all_of(ice_vars)) %>%
  pivot_longer(all_of(ice_vars), names_to = "ICE_type", values_to = "ICE_date") %>%
  filter(!is.na(ICE_date)) %>%
  mutate(days     = as.numeric(ICE_date - vax_date),
         ICE_type = factor(paste0("ICE", match(ICE_type, ice_vars)),
                           levels = paste0("ICE", seq_along(ice_vars))),
         Arm      = factor(Treat, 0:1, c("Control", "Vaccine")))
fig4 <- ggplot(ice_long, aes(x = days, fill = Arm)) +
  geom_histogram(binwidth = 10, position = "identity", alpha = .55) +
  # vertical scales differ across panels, since the types differ in frequency
  facet_wrap(~ ICE_type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#d95f02", Vaccine = "#1b9e77")) +
  # the thesis figure shows only the first year; a few ICE2 and ICE3 events
  # occur later and fall outside the panel
  coord_cartesian(xlim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 300, 100)) +
  labs(x = "Days from vaccination to intercurrent event", y = "Participants", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position    = "top",
        strip.text         = element_text(hjust = 0),   # panel labels left-aligned
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())           # horizontal gridlines only
ggsave("figures/ice_timing_by_type.pdf", fig4, width = 6.5, height = 8)

#  C.3  Kaplan-Meier survival (Figure 5)
# Same information as Figure 2, since F(t) = 1 - S(t) with no competing
# events. Participants experiencing an ICE are censored here.
km_data <- data.frame(Time = Y_surv, Event = Delta_surv, Treatment = A_surv)
km_surv <- survfit(Surv(Time, Event) ~ Treatment, data = km_data)
cat(sprintf("\nevent-free survival at end of follow-up: control %.3f | vaccine %.3f\n",
            min(summary(survfit(Surv(Time, Event) ~ 1, data = km_data[km_data$Treatment == 0, ]))$surv),
            min(summary(survfit(Surv(Time, Event) ~ 1, data = km_data[km_data$Treatment == 1, ]))$surv)))

vax_time <- km_data$Time[km_data$Treatment == 1]
for (d in sort(vax_time[km_data$Event[km_data$Treatment == 1] == 1]))
  cat(sprintf("vaccine-arm event at day %4.0f: risk set %5d, increment %.2f pp\n",
              d, sum(vax_time >= d), 100 / sum(vax_time >= d)))
fig5 <- ggsurvplot(
  km_surv, data = km_data, pval = TRUE, conf.int = TRUE, censor = FALSE, size = 1.2,
  palette = c("#d95f02", "#1b9e77"), ylim = c(0.94, 1), break.time.by = 365,
  xlab = "Days since randomisation", ylab = "Probability of remaining CIN2+ free",
  legend.title = "Group", legend.labs = c("Placebo (0)", "Vaccine (1)"),
  risk.table = "nrisk_cumcensor", ggtheme = theme_minimal())
pdf("figures/kaplan_meier.pdf", width = 8, height = 6)
print(fig5); dev.off()