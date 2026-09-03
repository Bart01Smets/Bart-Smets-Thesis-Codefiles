
# MARGINAL PLOT ANALYSIS 
#
# Checking every treated-control pair directly would require
# n_treat * n_cont comparisons, resulting in millions of rows.
# The calculation below avoids creating all these pairs while giving
# the same result as direct pairwise averaging.
#
# For each observed event in one arm, I check how many people in the
# opposite arm were still under follow-up at the time of that event.
# Only events that were observed and occurred before the analysis
# horizon t are considered.
#
# Because there are only 3 such events in the vaccine arm and 86 in
# the control arm, this greatly reduces the amount of computation.
# Follow-up times in the opposite arm are sorted once, so the number
# still under follow-up at each event time can be obtained efficiently.
#
# The resulting counts are summed and divided by n_treat * n_cont to
# obtain the proportion of treated-control pairs contributing to psi_1.
#
# Direct pairwise averaging can also be used as a check and gives the
# same result.

# how many values in the sorted vector sorted_v are >= tt 
count_ge <- function(sorted_v, tt) length(sorted_v) - findInterval(tt, sorted_v, left.open = TRUE)

# core PLOT calculation
# risk_treat = share of all possible treated-control pairs where the
#              treated person had a valid, in-window event
# risk_cont  = same thing, for control-person events
# risk_treat / risk_cont is the treatment-effect ratio we care about.
plot_pair_risks <- function(treat_df, cont_df, t_eval = Inf) {
  n_t <- nrow(treat_df); n_c <- nrow(cont_df)
  if (n_t == 0 || n_c == 0) return(c(risk_treat = NA_real_, risk_cont = NA_real_))
  
  # sort each arm's follow-up cutoff once, so I can look up "still at risk?" fast
  s_t <- sort(treat_df$t_ice_capped)
  s_c <- sort(cont_df$t_ice_capped)
  
  # keep only event times that are valid: event happened, before that
  # person's own cutoff, and within the evaluation horizon
  e_t <- treat_df$t_pos[treat_df$Event == 1 & treat_df$t_pos <= treat_df$t_ice_capped & treat_df$t_pos <= t_eval]
  e_c <- cont_df$t_pos[cont_df$Event == 1 & cont_df$t_pos <= cont_df$t_ice_capped & cont_df$t_pos <= t_eval]
  
  # for each treated event, count how many control partners were still
  # being followed at that time (and vice versa); divide by all possible
  # pairs to get the average event rate per arm
  c(risk_treat = if (length(e_t)) sum(count_ge(s_c, e_t)) / (n_t * n_c) else 0,
    risk_cont  = if (length(e_c)) sum(count_ge(s_t, e_c)) / (n_t * n_c) else 0)
}

# bootstrap: repeat the calculation on resampled data to get a CI 
# Resampling is done separately within each arm (with replacement) because
# the treated/control split itself is fixed by randomisation, not random.
boot_plot <- function(treat_df, cont_df, B = 2000, t_eval = Inf, seed = 2026) {
  set.seed(seed)
  n_t <- nrow(treat_df); n_c <- nrow(cont_df)
  
  # pull out plain vectors up front
  Tt <- treat_df$t_ice_capped; Pt <- treat_df$t_pos; Et <- treat_df$Event
  Tc <- cont_df$t_ice_capped;  Pc <- cont_df$t_pos;  Ec <- cont_df$Event
  
  rt <- rc <- numeric(B)   # store risk_treat / risk_cont from each replicate
  for (b in seq_len(B)) {
    # draw a random resample (with replacement) of each arm
    i <- sample.int(n_t, n_t, TRUE); j <- sample.int(n_c, n_c, TRUE)
    
    # same calculation as plot_pair_risks(), just inlined for speed
    s_t <- sort(Tt[i]); s_c <- sort(Tc[j])
    e_t <- Pt[i][Et[i] == 1 & Pt[i] <= Tt[i] & Pt[i] <= t_eval]
    e_c <- Pc[j][Ec[j] == 1 & Pc[j] <= Tc[j] & Pc[j] <= t_eval]
    rt[b] <- if (length(e_t)) sum(count_ge(s_c, e_t)) / (n_t * n_c) else 0
    rc[b] <- if (length(e_c)) sum(count_ge(s_t, e_c)) / (n_t * n_c) else 0
    
    if (b %% 500 == 0) cat("  replicate", b, "/", B, "\n")  # progress ping
  }
  # the spread of these B estimates is what becomes the bootstrap CI later
  list(risk_treat = rt, risk_cont = rc)
}

#####
# ANALYSIS
#####

# Horizon is month 50, meaning an event anywhere in month 50 still counts.
# Months are floor(t_pos / 30.44), so "month <= 50" runs to the end of month
# 50, which is 51 * 30.44 days.

T_EVAL_DAYS <- 51 * 30.44 - 1e-9   # 
# the -1e-9 turns the <= checks elsewhere in the script into a strict <

treat_df <- analysis_df[analysis_df$Treat == 1, ]
cont_df  <- analysis_df[analysis_df$Treat == 0, ]
cat("Vaccine participants:", nrow(treat_df), "\n")
cat("Control participants:", nrow(cont_df), "\n")

# Original-sample PLOT estimate
risks <- plot_pair_risks(treat_df, cont_df, t_eval = T_EVAL_DAYS)
risk_treat <- unname(risks["risk_treat"]); risk_cont <- unname(risks["risk_cont"])

plot_ratio    <- risk_treat / risk_cont
plot_logratio <- log(plot_ratio)
plot_ve       <- (1 - plot_ratio) * 100

cat("\n=========================================================\n")
cat(" ORIGINAL MARGINAL PLOT RESULTS\n")
cat("=========================================================\n")
cat(sprintf("Marginal PLOT vaccine component: %.6f\n", risk_treat))
cat(sprintf("Marginal PLOT control component: %.6f\n", risk_cont))
cat(sprintf("PLOT Log-Ratio: %.6f\n", plot_logratio))
cat(sprintf("PLOT Ratio: %.6f\n", plot_ratio))
cat(sprintf("PLOT Vaccine Efficacy: %.2f%%\n", plot_ve))
cat("=========================================================\n")

# Bootstrap (percentile CI) 
B  <- 2000
t0 <- Sys.time()
bt <- boot_plot(treat_df, cont_df, B = B, t_eval = T_EVAL_DAYS)
boot_risk_treat <- bt$risk_treat
boot_risk_cont  <- bt$risk_cont
cat(sprintf("\nBootstrap completed in %.1f seconds.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Bootstrap diagnostics
# log-ratio only defined (finite) when both components are > 0
boot_logratio <- ifelse(boot_risk_treat > 0 & boot_risk_cont > 0, log(boot_risk_treat / boot_risk_cont), NA_real_)
n_finite    <- sum(is.finite(boot_logratio))
n_nonfinite <- sum(!is.finite(boot_logratio))
####################
#Appendix
cat("\n=========================================================\n")
cat(" BOOTSTRAP STABILITY DIAGNOSTICS\n")
cat("=========================================================\n")
cat("Bootstrap replicates:", B, "\n")
cat("Finite log-ratio replicates:", n_finite, "\n")
cat("Non-finite log-ratio replicates:", n_nonfinite, "\n")
cat("Percentage non-finite:", round(100 * n_nonfinite / B, 2), "%\n")
cat("Replicates with vaccine component = 0:", sum(boot_risk_treat == 0, na.rm = TRUE), "\n")
cat("Replicates with control component = 0:", sum(boot_risk_cont == 0, na.rm = TRUE), "\n")
cat("=========================================================\n")

# Percentile bootstrap CI, computed on the ratio scale
boot_ratio  <- boot_risk_treat / boot_risk_cont
ci_ratio    <- quantile(boot_ratio, probs = c(0.025, 0.975), na.rm = TRUE)
ci_logratio <- log(ci_ratio)
ci_ve       <- c((1 - ci_ratio[2]) * 100, (1 - ci_ratio[1]) * 100)  # VE decreases in ratio, so bounds reverse

# Final results
cat("\n=========================================================\n")
cat(sprintf(" FINAL MARGINAL PLOT RESULTS (B=%d Bootstrap)\n", B))
cat("=========================================================\n")
cat(sprintf("PLOT Log-Ratio: %.4f (95%% CI: %.4f to %.4f)\n", plot_logratio, ci_logratio[1], ci_logratio[2]))
cat(sprintf("PLOT Ratio: %.4f (95%% CI: %.4f to %.4f)\n", plot_ratio, ci_ratio[1], ci_ratio[2]))
cat(sprintf("PLOT Vaccine Efficacy: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", plot_ve, ci_ve[1], ci_ve[2]))
cat(sprintf("Finite bootstrap replicates: %d of %d\n", n_finite, B))
cat("=========================================================\n")

# Support and cumulative events by month (thesis Table 27) 

months <- c(6, 12, 18, 24, 30, 36, 42, 48, 50)
mo_end <- (months + 1) * 30.44 - 1e-9
# Per arm: % still intercurrent-event-free at the start of month m, and events
# accumulated by the end of month m.
arm_cols <- function(a) {
  s  <- analysis_df[analysis_df$Treat == a, ]
  ev <- s$Event == 1 & s$t_pos <= s$t_ice_capped   # event inside that person's own cutoff
  list(pct = sapply(months, function(m) 100 * mean(s$t_ice_capped >= m * 30.44)),
       cum = sapply(mo_end, function(k) sum(ev & s$t_pos <= k)))
}
vac <- arm_cols(1); ctl <- arm_cols(0)
support_tab <- data.frame(Month = months,
                          Control_pct_observable = ctl$pct, Control_cum_events = ctl$cum,
                          Vaccine_pct_observable = vac$pct, Vaccine_cum_events = vac$cum)
cat("\n=========================================================\n")
cat(" SUPPORT AND CUMULATIVE EVENTS BY MONTH\n")
cat("=========================================================\n")
print(support_tab, row.names = FALSE, digits = 4)

# Pair-specific common observation time (thesis Table/Fig D.5.2)
tc     <- pmin(analysis_df$t_ice_capped, T_EVAL_DAYS)
grid   <- sort(unique(tc))
S_t    <- 1 - ecdf(tc[analysis_df$Treat == 1])(grid)
S_c    <- 1 - ecdf(tc[analysis_df$Treat == 0])(grid)
F_pair <- 1 - S_t * S_c                                  # CDF of the pair common time
med_pair <- grid[which(F_pair >= 0.5)[1]]
iqr_pair <- grid[c(which(F_pair >= 0.25)[1], which(F_pair >= 0.75)[1])]
p10_pair <- grid[which(F_pair >= 0.10)[1]]               # D.5.2 also reports the 10th percentile

cat("\n=========================================================\n")
cat(" PAIR-SPECIFIC COMMON OBSERVATION TIME\n")
cat("=========================================================\n")
cat(sprintf("Median pair common time: %.0f days (IQR %.0f to %.0f, 10th pct %.0f)\n",
            med_pair, iqr_pair[1], iqr_pair[2], p10_pair))
cat(sprintf("Individuals unobservable by day 1350: %.1f%%\n", 100 * mean(tc < 1350)))
# F_pair is a step function, not continuous -- look up the exact step at day
# 1350 (last grid point <= 1350) rather than linearly interpolating, which
# understates it slightly
pair_by_1350 <- 100 * F_pair[max(which(grid <= 1350))]
cat(sprintf("Pairs already compared by day 1350:    %.1f%%\n", pair_by_1350))

# Truncation by intercurrent event, by arm (thesis Table 28)

t_fu_cc_aligned <- t_fu_cc[Y_obs > 0]
stopifnot(length(t_fu_cc_aligned) == nrow(analysis_df))
truncated <- analysis_df$t_ice_capped < t_fu_cc_aligned - 1e-6

trunc_tab <- data.frame(
  Group        = c("Control", "Vaccine", "Overall"),
  Participants = c(sum(analysis_df$Treat == 0), sum(analysis_df$Treat == 1), nrow(analysis_df)))
trunc_tab$Truncated  <- c(sum(truncated & analysis_df$Treat == 0),
                          sum(truncated & analysis_df$Treat == 1), sum(truncated))
trunc_tab$Pct        <- 100 * trunc_tab$Truncated / trunc_tab$Participants
# Median and IQR of the truncation time itself, i.e. when the ICE stopped the
# clock, among those it stopped. Same three groups, so one list of masks.
tr_grp <- list(truncated & analysis_df$Treat == 0,
               truncated & analysis_df$Treat == 1,
               truncated)
trunc_tab$MedianDays <- sapply(tr_grp, function(k) median(analysis_df$t_ice_capped[k]))
trunc_tab$IQRDays    <- sapply(tr_grp, function(k)
  paste(round(quantile(analysis_df$t_ice_capped[k], c(0.25, 0.75))), collapse = "-"))

cat("\n=========================================================\n")
cat(" TRUNCATION BY INTERCURRENT EVENTS\n")
cat("=========================================================\n")
print(trunc_tab, row.names = FALSE, digits = 4)

# D.5.3 text: how much follow-up the truncated participants lost, and the
# median of the end-of-follow-up component for everyone else (for whom
# t_ice_capped is simply their administrative end).
cat(sprintf("Median follow-up lost by truncated participants: %.0f days\n",
            median(t_fu_cc_aligned[truncated] - analysis_df$t_ice_capped[truncated])))
cat(sprintf("Median end-of-follow-up component:               %.0f days\n",
            median(analysis_df$t_ice_capped[!truncated])))

# Assumption 1, condition (a): are truncation rates equal across arms?

x     <- c(sum(truncated & analysis_df$Treat == 0), sum(truncated & analysis_df$Treat == 1))
n_grp <- c(sum(analysis_df$Treat == 0), sum(analysis_df$Treat == 1))
pt    <- prop.test(x, n_grp, correct = FALSE)

cat("\n=========================================================\n")
cat(" ASSUMPTION 1, CONDITION (a)\n")
cat("=========================================================\n")
cat(sprintf("Control truncated: %.2f%% | Vaccine truncated: %.2f%%\n",
            100 * pt$estimate[1], 100 * pt$estimate[2]))
cat(sprintf("Risk difference: %.2f pp (95%% CI %.2f to %.2f), p = %.3g\n",
            100 * (pt$estimate[1] - pt$estimate[2]),
            100 * pt$conf.int[1], 100 * pt$conf.int[2], pt$p.value))

# Figures 13 and 14

dir.create("figures", showWarnings = FALSE)

# Figure 13: support curves by arm
pdf("figures/plot_support_by_arm.pdf", width = 7, height = 5)
mo <- seq(0, 51, by = 0.5)
pc <- sapply(c(0, 1), function(a) {
  s <- analysis_df$t_ice_capped[analysis_df$Treat == a]
  sapply(mo, function(m) 100 * mean(s >= m * 30.44))
})
matplot(mo, pc, type = "l", lty = c(1, 2), col = "black", ylim = c(0, 100),
        xlab = "Month", ylab = "Participants remaining observable (%)")
legend("bottomleft", lty = c(1, 2), legend = c("Control", "Vaccine"), bty = "n")
dev.off()

# Figure 14: distribution of pair-specific common observation times.
# F_pair is the CDF on `grid`, so its differences are the mass at each time.
pdf("figures/plot_pair_common_time.pdf", width = 7, height = 5)
plot(grid / 30.44, diff(c(0, F_pair)), type = "h",
     xlab = "Pair-specific common observation time (months)",
     ylab = "Proportion of treated-control pairs")
dev.off()
