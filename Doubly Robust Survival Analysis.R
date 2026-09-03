#  Doubly-robust survival analysis by Westling et al. (2024)
#
#  PART I   Main analysis (thesis 7.2.2: Table 11, Figure 3)
#  PART II  Appendix D.4, in the order of its subsections:
#  Evaluation grid and reporting interval   (D.4.1)
#  SuperLearner performance                 (D.4.2, Table 25)
#  Positivity and censoring support         (D.4.3)
#  Support and influence over follow-up     (Table 26, Fig 11-12)
#  Influence-function diagnostics           (D.4.4)
#  Sensitivity to truncation                (D.4.5)
#



G_FLOOR <- 0.01     # truncation level for the censoring survival
T0_Q    <- 0.05     # reporting interval, quantiles of the event times
T1_Q    <- 0.95
N_SIM   <- 2000     # multiplier-bootstrap draws for the uniform band
V       <- 5        # cross-fitting folds
SEED    <- 1
POISSON_VE <- 95.76 # Table 10, for the reference line in Figure 3
n <- length(Y_surv)
# grid: every distinct observed event time (D.4.1) 
event_times <- Y_surv[Delta_surv == 1 & is.finite(Y_surv)]
eval_grid   <- sort(unique(event_times))
cat(sprintf("grid: %d points, %.0f to %.0f days\n",
            length(eval_grid), min(eval_grid), max(eval_grid)))

# propensity score 
prop_df  <- cbind(data.frame(Treat = A_surv), W_surv)
prop_mod <- glm(Treat ~ ., data = prop_df, family = binomial("logit"))
manual_prop_scores <- predict(prop_mod, type = "response")
cat(sprintf("propensity range: %.4f to %.4f\n",
            min(manual_prop_scores), max(manual_prop_scores)))

# outcome nuisance: marginal Kaplan-Meier per arm 
km_0 <- survfit(Surv(Y_surv[A_surv == 0], Delta_surv[A_surv == 0]) ~ 1)
km_1 <- survfit(Surv(Y_surv[A_surv == 1], Delta_surv[A_surv == 1]) ~ 1)
surv_0 <- summary(km_0, times = eval_grid, extend = TRUE)$surv
surv_1 <- summary(km_1, times = eval_grid, extend = TRUE)$surv
manual_pred_0 <- matrix(rep(surv_0, n), nrow = n, byrow = TRUE)
manual_pred_1 <- matrix(rep(surv_1, n), nrow = n, byrow = TRUE)

# censoring nuisance: cross-fitted SuperLearner
set.seed(SEED)
SL_cens_lib <- c("survSL.expreg", "survSL.coxph", "survSL.km",
                 "survSL.gam", "survSL.rfsrc")
# Folds are formed within arm, so the arms stay balanced across folds.
fold_id <- integer(n)
for (a in c(0, 1)) {
  idx <- which(A_surv == a)
  fold_id[idx] <- sample(rep(1:V, length.out = length(idx)))
}
cens.pred.0 <- matrix(NA_real_, n, length(eval_grid))
cens.pred.1 <- matrix(NA_real_, n, length(eval_grid))
SL_diagnostics <- list()
for (v in 1:V) {
  train <- which(fold_id != v)
  test  <- which(fold_id == v)
  
  for (a in c(0, 1)) {
    tr_a <- train[A_surv[train] == a]
    
    fit_c <- survSuperLearner(
      time      = Y_surv[tr_a],
      event     = 1 - Delta_surv[tr_a],           # censoring is the event here
      X         = as.data.frame(W_surv[tr_a, , drop = FALSE]),
      newX      = as.data.frame(W_surv[test, , drop = FALSE]),
      new.times = eval_grid,
      event.SL.library = SL_cens_lib,
      cens.SL.library  = c("survSL.km"),          # clinical events censor here
      verbose = FALSE)
    
    if (a == 0) cens.pred.0[test, ] <- fit_c$event.SL.predict
    if (a == 1) cens.pred.1[test, ] <- fit_c$event.SL.predict
    
    SL_diagnostics[[paste0("arm", a, "_fold", v)]] <- data.frame(
      arm = a, fold = v,
      learner = apply(fit_c$event.libraryNames, 1, paste, collapse = "_"),
      weight  = as.numeric(fit_c$event.coef),
      cv_risk = as.numeric(fit_c$event.cvRisk))
  }
}
SL_diagnostics_df <- bind_rows(SL_diagnostics)

# D.4.5 compares the reported fit against one with the censoring survival
# floored at G_FLOOR, which is Westling et al.'s remedy for weak positivity.
cens.pred.0.tr <- pmax(cens.pred.0, G_FLOOR)
cens.pred.1.tr <- pmax(cens.pred.1, G_FLOOR)
n_below <- sum(cens.pred.0 < G_FLOOR) + sum(cens.pred.1 < G_FLOOR)
cat(sprintf("truncation at %.3f alters %d of %d entries (%.1f%%)\n",
            G_FLOOR, n_below, 2 * length(cens.pred.0),
            100 * n_below / (2 * length(cens.pred.0))))

# the two doubly-robust fits (!)

run_dr <- function(c0, c1) {
  CFsurvival(
    time = Y_surv, event = Delta_surv, treat = A_surv, confounders = W_surv,
    fit.times = eval_grid,
    contrasts = c("surv.diff", "surv.ratio"),
    nuisance.options = list(
      eval.times   = eval_grid,
      event.pred.0 = manual_pred_0, event.pred.1 = manual_pred_1,
      cens.pred.0  = c0,            cens.pred.1  = c1,
      prop.pred    = manual_prop_scores,
      cross.fit    = FALSE))
}
fit_cf    <- run_dr(cens.pred.0,    cens.pred.1)
fit_cf_tr <- run_dr(cens.pred.0.tr, cens.pred.1.tr)

# VE(t) with pointwise intervals 
# VE(t) = 1 - F1(t)/F0(t) at every grid point, with a confidence interval
# built from the estimator's influence functions. 
ve_table <- function(fit) {
  grid <- fit$fit.times
  n_if <- nrow(fit$IF.vals.1)
  stopifnot(ncol(fit$IF.vals.1) == length(grid), n_if == n)
  
  # the fitted survival curves, one row per (arm, time), lined up to the grid
  S1L <- fit$surv.df %>% filter(trt == 1) %>% select(time, surv) %>% distinct()
  S0L <- fit$surv.df %>% filter(trt == 0) %>% select(time, surv) %>% distinct()
  S1g <- S1L$surv[match(grid, S1L$time)]
  S0g <- S0L$surv[match(grid, S0L$time)]
  stopifnot(!anyNA(S1g), !anyNA(S0g))
  
  IF_log <- matrix(NA_real_, n_if, length(grid))   # filled in as we go
  out <- lapply(seq_along(grid), function(j) {
    F1 <- 1 - S1g[j]; F0 <- 1 - S0g[j]             # cumulative incidence
    # theta = F1/F0 is undefined on the log scale if either is zero; the
    # point estimate is still 100% when only the vaccine arm is empty.
    if (F0 <= 0 || F1 <= 0)
      return(data.frame(time = grid[j], F0 = F0, F1 = F1, theta = NA_real_,
                        se_log = NA_real_, VE = if (F0 > 0) 100 else NA_real_,
                        VE_lo = NA_real_, VE_hi = NA_real_))
    theta <- F1 / F0                               # cumulative-incidence ratio
    # Delta method on the log scale. F_a = 1 - S_a, so the IF for F_a is minus
    # the IF for S_a; dividing the IF for theta by theta gives log(theta)'s.
    ifl <- ((-1 / F0) * fit$IF.vals.1[, j] +
              (F1 / F0^2) * fit$IF.vals.0[, j]) / theta
    IF_log[, j] <<- ifl                            # keep for the band
    se <- sqrt(mean(ifl^2) / n_if)                 # SE of log(theta)
    # interval built on the log scale, then mapped to VE = 1 - theta, which
    # flips the endpoints since VE decreases in theta
    data.frame(time = grid[j], F0 = F0, F1 = F1, theta = theta, se_log = se,
               VE    = (1 - theta) * 100,
               VE_lo = (1 - theta * exp(qnorm(.975) * se)) * 100,
               VE_hi = (1 - theta * exp(-qnorm(.975) * se)) * 100)
  }) %>% bind_rows()
  
  list(tab = out, IF_log = IF_log, grid = grid)
}
ve_main <- ve_table(fit_cf)      # reported analysis
ve_tr   <- ve_table(fit_cf_tr)   # truncated, for the D.4.5 comparison only
#  reporting interval (D.4.1) 
t0 <- as.numeric(quantile(event_times, T0_Q))
t1 <- as.numeric(quantile(event_times, T1_Q))
cat(sprintf("\nreporting interval: day %.0f to %.0f (quantiles %.2f-%.2f of event times)\n",
            t0, t1, T0_Q, T1_Q))
cat(sprintf("full grid runs to day %.0f; %d of %d grid points lie inside\n",
            max(eval_grid), sum(eval_grid >= t0 & eval_grid <= t1), length(eval_grid)))

##  simultaneous confidence band 
# A pointwise interval is only valid at one time. Reading a whole curve needs
# a band valid at every time at once, so the critical value comes from the
# distribution of the largest standardised deviation over the interval rather
# than from a normal quantile. The multiplier bootstrap simulates that. Give
# every participant an independent random normal, multiply their influence-
# function curve by it, and add the curves up. That sum is one plausible draw
# of the estimation error across the whole interval, and its largest absolute
# value is one draw of the maximum deviation. Repeat 2000 times.
uniform_band <- function(vt, t0, t1, n_sim = N_SIM, seed = SEED) {
  # only times inside the reporting interval, and only where VE is defined
  keep <- which(vt$grid >= t0 & vt$grid <= t1 & !is.na(vt$tab$theta))
  PHI <- vt$IF_log[, keep, drop = FALSE]   # participants x retained times
  sig <- sqrt(colMeans(PHI^2))             # pointwise SD at each time
  Zs  <- sweep(PHI, 2, sig, "/") / sqrt(nrow(PHI))   # standardise each column
  
  set.seed(seed)
  sup <- replicate(n_sim, max(abs(colSums(rnorm(nrow(PHI)) * Zs))))
  cc  <- as.numeric(quantile(sup, 0.95))   # 95th percentile of those maxima
  
  # same half-width construction as the pointwise interval, but with cc in
  # place of 1.96, then mapped from log(theta) back to the VE scale
  lr  <- log(vt$tab$theta[keep])
  hw  <- cc * sig / sqrt(nrow(PHI))
  data.frame(time = vt$grid[keep],
             VE    = (1 - exp(lr))      * 100,
             VE_lo = (1 - exp(lr + hw)) * 100,
             VE_hi = (1 - exp(lr - hw)) * 100,
             crit  = cc)
}
band_main <- uniform_band(ve_main, t0, t1)
# the band is wider than the pointwise interval by exactly this ratio
cat(sprintf("uniform-band critical value: %.3f (pointwise would be %.3f)\n",
            band_main$crit[1], qnorm(.975)))
rep_tab <- ve_main$tab %>%
  mutate(inside = time >= t0 & time <= t1) %>%
  left_join(band_main %>% select(time, band_lo = VE_lo, band_hi = VE_hi),
            by = "time")

# 7.2.2 quotes a single-event dip that falls between Table 11's rows
print(rep_tab %>% filter(time %in% c(1237, 1242)) %>%
        transmute(day = round(time), VE = round(VE, 2)), row.names = FALSE)
# Table 11
show_days <- c(562, 698, 784, 895, 970, 1075, 1176, 1301, 1428, 1501, 1581,
               max(eval_grid))
tab_out <- rep_tab %>%
  filter(time %in% show_days) %>%
  transmute(day = round(time), S1 = round(1 - F1, 4), S0 = round(1 - F0, 4),
            VE = round(VE, 2), pw_lo = round(VE_lo, 2), pw_hi = round(VE_hi, 2),
            band_lo = round(band_lo, 2), band_hi = round(band_hi, 2), inside)
cat("\n===== Table 11 =====\n")
print(tab_out, row.names = FALSE)
write.csv(tab_out, "tables/dr_results_table.csv", row.names = FALSE)

# incidence-rate summary at t1 (7.2.2, "Scale") 
# Each arm's survival curve is integrated to t1 by a left-Riemann sum over
# the step curve, giving restricted mean event-free survival times.
j1  <- which.min(abs(ve_main$grid - t1))
r1  <- ve_main$tab[j1, ]
tt  <- c(0, ve_main$grid[1:j1])
s1v <- c(1, 1 - ve_main$tab$F1[1:j1])
s0v <- c(1, 1 - ve_main$tab$F0[1:j1])
dt  <- diff(tt)
mu1 <- sum(dt * head(s1v, -1))
mu0 <- sum(dt * head(s0v, -1))
theta_IR <- (r1$F1 / mu1) / (r1$F0 / mu0)
cat(sprintf("\nrestricted mean event-free survival to day %.0f: mu0 %.1f | mu1 %.1f days\n",
            r1$time, mu0, mu1))
cat(sprintf("incidence-rate ratio %.4f -> VE %.2f%% (cumulative-incidence VE %.2f%%)\n",
            theta_IR, 100 * (1 - theta_IR), r1$VE))
cat("\n===== headline, end of the reporting interval =====\n")
cat(sprintf("day %.0f | S1 %.4f | S0 %.4f | F0 %.5f | F1 %.5f\n",
            r1$time, 1 - r1$F1, 1 - r1$F0, r1$F0, r1$F1))
cat(sprintf("VE %.2f%% | pointwise %.2f to %.2f | SE(log) %.4f\n",
            r1$VE, r1$VE_lo, r1$VE_hi, r1$se_log))
b1 <- band_main[which.min(abs(band_main$time - t1)), ]
cat(sprintf("uniform band %.2f to %.2f (critical value %.3f)\n",
            b1$VE_lo, b1$VE_hi, b1$crit))

# Figure 3: counterfactual survival and VE(t) 
col_vac  <- "#1B6D74"; col_ctrl <- "#C46B3E"
col_ve   <- "#2A2A33"; col_ref  <- "#B0403A"
# SE of a survival estimate at each time, from its own influence functions
se_surv <- function(IF) sqrt(colMeans(IF^2) / nrow(IF))
lgt  <- function(p) log(p / (1 - p))       # logit
ilgt <- function(x) 1 / (1 + exp(-x))      # inverse logit
z <- qnorm(.975)
s1 <- 1 - ve_main$tab$F1                   # counterfactual survival, vaccine
s0 <- 1 - ve_main$tab$F0                   # counterfactual survival, control
e1 <- se_surv(fit_cf$IF.vals.1)            # S = 1 - F, so same SE as F
e0 <- se_surv(fit_cf$IF.vals.0)
# Intervals are formed on the logit scale so they stay inside (0,1)
surv_long <- rbind(
  data.frame(time = ve_main$grid, surv = s1, arm = "Vaccine",
             lo = ilgt(lgt(s1) - z * e1 / (s1 * (1 - s1))),
             hi = ilgt(lgt(s1) + z * e1 / (s1 * (1 - s1)))),
  data.frame(time = ve_main$grid, surv = s0, arm = "Control",
             lo = ilgt(lgt(s0) - z * e0 / (s0 * (1 - s0))),
             hi = ilgt(lgt(s0) + z * e0 / (s0 * (1 - s0)))))
surv_long$arm <- factor(surv_long$arm, levels = c("Vaccine", "Control"))
# Figure 3 shows the reporting interval only; the rest is in Table 11
surv_long <- surv_long[surv_long$time >= t0 & surv_long$time <= t1, ]
y_lo <- min(surv_long$lo, na.rm = TRUE)    # axis floor, from the ribbons
pA <- ggplot(surv_long, aes(time, surv, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_step(linewidth = .9, direction = "hv") +
  scale_colour_manual(values = c(Vaccine = col_vac, Control = col_ctrl), name = NULL) +
  scale_fill_manual(values = c(Vaccine = col_vac, Control = col_ctrl), name = NULL) +
  scale_y_continuous("Event-free survival") +
  coord_cartesian(xlim = c(t0, t1), ylim = c(y_lo, 1)) +
  labs(tag = "A") +
  theme_minimal(base_size = 11) +
  theme(legend.position = c(.15, .25),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.height = unit(10, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        plot.tag = element_text(face = "bold", size = 13),
        plot.margin = margin(6, 10, 2, 6))
ve_in <- rep_tab %>% filter(time >= t0, time <= t1, !is.na(VE))
pB <- ggplot(ve_in, aes(time, VE)) +
  geom_hline(yintercept = POISSON_VE, colour = col_ref,
             linetype = "22", linewidth = .6) +
  annotate("text", x = min(ve_in$time), y = POISSON_VE - 7,
           label = sprintf("Poisson  %.2f%%", POISSON_VE),
           hjust = 0, size = 3.1, colour = col_ref) +
  geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = col_ve, alpha = .10) +
  geom_ribbon(aes(ymin = VE_lo,   ymax = VE_hi),   fill = col_ve, alpha = .18) +
  geom_line(linewidth = .9, colour = col_ve) +
  geom_point(size = 1.6, colour = col_ve) +
  scale_y_continuous("Vaccine efficacy (%)", breaks = seq(-50, 100, 25)) +
  scale_x_continuous("Time since randomisation (days)") +
  coord_cartesian(xlim = c(t0, t1), ylim = c(-50, 100)) +
  labs(tag = "B") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        plot.tag = element_text(face = "bold", size = 13),
        plot.margin = margin(2, 10, 6, 6))
fig <- pA / pB + plot_layout(heights = c(1, 1.15))
ggsave("figures/cfsurvival_ve_figure.pdf", fig, width = 6.8, height = 6.6)
ggsave("figures/cfsurvival_ve_figure.png", fig, width = 6.8, height = 6.6, dpi = 300)

# PART II: APPENDIX D.4

# Evaluation grid and reporting interval (D.4.1)
grid   <- fit_cf$fit.times          
inside <- grid >= t0 & grid <= t1   
# SuperLearner performance (D.4.2, Table 25) 
SL_tab <- SL_diagnostics_df %>%
  group_by(arm, learner) %>%
  summarise(mean_weight  = mean(weight), min_weight = min(weight),
            max_weight   = max(weight),  mean_cv_risk = mean(cv_risk),
            .groups = "drop") %>%
  arrange(arm, learner)
cat("\n===== Table 25: SuperLearner, censoring nuisance =====\n")
print(as.data.frame(SL_tab), row.names = FALSE, digits = 4)

# Positivity and censoring support (D.4.3)
pi1 <- as.numeric(manual_prop_scores)
cat(sprintf("\nestimated P(A=1|W): %.4f to %.4f (mean %.4f)\n",
            min(pi1), max(pi1), mean(pi1)))
#  Support and influence over follow-up (Table 26)
# survival = 1 - cumulative incidence, already worked out per grid point
S1g <- 1 - ve_main$tab$F1
S0g <- 1 - ve_main$tab$F0
zmax <- function(x) max(abs(x / sd(x, na.rm = TRUE)), na.rm = TRUE)
# One row per grid point: support and influence both change over follow-up
dr_tab <- lapply(seq_along(grid), function(j) {
  t   <- grid[j]
  IF0 <- fit_cf$IF.vals.0[, j]; IF1 <- fit_cf$IF.vals.1[, j]   # one per person
  F0  <- 1 - S0g[j]; F1 <- 1 - S1g[j]                          # risk in each arm
  ok  <- F0 > 0 && F1 > 0        # false early on as no vaccine-arm event yet
  # per-person pull on the log ratio; NA where the ratio is undefined
  IF_logRR <- if (ok) -IF1 / F1 + IF0 / F0 else rep(NA_real_, length(IF0))
  data.frame(day = round(t), inside = t >= t0 && t <= t1,
             minG0  = min(cens.pred.0[, j]),   # weakest censoring support here
             minG1  = min(cens.pred.1[, j]),
             nrisk0 = sum(Y_surv[A_surv == 0] >= t),   # still under observation
             nrisk1 = sum(Y_surv[A_surv == 1] >= t),
             VE = if (ok) 100 * (1 - F1 / F0) else NA_real_,
             maxz_logRR    = if (ok) zmax(IF_logRR) else NA_real_,  # in SDs
             mean_IF_logRR = if (ok) mean(IF_logRR) else NA_real_)  # should be ~0
}) %>% bind_rows()   # one row per time point, stacked into the table
write.csv(dr_tab, "tables/dr_support_eif.csv", row.names = FALSE)
# reported rows, plus everything outside the interval
show <- dr_tab$day %in% c(562, 784, 970, 1176, 1301, 1428, 1501, 1581) |
  !dr_tab$inside
cat("\n===== Table 26: censoring support and EIF concentration =====\n")
print(dr_tab[show, c("day", "minG0", "minG1", "nrisk0", "nrisk1",
                     "maxz_logRR", "inside")],
      row.names = FALSE, digits = 4)
cat(sprintf("largest |mean| of the log-ratio influence function: %.3g\n",
            max(abs(dr_tab$mean_IF_logRR), na.rm = TRUE)))

# Figure 11: EIF concentration over follow-up, reporting interval shaded.
pdf("figures/EIF_over_followup.pdf", width = 7, height = 5)
plot(dr_tab$day, dr_tab$maxz_logRR, type = "n",
     xlab = "Follow-up time (days)",
     ylab = "Maximum absolute standardised EIF")
rect(t0, par("usr")[3], t1, par("usr")[4], col = "grey93", border = NA)
lines(dr_tab$day, dr_tab$maxz_logRR, type = "b", pch = 20)
box(); dev.off()
# Figure 12: censoring support by arm, log scale, interval shaded.
pdf("figures/censoring_support.pdf", width = 7, height = 5)
plot(dr_tab$day, dr_tab$minG0, type = "n", log = "y",
     ylim = range(c(dr_tab$minG0, dr_tab$minG1)),
     xlab = "Follow-up time (days)",
     ylab = "Smallest estimated P(uncensored) (log scale)")
rect(t0, 10^par("usr")[3], t1, 10^par("usr")[4], col = "grey93", border = NA)
lines(dr_tab$day, dr_tab$minG0, type = "b", pch = 20)
lines(dr_tab$day, dr_tab$minG1, type = "b", pch = 1, lty = 2)
legend("bottomleft", c("Control", "Vaccine"), pch = c(20, 1),
       lty = c(1, 2), bty = "n")
box(); dev.off()

# Influence-function diagnostics (D.4.4) 
# At one time point see how much of the variance sits on the few most
# influential participants, and whether those few are the vaccine-arm events.
report_at <- function(at, label) {
  j <- which.min(abs(grid - at))               # nearest grid point
  F0t <- 1 - S0g[j]; F1t <- 1 - S1g[j]
  ifl <- -fit_cf$IF.vals.1[, j] / F1t + fit_cf$IF.vals.0[, j] / F0t
  sq  <- (ifl - mean(ifl))^2                   # each person's variance share
  ord <- order(sq, decreasing = TRUE)          # most influential first
  cat(sprintf("\n--- %s (day %.0f) ---\n", label, grid[j]))
  cat(sprintf("variance share, top 1 / top 5: %.1f%% / %.1f%%\n",
              100 * sq[ord[1]] / sum(sq), 100 * sum(sq[ord[1:5]]) / sum(sq)))
  cat(sprintf("of the 5 most influential, %d are vaccine-arm events\n",
              sum(A_surv[ord[1:5]] == 1 & Delta_surv[ord[1:5]] == 1)))
}
report_at(t1,        "end of the reporting interval")
report_at(max(grid), "last grid point, outside the interval")

# Sensitivity to truncation (D.4.5) 
# Flooring the censoring probabilities is the remedy Westling et al. suggest
# when positivity is violated.
for (at in c(t1, max(grid))) {
  j <- which.min(abs(ve_main$grid - at))
  cat(sprintf("day %.0f (%s): untruncated %.2f%% | G >= %.2f %.2f%% | diff %.2f pp\n",
              ve_main$grid[j], if (at == t1) "inside" else "outside",
              ve_main$tab$VE[j], G_FLOOR, ve_tr$tab$VE[j],
              ve_main$tab$VE[j] - ve_tr$tab$VE[j]))
}
cat("\nfigures/cfsurvival_ve_figure.{pdf,png}\n")
cat("figures/EIF_over_followup.pdf\n")
cat("figures/censoring_support.pdf\n")
cat("tables/dr_results_table.csv\n")
cat("tables/dr_support_eif.csv\n")