# CPLOT analysis code. Part I: Main analysis, Part II: Appendix
T_EVAL <- 50
## Main analysis

# Discretize to months and cap at the horizon (month 50).
analysis_df <- analysis_df %>%
  mutate(
    time_month_raw  = floor(PersonTime / 30.44),                 # month of event or censoring
    t_lastobs_month = pmin(floor(t_ice_capped / 30.44), T_EVAL), # ICE-free window, capped
    time_month_out  = pmin(time_month_raw, T_EVAL),              # outcome grid, capped
    # An event after month 50 did not occur within the horizon, so that person
    # contributes Y(s) = 0 at every s <= 50 
    Event_t         = as.integer(Event == 1 & time_month_raw <= T_EVAL)
  )

# One row per participant, wide Y0..Y50: each person's event flag is placed
# in the column for the month it happened, every other month gets a 0.
analysis_widerr <- analysis_df %>%
  mutate(patient_id      = row_number(),
         time_month      = time_month_out,      # month the column name comes from
         max_time_month  = t_lastobs_month) %>% # T_i, the ICE-free window
  pivot_wider(
    id_cols      = c(patient_id, max_time_month, Treat, any_of(colnames(W_covars_clean))),
    names_from   = time_month, values_from = Event_t,            # month -> column, flag -> value
    names_prefix = "Y", names_sort = TRUE,                       # Y0, Y1, ... in order
    values_fn    = sum, values_fill = list(Event_t = 0)          # months w/o event -> 0
  )

# Enforce a complete Y0..Y50 grid
# pivot_wider only creates a column for months that actually occur in the
# data, so a month in which nobody was censored or had an event would be
# missing. Fill those with 0 and fix the column order.
ycols          <- paste0("Y", 0:T_EVAL)                          # the grid CPLOT expects
missing_months <- setdiff(ycols, names(analysis_widerr))         # months the pivot skipped
if (length(missing_months))
  message("no observations in month(s): ", paste(missing_months, collapse = ", "))
analysis_widerr[missing_months] <- 0                             # nobody at risk there -> 0
analysis_widerr <- analysis_widerr %>%
  select(patient_id, max_time_month, Treat,                      # keep the id columns first,
         any_of(colnames(W_covars_clean)), all_of(ycols))        # then Y0..Y50 in order

# Cumulative outcome Y(s): carry the first event forward 
# The pivot puts a 1 only in the participant's event month. CPLOT needs
# "has the event occurred by month s", so the first 1 is carried forward.
analysis_widerr[ycols] <- t(apply(as.matrix(analysis_widerr[, ycols]), 1, cummax))
cat("Events by treatment arm within the horizon:\n")
print(table(Event = as.integer(rowSums(analysis_widerr[, ycols]) > 0),
            Treat = analysis_widerr$Treat))

# Fit CPLOT nuisance models 
source("CPLOTcode.R")
family    <- "poisson"
SL.library <- c("SL.mean", "SL.glmnet")                          # restricted lib because of sparse vaccine arm events
fit_nuisance <- CFnuisance(
  as.data.frame(analysis_widerr), t = T_EVAL, nfolds = 1,
  treatment = "Treat", colnames(W_covars_clean), lotp = "max_time_month",
  family, SL.library, RCT = TRUE                                 # no cross-fitting
)

# Estimate CPLOT at month 50 
cplot <- CPLOT(as.data.frame(analysis_widerr), T_EVAL, treatment = "Treat", "max_time_month", fit_nuisance)
str(cplot)
# Log-ratio -> VE-like summary 
lr    <- cplot$logratio$effect
se_lr <- cplot$logratio$se
lr_low  <- lr - qnorm(.975) * se_lr
lr_high <- lr + qnorm(.975) * se_lr
p_lr    <- 2 * pnorm(-abs(lr / se_lr))
VE      <- (1 - exp(lr))      * 100
VE_low  <- (1 - exp(lr_high)) * 100                              # VE is decreasing in lr -> bounds reverse
VE_high <- (1 - exp(lr_low))  * 100
# Print results 
cat(sprintf("\nCPLOT log-ratio at month %d: %.4f (95%% CI: %.4f to %.4f), p = %.4g\n",
            T_EVAL, lr, lr_low, lr_high, p_lr))
cat(sprintf("CPLOT VE-like summary at month %d: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n",
            T_EVAL, VE, VE_low, VE_high))
###############
#  CPLOT APPENDIX DIAGNOSTICS (D.6):
#  Influence-function diagnostics (Table 29)
#  Sensitivity to the evaluation horizon (D.6.2, Figure 15)
#  Inference without cross-fitting (D.6.3)
#  Sensitivity to the nuisance library (Table 30)

Z      <- qnorm(0.975)
cp_dat <- as.data.frame(analysis_widerr)
event <- as.integer(rowSums(analysis_widerr[, ycols]) > 0)
vax_event <- analysis_widerr$Treat == 1 & event == 1   # the three vaccine-arm events

# Influence-function diagnostics (Table 29) 
# Inference is on log(psi1/psi0) with the delta method
psi1 <- as.numeric(cplot$psi1)[1]
psi0 <- as.numeric(cplot$psi0)[1]
stopifnot(psi1 > 0, psi0 > 0)   # a negative arm functional is the failure described in section 6.3.
eif  <- as.numeric(cplot$EIF1)/psi1 - as.numeric(cplot$EIF0)/psi0
cor10 <- cor(as.numeric(cplot$EIF1), as.numeric(cplot$EIF0))
# The delta-method SE should reproduce the one the estimation routine returns.
cat(sprintf("SE from delta-method EIF: %.4f | SE reported by CPLOT(): %.4f\n",
            sd(eif)/sqrt(length(eif)), se_lr))
z   <- (eif - mean(eif)) / sd(eif)   # each participant's pull, in SDs
sq  <- (eif - mean(eif))^2           # each one's contribution to the variance
ord <- order(sq, decreasing = TRUE)  # participants, most influential first
shr <- function(k) 100 * sum(sq[ord[seq_len(k)]]) / sum(sq)
print(data.frame(
  Diagnostic = c("Log-ratio", "Standard error", "p-value",
                 "Correlation between arm-specific EIFs",
                 "Max |standardised EIF|", "N with |z| > 3", "Prop with |z| > 3 (%)",
                 "Vaccine-arm events among five most influential",
                 "Variance share, top 1 (%)", "Variance share, top 5 (%)",
                 "Variance share, top 10 (%)"),
  Value = c(lr, se_lr, p_lr, cor10,
            max(abs(z)), sum(abs(z) > 3), 100 * mean(abs(z) > 3),
            sum(vax_event[ord[1:5]]),
            shr(1), shr(5), shr(10))),
  row.names = FALSE, digits = 4)

# Sensitivity to the evaluation horizon (D.6.2, Figure 15)
# The same nuisance fit is reused at every horizon. CFnuisance is not
# re-run (so not multiple independent analyses)
grid <- 0:T_EVAL
res  <- do.call(rbind, lapply(grid, function(s) {
  # one CPLOT call per horizon s, on the nuisances already fitted at t = 50
  f  <- CPLOT(cp_dat, s, treatment = "Treat", "max_time_month", fit_nuisance)
  p1 <- f$psi1; p0 <- f$psi0
  # before month 13 the vaccine arm has no events, so psi1 = 0 and the
  # log-ratio is -Inf: those horizons are dropped rather than plotted
  if (!isTRUE(p1 > 0 && p0 > 0)) return(NULL)
  e <- f$EIF1/p1 - f$EIF0/p0        # delta-method EIF, as in Table 29
  data.frame(s = s, lr = log(p1/p0), se = sd(e)/sqrt(length(e)),
             # cumulative events per arm at s: used to drop horizons with an
             # empty arm and to locate the month the third vaccine event enters
             ev_vax = sum(cp_dat[[paste0("Y", s)]][cp_dat$Treat == 1]),
             ev_ctl = sum(cp_dat[[paste0("Y", s)]][cp_dat$Treat == 0]))
}))

# Marginal distribution of common observation times for a randomly selected
# treated-control pair.
s1 <- cp_dat$max_time_month[cp_dat$Treat == 1]        # T_i, vaccine arm
s0 <- cp_dat$max_time_month[cp_dat$Treat == 0]        # T_i, control arm
# D_s = P(both still ICE-free at s) = S1(s)*S0(s). The two members of a pair
# are independent draws, so this is the marginal version, as the caption says.
D_full    <- vapply(grid, function(s) mean(s1 >= s) * mean(s0 >= s), numeric(1))
# P(M = s) = D_s - D_{s+1}; the last entry is the tail P(M >= 50), equal to
# P(M = 50) because max_time_month is capped at the horizon. Sums to 1.
dens_full <- c(-diff(D_full), D_full[length(D_full)])
res$dens  <- dens_full[match(res$s, grid)]      # keep only the plotted horizons
# Second check on the event counts themselves: the psi > 0 guard in the loop
# has already dropped the horizons with an empty arm.
res <- res[res$ev_vax >= 1 & res$ev_ctl >= 1, ]
# Inference is on the log-ratio; VE is a display transformation of it.
res$VE    <- (1 - exp(res$lr)) * 100
res$VE_lo <- (1 - exp(res$lr + Z * res$se)) * 100   # VE decreases in lr,
res$VE_hi <- (1 - exp(res$lr - Z * res$se)) * 100   # so the bounds reverse
# Months 0-12 have no vaccine-arm event, so the log-ratio is undefined there:
# the sweep starts at 13.
cat(sprintf("\n%d evaluation horizons, month %d to %d\n",
            nrow(res), min(res$s), max(res$s)))
# the four horizons quoted in the text, with the SE falling over follow-up
for (s in c(13, 18, 30, 50))
  cat(sprintf("s = %2d: VE %6.2f%% | SE %.3f\n", s,
              res$VE[res$s == s], res$se[res$s == s]))
# First month at which the cumulative vaccine-arm count reaches its maximum,
# i.e. when the third event enters. After that no new treated-arm information
# arrives, so the remaining movement is horizon effect alone.
m_last <- min(res$s[res$ev_vax == max(res$ev_vax)])
cat(sprintf("last vaccine-arm event enters at month %d; VE from there on: %.1f%% to %.1f%%\n",
            m_last, min(res$VE[res$s >= m_last]), max(res$VE[res$s >= m_last])))
# how much of the estimand rests on the end-of-trial block rather than on
# clinically determined observation times
cat(sprintf("months 45-48 carry %.1f%% of the estimand weight\n",
            100 * sum(dens_full[match(45:48, grid)])))

# Figure 15. VE is unbounded below, so the axis is fixed to [0, 100]:
# limits outside it are clipped, and point estimates outside it are drawn as
# open downward triangles on the axis rather than silently dropped by ggplot.
VE_YLIM <- c(0, 100)
d <- res
d$trunc <- d$VE < VE_YLIM[1] | d$VE > VE_YLIM[2]
d$VE    <- pmin(pmax(d$VE, VE_YLIM[1]), VE_YLIM[2])
d$VE_lo <- pmax(d$VE_lo, VE_YLIM[1])
d$VE_hi <- pmin(d$VE_hi, VE_YLIM[2])
cat(sprintf("clipped: %d of %d confidence limits and %d of %d point estimates (s = %s)\n",
            sum(res$VE_lo < VE_YLIM[1]) + sum(res$VE_hi > VE_YLIM[2]), 2 * nrow(res),
            sum(d$trunc), nrow(d), paste(d$s[d$trunc], collapse = ", ")))
off <- VE_YLIM[1]; scl <- diff(VE_YLIM) / max(d$dens)
p15 <- ggplot(d, aes(x = s)) +
  # ribbon rather than area: geom_area fills from y = 0 and inverts on an
  # axis that does not contain zero
  geom_ribbon(aes(ymin = off, ymax = dens * scl + off), fill = "#9ecae1",
              alpha = .45, colour = "#6baed6", linewidth = .3) +
  geom_linerange(aes(ymin = VE_lo, ymax = VE_hi), colour = "grey40", linewidth = .35) +
  geom_point(data = d[!d$trunc, ], aes(y = VE), size = 1.6) +
  geom_point(data = d[d$trunc, ], aes(y = VE), shape = 25, size = 2,
             fill = "white", colour = "black") +
  scale_y_continuous(name = "CPLOT VE-like estimate at horizon s (%)", limits = VE_YLIM,
                     sec.axis = sec_axis(~ (. - off) / scl,
                                         name = "Density of first ICE time in a matched pair",
                                         labels = scales::percent_format(accuracy = 0.1))) +
  labs(x = "Evaluation horizon s (months)") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("figures/cplot_horizon_VE.pdf", p15, width = 7, height = 4.2)
# Note: this is Phi_s, the CPLOT estimand at horizon s, not the de-attenuated
# psi_s of Baklicharov et al. 

# Inference without cross-fitting (D.6.3) 
# Cross-fitting was attempted with K = 2 and could not be used: the
# vaccine-arm functional comes back negative, i.e. outside [0, 1].
nuis_cf2 <- CFnuisance(cp_dat, t = T_EVAL, nfolds = 2, treatment = "Treat",
                       colnames(W_covars_clean), lotp = "max_time_month",
                       "poisson", c("SL.mean", "SL.glmnet"), RCT = TRUE)
cf2 <- CPLOT(cp_dat, T_EVAL, treatment = "Treat", "max_time_month", nuis_cf2)
cat(sprintf("\npsi1: %.5f with K = 2 cross-fitting against %.5f without\n",
            as.numeric(cf2$psi1)[1], psi1))
cat(sprintf("psi0: %.5f against %.5f (%.2f%% change)\n",
            as.numeric(cf2$psi0)[1], psi0,
            100 * abs(as.numeric(cf2$psi0)[1] - psi0) / psi0))
#Sensitivity to the nuisance library (Table 30)
cat(sprintf("restricted library: 3 vaccine-arm events hold %.1f%% of the EIF variance\n",
            100 * sum(z[vax_event]^2) / sum(z^2)))
# Row 2: same horizon, same capped clock, flexible outcome library.
nuis_flex <- CFnuisance(cp_dat, t = T_EVAL, nfolds = 1, treatment = "Treat",
                        colnames(W_covars_clean), lotp = "max_time_month",
                        "poisson", c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"),
                        RCT = TRUE)
flex <- CPLOT(cp_dat, T_EVAL, treatment = "Treat", "max_time_month", nuis_flex)
p1_f  <- flex$psi1; p0_f <- flex$psi0
eif_f <- flex$EIF1/p1_f - flex$EIF0/p0_f
lr_f  <- log(p1_f/p0_f)
se_f  <- sd(eif_f)/sqrt(length(eif_f))
q_f   <- (eif_f - mean(eif_f))^2          # contribution to the variance
top_f <- which.max(q_f)                   # single most influential participant
cat(sprintf("flexible library: VE %.1f%% (%.1f to %.1f) | SE %.3f\n",
            (1 - exp(lr_f)) * 100,
            (1 - exp(lr_f + Z * se_f)) * 100,   # VE decreases in lr
            (1 - exp(lr_f - Z * se_f)) * 100, se_f))
cat(sprintf("flexible library: 3 vaccine-arm events hold %.1f%% of the EIF variance | leading contributor is a %s (%.1f%%)\n",
            100 * sum(q_f[vax_event]) / sum(q_f),
            if (vax_event[top_f]) "vaccine-arm event" else "non-event",
            100 * q_f[top_f] / sum(q_f)))