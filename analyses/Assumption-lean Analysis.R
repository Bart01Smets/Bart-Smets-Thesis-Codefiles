
################################
######Assumption-lean analysis##
################################
# Part I:Primary analysis, Part II: Appendix
library(dplyr)

# Move log(time) to the last column
df <- df %>%
  mutate(Y_time = log(Y_time)) %>%      
  relocate(Y_time, .after = last_col()) 

#The fit
SL.library <-c("SL.mean", "SL.glmnet")
fit_alglm <- alglm(
  data = df,
  family.a = "binomial",
  family.y = "poisson",
  SL.library.y = SL.library,
  SL.library.a = SL.library,
  k = 1
)

# extract irr
irr_alglm <- exp(fit_alglm$estimates$l13)
irr_alglm
# Extract CI's

se_alglm <- 1/sqrt(nrow(df))*sd(fit_alglm$estimates$eic.l13)
se_alglm
ci_alglm_upper <- exp(fit_alglm$estimates$l13+1.96*se_alglm)
ci_alglm_lower <- exp(fit_alglm$estimates$l13-1.96*se_alglm)

ci_alglm<-cbind(ci_alglm_lower, ci_alglm_upper)

# Vaccine efficacy
ve_alglm <- (1 - irr_alglm) * 100
ve_lower_alglm <- (1 - ci_alglm[2]) * 100 
ve_upper_alglm <- (1 - ci_alglm[1]) * 100 

# Print
cat(sprintf("IRR: %.4f (95%% CI: %.4f to %.4f)\n", irr_alglm, ci_alglm[1], ci_alglm[2]))
cat(sprintf("Vaccine Efficacy: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", ve_alglm, ve_lower_alglm, ve_upper_alglm))
#############################
#  APPENDIX: Section D.2.

#  Content:
#  SuperLearner performance (Table 17)
#  Fitted nuisance functions (D.2.1 text)
#  The convergence warning + SL.mean-only sensitivity (D.2.1 text)
#  Influence-function diagnostics (Table 18, Figure 7)
#  Overlap weights (Table 19)
#  Sensitivity to the nuisance-learning specification (Table 20)
#


Z <- qnorm(0.975)
n <- nrow(df)

# SuperLearner performance (Table 17) 

sl_print <- function(m, label) {
  fit <- m[[1]]
  cat("\n", label, "\n", sep = "")
  print(data.frame(Risk = fit$cvRisk, Weight = fit$coef), digits = 5)
}
sl_print(fit_alglm$nuis$my0.model, "Control outcome")
sl_print(fit_alglm$nuis$my1.model, "Vaccine outcome")
sl_print(fit_alglm$nuis$ma.model,  "Treatment")

#  Fitted nuisance functions 
# Min, max and mean of each fitted nuisance.
print(sapply(list("E[Y|A=1,W]" = fit_alglm$nuis$my1,
                  "E[Y|A=0,W]" = fit_alglm$nuis$my0,
                  "P(A=1|W)"   = fit_alglm$nuis$pia),
             summary), digits = 4)
cat(sprintf("observed P(A=1): %.4f\n", mean(df$A_treat)))
#  The convergence warning + SL.mean-only sensitivity 
# (because SL.mean alone contains no glm fit. )
fit_mean <- alglm(data = df, family.a = "binomial", family.y = "poisson",
                  SL.library.y = "SL.mean", SL.library.a = "SL.mean", k = 1)
cat(sprintf("SL.mean-only fit: l13 %.5f | VE %.2f%%\n",
            fit_mean$estimates$l13,
            (1 - exp(fit_mean$estimates$l13)) * 100))

# Influence-function diagnostics (Table 18, Figure 7)
eic  <- fit_alglm$estimates$eic.l13
beta <- fit_alglm$estimates$l13
stopifnot(length(eic) == n)
se_al <- sd(eic) / sqrt(n)
ve_al <- (1 - exp(beta)) * 100
# VE = 1 - exp(beta) decreases in beta, so the upper limit of beta gives the
# lower limit of VE: hence c(1, -1) here, not the usual c(-1, 1).
ci_al <- (1 - exp(beta + c(1, -1) * Z * se_al)) * 100
cat(sprintf("l13 %.5f | SE %.5f | VE %.2f%% (95%% CI %.2f%% to %.2f%%)\n",
            beta, se_al, ve_al, ci_al[1], ci_al[2]))

# The EIF is mean-zero in theory, so standardise by the SD without centring.
zs  <- eic / sd(eic)                # each participant's pull, in SDs
sq  <- (eic - mean(eic))^2          # each one's contribution to the variance
ord <- order(sq, decreasing = TRUE) # participants, most influential first
top <- ord[1:5]                     # the five most influential
shr <- function(k) 100 * sum(sq[ord[seq_len(k)]]) / sum(sq)  # % of variance from the top k

print(data.frame(
  Diagnostic = c("SD of EIF", "Minimum EIF", "Maximum EIF",
                 "Max |standardised EIF|", "N with |z| > 3",
                 "Variance share, top 1", "Variance share, top 5",
                 "Vaccine-arm events in top 5"),
  Value = c(sd(eic), min(eic), max(eic), max(abs(zs)),
            sum(abs(zs) > 3), shr(1), shr(5),
            sum(df$A_treat[top] == 1 & df$Delta_event[top] == 1))),
  row.names = FALSE, digits = 6)

# Plot EIF
pdf("figures/AL_EIF_influence.pdf", width = 7, height = 5)
plot(seq_along(eic), eic, pch = 20, cex = .4, col = "grey40",
     xlab = "Participant", ylab = "Estimated influence-function value")
abline(h = c(-3, 3) * sd(eic), lty = 3)
ev <- which(df$A_treat == 1 & df$Delta_event == 1)
points(ev, eic[ev], pch = 19, cex = .9)
legend("topleft", pch = 19, legend = "vaccine-arm event", bty = "n")
dev.off()

# Overlap weights (Table 19) 

pi_hat <- as.numeric(fit_alglm$nuis$pia)
w      <- pi_hat * (1 - pi_hat)
qs     <- c(0, .05, .25, .5, .75, .95, 1)
ow <- rbind(pi = quantile(pi_hat, qs), w = quantile(w, qs))
print(round(ow, 4))
write.csv(ow, "tables/al_overlap_weights.csv")

cat(sprintf("\nweight range: %.4f to %.4f (factor %.3f)\n",
            min(w), max(w), max(w) / min(w)))
cat(sprintf("outside [0.4, 0.6]: %d participants\n",
            sum(pi_hat < 0.4 | pi_hat > 0.6)))
cat(sprintf("share of total weight held by the most-weighted quarter: %.2f%%\n",
            100 * sum(sort(w, decreasing = TRUE)[seq_len(floor(n/4))]) / sum(w)))

# Sensitivity to the nuisance-learning specification (Table 20) 

# Run the main analysis with the following superlearner libraries:

LIBS <- list(
  "Regularised"                   = c("SL.mean", "SL.glmnet", "SL.ridge", "SL.glmnet.0.5"),
  "Bayesian"                      = c("SL.mean", "SL.glmnet", "SL.bayesglm"),
  "Expanded regularised/Bayesian" = c("SL.mean", "SL.glmnet", "SL.ridge",
                                      "SL.glmnet.0.5", "SL.bayesglm"),
  "Full flexible"                 = c("SL.glm", "SL.mean", "SL.glmnet",
                                      "SL.ranger", "SL.gam"),
  "Full without SL.ranger"        = c("SL.glm", "SL.mean", "SL.glmnet", "SL.gam"))
