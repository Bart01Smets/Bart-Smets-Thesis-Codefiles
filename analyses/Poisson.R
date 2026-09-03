# CONVENTIONAL POISSON ANALYSIS + DIAGNOSTICS (thesis Appendix D.1)

# PRIMARY ANALYSIS 
cat("PRIMARY ANALYSIS: Poisson\n")

formula_exact <- "Event ~ Treat + AGE + CTNNorthern.Europe + CTNWestern.Europe + CTNNorthern.America + CTNSouthern.America + CTNSoutheastern.Asia + CTNEastern.Asia + RACE2 + RACE3 + RACE4 + HORM_FLGY + CYTOD0P + NB_PARLess.than.3 + SMKY + NB_PAR_I1 + offset(log(PersonTime))"

mod <- glm(as.formula(formula_exact), family = poisson(link = "log"), data = analysis_df)

irr        <- exp(coef(mod)["Treat"])
robust_cov <- vcovHC(mod, type = "HC0")               # sandwich SEs
ci         <- exp(coefci(mod, vcov = robust_cov)["Treat", ])
ve         <- (1 - irr) * 100
ve_lower   <- (1 - ci[2]) * 100
ve_upper   <- (1 - ci[1]) * 100

cat(sprintf("IRR: %.4f (95%% CI: %.4f to %.4f)\n", irr, ci[1], ci[2]))
cat(sprintf("Vaccine Efficacy: %.1f%% (95%% CI: %.1f%% to %.1f%%)\n", ve, ve_lower, ve_upper))
##############################
# Appendix D.1.1
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

# Influence: Cook's distance + DFBETAS for the treatment coefficient
cook <- cooks.distance(mod)
dfb  <- as.matrix(dfbetas(mod))
dft  <- dfb[, "Treat"]


cat(sprintf("max Cook's distance: %.4g\nmax |DFBETA| for Treat: %.4g\n",
            max(cook, na.rm = TRUE), max(abs(dft), na.rm = TRUE)))

# Table 15: five observations with the largest |DFBETA|
idx  <- head(order(abs(dft), decreasing = TRUE), 5)
infl <- data.frame(Observation = idx,
                   Event = analysis_df$Event[idx],
                   Treat = analysis_df$Treat[idx],
                   Cook = cook[idx], DFBETA_Treat = dft[idx], row.names = NULL)
print(infl, row.names = FALSE, digits = 4)
write.csv(infl, "tables/poisson_influential_observations.csv", row.names = FALSE)

# Figure 6: DFBETA against participant index 
pdf("figures/poisson_dfbetas_treat.pdf", width = 7, height = 5)
plot(dft, type = "h", xlab = "Participant", ylab = "DFBETA for treatment coefficient")
abline(h = 0, lty = 2)
dev.off()

# Functional form of age: linear vs. natural spline (df = 3)
mod_spline <- update(mod, . ~ . - AGE + ns(AGE, df = 3))
lrt <- anova(mod, mod_spline, test = "LRT")


print(lrt)
cat(sprintf("AIC: linear = %.1f, spline = %.1f\n", AIC(mod), AIC(mod_spline)))

ve_linear <- (1 - exp(coef(mod)["Treat"])) * 100
ve_spline <- (1 - exp(coef(mod_spline)["Treat"])) * 100
cat(sprintf("VE: linear = %.3f%%, spline = %.3f%%\n", ve_linear, ve_spline))
