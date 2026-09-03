
#  Cox models as comparison for the poisson analysis (thesis D.1.2)

form_cox <- as.formula(paste("Surv(PersonTime, Event) ~ Treat +", covariate_names))

# Robust standard errors in both: in the weighted fit because the weights
# are estimated, and in the unweighted fit so the two are on the same footing.
cox_unw <- coxph(form_cox, data = analysis_ipcw, ties = "efron", robust = TRUE)
cox_w   <- coxph(form_cox, data = analysis_ipcw, ties = "efron", robust = TRUE,
                 weights = analysis_ipcw$ipcw_weight)

# Table 16
cox_row <- function(fit, label) {
  s  <- summary(fit)
  hr <- s$conf.int["Treat", c("exp(coef)", "lower .95", "upper .95")]
  data.frame(Model = label,
             HR = round(hr[1], 4),
             CI = sprintf("%.4f-%.4f", hr[2], hr[3]),
             VE = round(100 * (1 - hr[1]), 2),
             VE_CI = sprintf("%.2f-%.2f", 100 * (1 - hr[3]), 100 * (1 - hr[2])),
             C = round(s$concordance[1], 3))
}
cat("\nTable 16\n")
print(rbind(cox_row(cox_unw, "Unweighted"), cox_row(cox_w, "IPCW-weighted")),
      row.names = FALSE)
cat(sprintf("n = %d | events %d (vaccine %d, control %d)\n",
            cox_unw$n, cox_unw$nevent,
            sum(analysis_ipcw$Event == 1 & analysis_ipcw$Treat == 1),
            sum(analysis_ipcw$Event == 1 & analysis_ipcw$Treat == 0)))

# Weighting moves the estimate very little
ve_unw <- 100 * (1 - exp(coef(cox_unw)["Treat"]))
ve_w   <- 100 * (1 - exp(coef(cox_w)["Treat"]))
cat(sprintf("weighting moves VE by %.2f pp; unweighted Cox sits %.2f pp from the Poisson estimate\n",
            abs(ve_w - ve_unw), abs(ve_unw - ve)))

# Proportional hazards in the outcome model. 
ph_out <- cox.zph(cox_unw)
g <- ph_out$table["GLOBAL", ]
cat(sprintf("global Schoenfeld: chi2 = %.2f on %.0f df, p = %.2f | terms failing at 5%%: %d of %d\n",
            g["chisq"], g["df"], g["p"],
            sum(ph_out$table[rownames(ph_out$table) != "GLOBAL", "p"] < 0.05),
            nrow(ph_out$table) - 1))