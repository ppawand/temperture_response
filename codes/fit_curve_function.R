# Fit one A/Ci curve with plantecophys::fitaci

#   1. Tleaf must be the MEASURED leaf temperature of each point, not the
#      set point. fitaci uses it to calculate Km and Gamma* (Bernacchi et al.
#      2001), which change ~9% per °C at 45 °C, so a wrong temperature biases
#      Vcmax and Jmax.
#   2. Atmospheric pressure from the LI-6800 (column Pa, kPa) is passed in;
#      fitaci otherwise assumes 100 kPa.
#   3. The output also records the parameter SEs, fit RMSE, the curve's mean
#      measured leaf temperature and its set point, so the temperature-response
#      fits downstream use the real leaf temperature.

fit_curve <- function(data = final_aci_data, id_name = "id", fitmethod = "default",
                        citrans = NULL) {
  d <- data %>% filter(id == id_name)
  
  patm <- if ("Pa" %in% names(d)) mean(d$Pa, na.rm = TRUE) else 100
  
  curve <- fitaci(d,
                  varnames = list(ALEAF = "A",
                                  Tleaf = "Tleaf",      # measured, per point
                                  Ci    = "Ci",
                                  PPFD  = "Qin"),
                  Patm = patm,
                  fitmethod = "default",
                  citransition = citrans,
                  fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
  
  se <- function(p) if (p %in% rownames(curve$pars)) curve$pars[p, 2] else NA
  
  coef_data <- data.frame(id = id_name, t(coef(curve)),
                          Vcmax_SE   = se("Vcmax"),
                          Jmax_SE    = se("Jmax"),
                          RMSE       = if (!is.null(curve$RMSE)) curve$RMSE else NA,
                          n_points   = nrow(d),
                          Tleaf_mean = mean(d$Tleaf),
                          Tleaf_sd   = sd(d$Tleaf),
                          Tleaf_set  = if ("temp_setpoint" %in% names(d))
                            d$temp_setpoint[1] else NA,
                          citrans    = if (is.null(citrans)) NA else citrans)
  return(list(curve, coef_data))
}

