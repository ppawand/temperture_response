###############################################################################
# Gas exchange at ambient CO2: Asat, gsw, Ci/Ca
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(lme4)
library(car)
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(egg)

###############################################################################
# prepare the variables
###############################################################################
# Leaf-temperature set point that matches each growth regime
# ( 34 and 40 °C as "growth temperature")
growth_T <- c(Ambient = 34, Heat = 40)

# Reference CO2 window defining "ambient" (µmol mol-1)
co2_window <- c(419.5, 420.5)

# re-measurements that replaced the original curve...
use_remeasure_pima <- c("4_16_re", "4_22_re", "12_46_re", "20_16re")
use_remeasure_upland <- c("47_28re", "52_16re")



###############################################################################
# 1. Raw data, same curves as the A/Ci fits
###############################################################################
pima_aci_data <- read.csv("Final_data/pima_final.csv")

pima_raw <- pima_aci_data |> 
  dplyr::rename(plant_id = plants)

pima_raw <- pima_raw %>%
  filter(!grepl("re$", id) | id %in% use_remeasure_pima) |> 
  filter(
    # Remove the original curves corresponding to the kept remeasured curves
    !id %in% sub("_re$", "", use_remeasure_pima)
  )
pima_raw <- pima_raw %>%
  mutate(id = sub("_?re$", "", id))

pima_raw <- pima_raw %>%
  group_by(id) %>%
  mutate(row_id = row_number()) %>%
  filter(
    (id == "1_46"  & row_id <= 96) |
      (id == "20_28" & row_id > 96) |
      (!id %in% c("1_46", "20_28"))
  ) %>%
  ungroup() %>%
  select(-row_id)

pima_raw <- meta_pima |> 
  left_join(pima_raw, by = "plant_id")

# clean upland

upland_aci_data <- read.csv("Final_data/upland_final.csv")

upland_raw <- upland_aci_data |> 
  dplyr::rename(plant_id = plants)

upland_raw <- upland_raw %>%
  filter(!grepl("re$", id) | id %in% use_remeasure_upland) |> 
  filter(
    # Remove the original curves corresponding to the kept remeasured curves
    !id %in% sub("re$", "", use_remeasure_upland)
  )


upland_raw <- upland_raw %>%
  mutate(id = sub("_?re$", "", id))

upland_raw <- meta_upland |> 
  left_join(upland_raw, by = "plant_id")

# combine both

raw_all <- bind_rows(pima_raw, upland_raw)

raw_all[, c("plant_date", "Species", "irrigation_trt", "growth_temp", "treatment", "plant_id", "temp_setpoint")] <-
  lapply(raw_all[, c("plant_date", "Species", "irrigation_trt", "growth_temp", "treatment", "plant_id", "temp_setpoint")], factor)

str(raw_all)

###############################################################################
# Values at ambient CO2: first contiguous block of points at ~420 ppm
#   (the steady-state hold before the CO2 ramp, not points the ramp
#   passes through later)
###############################################################################

ambient <- raw_all %>%
  group_by(Species, plant_id, id, temp_setpoint, treatment,
           growth_temp, irrigation_trt, plant_date) %>%
  mutate(in_win = CO2_r > co2_window[1] & CO2_r < co2_window[2],
         block  = cumsum(in_win != lag(in_win, default = FALSE))) %>%
  filter(in_win) %>%
  filter(block == min(block)) %>%
  summarise(n_pts     = n(),
            first_row = min(row_number()),
            A       = mean(A),
            gsw     = mean(gsw),
            Ci      = mean(Ci),
            Ca      = mean(CO2_s),
            VPDleaf = mean(VPDleaf),
            Qin     = mean(Qin),
            Tleaf   = mean(Tleaf),
            Patm    = mean(Pa),
            .groups = "drop") %>%
  mutate(CiCa = Ci / Ca,
         iWUE = A / gsw)        # intrinsic WUE, µmol CO2 mol-1 H2O


###############################################################################
# 4. Labels
###############################################################################

label_factors <- function(df) {
  df %>% mutate(
    growth_temp = factor(dplyr::recode(growth_temp,
                                       "Ambient" = "Low growth temperature",
                                       "Heat"    = "High growth temperature"),
                         levels = c("Low growth temperature", "High growth temperature")),
    irrigation_trt = factor(dplyr::recode(irrigation_trt,
                                          "Irrigated" = "High water availability",
                                          "Drought"   = "Low water availability"),
                            levels = c("High water availability", "Low water availability")),
    Species    = factor(dplyr::recode(Species, "pima" = "Pima Cotton",
                                      "upland" = "Upland Cotton")),
    plant_date = factor(plant_date),
    plant_id  = factor(plant_id))
}

# Growth temperature of each plant, as a leaf-temperature set point
ambient <- ambient %>%
  mutate(growth_set = unname(growth_T[as.character(growth_temp)]))

amb <- label_factors(ambient)

###############################################################################
# A at growth temperature (one value per plant)
###############################################################################

at_growth <- amb %>% filter(temp_setpoint == growth_set)

fit_cell_model <- function(df, y, log_y = FALSE) {
  lhs <- if (log_y) paste0("log(", y, ")") else y
  f <- as.formula(paste(lhs, "~ growth_temp * irrigation_trt * Species + (1 | plant_date)"))
  m <- lmer(f, data = df)
  attr(m, "data") <- df
  attr(m, "y") <- y
  m
}

growth_models <- list(
  A_growth    = fit_cell_model(at_growth, "A"),
  gsw_growth  = fit_cell_model(at_growth, "gsw"),
  CiCa_growth = fit_cell_model(at_growth, "CiCa"),
  iWUE_growth = fit_cell_model(at_growth, "iWUE")
)


# Diagnostics: one DHARMa panel per model
for (nm in names(growth_models)) {
  pdf(file.path("Results/Dharma", paste0("diagnostics_", nm, ".pdf")), width = 9, height = 5)
  plot(simulateResiduals(growth_models[[nm]]), title = nm)
  dev.off()
}


###############################################################################
#    A-T curve shift per plant (June et al. 2004)
###############################################################################

fit_june <- function(d) {
  d <- d[is.finite(d$A), ]
  out <- tibble(Topt_A = NA_real_, Topt_A_se = NA_real_, Aopt = NA_real_,
                Omega = NA_real_, converged = FALSE,
                T_min = min(d$Tleaf), T_max = max(d$Tleaf), n_T = nrow(d))
  if (nrow(d) < 4) return(out)
  starts <- expand.grid(Topt = c(d$Tleaf[which.max(d$A)], 30, 36), Omega = c(10, 18))
  fits <- lapply(seq_len(nrow(starts)), function(i) tryCatch(
    nlsLM(A ~ Aopt * exp(-((Tleaf - Topt) / Omega)^2), data = d,
          start = list(Aopt = max(d$A), Topt = starts$Topt[i], Omega = starts$Omega[i]),
          lower = c(Aopt = 0, Topt = 5, Omega = 2),
          upper = c(Aopt = 100, Topt = 60, Omega = 80),
          control = nls.lm.control(maxiter = 500)),
    error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) return(out)
  fit <- fits[[which.min(sapply(fits, deviance))]]
  co  <- summary(fit)$coefficients
  out$Topt_A <- co["Topt", 1]; out$Topt_A_se <- co["Topt", 2]
  out$Aopt   <- co["Aopt", 1]; out$Omega     <- co["Omega", 1]
  out$converged <- TRUE
  out
}

topt_A <- amb %>%
  group_by(Species, plant_id, growth_temp, irrigation_trt, plant_date) %>%
  group_modify(~ fit_june(.x)) %>%
  ungroup() %>%
  mutate(topt_in_range = converged & Topt_A >= T_min & Topt_A <= T_max)

topt_models <- list(
  Topt_A = fit_cell_model(filter(topt_A, topt_in_range), "Topt_A"),
  Aopt   = fit_cell_model(filter(topt_A, topt_in_range), "Aopt"),
  Omega  = fit_cell_model(filter(topt_A, topt_in_range), "Omega")
)

# Diagnostics: one DHARMa panel per model
for (nm in names(topt_models)) {
  pdf(file.path("Results/Dharma", paste0("diagnostics_", nm, ".pdf")), width = 9, height = 5)
  plot(simulateResiduals(topt_models[[nm]]), title = nm)
  dev.off()
}


###############################################################################
#   all measurement temperatures (supplementary)
#    plant_id random effect: each plant is measured at six temperatures
###############################################################################

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
all_T_models <- list(
  A_allT   = glmmTMB(A ~ growth_temp * irrigation_trt * Species * temp_setpoint +
                       (1 | plant_id) + (1 | plant_date),
                     dispformula = ~ temp_setpoint + growth_temp, family = gaussian, data = amb,
                     control = ctrl),
  gsw_allT = glmmTMB(log(gsw) ~ growth_temp * irrigation_trt * Species * temp_setpoint +
                       (1 | plant_id) + (1 | plant_date),
                     dispformula = ~ temp_setpoint, family = gaussian, data = amb,
                     control = ctrl)
)

# Diagnostics: one DHARMa panel per model
for (nm in names(all_T_models)) {
  pdf(file.path("Results/Dharma", paste0("diagnostics_", nm, ".pdf")), width = 9, height = 5)
  plot(simulateResiduals(all_T_models[[nm]]), title = nm)
  dev.off()
}

###############################################################################
# output tables
###############################################################################

all_models <- c(growth_models, topt_models, all_T_models)

anova_table <- imap_dfr(all_models, function(m, nm) {
  a <- as.data.frame(Anova(m))
  tibble(model = nm, term = rownames(a), chisq = a$Chisq, df = a$Df,
         p = a$`Pr(>Chisq)`, n = nobs(m))
}) %>%
  mutate(sig = cut(p, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   c("***", "**", "*", ".", "")))
print(anova_table, n = Inf)

###############################################################################
# Figures
###############################################################################

water_cols <- c("High water availability" = "blue", "Low water availability" = "brown")
wrap_lab <- scales::label_wrap(12)

plot_cell <- function(m, ylab) {
  df  <- attr(m, "data")
  emm <- as.data.frame(emmeans(m, ~ growth_temp * irrigation_trt * Species, type = "response"))
  est <- intersect(c("emmean", "response"), names(emm))[1]
  emm$lo <- emm[[est]] - emm$SE
  emm$hi <- emm[[est]] + emm$SE
  pd <- position_dodge(width = 0.6)
  ggplot(emm, aes(growth_temp, .data[[est]], colour = irrigation_trt)) +
    geom_point(data = df, aes(y = .data[[attr(m, "y")]]), alpha = 0.35, size = 2,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.8, position = pd) +
    geom_point(shape = 15, size = 4, position = pd) +
    facet_wrap(~ Species) +
    scale_x_discrete(labels = wrap_lab) +
    scale_colour_manual(values = water_cols, labels = wrap_lab) +
    labs(x = NULL, y = ylab, colour = NULL) +
    theme_bw(base_size = 14) + theme(legend.position = "top")
}

p_Agrowth <-  plot_cell(growth_models$A_growth,
                        expression(A[sat] ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1}))) +
  theme(axis.text.x = element_blank())

p_Agrowth <-
  tag_facet(
  p_Agrowth,
  open = "", close = "",
  tag_pool = letters,
  fontface = 2,
  hjust = -0.5, vjust = 1.5
) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )

p_gsw <- plot_cell(growth_models$gsw_growth,
                   expression(g[sw] ~ (mol ~ m^{-2} ~ s^{-1}))) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank())
p_gsw <-
  tag_facet(
    p_gsw,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


p_cica <- plot_cell(growth_models$CiCa_growth, expression(C[i] / C[a])) +
  scale_y_continuous(limits = c(0.6, 1.0)) +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
p_cica <-
  tag_facet(
    p_cica,
    open = "", close = "",
    tag_pool = letters[5:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )

p_growth <- p_Agrowth / p_gsw / p_cica +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(file.path("Results/gas_exchange_at_growth_T.pdf"), p_growth, width = 7, height = 8)


p_toptA <- plot_cell(topt_models$Topt_A, expression(T[opt] ~ of ~ A ~ (degree*C))) +
  theme(axis.text.x = element_blank())

p_toptA <-
  tag_facet(
    p_toptA,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )
  
p_Aopt <- plot_cell(topt_models$Aopt,   expression(A[opt] ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1}))) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank())
p_Aopt <-
  tag_facet(
    p_Aopt,
    open = "", close = "",
    tag_pool = letters[5:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


p_omega <- plot_cell(topt_models$Omega,  expression(Omega ~ (degree*C))) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) 

p_omega <-
  tag_facet(
    p_omega,
    open = "", close = "",
    tag_pool = letters[7:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


p_topt <- p_toptA / p_Aopt/ p_omega +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

ggsave(file.path("Results/Topt_A.pdf"), p_topt, width = 6, height = 8)


# A-T curves: fitted June curve from the mean parameters of each treatment,
# over plant-level means ± SE at each set point (same plants as the models)
june_curves <- topt_A %>%
  filter(topt_in_range) %>%
  group_by(Species, growth_temp, irrigation_trt) %>%
  summarise(across(c(Topt_A, Aopt, Omega), mean), .groups = "drop") %>%
  tidyr::crossing(Tleaf = seq(14, 48, 0.5)) %>%
  mutate(A = Aopt * exp(-((Tleaf - Topt_A) / Omega)^2))

A_obs <- amb %>%
  filter(plant_id %in% topt_A$plant_id[topt_A$topt_in_range]) %>%
  group_by(Species, growth_temp, irrigation_trt, temp_setpoint) %>%
  summarise(mean = mean(A), se = sd(A) / sqrt(n()), .groups = "drop")
A_obs$temp_setpoint<- as.numeric(as.character(A_obs$temp_setpoint))

p_Acurve <- ggplot(june_curves, aes(Tleaf, A, colour = irrigation_trt, linetype = growth_temp)) +
  geom_line(linewidth = 0.9) +
  geom_errorbar(data = A_obs, aes(x = temp_setpoint, y = mean, ymin = mean - se, ymax = mean + se),
                width = 0.6, linetype = "solid") +
  geom_point(data = A_obs, aes(x = temp_setpoint, y = mean, shape = growth_temp), size = 2.5) +
  facet_wrap(~ Species) +
  scale_colour_manual(values = water_cols, labels = wrap_lab) +
  scale_linetype_manual(values = c("dashed", "solid"), labels = wrap_lab) +
  scale_shape_manual(values = c(17, 16), labels = wrap_lab) +
  labs(x = expression(Leaf ~ temperature ~ (degree*C)),
       y = expression(A[sat] ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1})),
       colour = NULL, linetype = NULL, shape = NULL,
       tag = "a") +
  theme_bw(base_size = 14) +
  scale_x_continuous(
    breaks = c(16, 22, 28, 34, 40, 46),
    labels = c(16, 22, 28, 34, 40, 46)
  ) +
  theme(legend.position = "top", legend.key.width = unit(1.2, "cm"))

p_Acurve <-
  tag_facet(
    p_Acurve,
    open = "", close = "",
    tag_pool = letters,
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )

ggsave(file.path("Results/A_temperature_response.pdf"), p_Acurve, width = 11, height = 5)


p_Aresponse <- p_Acurve/p_topt +
  plot_layout(heights = c(1, 2.5))
  
ggsave("~/Desktop/Photosynthesis_paper/Results/tem_response_A.pdf", p_Aresponse,
       height = 12,
       width = 8,
       dpi = 300)



# A and gsw vs leaf temperature: plant-level means ± SE at each set point
curve_means <- function(y, ylab) {
  amb %>%
    group_by(Species, growth_temp, irrigation_trt, temp_setpoint) %>%
    summarise(mean = mean(.data[[y]]), se = sd(.data[[y]]) / sqrt(n()), .groups = "drop") %>%
    ggplot(aes(as.numeric(as.character(temp_setpoint)), mean, colour = irrigation_trt, linetype = growth_temp,
               shape = growth_temp)) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linetype = "solid") +
    geom_point(size = 2.5) +
    facet_wrap(~ Species) +
    scale_colour_manual(values = water_cols, labels = wrap_lab) +
    scale_linetype_manual(values = c("dashed", "solid"), labels = wrap_lab) +
    scale_shape_manual(values = c(17, 16), labels = wrap_lab) +
    labs(x = expression(Leaf ~ temperature ~ (degree*C)), y = ylab,
         colour = NULL, linetype = NULL, shape = NULL) +
    theme_bw(base_size = 14)
}

p_allT <- curve_means("A", expression(A[sat] ~ (mu*mol ~ CO[2] ~ m^{-2} ~ s^{-1}))) /
  curve_means("gsw", expression(g[sw] ~ (mol ~ m^{-2} ~ s^{-1}))) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", legend.key.width = unit(1.2, "cm"))

ggsave(file.path("results/A_gsw_vs_Tleaf.pdf"), p_allT, width = 11, height = 8)


