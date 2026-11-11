###############################################################################
# Thermal acclimation of Vcmax and Jmax in Pima and Upland cotton
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(minpack.lm)
library(lme4)
library(car)
library(emmeans)
library(stringr)
library(DHARMa)

###############################################################################
# set up data and variables

setwd("~/Desktop/Photosynthesis_paper/")

R_gas <- 8.314     # J mol-1 K-1
Hd    <- 200000    # J mol-1, fixed (Medlyn et al. 2002)
T_ref <- 298.15    # K, reference temperature for Vcmax25 / Jmax25
setpoints <- c(16, 22, 28, 34, 40, 46)

# Leaf-temperature set points (°C). Curves are FITTED on the measured mean leaf
# temperature; set points are used only to group points when
# plotting observed means. 


###############################################################################
# Data
###############################################################################

pima_parameters <- read.csv("Final_data/aci_parameters_pima.csv")
pima_parameters <- pima_parameters |> 
  mutate(plant_id = str_replace(id, "(^\\d+)_.*", "P\\1")) |> 
  dplyr::select(plant_id, Vcmax, Jmax, Tleaf_mean, Tleaf_set)

meta_pima <- read.csv("Final_data/metadata_pima.csv")

pima_parameters <- meta_pima |> 
  left_join(pima_parameters, by = "plant_id")


upland_parameters <- read.csv("Final_data/aci_parameters_upland.csv")
upland_parameters <- upland_parameters |> 
  mutate(plant_id = str_replace(id, "(^\\d+)_.*", "U\\1")) |> 
  dplyr::select(plant_id, Vcmax, Jmax, Tleaf_mean, Tleaf_set)

meta_upland <- read.csv("Final_data/metadata_upland.csv")
upland_parameters <- meta_upland |> 
  left_join(upland_parameters, by = "plant_id")

aci_parameter <- bind_rows(pima_parameters, upland_parameters) |> 
  mutate(Tleaf_K   = Tleaf_mean + 273.15)

###############################################################################
# fit the peaked Arrhenius function (Medlyn, 2002) to each curve
###############################################################################

peaked_medlyn <- function(Tk, kopt, Ha, Topt) {
  kopt * (Hd * exp(Ha * (Tk - Topt) / (Tk * R_gas * Topt))) /
    (Hd - Ha * (1 - exp(Hd * (Tk - Topt) / (Tk * R_gas * Topt))))
}

# Fits from several starting values; keeps the converged fit with lowest RSS
fit_curve <- function(df, response) {
  df <- df[is.finite(df[[response]]), ]
  if (nrow(df) < 4) return(NULL)
  f <- as.formula(paste(response, "~ peaked_medlyn(Tleaf_K, kopt, Ha, Topt)"))
  starts <- expand.grid(Ha = c(115000, 60000), Topt = c(305, 310, 315))
  fits <- map(seq_len(nrow(starts)), function(i) tryCatch(
    nlsLM(f, data = df,
          start = list(kopt = max(df[[response]]),
                       Ha = starts$Ha[i], Topt = starts$Topt[i]),
          lower = c(kopt = 0,    Ha = 1e4, Topt = 290),
          upper = c(kopt = 1000, Ha = 3e5, Topt = 330),
          control = nls.lm.control(maxiter = 1000)),
    error = function(e) NULL))
  fits <- compact(fits)
  if (length(fits) == 0) return(NULL)
  fits[[which.min(map_dbl(fits, deviance))]]
}

extract_fit <- function(fit, df) {
  if (is.null(fit)) return(tibble(converged = FALSE))
  co <- summary(fit)$coefficients
  tibble(converged = TRUE,
         kopt = co["kopt", 1], kopt_se = co["kopt", 2],
         Ha   = co["Ha", 1],   Ha_se   = co["Ha", 2],
         Topt = co["Topt", 1], Topt_se = co["Topt", 2],
         k25  = peaked_medlyn(T_ref, co["kopt", 1], co["Ha", 1], co["Topt", 1]),
         n_points = nrow(df),
         Tleaf_min = min(df$Tleaf_K), Tleaf_max = max(df$Tleaf_K))
}

fit_all <- function(response) {
  aci_parameter %>%
    group_by(plant_id, Species, growth_temp, irrigation_trt,
             treatment, plant_date) %>%
    nest() %>%
    ungroup() %>%
    mutate(fit  = map(data, fit_curve, response = response),
           pars = map2(fit, data, extract_fit)) %>%
    unnest(pars) %>%
    mutate(response = !!response,
           # peak not observed if Topt lies outside the measured range
           topt_in_range = converged & Topt >= Tleaf_min & Topt <= Tleaf_max)
}

fits_v <- fit_all("Vcmax")
fits_j <- fit_all("Jmax")

###############################################################################
# parameter table for modelling
###############################################################################

label_factors <- function(df) {
  df %>% mutate(
    growth_temp = factor(dplyr::recode(growth_temp,
                                       "Ambient" = "Low growth temperature",
                                       "Heat"    = "High growth temperature"),
                         levels = c("Low growth temperature",
                                    "High growth temperature")),
    irrigation_trt = factor(dplyr::recode(irrigation_trt,
                                          "Irrigated" = "High water availability",
                                          "Drought"   = "Low water availability"),
                            levels = c("High water availability",
                                       "Low water availability")),
    Species = factor(dplyr::recode(Species,
                                   "Pima" = "Pima Cotton",
                                   "Upland" = "Upland Cotton")),
    plant_date = factor(plant_date))
}


pars_v <- label_factors(fits_v)
pars_j <- label_factors(fits_j)

# Jmax:Vcmax at 25 °C, for plants with both fits retained
ratio25 <- inner_join(
  pars_v %>% dplyr::select(plant_id, growth_temp, irrigation_trt, Species,
                           plant_date, V25 = k25),
  pars_j %>% dplyr::select(plant_id, J25 = k25),
  by = "plant_id") %>%
  mutate(JV25 = J25 / V25)

###############################################################################
# mixed models
###############################################################################

# Fits y ~ growth_temp * irrigation_trt * Species + (1|plant_date),

fit_param_model <- function(df, y, log_y = FALSE) {
  df$w <- 1
  lhs <- if (log_y) paste0("log(", y, ")") else y
  f <- as.formula(paste(lhs, "~ growth_temp * irrigation_trt * Species + (1 | plant_date)"))
  m <- lmer(f, data = df)
  attr(m, "data") <- df
  attr(m, "y")    <- y
  m
}

models <- list(
  Vcmax_kopt = fit_param_model(pars_v, "kopt"),
  Vcmax_Topt = fit_param_model(pars_v, "Topt"),
  Vcmax_Ha   = fit_param_model(pars_v, "Ha"),
  Vcmax_25   = fit_param_model(pars_v, "k25"),
  Jmax_kopt  = fit_param_model(pars_j, "kopt"),
  Jmax_Topt  = fit_param_model(pars_j, "Topt"),
  Jmax_Ha    = fit_param_model(pars_j, "Ha"),
  Jmax_25    = fit_param_model(pars_j, "k25"),
  JV_ratio25 = fit_param_model(ratio25, "JV25")
)

# Diagnostics: one DHARMa panel per model
for (nm in names(models)) {
  pdf(file.path("Results/Dharma", paste0("diagnostics_", nm, ".pdf")), width = 9, height = 5)
  plot(simulateResiduals(models[[nm]]), title = nm)
  dev.off()
}

# all model passed the DHARMa diagnostic tests

# Wald chi-square tests (Type II), one table for the paper
anova_table <- imap_dfr(models, function(m, nm) {
  a <- as.data.frame(Anova(m))
  tibble(model = nm, term = rownames(a), chisq = a$Chisq, df = a$Df,
         p = a$`Pr(>Chisq)`, n = nrow(attr(m, "data")))
}) %>%
  mutate(sig = cut(p, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                   c("***", "**", "*", ".", "")))
print(anova_table, n = Inf)

write.csv(anova_table, file.path("Results/", "anova_parameters.csv"), row.names = FALSE)

# Post-hoc contrasts: growth-temperature effect within each water x species
# cell (tests H1, H2 and H3 directly)
contrasts_gt <- imap_dfr(models, function(m, nm) {
  as.data.frame(pairs(emmeans(m, ~ growth_temp | irrigation_trt * Species,
                              type = "response"))) %>%
    mutate(model = nm, .before = 1)
})

write.csv(contrasts_gt, file.path("Results/", "contrasts_growth_temp.csv"),
          row.names = FALSE)

# Water effect within each growth temperature x species cell
contrasts_water <- imap_dfr(models, function(m, nm) {
  as.data.frame(pairs(emmeans(m, ~ irrigation_trt | growth_temp * Species,
                              type = "response"))) %>%
    mutate(model = nm, .before = 1)
})

write.csv(contrasts_water, file.path("Results/", "contrasts_water.csv"),
          row.names = FALSE)

###############################################################################
# Figures
###############################################################################

water_cols <- c("High water availability" = "blue",
                "Low water availability"  = "brown")
wrap_lab <- scales::label_wrap(12)   # two-line axis / legend labels

#  Parameter panels: raw points + model means ± SE, faceted by species
plot_param <- function(m, ylab, transform = identity) {
  df  <- attr(m, "data")
  emm <- as.data.frame(emmeans(m, ~ growth_temp * irrigation_trt * Species,
                               type = "response"))
  est <- intersect(c("emmean", "response"), names(emm))[1]
  df$y_plot <- transform(df[[attr(m, "y")]])
  # mean ± SE: add/subtract SE on the model scale FIRST, then convert units
  # (K -> °C, J -> kJ), so the SE itself is never shifted by 273.15
  emm$est <- transform(emm[[est]])
  emm$lo  <- transform(emm[[est]] - emm$SE)
  emm$hi  <- transform(emm[[est]] + emm$SE)
  pd <- position_dodge(width = 0.6)
  ggplot(emm, aes(growth_temp, est, colour = irrigation_trt)) +
    geom_point(data = df, aes(y = y_plot), alpha = 0.35, size = 2,
               position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.6)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.8,
                  position = pd) +
    geom_point(shape = 15, size = 4, position = pd) +
    facet_wrap(~ Species) +
    scale_x_discrete(labels = wrap_lab) +
    scale_colour_manual(values = water_cols, labels = wrap_lab) +
    labs(x = NULL, y = ylab, colour = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
}

to_C  <- function(x) x - 273.15
to_kJ <- function(x) x / 1000

lab_v25  <- expression(V[cmax25] ~ (mu*mol ~ m^{-2} ~ s^{-1}))
lab_j25  <- expression(J[max25] ~ (mu*mol ~ m^{-2} ~ s^{-1}))
lab_vopt <- expression(V[cmax] ~ "at" ~ T[opt] ~ (mu*mol ~ m^{-2} ~ s^{-1}))
lab_jopt <- expression(J[max] ~ "at" ~ T[opt] ~ (mu*mol ~ m^{-2} ~ s^{-1}))



p_params <-
  (plot_param(models$Vcmax_Topt, expression(T[opt] ~ of ~ V[cmax] ~ (degree*C)), to_C) +
     plot_param(models$Jmax_Topt,  expression(T[opt] ~ of ~ J[max] ~ (degree*C)),  to_C)) /
  (plot_param(models$Vcmax_25, lab_v25) +
     plot_param(models$Jmax_25,  lab_j25)) /
  (plot_param(models$Vcmax_Ha, expression(H[a] ~ of ~ V[cmax] ~ (kJ ~ mol^{-1})), to_kJ) +
     plot_param(models$Jmax_Ha,  expression(H[a] ~ of ~ J[max] ~ (kJ ~ mol^{-1})),  to_kJ)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top")

ggsave(file.path("Results/", "acclimation_parameters.pdf"), p_params,
       width = 11, height = 12)

p_supp <- (plot_param(models$Vcmax_kopt, lab_vopt) +
             plot_param(models$Jmax_kopt,  lab_jopt)) /
  plot_param(models$JV_ratio25, expression(J[max] : V[cmax] ~ "at 25" ~ degree*C)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top")

ggsave(file.path("Results/", "supp_kopt_ratio.pdf"), p_supp, width = 11, height = 8)

# Temperature response curves: observed means ± SE at each set point
#     plus the curve drawn from the mean fitted parameters of each treatment.
#     Uses exactly the same plants as the models.
curve_panel <- function(pars, response, ylab) {
  ids <- pars$plant_id
  
  obs <- aci_parameter %>%
    filter(plant_id %in% ids) %>%
    label_factors() %>%
    group_by(Species, growth_temp, irrigation_trt, Tleaf_set) %>%
    summarise(mean = mean(.data[[response]], na.rm = TRUE),
              se   = sd(.data[[response]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response]]))),
              .groups = "drop")
  obs$Tleaf_set <- as.numeric(as.character(obs$Tleaf_set))
  
  curves <- pars %>%
    group_by(Species, growth_temp, irrigation_trt) %>%
    summarise(across(c(kopt, Ha, Topt), mean), .groups = "drop") %>%
    crossing(Tk = seq(285.15, 320.15, 0.5)) %>%
    mutate(value = peaked_medlyn(Tk, kopt, Ha, Topt), Tleaf = Tk - 273.15)
  
  ggplot(curves, aes(Tleaf, value, colour = irrigation_trt, linetype = growth_temp)) +
    geom_line(linewidth = 0.9) +
    geom_errorbar(data = obs, aes(x = Tleaf_set, y = mean,
                                  ymin = mean - se, ymax = mean + se),
                  width = 0.6, linetype = "solid") +
    geom_point(data = obs, aes(x = Tleaf_set, y = mean, shape = growth_temp),
               size = 2.5) +
    facet_wrap(~ Species) +
    scale_colour_manual(values = water_cols, labels = wrap_lab) +
    scale_linetype_manual(values = c("dashed", "solid"), labels = wrap_lab) +
    scale_shape_manual(values = c(17, 16), labels = wrap_lab) +
    scale_x_continuous(
      breaks = c(16, 22, 28, 34, 40, 46),
      labels = c(16, 22, 28, 34, 40, 46)) +
    labs(x = expression(Leaf ~ temperature ~ (degree*C)), y = ylab,
         colour = NULL, linetype = NULL, shape = NULL) +
    theme_bw(base_size = 14)
}

p_curves <- curve_panel(pars_v, "Vcmax", expression(V[cmax] ~ (mu*mol ~ m^{-2} ~ s^{-1}))) /
  curve_panel(pars_j, "Jmax",  expression(J[max] ~ (mu*mol ~ m^{-2} ~ s^{-1}))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", legend.key.width = unit(1.2, "cm"))

ggsave(file.path("Results/", "temperature_response_curves.pdf"), p_curves,
       width = 11, height = 9)

###############################################################################
# Final figures
###############################################################################
Tresponse_Vcmax <- curve_panel(pars_v, "Vcmax", expression(V[cmax] ~ (mu*mol ~ m^{-2} ~ s^{-1}))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

Tresponse_Vcmax <-
  tag_facet(
    Tresponse_Vcmax,
    open = "", close = "",
    tag_pool = letters,
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect() ,
    legend.position = "top"
  )

Tresponse_Jamx <-  curve_panel(pars_j, "Jmax",  expression(J[max] ~ (mu*mol ~ m^{-2} ~ s^{-1}))) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

Tresponse_Jamx <-
  tag_facet(
    Tresponse_Jamx,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )

t_response <- Tresponse_Vcmax/Tresponse_Jamx
  

p_Topt_Vcmax <- plot_param(models$Vcmax_Topt, expression(T[opt] ~ of ~ V[cmax] ~ (degree*C)), to_C) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

p_Topt_Vcmax <-
  tag_facet(
    p_Topt_Vcmax,
    open = "", close = "",
    tag_pool = letters[5:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )
  
p_Topt_Jmax <- plot_param(models$Jmax_Topt,  expression(T[opt] ~ of ~ J[max] ~ (degree*C)),  to_C) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")
p_Topt_Jmax <-
  tag_facet(
    p_Topt_Jmax,
    open = "", close = "",
    tag_pool = letters[7:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


Topt_plot <- p_Topt_Vcmax/p_Topt_Jmax

t_response_final <- (t_response / Topt_plot) +
  plot_layout(heights = c(1, 1, 2))

ggsave("Results/Tresponse.pdf", t_response_final, height = 10, width = 8, dpi = 300)

# Supplemental fig (Ha and Kopt)

Vcamx_Ha <- plot_param(models$Vcmax_Ha, expression(H[a] ~ of ~ V[cmax] ~ (kJ ~ mol^{-1})), to_kJ) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")
Vcamx_Ha <-
  tag_facet(
    Vcamx_Ha,
    open = "", close = "",
    tag_pool = letters,
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )             


Jmax_Ha <- plot_param(models$Jmax_Ha,  expression(H[a] ~ of ~ J[max] ~ (kJ ~ mol^{-1})),  to_kJ) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")  

Jmax_Ha <-
  tag_facet(
    Jmax_Ha,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )

Kopt_Vcmax <- plot_param(models$Vcmax_kopt, lab_vopt) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")  

Kopt_Vcmax <-
  tag_facet(
    Kopt_Vcmax,
    open = "", close = "",
    tag_pool = letters[5:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )

Kopt_Jmax <-plot_param(models$Jmax_kopt,  lab_jopt) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") 
Kopt_Jmax <-
  tag_facet(
    Kopt_Jmax,
    open = "", close = "",
    tag_pool = letters[7:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) 

tresponse_para_supp <- (Vcamx_Ha / Jmax_Ha / Kopt_Vcmax/ Kopt_Jmax)

ggsave("Results/Tresponse_para_supp.pdf", tresponse_para_supp, height = 12, width = 8, dpi = 300)


# photosynthetic capacites

p_Vcmax25 <- plot_param(models$Vcmax_25, lab_v25) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top")
p_Vcmax25 <-
  tag_facet(
    p_Vcmax25,
    open = "", close = "",
    tag_pool = letters,
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )
  

p_Jmax25 <- plot_param(models$Jmax_25,  lab_j25) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

p_Jmax25 <-
  tag_facet(
    p_Jmax25,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


p_ratio <- plot_param(models$JV_ratio25, expression(J[max] : V[cmax] ~ "at 25" ~ degree*C)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")
p_ratio <-
  tag_facet(
    p_ratio,
    open = "", close = "",
    tag_pool = letters[5:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


rubisco_conc_graph

p_photo_capacity <- (p_Vcmax25 / p_Jmax25 / p_ratio / rubisco_conc_graph)

ggsave("Results/photo_capacities.pdf", p_photo_capacity, height = 12, width = 8, dpi = 300)





