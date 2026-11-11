
###############################################################################
# Biomass, root:shoot
###############################################################################

library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(lme4)
library(car)
library(emmeans)
library(DHARMa)

biomass_data <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/biomass_data.xlsx") 
metadata <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/metadata.xlsx")
biomass_data <- biomass_data|> 
  left_join(metadata, by = "plant_id") |> 
  mutate(
    # factor by NAME, not column position; clean labels (line breaks only in plots)
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
    plant_date = factor(plant_date))      # CHECK: column name of planting date


fit_biomass <- function(y, log_y = FALSE) {
  lhs <- if (log_y) paste0("log(", y, ")") else y
  m <- lmer(as.formula(paste(lhs, "~ growth_temp * irrigation_trt * Species + (1 | plant_date)")),
            data = biomass_data)
  attr(m, "data") <- biomass_data
  attr(m, "y")    <- y
  m
}

models <- list(
  aboveground = fit_biomass("aboveground_biomass"),
  belowground = fit_biomass("root_biomass"),
  root_shoot  = fit_biomass("root_shoot_ratio")
)

# Diagnostics 
for (nm in names(models)) {
  plot(simulateResiduals(models[[nm]]), title = nm)
}

# Wald chi-square tests, one table
anova_biomass <- bind_rows(lapply(names(models), function(nm) {
  a <- as.data.frame(Anova(models[[nm]]))
  data.frame(model = nm, term = rownames(a), chisq = a$Chisq, df = a$Df,
             p = a$`Pr(>Chisq)`, n = nobs(models[[nm]]))
}))
print(anova_biomass, digits = 3)

write.csv(anova_biomass, "~/Desktop/Photosynthesis_paper/Results/anova_biomass.csv", row.names = FALSE)


# Contrasts with % change 

contrast_pct <- function(m, factor, by = NULL) {
  f   <- if (is.null(by)) as.formula(paste("~", factor)) else
    as.formula(paste("~", factor, "|", paste(by, collapse = " * ")))
  emm <- emmeans(m, f, type = "response")
  as.data.frame(pairs(emm, reverse = TRUE)) %>%
    mutate(pct = if ("ratio" %in% names(.)) 100 * (ratio - 1) else NA_real_)
}


# Figures


water_cols <- c("High water availability" = "blue", "Low water availability" = "brown")
wrap_lab   <- scales::label_wrap(12)

plot_biomass <- function(m, ylab) {
  df  <- attr(m, "data")
  emm <- as.data.frame(emmeans(m, ~ growth_temp * irrigation_trt * Species, type = "response"))
  est <- intersect(c("response", "emmean"), names(emm))[1]
  pd  <- position_dodge(width = 0.6)
  ggplot(emm, aes(growth_temp, .data[[est]], colour = irrigation_trt)) +
    geom_point(data = df, aes(y = .data[[attr(m, "y")]]), alpha = 0.35, size = 2,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
    geom_errorbar(aes(ymin = .data[[est]] - SE, ymax = .data[[est]] + SE),
                  width = 0.15, linewidth = 0.8, position = pd) +
    geom_point(shape = 15, size = 4, position = pd) +
    facet_wrap(~ Species) +
    scale_x_discrete(labels = wrap_lab) +
    scale_colour_manual(values = water_cols, labels = wrap_lab) +
    labs(x = NULL, y = ylab, colour = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
}

p_abg <- plot_biomass(models$aboveground, "Aboveground biomass (g)") +
  theme(axis.text.x = element_blank())

p_bgb <- plot_biomass(models$belowground, "Belowground biomass (g)") +
  theme(axis.text.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())

p_rsratio <- plot_biomass(models$root_shoot, "Root : shoot ratio") +
  theme(strip.text = element_blank(),
        strip.background = element_blank())

p_biomass <-
  p_abg /
  p_bgb /
  p_rsratio +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top")

p_biomass

ggsave("~/Desktop/Photosynthesis_paper/Results/plant_biomass.pdf", p_biomass, 
       height = 9,
       width = 6,
       dpi = 300)
