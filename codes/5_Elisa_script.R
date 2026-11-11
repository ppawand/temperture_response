library(drc)
library(readxl)
library(ggplot2)
library(car)
library(emmeans)
library(lme4)
library(multcomp)

sc_data <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/Elisa_final.xlsx",
                   sheet = "standard_curve")
# fit 4PL model

model_4pl <- drm(OD ~ conc, data = sc_data, fct = LL.4())
summary(model_4pl)

# plot
plot(model_4pl, log = "x", type = "all",
     main = "4PL Standard Curve", 
     xlab = "Concentration", ylab = "Optical Density")

# prediction
elisa_data <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/Elisa_final.xlsx",
                         sheet = "average_OD")

predicted_conc <- as.data.frame(ED(model_4pl, respLev = elisa_data$OD,
                     type = "absolute"))

elisa_data$conc <- predicted_conc$Estimate

metadata <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/metadata.xlsx")
metadata <- metadata |> 
  dplyr::select(plant_id, plant_date)

elisa_data <- left_join(elisa_data, metadata, by = "plant_id")

elisa_data <- elisa_data %>%
  mutate(
    growth_temp = relevel(factor(growth_temp), ref = "Ambient"),      
    water_trt   = relevel(factor(water_trt),   ref = "Irrigated"),
    species     = factor(species),
    plant_date = factor(plant_date)
  )

# 
ggplot(elisa_data, aes(trt, conc)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(~species, scale = "free")

# 
conc_model <- lmer(log(conc) ~ growth_temp * water_trt * species + (1|plant_date),
       data = elisa_data)

plot(simulateResiduals(conc_model))

Anova(conc_model, type = 2)
summary(conc_model)

rubisco_conc <- as.data.frame(emmeans(conc_model, ~ growth_temp:water_trt:species, type = "response"))

emmeans(conc_model, ~ growth_temp:water_trt) # only significant interaction


# plot three way interaction means


(rubisco_conc_graph <- 
    ggplot(rubisco_conc, aes(growth_temp, response, color = water_trt)) +
    geom_point(data = elisa_data,
               aes(x = growth_temp, y = conc), alpha = 0.35, size = 2,
               position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.6)) +
    geom_point(shape = 15, size = 4, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                  position = position_dodge(width = 0.6), width = 0.15, linewidth = 0.8) +
    scale_color_manual(values = c("Drought" = "brown", "Irrigated" = "blue"),
                      labels = c("Drought"   = "Low water\navailability",
                                 "Irrigated" = "High water\navailability"),
                      breaks = c("Irrigated", "Drought"))  +
    facet_wrap(~ species,
               labeller = labeller(species = c("Pima"   = "Pima Cotton",
                                               "Upland" = "Upland Cotton"))) +
    labs(x = "Growth Temperatures",
         y =  "Rubisco concentration (µg/mL)",
         color = NULL) +
    scale_y_continuous(limits = c(0, 420)) +
    scale_x_discrete(labels = c("Low growth\ntemperature",
                                "High growth\ntemperature")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_blank()))

rubisco_conc_graph <-
  tag_facet(
    rubisco_conc_graph,
    open = "", close = "",
    tag_pool = letters[7:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )


# plot growth temp water interaction effects only
ggplot(as.data.frame(emmeans(conc_model, ~ growth_temp:water_trt, type = "response")), aes(x = growth_temp, y = response, fill = water_trt)) +
  geom_col(position = position_dodge(0.9), color = "black", width = 0.75, alpha = 0.85) +
  geom_errorbar(aes(ymin = response,                
                    ymax = response + SE),        
                position = position_dodge(0.9), width = 0.18, linewidth = 0.5) +
  scale_fill_manual(values = c("Drought" = "#C0392B", "Irrigated" = "#2C6FBB"),
                    labels = c("LW", "HW")) +
  labs(x = "Growth Temperatures",
       y =  "RuBisco concentration (µg/mL)",
       fill = NULL) +
  scale_x_discrete(labels = c("LT",
                              "HT")) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
        panel.grid = element_blank())


# Vcmax Topt (K difference = °C difference)
print(contrast_pct(models$Vcmax_Topt, "growth_temp", "irrigation_trt"), width = Inf)

# main effects of growth temperature (% change and ratio)
for (nm in c("Vcmax_25", "Jmax_25", "Vcmax_kopt", "Jmax_kopt", "JV_ratio25")) {
  cat("\n", nm, "\n")
  print(emmeans(models[[nm]], ~ growth_temp))
}

# Rubisco: both directions of the interaction
print(contrast_pct(conc_model, "growth_temp", "water_trt"), width = Inf)
print(contrast_pct(conc_model, "water_trt", "growth_temp"), width = Inf)

