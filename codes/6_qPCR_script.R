library(readxl)
library(dplyr)
library(ggplot2)
library(car)
library(emmeans)
library(lme4)
library(patchwork)

qpcr_data <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/qPCR_data.xlsx",
                        sheet = 1) |> 
  filter(!is.na(Ct_alpha), !is.na(Ct_beta), !is.na(Ct_UBQ7))

metadata <- read_excel("~/Desktop/Photosynthesis_paper/Final_data/metadata.xlsx")
metadata <- metadata |> 
  dplyr::select(plant_id, plant_date)

qpcr_data <- left_join(qpcr_data, metadata, by = "plant_id")


# calculate delta-Ct = target - reference gene (UBQ7). 
#    We model -delta-Ct: HIGHER = MORE expression (easier to read).

qpcr_data <- qpcr_data %>%
  mutate(
    neg_dCt_alpha = -(Ct_alpha - Ct_UBQ7),          # log2 relative expression
    neg_dCt_beta  = -(Ct_beta  - Ct_UBQ7),
    log2_beta_alpha = (Ct_alpha - Ct_UBQ7) - (Ct_beta - Ct_UBQ7),  # UBQ7 cancels
    growth_temp = relevel(factor(growth_temp), ref = "Ambient"),
    water_trt   = relevel(factor(water_trt),   ref = "Irrigated"),
    Species     = factor(Species),
    plant_date  = factor(plant_date)
  )


# Model on the delta-Ct scale

fit_qpcr <- function(y) {
  df <- qpcr_data[is.finite(qpcr_data[[y]]), ]
  m  <- lmer(as.formula(paste(y, "~ growth_temp * water_trt * Species + (1 | plant_date)")),
             data = df)
  attr(m, "data") <- df
  m
}

qpcr_models <- list(
  RCA_alpha  = fit_qpcr("neg_dCt_alpha"),
  RCA_beta   = fit_qpcr("neg_dCt_beta"),
  beta_alpha = fit_qpcr("log2_beta_alpha")
)

for (nm in names(qpcr_models)) {
  plot(simulateResiduals(qpcr_models[[nm]]), title = nm)
}

anova_qpcr <- bind_rows(lapply(names(qpcr_models), function(nm) {
  a <- as.data.frame(Anova(qpcr_models[[nm]], type = 2))
  data.frame(model = nm, term = rownames(a), chisq = a$Chisq, df = a$Df,
             p = a$`Pr(>Chisq)`, n = nobs(qpcr_models[[nm]]))
}))
print(anova_qpcr, digits = 3)

write.csv(anova_qpcr, "~/Desktop/Photosynthesis_paper/Results/anova_qpcr.csv", 
          row.names = FALSE)

d <- fold_vs_control(qpcr_models$RCA_alpha); print(d)
sapply(qpcr_data[c("growth_temp", "water_trt", "Species")], levels)

# Fold change vs. control (Ambient + Irrigated) of the same species

fold_vs_control <- function(m) {
  emm <- emmeans(m, ~ growth_temp * water_trt | Species)        # control = Ambient Irrigated
  ctr <- as.data.frame(confint(contrast(emm, method = "trt.vs.ctrl1")))
  trt <- ctr |>
    tidyr::separate(contrast, into = c("growth_temp", "water_trt"),
                    sep = " ", extra = "drop") |>                 # "Heat Irrigated - Ambient Irrigated"
    dplyr::transmute(growth_temp, water_trt, Species,
                     log2fc = estimate, se = SE, lo = lower.CL, hi = upper.CL)
  ctrl <- expand.grid(growth_temp = "Ambient", water_trt = "Irrigated",
                      Species = levels(qpcr_data$Species), stringsAsFactors = FALSE) |>
    dplyr::mutate(log2fc = 0, se = NA_real_, lo = NA_real_, hi = NA_real_)
  dplyr::bind_rows(trt, ctrl) |>
    dplyr::mutate(growth_temp = factor(growth_temp, levels = c("Ambient", "Heat")),
                  water_trt   = factor(water_trt,   levels = c("Irrigated", "Drought")),
                  Species     = factor(Species),
                  fold        = 2^log2fc)
}


print(fold_vs_control(qpcr_models$RCA_alpha))

plot_fold <- function(m, ylab) {
  d  <- fold_vs_control(m)
  pd <- position_dodge(0.9)
  ggplot(d, aes(growth_temp, log2fc, fill = water_trt)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_col(position = pd, colour = "black", alpha = 0.85, width = 0.75) +
    geom_errorbar(aes(ymin = log2fc - se, ymax = log2fc + se),
                  position = pd, width = 0.18, linewidth = 0.5, na.rm = TRUE) +
    # log2 axis labelled as fold change: 2^-2 = 0.25, not "-4"
    scale_y_continuous(
      breaks = function(l) seq(2 * floor(l[1] / 2), 2 * ceiling(l[2] / 2), 2),
      labels = function(b) ifelse(b >= 0, as.character(2^b), paste0("1/", 2^-b))
    )+
    scale_fill_manual(values = c("Irrigated" = "blue", "Drought" = "brown"),
                      labels = c("Irrigated" = "High water\navailability",
                                 "Drought"   = "Low water\navailability")) +
    scale_x_discrete(labels = c("Ambient" = "Low growth\ntemperature",
                                "Heat"    = "High growth\ntemperature")) +
    facet_wrap(~ Species, labeller = labeller(Species = c("Pima"   = "Pima Cotton",
                                                          "Upland" = "Upland Cotton"))) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
}

p_alpha <- plot_fold(qpcr_models$RCA_alpha,  expression("RCA-" * alpha * " (fold change)")) +
  theme(axis.text.x = element_blank(),
        legend.position = "top")

p_alpha <-
  tag_facet(
    p_alpha,
    open = "", close = "",
    tag_pool = letters,
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 6, b = 6, l = 12, r = 12)),             
    strip.background = element_rect()   
  )

p_beta <- plot_fold(qpcr_models$RCA_beta,   expression("RCA-" * beta  * " (fold change)")) +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

p_beta <-
  tag_facet(
    p_beta,
    open = "", close = "",
    tag_pool = letters[3:26],
    fontface = 2,
    hjust = -0.5, vjust = 1.5
  )

p_beta_alpha <-   plot_fold(qpcr_models$beta_alpha, expression("RCA-" * beta * ":" * alpha * " (fold change)")) 

p_rca <- p_alpha / p_beta

p_rca

ggsave("~/Desktop/Photosynthesis_paper/Results/RCA_plot.pdf", p_rca, height = 7,
       width = 6, dpi = 300)


# log2 difference -> signed fold change
signed_fold <- function(log2fc) ifelse(log2fc >= 0, 2^log2fc, -(2^-log2fc))

## Same numbers as the figure bars: each treatment vs Ambient + Irrigated, within species
fold_table <- function(m) {
  emm <- emmeans(m, ~ growth_temp * water_trt | Species)
  as.data.frame(contrast(emm, method = "trt.vs.ctrl1")) |>
    mutate(fold = round(signed_fold(estimate), 2)) |>
    dplyr::select(Species, contrast, fold, p.value)
}

fold_table(qpcr_models$RCA_alpha)
fold_table(qpcr_models$RCA_beta)

## Drought vs Irrigated within each species (averaged over growth temperature)
as.data.frame(pairs(emmeans(qpcr_models$RCA_alpha, ~ water_trt | Species), reverse = TRUE)) |>
  mutate(fold = round(signed_fold(estimate), 2)) |> dplyr::select(Species, contrast, fold, p.value)

as.data.frame(pairs(emmeans(qpcr_models$RCA_beta, ~ water_trt | Species), reverse = TRUE)) |>
  mutate(fold = round(signed_fold(estimate), 2)) |> dplyr::select(Species, contrast, fold, p.value)

## Heat vs Ambient for RCA-alpha (averaged over water and species)
as.data.frame(pairs(emmeans(qpcr_models$RCA_alpha, ~ growth_temp), reverse = TRUE)) |>
  mutate(fold = round(signed_fold(estimate), 2)) |> dplyr::select(contrast, fold, p.value)

