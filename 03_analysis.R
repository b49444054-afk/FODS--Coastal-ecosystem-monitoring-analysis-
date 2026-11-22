# =============================================================================
# 03_analysis.R
# Coastal Ecosystem Monitoring Analysis - Step 3: relationships & conclusions
# -----------------------------------------------------------------------------
# Questions we ask of the processed data
#   Q1. How are the variables related to each other?            (correlation)
#   Q2. Do sites and seasons differ in water quality?            (group tests)
#   Q3. What drives dissolved oxygen and species richness?       (regression)
#   Q4. Can we predict "Good" ecosystem health?                  (logistic reg.)
#   Q5. Is anything changing over the years?                     (trend checks)
#
# Run after 02_clean_data.R:  source("03_analysis.R")
# Output: figures in outputs/figures, tables in outputs/tables
# =============================================================================

library(tidyverse)
library(lubridate)
library(broom)      # turns model results into tidy tables (installed with tidyverse)

theme_set(theme_minimal(base_size = 12))
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

save_plot  <- function(p, name, w = 8, h = 5) {
  ggsave(file.path("outputs/figures", name), p, width = w, height = h, dpi = 300)
}
save_table <- function(tbl, name) write_csv(tbl, file.path("outputs/tables", name))

# ---- Load processed data -----------------------------------------------------
# Re-declare factor types: CSV files do not remember them.
df <- read_csv("data/processed/coastal_monitoring_clean.csv", show_col_types = FALSE) |>
  mutate(
    habitat      = factor(habitat),
    season       = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    health_class = factor(health_class, levels = c("Poor", "Fair", "Good"), ordered = TRUE),
    date_dec     = decimal_date(sample_date)          # e.g. 2021.37, for trend models
  )

# =============================================================================
# 1. DESCRIBE BEFORE YOU MODEL
# WHY: summaries and histograms reveal skew, scale and oddities that decide
#      which tests and transformations are appropriate.
# =============================================================================
summary_by_habitat <- df |>
  group_by(habitat) |>
  summarise(
    n = n(),
    across(c(water_temp_c, dissolved_oxygen_mgl, ph, salinity_psu,
             turbidity_ntu, nitrate_mgl, wqi),
           list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
print(summary_by_habitat)
save_table(summary_by_habitat, "summary_by_habitat.csv")

p_dist <- df |>
  select(water_temp_c, dissolved_oxygen_mgl, ph, salinity_psu,
         turbidity_ntu, nitrate_mgl, chlorophyll_ugl, wqi) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(value)) +
  geom_histogram(bins = 30, fill = "#2C7FB8", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Distribution of key variables",
       subtitle = "Turbidity, nitrate and chlorophyll are right-skewed -> we use logs / rank methods",
       x = NULL, y = "Count")
save_plot(p_dist, "01_distributions.png", 10, 6)

# =============================================================================
# 2. CORRELATIONS (Q1)
# WHY: correlation gives a fast map of which variables move together.
#      We use SPEARMAN (rank-based) because several variables are skewed and
#      Spearman does not need a straight-line / normal-data assumption.
# REMINDER: correlation is not causation, and wqi is built FROM some of these
#      variables, so those correlations are partly by construction.
# =============================================================================
cor_vars <- c("water_temp_c", "salinity_psu", "dissolved_oxygen_mgl", "ph",
              "turbidity_ntu", "nitrate_mgl", "phosphate_mgl",
              "chlorophyll_ugl", "species_richness", "wqi")

cor_mat  <- cor(df[cor_vars], use = "pairwise.complete.obs", method = "spearman")
cor_long <- as.data.frame(cor_mat) |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "rho")
save_table(round(as.data.frame(cor_mat), 3) |> rownames_to_column("variable"),
           "spearman_correlations.csv")

p_cor <- ggplot(cor_long, aes(var1, var2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(rho, 2)), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-1, 1)) +
  labs(title = "Spearman correlation matrix", x = NULL, y = NULL, fill = "rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_cor, "02_correlation_heatmap.png", 9, 7)

# A formal test for one headline relationship: temperature vs dissolved oxygen.
# The p-value asks "could a correlation this strong appear by chance?"
temp_do_test <- tidy(cor.test(~ water_temp_c + dissolved_oxygen_mgl, data = df,
                              method = "spearman", exact = FALSE))
print(temp_do_test)

p_temp_do <- ggplot(df, aes(water_temp_c, dissolved_oxygen_mgl, colour = habitat)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  labs(title = "Warmer water holds less dissolved oxygen",
       x = "Water temperature (C)", y = "Dissolved oxygen (mg/L)", colour = "Habitat")
save_plot(p_temp_do, "03_temp_vs_do.png")

# =============================================================================
# 3. GROUP COMPARISONS (Q2)
# WHY: are differences between sites/seasons real or just noise? Kruskal-Wallis
#      is the rank-based version of ANOVA (no normality assumption). If it is
#      significant, pairwise Wilcoxon tests with Holm correction show WHICH
#      groups differ; the correction stops "many tests" producing false alarms.
# =============================================================================
compare_groups <- function(data, response, group) {
  kw <- kruskal.test(reformulate(group, response), data = data)
  pw <- pairwise.wilcox.test(data[[response]], data[[group]],
                             p.adjust.method = "holm", exact = FALSE)
  list(kruskal = tidy(kw), pairwise_p = round(pw$p.value, 4))
}

wqi_by_habitat  <- compare_groups(df, "wqi", "habitat")
do_by_season    <- compare_groups(df, "dissolved_oxygen_mgl", "season")
print(wqi_by_habitat)
print(do_by_season)
save_table(wqi_by_habitat$kruskal, "kruskal_wqi_by_habitat.csv")
save_table(do_by_season$kruskal,   "kruskal_do_by_season.csv")

p_wqi <- ggplot(df, aes(reorder(habitat, wqi, median, na.rm = TRUE), wqi, fill = habitat)) +
  geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +
  labs(title = "Water Quality Index differs strongly between habitats",
       x = "Habitat (ordered by median WQI)", y = "WQI (0-100, higher = better)")
save_plot(p_wqi, "04_wqi_by_habitat.png")

p_do_season <- ggplot(df, aes(season, dissolved_oxygen_mgl, fill = season)) +
  geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +
  labs(title = "Dissolved oxygen is lowest in summer",
       x = "Season", y = "Dissolved oxygen (mg/L)")
save_plot(p_do_season, "05_do_by_season.png")

# =============================================================================
# 4. LINEAR REGRESSION (Q3a): what predicts dissolved oxygen?
# WHY: regression measures each predictor's effect while holding the others
#      fixed, which plain correlation cannot do. We log-transform skewed
#      predictors so the relationship is closer to linear.
# Read the table: 'estimate' = change in DO per one-unit change in the
#      predictor; the confidence interval shows the uncertainty.
# =============================================================================
m_do <- lm(dissolved_oxygen_mgl ~ water_temp_c + log(chlorophyll_ugl) +
             log(nitrate_mgl) + salinity_psu + habitat, data = df)

do_coefs <- tidy(m_do, conf.int = TRUE) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
print(do_coefs)
print(glance(m_do) |> select(r.squared, adj.r.squared, sigma, nobs))
save_table(do_coefs, "lm_dissolved_oxygen.csv")

# Diagnostics: a model is only trustworthy if its assumptions roughly hold.
aug <- augment(m_do)
p_resid <- ggplot(aug, aes(.fitted, .resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, colour = "red", method = "loess", formula = y ~ x) +
  labs(title = "Residuals vs fitted (want a flat, patternless cloud)",
       x = "Fitted DO", y = "Residual")
save_plot(p_resid, "06_lm_residuals.png")

p_qq <- ggplot(aug, aes(sample = .std.resid)) +
  stat_qq(alpha = 0.5) + stat_qq_line(colour = "red") +
  labs(title = "Normal Q-Q plot of residuals (points near the line = good)",
       x = "Theoretical quantiles", y = "Standardised residuals")
save_plot(p_qq, "07_lm_qq.png")

# =============================================================================
# 5. POISSON REGRESSION (Q3b): what predicts species richness?
# WHY: species richness is a COUNT (0, 1, 2, ...), so an ordinary straight-line
#      model is a poor fit; Poisson regression is built for counts.
#      Exponentiated coefficients are rate ratios (1.10 = 10% more species).
#      Dispersion should be near 1; much larger means the model is
#      overconfident and a quasi-Poisson / negative binomial model is safer.
# =============================================================================
m_sp <- glm(species_richness ~ dissolved_oxygen_mgl + log(nitrate_mgl) + habitat,
            family = poisson, data = df)
dispersion <- sum(residuals(m_sp, type = "pearson")^2) / df.residual(m_sp)
message("Poisson dispersion = ", round(dispersion, 2), " (about 1 is fine)")

sp_coefs <- tidy(m_sp, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(sp_coefs)
save_table(sp_coefs, "poisson_species_richness.csv")

# =============================================================================
# 6. LOGISTIC REGRESSION (Q4): predicting "Good" ecosystem health
# WHY: the target is yes/no, so we model the PROBABILITY of "Good". To test
#      honestly we train on 70% of the data and evaluate on the 30% the model
#      has never seen; scoring a model on its own training data flatters it.
# Read the table: exponentiated coefficients are ODDS RATIOS (above 1 raises
#      the odds of Good health, below 1 lowers them).
# CAUTION: the health class and our predictors are physically linked, so high
#      accuracy partly reflects that, not a discovery.
# =============================================================================
set.seed(42)
model_df <- df |>
  drop_na(good_health, water_temp_c, dissolved_oxygen_mgl, nitrate_mgl, turbidity_ntu)

train_idx <- sample(nrow(model_df), size = round(0.7 * nrow(model_df)))
train <- model_df[train_idx, ]
test  <- model_df[-train_idx, ]

m_log <- glm(good_health ~ dissolved_oxygen_mgl + log(nitrate_mgl) +
               turbidity_ntu + water_temp_c, family = binomial, data = train)

or_table <- tidy(m_log, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(or_table)
save_table(or_table, "logistic_odds_ratios.csv")

test$p_good <- predict(m_log, newdata = test, type = "response")
test$pred   <- as.integer(test$p_good >= 0.5)

conf_matrix <- table(Actual = test$good_health, Predicted = test$pred)
accuracy    <- mean(test$pred == test$good_health)
baseline    <- max(mean(test$good_health), 1 - mean(test$good_health))  # always guess majority

# AUC via the rank formula: probability a random "Good" site scores higher than
# a random non-"Good" site (0.5 = coin flip, 1 = perfect).
auc <- function(score, y) {
  r  <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
test_auc <- auc(test$p_good, test$good_health)

print(conf_matrix)
message("Test accuracy: ", round(accuracy, 3),
        " | majority-class baseline: ", round(baseline, 3),
        " | AUC: ", round(test_auc, 3))

p_pred <- ggplot(test, aes(p_good, fill = factor(good_health, labels = c("Not Good", "Good")))) +
  geom_density(alpha = 0.5) +
  labs(title = "Predicted probability of 'Good' health on unseen data",
       subtitle = "Clear separation between the two groups = useful model",
       x = "Predicted P(Good)", y = "Density", fill = "Actual class")
save_plot(p_pred, "08_logistic_predictions.png")

# =============================================================================
# 7. TREND CHECKS (Q5)
# WHY: monitoring programmes exist to detect change. We include season in the
#      model so the strong yearly cycle is not mistaken for a trend.
# LIMITATION: successive monthly readings at one site are not independent
#      (autocorrelation), so these p-values are somewhat too optimistic.
# =============================================================================
temp_trend <- lm(water_temp_c ~ date_dec + season, data = df)
temp_slope <- tidy(temp_trend, conf.int = TRUE) |>
  filter(term == "date_dec") |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(temp_slope)   # estimate = change in degrees C per YEAR

wqi_trends <- df |>
  drop_na(wqi) |>
  group_by(habitat) |>
  group_modify(~ tidy(lm(wqi ~ date_dec + season, data = .x), conf.int = TRUE)) |>
  filter(term == "date_dec") |>
  ungroup() |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(wqi_trends)   # estimate = WQI points per YEAR at each habitat
save_table(temp_slope, "trend_temperature.csv")
save_table(wqi_trends, "trend_wqi_by_habitat.csv")

p_trend <- ggplot(df, aes(sample_date, water_temp_c)) +
  geom_point(alpha = 0.25, colour = "grey40") +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, colour = "#D95F0E", formula = y ~ x) +
  labs(title = "Water temperature: strong seasonal cycle with a slight upward drift",
       x = "Date", y = "Water temperature (C)")
save_plot(p_trend, "09_temperature_timeseries.png")

message("\nAnalysis complete. See outputs/figures and outputs/tables.")
