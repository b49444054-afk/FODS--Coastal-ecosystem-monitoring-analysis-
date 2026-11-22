# =============================================================================
# 02_clean_data.R
# Coastal Ecosystem Monitoring Analysis - Step 2: RAW -> PROCESSED
# -----------------------------------------------------------------------------
# Each cleaning step is its own small function so you can read, test, and
# explain them one at a time. The pipeline at the bottom chains them together.
#
# Golden rules followed here
#   1. Never edit the raw file; always write a new processed file.
#   2. Every step reports what it changed (an audit trail).
#   3. Prefer FLAGGING over silently deleting when a value is merely unusual.
#
# Run from the project root:  source("02_clean_data.R")
# Input : data/raw/coastal_monitoring_raw.csv
# Output: data/processed/coastal_monitoring_clean.csv
# =============================================================================

library(tidyverse)
library(lubridate)

raw_path <- "data/raw/coastal_monitoring_raw.csv"
out_path <- "data/processed/coastal_monitoring_clean.csv"

numeric_cols <- c("water_temp_c", "dissolved_oxygen_mgl", "ph", "salinity_psu",
                  "turbidity_ntu", "nitrate_mgl", "phosphate_mgl",
                  "chlorophyll_ugl", "species_richness")

# Variables we are willing to impute (species_richness is a count/outcome, so
# we do not invent values for it).
impute_vars <- c("water_temp_c", "dissolved_oxygen_mgl", "ph", "salinity_psu",
                 "turbidity_ntu", "nitrate_mgl", "phosphate_mgl", "chlorophyll_ugl")

# Physically plausible ranges. A value outside these is not "unusual", it is
# impossible (e.g. pH 14.2 in seawater), so we treat it as an error.
plausible_limits <- tribble(
  ~variable,              ~lower, ~upper,
  "water_temp_c",            0,     38,
  "dissolved_oxygen_mgl",    0,     20,
  "ph",                      6.5,   9.5,
  "salinity_psu",            0,     45,
  "turbidity_ntu",           0,    500,
  "nitrate_mgl",             0,     10,
  "phosphate_mgl",           0,      2,
  "chlorophyll_ugl",         0,    200,
  "species_richness",        0,    200
)

# -----------------------------------------------------------------------------
# STEP 1: Load everything as text
# WHY: if R guesses types, a stray "ND" or "-999" can silently turn a whole
#      column into text or hide errors. Reading as text lets us convert on
#      purpose, and we can count exactly what failed to convert.
# -----------------------------------------------------------------------------
load_raw <- function(path) {
  read_csv(
    path,
    col_types = cols(.default = col_character()),
    na        = c("", "NA", "N/A", "n/a", "ND"),   # "ND" = not detected -> NA
    trim_ws   = TRUE,
    show_col_types = FALSE
  )
}

# -----------------------------------------------------------------------------
# STEP 2: Inspect structure
# WHY: you cannot clean what you have not looked at. Structure, missingness and
#      duplicates drive every later decision.
# -----------------------------------------------------------------------------
inspect_data <- function(df, title) {
  message("\n==== ", title, " ====")
  message("Rows: ", nrow(df), " | Columns: ", ncol(df))
  glimpse(df)
  message("Missing values per column:")
  print(colSums(is.na(df)))
  message("Exact duplicate rows: ", sum(duplicated(df)))
  invisible(df)
}

# -----------------------------------------------------------------------------
# STEP 3: Fix data types
# WHY: arithmetic, correlations and models need numbers. Anything that cannot be
#      converted becomes NA; we report how many so nothing disappears silently.
# -----------------------------------------------------------------------------
convert_types <- function(df) {
  before <- sum(is.na(df[numeric_cols]))
  df <- df |>
    mutate(across(all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x))))
  after <- sum(is.na(df[numeric_cols]))
  message("convert_types: ", after - before, " value(s) could not be read as numbers -> NA")
  df
}

# -----------------------------------------------------------------------------
# STEP 4: Standardise categorical labels
# WHY: "st01", " ST-01 " and "ST-01" are one station to a human but three
#      groups to R. Inconsistent labels would split group comparisons.
# -----------------------------------------------------------------------------
standardise_labels <- function(df) {
  df |>
    mutate(
      site_id      = sprintf("ST-%02d", as.integer(str_extract(site_id, "\\d+"))),
      health_class = str_to_title(str_trim(health_class)),
      health_class = factor(health_class, levels = c("Poor", "Fair", "Good"), ordered = TRUE)
    )
}

# -----------------------------------------------------------------------------
# STEP 5: Parse dates
# WHY: dates stored as text cannot be sorted, grouped by season, or used in
#      trend models. Two formats are mixed, and "03/04/2021" is ambiguous, so we
#      decide the format from the text pattern instead of letting R guess.
# -----------------------------------------------------------------------------
parse_dates <- function(df) {
  df <- df |>
    mutate(
      iso_txt = if_else(str_detect(sample_date, "^\\d{4}-\\d{2}-\\d{2}$"), sample_date, NA_character_),
      dmy_txt = if_else(str_detect(sample_date, "^\\d{2}/\\d{2}/\\d{4}$"), sample_date, NA_character_),
      sample_date = coalesce(ymd(iso_txt, quiet = TRUE), dmy(dmy_txt, quiet = TRUE))
    ) |>
    select(-iso_txt, -dmy_txt)
  bad <- sum(is.na(df$sample_date))
  message("parse_dates: ", bad, " unparseable date(s) removed")
  filter(df, !is.na(sample_date))
}

# -----------------------------------------------------------------------------
# STEP 6: Remove duplicates
# WHY: duplicates double-count observations, which makes results look more
#      certain than they are. Done AFTER standardising so "st01" == "ST-01".
# -----------------------------------------------------------------------------
remove_duplicates <- function(df) {
  n0 <- nrow(df)
  df <- distinct(df)
  message("remove_duplicates: ", n0 - nrow(df), " duplicate row(s) removed")
  dup_keys <- sum(duplicated(df[c("site_id", "sample_date")]))
  if (dup_keys > 0) message("  WARNING: ", dup_keys, " site+date pair(s) still repeated - investigate")
  df
}

# -----------------------------------------------------------------------------
# STEP 7: Sentinel codes -> NA
# WHY: -999 is a placeholder, not a measurement. Left in, it would drag means
#      and correlations to nonsense.
# -----------------------------------------------------------------------------
replace_sentinels <- function(df, sentinel = -999) {
  n <- sum(df[numeric_cols] == sentinel, na.rm = TRUE)
  message("replace_sentinels: ", n, " sentinel value(s) set to NA")
  df |> mutate(across(all_of(numeric_cols), ~ if_else(.x == sentinel, NA_real_, .x)))
}

# -----------------------------------------------------------------------------
# STEP 8a: Outliers, part 1 - impossible values (domain knowledge)
# WHY: some values cannot be real. Removing them is justified by science,
#      not by the shape of the data.
# -----------------------------------------------------------------------------
apply_plausibility_limits <- function(df, limits = plausible_limits) {
  for (i in seq_len(nrow(limits))) {
    v   <- limits$variable[i]
    bad <- !is.na(df[[v]]) & (df[[v]] < limits$lower[i] | df[[v]] > limits$upper[i])
    if (any(bad)) message("  ", v, ": ", sum(bad), " impossible value(s) set to NA")
    df[[v]][bad] <- NA
  }
  df
}

# -----------------------------------------------------------------------------
# STEP 8b: Outliers, part 2 - statistically extreme but possible values (flag)
# WHY: a storm can truly produce very turbid water. Deleting real extremes
#      would hide the very events monitoring exists to detect. So we only FLAG
#      them (3 x IQR rule, within each station) and decide later.
#      Skewed variables are judged on the log scale, otherwise their long
#      natural tail would be flagged constantly.
# -----------------------------------------------------------------------------
flag_extremes <- function(df, vars, transform = identity, k = 3) {
  df |>
    group_by(site_id) |>
    mutate(across(all_of(vars), ~ {
      x      <- transform(.x)
      q      <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
      spread <- q[[2]] - q[[1]]
      !is.na(x) & (x < q[[1]] - k * spread | x > q[[2]] + k * spread)
    }, .names = "{.col}_extreme")) |>
    ungroup()
}

# -----------------------------------------------------------------------------
# STEP 9: Time features
# WHY: seasonality is the biggest driver in coastal data. Season/year let us
#      compare like with like and separate trends from seasonal cycles.
#      (Northern-hemisphere meteorological seasons; adapt for monsoon regions.)
# -----------------------------------------------------------------------------
derive_time_features <- function(df) {
  df |>
    mutate(
      year      = year(sample_date),
      month_num = month(sample_date),
      season    = case_when(
        month_num %in% c(12, 1, 2) ~ "Winter",
        month_num %in% 3:5         ~ "Spring",
        month_num %in% 6:8         ~ "Summer",
        TRUE                       ~ "Autumn"
      ),
      season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))
    )
}

# -----------------------------------------------------------------------------
# STEP 10: Missing values
# Strategy:
#   (a) Rows with no health_class are dropped: it is the outcome we want to
#       study, and inventing an outcome would be circular.
#   (b) Missing predictors are filled with the median of the same station and
#       season, because conditions differ strongly by both. Medians are robust
#       to the outliers we flagged. We record HOW MANY values were imputed per
#       row so you can check results do not depend on imputation.
# Caveat: imputation shrinks variability and assumes the data are missing at
#         random. The hot-water DO dropout in this dataset breaks that
#         assumption a little - mention it in your limitations.
# -----------------------------------------------------------------------------
drop_missing_outcome <- function(df) {
  n0 <- nrow(df)
  df <- filter(df, !is.na(health_class))
  message("drop_missing_outcome: ", n0 - nrow(df), " row(s) without health_class removed")
  df
}

impute_site_season_median <- function(df, vars) {
  df <- df |> mutate(n_imputed = rowSums(is.na(pick(all_of(vars)))))
  message("impute: ", sum(df$n_imputed), " value(s) to fill in ", sum(df$n_imputed > 0), " row(s)")
  df |>
    group_by(site_id, season) |>
    mutate(across(all_of(vars), ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) |>
    ungroup()
}

# -----------------------------------------------------------------------------
# STEP 11: Composite Water Quality Index (WQI, 0-100, higher = healthier)
# WHY: six noisy variables are hard to compare across sites. An index puts them
#      on one scale. Each variable is scored 0-100 between a "bad" and a "good"
#      threshold, then combined with weights. THESE THRESHOLDS AND WEIGHTS ARE
#      OUR OWN ASSUMPTIONS, not an official standard - say so in your report.
#      Temperature and salinity are left out because they vary naturally with
#      habitat and season rather than signalling "good" or "bad".
# -----------------------------------------------------------------------------
score_linear <- function(x, bad, good) {
  pmin(pmax((x - bad) / (good - bad), 0), 1) * 100   # works whether high or low is good
}

add_wqi <- function(df) {
  df |>
    mutate(
      score_do   = score_linear(dissolved_oxygen_mgl, bad = 3,    good = 7),
      score_ph   = pmax(0, 100 * (1 - abs(ph - 8.1) / 1)),         # best near 8.1
      score_turb = score_linear(turbidity_ntu,        bad = 25,   good = 3),
      score_n    = score_linear(nitrate_mgl,          bad = 0.8,  good = 0.1),
      score_p    = score_linear(phosphate_mgl,        bad = 0.12, good = 0.02),
      score_chl  = score_linear(chlorophyll_ugl,      bad = 12,   good = 2),
      wqi = 0.25 * score_do + 0.10 * score_ph + 0.15 * score_turb +
            0.20 * score_n  + 0.10 * score_p  + 0.20 * score_chl,   # weights sum to 1
      wqi_band = cut(wqi, breaks = c(-Inf, 50, 70, Inf),
                     labels = c("Low", "Moderate", "High"),
                     right = FALSE, ordered_result = TRUE)
    )
}

# -----------------------------------------------------------------------------
# STEP 12: Other derived features
# WHY: models and plots often need re-expressed variables. A binary target is
#      required for logistic regression; N:P ratio is a standard ecology metric.
# -----------------------------------------------------------------------------
add_features <- function(df) {
  df |>
    mutate(
      good_health = as.integer(health_class == "Good"),    # 1 = Good, 0 = Fair/Poor
      np_ratio    = nitrate_mgl / phosphate_mgl,
      bloom_risk  = chlorophyll_ugl > 10                    # simple algal-bloom flag
    )
}

# -----------------------------------------------------------------------------
# STEP 13: Validate before saving
# WHY: automatic checks catch mistakes you would otherwise only notice later.
# -----------------------------------------------------------------------------
validate_clean <- function(df) {
  stopifnot(
    "duplicate site+date rows remain" = !any(duplicated(df[c("site_id", "sample_date")])),
    "dates outside 2018-2023"         = all(df$sample_date >= ymd("2018-01-01") &
                                            df$sample_date <= ymd("2023-12-31")),
    "site_id format wrong"            = all(str_detect(df$site_id, "^ST-\\d{2}$")),
    "pH out of range"                 = all(df$ph >= 6.5 & df$ph <= 9.5, na.rm = TRUE),
    "negative DO"                     = all(df$dissolved_oxygen_mgl >= 0, na.rm = TRUE),
    "WQI outside 0-100"               = all(df$wqi >= 0 & df$wqi <= 100, na.rm = TRUE)
  )
  message("validate_clean: all checks passed")
  invisible(df)
}

# =============================================================================
# THE PIPELINE: raw -> processed
# =============================================================================
clean_pipeline <- function(raw_path, out_path) {
  raw <- load_raw(raw_path)
  inspect_data(raw, "RAW DATA (as loaded)")

  clean <- raw |>
    convert_types() |>
    standardise_labels() |>
    parse_dates() |>
    remove_duplicates() |>
    replace_sentinels() |>
    apply_plausibility_limits() |>
    derive_time_features() |>
    flag_extremes(c("water_temp_c", "dissolved_oxygen_mgl", "ph", "salinity_psu")) |>
    flag_extremes(c("turbidity_ntu", "nitrate_mgl", "phosphate_mgl", "chlorophyll_ugl"),
                  transform = log) |>
    drop_missing_outcome() |>
    impute_site_season_median(impute_vars) |>
    add_wqi() |>
    add_features() |>
    arrange(site_id, sample_date)

  validate_clean(clean)
  inspect_data(clean, "PROCESSED DATA")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(clean, out_path)
  message("\nSaved ", out_path, " (", nrow(clean), " rows x ", ncol(clean), " columns)")
  invisible(clean)
}

clean_data <- clean_pipeline(raw_path, out_path)
