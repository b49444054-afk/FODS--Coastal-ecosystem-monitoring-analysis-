# =============================================================================
# 01_generate_data.R
# Coastal Ecosystem Monitoring Analysis - Step 1: create the RAW dataset
# -----------------------------------------------------------------------------
# What this does
#   Builds a realistic *synthetic* monitoring dataset (6 stations x 72 months)
#   whose variables are related to each other the way they are in real coastal
#   water (warm water holds less oxygen, nutrients drive algae, etc.), then
#   deliberately "dirties" it the way real field data is dirty:
#     - missing values (random, and one non-random pattern)
#     - impossible sensor readings (outliers)
#     - sentinel codes (-999) and text codes ("ND" = not detected)
#     - inconsistent site IDs, labels and date formats
#     - duplicated rows
#
# Why synthetic?  You get a dataset where the "truth" is known, so you can check
# whether your cleaning/analysis recovers it. For a real-data version see README.
#
# Run from the project root:  source("01_generate_data.R")
# Output: data/raw/coastal_monitoring_raw.csv
# Requires R >= 4.1 and the tidyverse (install.packages("tidyverse")).
# =============================================================================

library(tidyverse)
library(lubridate)

set.seed(2026)  # makes every "random" step reproducible - essential for science

# ---- 1. Station metadata -----------------------------------------------------
# Each habitat has different baseline conditions. These hidden parameters are
# what make site comparisons in the analysis meaningful.
stations <- tibble(
  site_id       = sprintf("ST-%02d", 1:6),
  site_name     = c("Mangrove Creek", "River Estuary", "Reef Point",
                    "Tidal Lagoon", "Seagrass Bay", "Harbour Basin"),
  habitat       = c("Mangrove", "Estuary", "Reef", "Lagoon", "Seagrass", "Harbour"),
  base_salinity = c(28, 22, 35, 31, 34, 33),          # PSU
  freshwater    = c(0.8, 1.5, 0.2, 0.5, 0.2, 0.4),    # river/rain influence
  nutrient_load = c(1.0, 1.8, 0.4, 0.9, 0.5, 2.2),    # runoff / sewage pressure
  impact        = c(0.3, 0.8, 0.1, 0.4, 0.2, 1.0),    # general human pressure
  base_richness = c(30, 22, 45, 28, 38, 15)           # typical species count
)

# ---- 2. Sampling grid: every station, monthly, 2018-2023 ---------------------
months <- seq(ymd("2018-01-01"), ymd("2023-12-01"), by = "month")

true_data <- crossing(stations, month_start = months) |>
  mutate(
    sample_date = month_start + days(sample(0:27, n(), replace = TRUE)),
    doy         = yday(sample_date),
    years_since = decimal_date(sample_date) - 2018,
    heat_wave   = sin(2 * pi * (doy - 110) / 365),         # peaks in July
    wetness     = pmax(0, cos(2 * pi * (doy - 20) / 365))  # wet season peaks Jan
  ) |>
  # ---- 3. Simulate the measurements (order matters: later ones use earlier) --
  mutate(
    # Temperature: seasonal cycle + slow warming (0.2 C/yr) + warmer harbour
    water_temp_c  = 22 + 6 * heat_wave + 0.2 * years_since +
                    1.0 * (habitat == "Harbour") + rnorm(n(), 0, 1.2),
    # Salinity: dips in the wet season, most at river-influenced stations
    salinity_psu  = base_salinity - 5 * freshwater * wetness + rnorm(n(), 0, 1),
    # Turbidity, nutrients, chlorophyll: positive and right-skewed -> lognormal
    turbidity_ntu = exp(rnorm(n(), log(2 + 4 * freshwater * (1 + wetness) + 2 * impact), 0.45)),
    nitrate_mgl   = exp(rnorm(n(), log(0.15 * nutrient_load * (1 + 0.8 * wetness)), 0.4)),
    phosphate_mgl = exp(rnorm(n(), log(0.03 * nutrient_load * (1 + 0.5 * wetness)), 0.4)),
    chlorophyll_ugl = exp(rnorm(n(),
                          log(pmax(0.5, 1 + 10 * nitrate_mgl + 0.15 * (water_temp_c - 18))),
                          0.35)),
    # Dissolved oxygen: falls with temperature and algae load
    dissolved_oxygen_mgl = pmax(1, 11.5 - 0.17 * water_temp_c - 0.12 * chlorophyll_ugl -
                                   0.3 * impact + rnorm(n(), 0, 0.5)),
    ph = 8.1 + 0.02 * (chlorophyll_ugl - 3) - 0.1 * impact -
         0.1 * freshwater * wetness + rnorm(n(), 0, 0.08),
    # Biodiversity indicator: a COUNT, so Poisson noise is appropriate
    species_richness = rpois(n(), exp(log(base_richness) +
                                      0.10 * (dissolved_oxygen_mgl - 7) -
                                      0.8  * nitrate_mgl +
                                      0.4  * (ph - 8.1))),
    # Hidden "true" condition -> observer-assigned health class (noisy)
    health_latent = 0.6 * (dissolved_oxygen_mgl - 7) - 2 * nitrate_mgl -
                    0.04 * turbidity_ntu - 0.1 * (chlorophyll_ugl - 3) +
                    rnorm(n(), 0, 0.4)
  )

cuts <- quantile(true_data$health_latent, c(0.25, 0.65))
true_data <- true_data |>
  mutate(health_class = as.character(
    cut(health_latent, breaks = c(-Inf, cuts, Inf), labels = c("Poor", "Fair", "Good"))
  ))

# ---- 4. Keep only what a field team would actually record --------------------
# (heat_wave, wetness, health_latent etc. are hidden simulation internals)
observed <- true_data |>
  transmute(
    site_id, site_name, habitat, sample_date,
    water_temp_c         = round(water_temp_c, 1),
    dissolved_oxygen_mgl = round(dissolved_oxygen_mgl, 2),
    ph                   = round(ph, 2),
    salinity_psu         = round(salinity_psu, 1),
    turbidity_ntu        = round(turbidity_ntu, 1),
    nitrate_mgl          = round(nitrate_mgl, 3),
    phosphate_mgl        = round(phosphate_mgl, 4),
    chlorophyll_ugl      = round(chlorophyll_ugl, 2),
    species_richness,
    health_class
  )

# ---- 5. Make it messy, like real data ---------------------------------------
n_rows <- nrow(observed)
messy  <- observed

# Helper: random row positions covering a proportion p of the data
pick <- function(p) sample(n_rows, round(p * n_rows))

# 5a. Missing values at random (sensor dropout, lost samples)
miss_rates <- c(water_temp_c = .03, dissolved_oxygen_mgl = .05, ph = .04,
                salinity_psu = .03, turbidity_ntu = .06, nitrate_mgl = .08,
                phosphate_mgl = .10, chlorophyll_ugl = .07,
                species_richness = .02, health_class = .02)
for (col in names(miss_rates)) {
  messy[[col]][pick(miss_rates[[col]])] <- NA
}

# 5b. Missingness NOT at random: the DO sensor fails more often in hot water.
#     (A key concept: *why* data is missing affects how safely you can fill it.)
hot <- which(messy$water_temp_c > 27)
messy$dissolved_oxygen_mgl[sample(hot, round(0.15 * length(hot)))] <- NA

# 5c. Impossible readings = sensor faults / typing errors (outliers)
messy$water_temp_c[sample(n_rows, 3)]         <- c(45.0, 99.9, -5.0)
messy$dissolved_oxygen_mgl[sample(n_rows, 3)] <- c(25.4, -3.1, 19.9)
messy$ph[sample(n_rows, 2)]                   <- c(14.2, 2.8)
messy$salinity_psu[sample(n_rows, 2)]         <- c(0.0, 78.0)
messy$turbidity_ntu[sample(n_rows, 2)]        <- c(850, 420)   # 420 = plausible storm spike

# 5d. Sentinel codes: -999 means "no data" in many legacy instruments
messy$dissolved_oxygen_mgl[sample(n_rows, 5)] <- -999
messy$ph[sample(n_rows, 4)]                   <- -999

# 5e. Text codes inside numeric columns: "ND" = not detected (below lab limit)
messy <- messy |>
  mutate(nitrate_mgl = as.character(nitrate_mgl),
         phosphate_mgl = as.character(phosphate_mgl))
messy$nitrate_mgl[sample(n_rows, 6)]    <- "ND"
messy$phosphate_mgl[sample(n_rows, 10)] <- "ND"
messy$phosphate_mgl[sample(n_rows, 3)]  <- "-999"

# 5f. Mixed date formats (ISO vs day/month/year typed by hand)
messy <- messy |> mutate(sample_date = as.character(sample_date))
idx <- pick(0.12)
messy$sample_date[idx] <- format(observed$sample_date[idx], "%d/%m/%Y")

# 5g. Inconsistent IDs and labels. Only touch non-missing labels, because
#     paste0(NA, " ") would silently turn a real NA into the text "NA ".
idx <- pick(0.08)
messy$site_id[idx] <- tolower(gsub("-", "", messy$site_id[idx]))           # "st01"
idx <- pick(0.05)
messy$site_id[idx] <- paste0(" ", messy$site_id[idx], " ")                 # " ST-03 "

idx <- intersect(pick(0.10), which(!is.na(messy$health_class)))
messy$health_class[idx] <- toupper(messy$health_class[idx])                # "GOOD"
idx <- intersect(pick(0.05), which(!is.na(messy$health_class)))
messy$health_class[idx] <- paste0(messy$health_class[idx], " ")            # "Fair "

# 5h. Duplicate records (e.g. a file uploaded twice)
messy <- bind_rows(messy, slice_sample(messy, n = 6))

# ---- 6. Save -----------------------------------------------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
write_csv(messy, "data/raw/coastal_monitoring_raw.csv")

message("Saved data/raw/coastal_monitoring_raw.csv  (",
        nrow(messy), " rows x ", ncol(messy), " columns)")
