library(dplyr)
library(tidyr)
library(readr)

df <- read_csv("CVD_Data.csv", skip = 1)

# ---- Aggregate to long format: deaths summed by county and day ----
df_long <- df %>%
  group_by(FAC_COUNTY, DOD_STRING) %>%
  summarise(
    heart_disease_deaths   = sum(heart_disease_deaths, na.rm = TRUE),
    stroke_deaths          = sum(stroke_deaths, na.rm = TRUE),
    hypertension_deaths    = sum(hypertension_deaths, na.rm = TRUE),
    other_vascular_deaths  = sum(other_vascular_deaths, na.rm = TRUE),
    total_cvd_deaths       = sum(total_cvd_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(heart_disease_deaths, stroke_deaths, hypertension_deaths,
             other_vascular_deaths, total_cvd_deaths),
    names_to = "cvd_type",
    values_to = "deaths"
  ) %>%
  arrange(FAC_COUNTY, DOD_STRING)

# ---- Daily death rate per county (total_cvd_deaths only) ----
# Full study period: 2005-01-01 through 2025-12-31
n_days_total <- as.integer(as.Date("2025-12-31") - as.Date("2005-01-01")) + 1

county_daily_deaths <- df_long %>%
  filter(cvd_type == "total_cvd_deaths") %>%
  group_by(FAC_COUNTY) %>%
  summarise(
    total_cvd_deaths = sum(deaths),
    n_days           = n_days_total,
    daily_death_rate = total_cvd_deaths / n_days_total,
    .groups = "drop"
  ) %>%
  arrange(desc(daily_death_rate))

county_daily_deaths

# ---- Clean county names & combine Sarpy + Douglas into "Omaha" ----
df_long_clean <- df_long %>%
  mutate(FAC_COUNTY = case_when(
    FAC_COUNTY %in% c("butler", "Butler")     ~ "Butler",
    FAC_COUNTY %in% c("cherry", "Cherry")     ~ "Cherry",
    FAC_COUNTY %in% c("hall", "Hall")         ~ "Hall",
    FAC_COUNTY %in% c("keith", "Keith")       ~ "Keith",
    FAC_COUNTY %in% c("madison", "Madison")   ~ "Madison",
    FAC_COUNTY %in% c("platte", "Platte")     ~ "Platte",
    FAC_COUNTY %in% c("thurston", "Thurston") ~ "Thurston",
    FAC_COUNTY %in% c("Sarpy", "Douglas")     ~ "Omaha",
    TRUE ~ FAC_COUNTY
  )) %>%
  filter(!FAC_COUNTY %in% c("NULL", "United States", "Unknown")) %>%
  filter(!is.na(FAC_COUNTY)) %>%
  group_by(FAC_COUNTY, DOD_STRING, cvd_type) %>%
  summarise(deaths = sum(deaths), .groups = "drop") %>%
  arrange(FAC_COUNTY, DOD_STRING)

df_long_clean

# ---- New daily death rate table: manually combined names, nulls removed ----
new_county_daily_deaths <- df_long %>%
  filter(cvd_type == "total_cvd_deaths") %>%
  mutate(FAC_COUNTY = case_when(
    FAC_COUNTY %in% c("butler", "Butler")     ~ "Butler",
    FAC_COUNTY %in% c("cherry", "Cherry")     ~ "Cherry",
    FAC_COUNTY %in% c("hall", "Hall")         ~ "Hall",
    FAC_COUNTY %in% c("keith", "Keith")       ~ "Keith",
    FAC_COUNTY %in% c("madison", "Madison")   ~ "Madison",
    FAC_COUNTY %in% c("platte", "Platte")     ~ "Platte",
    FAC_COUNTY %in% c("thurston", "Thurston") ~ "Thurston",
    FAC_COUNTY %in% c("Sarpy", "Douglas")     ~ "Omaha",
    TRUE ~ FAC_COUNTY
  )) %>%
  filter(!FAC_COUNTY %in% c("NULL", "United States", "Unknown")) %>%
  filter(!is.na(FAC_COUNTY)) %>%
  group_by(FAC_COUNTY) %>%
  summarise(
    total_cvd_deaths = sum(deaths),
    n_days           = n_days_total,
    daily_death_rate = total_cvd_deaths / n_days_total,
    .groups = "drop"
  ) %>%
  arrange(desc(daily_death_rate))

new_county_daily_deaths