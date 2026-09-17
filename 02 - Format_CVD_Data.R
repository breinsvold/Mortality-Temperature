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

# ---- Clean names, combine into six groups, update the dataframe ----
# Three city groups (pulled out of the regions)
omaha_counties     <- c("Douglas", "Sarpy")
lincoln_counties   <- c("Lancaster")              # Lincoln the CITY = Lancaster County
tricities_counties <- c("Hall", "Buffalo", "Adams", "Kearney")

# Three regions (remaining counties; all 93 assigned exactly once)
west_counties <- c(
  "Dundy","Hitchcock","Red Willow","Chase","Hayes","Frontier","Perkins",
  "Kimball","Cheyenne","Deuel","Keith","Lincoln","Banner","Morrill","Garden",
  "Scotts Bluff","Sioux","Dawes","Box Butte","Sheridan","Grant","Hooker",
  "Thomas","Arthur","McPherson","Logan","Blaine","Loup","Cherry","Brown",
  "Rock","Keya Paha", "Custer"
)

southeast_counties <- c(
  "Sherman","Howard","Merrick","Polk","Butler","Saunders","Dawson",
  "Hamilton","York","Seward","Cass","Gosper","Phelps","Clay","Fillmore",
  "Saline","Otoe","Harlan","Franklin","Webster","Nuckolls","Thayer",
  "Jefferson","Gage","Johnson","Nemaha","Pawnee","Richardson","Furnas"
)

northeast_counties <- c(
  "Garfield","Wheeler","Valley","Greeley","Boone","Nance","Antelope","Holt",
  "Madison","Platte","Colfax","Stanton","Pierce","Wayne","Cedar","Knox","Boyd",
  "Dixon","Dakota","Thurston","Cuming","Burt","Dodge","Washington"
)

CVD_data <- df_long %>%
  mutate(FAC_COUNTY = case_when(
    FAC_COUNTY %in% c("butler", "Butler")     ~ "Butler",
    FAC_COUNTY %in% c("cherry", "Cherry")     ~ "Cherry",
    FAC_COUNTY %in% c("hall", "Hall")         ~ "Hall",
    FAC_COUNTY %in% c("keith", "Keith")       ~ "Keith",
    FAC_COUNTY %in% c("madison", "Madison")   ~ "Madison",
    FAC_COUNTY %in% c("platte", "Platte")     ~ "Platte",
    FAC_COUNTY %in% c("thurston", "Thurston") ~ "Thurston",
    TRUE ~ FAC_COUNTY
  )) %>%
  filter(!FAC_COUNTY %in% c("NULL", "United States", "Unknown")) %>%
  filter(!is.na(FAC_COUNTY)) %>%
  mutate(GROUP = case_when(
    FAC_COUNTY %in% omaha_counties     ~ "Omaha",
    FAC_COUNTY %in% lincoln_counties   ~ "Lincoln",
    FAC_COUNTY %in% tricities_counties ~ "Tri-Cities",
    FAC_COUNTY %in% west_counties      ~ "West",
    FAC_COUNTY %in% southeast_counties ~ "Southeast",
    FAC_COUNTY %in% northeast_counties ~ "Northeast",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(GROUP)) %>%
  group_by(GROUP, DOD_STRING, cvd_type) %>%
  summarise(deaths = sum(deaths), .groups = "drop") %>%
  arrange(GROUP, DOD_STRING)

# ---- Group daily death rate table (rates only; combining already done) ----
group_daily_deaths <- CVD_data %>%
  filter(cvd_type == "total_cvd_deaths") %>%
  group_by(GROUP) %>%
  summarise(
    total_cvd_deaths = sum(deaths),
    n_days           = n_days_total,
    daily_death_rate = total_cvd_deaths / n_days_total,
    .groups = "drop"
  ) %>%
  arrange(desc(daily_death_rate))

group_daily_deaths