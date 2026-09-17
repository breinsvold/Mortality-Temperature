# ============================================================
#  UTCI -> Nebraska region daily temperature attachment
#
#  Region geometry comes from CVD_regions.geojson (exported from
#  ArcGIS, WGS84 / EPSG:4326).  No sf/terra/GDAL needed - jsonlite
#  reads the polygons and point-in-polygon is done in plain R.
#  Only ncdf4 is compiled.
#
#  ncdf4 note: the custom build in ~/R/x86_64-module-4.5/library has
#  netcdf + hdf5 baked into its RPATH, so NO Lmod modules and NO
#  LD_LIBRARY_PATH are needed (the RStudio session has neither).
#
#  Requires CVD_data in the environment - run the CVD prep script first.
# ============================================================

.libPaths(c(path.expand("~/R/x86_64-module-4.5/library"), .libPaths()))

library(dplyr); library(tidyr); library(readr)
library(stringr); library(purrr); library(lubridate)
library(jsonlite); library(ncdf4)

# --- Config ---------------------------------------------------
geojson_path <- "../CVD_regions.geojson"   # adjust if stored elsewhere
utci_root    <- ".."                       # year folders live in the PARENT dir
year_dirs    <- paste0("utci_nebraska_", 2005:2025)
region_field <- "Region"                   # region-name property in the GeoJSON

# Both daily statistics are carried through to the output.
utci_vars <- c(utci_min = "utci_daily_min",
               utci_max = "utci_daily_max")

kelvin_offset <- 273.15   # source data is Kelvin (verified: range ~259-294 K)

# Cities take the single grid cell nearest the region centroid (urban
# exposure).  Rural regions average every grid cell inside the outline.
city_groups  <- c("Omaha", "Lincoln", "Tri-Cities")
rural_groups <- c("West", "Southeast", "Northeast")

stopifnot(
  "CVD_data not found - run the CVD prep script first." = exists("CVD_data"),
  "GeoJSON not found at geojson_path."                  = file.exists(geojson_path)
)

# --- Read region polygons from GeoJSON ------------------------
# Each region becomes a list of rings; a ring is an n x 2 matrix of
# lon/lat.  Ring 1 is the exterior, any further rings are holes (the
# Southeast polygon has Lancaster/Lincoln punched out this way).
read_regions <- function(path, field) {
  gj <- fromJSON(path, simplifyVector = FALSE)
  if (!identical(gj$type, "FeatureCollection"))
    stop("Expected a GeoJSON FeatureCollection, got: ", gj$type)
  
  out <- map(gj$features, function(f) {
    nm <- f$properties[[field]]
    if (is.null(nm) || is.na(nm) || !nzchar(nm))
      stop("A feature has an empty '", field, "' property. ",
           "Populate it in ArcGIS and re-export.")
    g <- f$geometry
    raw_rings <- switch(g$type,
                        "Polygon"      = g$coordinates,
                        "MultiPolygon" = unlist(g$coordinates, recursive = FALSE),
                        stop("Unsupported geometry type: ", g$type))
    rings <- map(raw_rings, function(r) matrix(unlist(r), ncol = 2, byrow = TRUE))
    list(name = as.character(nm), rings = rings)
  })
  names(out) <- map_chr(out, "name")
  out
}

regions <- read_regions(geojson_path, region_field)
message("Regions read: ", paste(names(regions), collapse = ", "))

expected  <- c(city_groups, rural_groups)
missing_r <- setdiff(expected, names(regions))
extra_r   <- setdiff(names(regions), expected)
if (length(missing_r))
  stop("GeoJSON is missing region(s): ", paste(missing_r, collapse = ", "),
       ". Check the spelling of the '", region_field, "' field.")
if (length(extra_r))
  warning("GeoJSON has unexpected region(s), ignored: ",
          paste(extra_r, collapse = ", "))

# Coordinates must be lon/lat degrees, not a projected system.
all_xy <- do.call(rbind, unlist(map(regions, "rings"), recursive = FALSE))
if (!(all(between(all_xy[, 1], -180, 180)) && all(between(all_xy[, 2], -90, 90))))
  stop("GeoJSON coordinates are not lon/lat degrees - re-export with ",
       "'Project to WGS_1984' checked.")

# --- Geometry helpers -----------------------------------------
# Even-odd ray casting across every ring at once.  A point inside a hole
# crosses that ring too, flipping parity back to outside, so holes need
# no special handling.
points_in_region <- function(px, py, rings) {
  inside <- logical(length(px))
  for (p in seq_along(px)) {
    crossings <- 0L
    for (r in rings) {
      x <- r[, 1]; y <- r[, 2]; n <- nrow(r)
      x1 <- x[-n]; y1 <- y[-n]; x2 <- x[-1]; y2 <- y[-1]
      straddles <- (y1 > py[p]) != (y2 > py[p])
      if (any(straddles)) {
        xi <- (x2[straddles] - x1[straddles]) *
          (py[p] - y1[straddles]) /
          (y2[straddles] - y1[straddles]) + x1[straddles]
        crossings <- crossings + sum(px[p] < xi)
      }
    }
    inside[p] <- (crossings %% 2L) == 1L
  }
  inside
}

# Area-weighted centroid of a region's exterior ring.
ring_centroid <- function(r) {
  x <- r[, 1]; y <- r[, 2]; n <- nrow(r)
  x1 <- x[-1]; y1 <- y[-1]; x0 <- x[-n]; y0 <- y[-n]
  cross <- x0 * y1 - x1 * y0
  A <- sum(cross) / 2
  c(lon = sum((x0 + x1) * cross) / (6 * A),
    lat = sum((y0 + y1) * cross) / (6 * A))
}

# --- List all daily NetCDF files, sorted by date --------------
# Filenames look like:
#   ECMWF_utci_daily_stats_20050101_v1.1_con.area-subset.<...>.nc
# Anchor the date to the _YYYYMMDD_ slot rather than grabbing the first
# 8-digit run, since the trailing bbox is also full of digits.
date_from_path <- function(path) {
  ymd(str_extract(basename(path), "(?<=_)\\d{8}(?=_)"), quiet = TRUE)
}

all_nc <- list.files(file.path(utci_root, year_dirs),
                     pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
if (!length(all_nc)) stop("No .nc files found under ", utci_root)

files <- tibble(path = all_nc) %>% mutate(date = date_from_path(path))
n_bad <- sum(is.na(files$date))
if (n_bad) warning(n_bad, " file(s) had no parseable date and were skipped.")

files <- files %>% filter(!is.na(date)) %>% arrange(date)
dup <- files %>% count(date) %>% filter(n > 1)
if (nrow(dup))
  warning(nrow(dup), " date(s) appear in more than one file, e.g. ",
          paste(format(head(dup$date, 3)), collapse = ", "))

message("NetCDF files to process: ", nrow(files),
        "  (", format(min(files$date)), " to ", format(max(files$date)), ")")

# --- Grid geometry: read once from the first file -------------
# The lon/lat grid is identical across every file, so the nearest-cell
# lookups and the point-in-polygon masks are computed a single time.
grid_ref <- local({
  nc <- nc_open(files$path[1]); on.exit(nc_close(nc), add = TRUE)
  lon_name <- intersect(c("longitude", "lon"), names(nc$dim))[1]
  lat_name <- intersect(c("latitude",  "lat"), names(nc$dim))[1]
  if (is.na(lon_name) || is.na(lat_name))
    stop("Could not find lon/lat dimensions. Found: ",
         paste(names(nc$dim), collapse = ", "))
  for (v in utci_vars)
    if (!v %in% names(nc$var))
      stop("Variable '", v, "' not in file. Available: ",
           paste(names(nc$var), collapse = ", "))
  lon <- as.numeric(ncvar_get(nc, lon_name))
  lat <- as.numeric(ncvar_get(nc, lat_name))
  list(lon = ifelse(lon > 180, lon - 360, lon), lat = lat)
})

lon <- grid_ref$lon; lat <- grid_ref$lat
n_lon <- length(lon); n_lat <- length(lat)
message("Grid: ", n_lon, " lon x ", n_lat, " lat  |  lon ",
        sprintf("%.2f..%.2f", min(lon), max(lon)), "  lat ",
        sprintf("%.2f..%.2f", min(lat), max(lat)))

# Flattened cell coordinates, matching as.vector() of a [lon, lat] matrix.
cells <- expand_grid(j = seq_len(n_lat), i = seq_len(n_lon)) %>%
  arrange(j, i) %>%
  mutate(lon = lon[i], lat = lat[j], flat = i + (j - 1L) * n_lon)

# Warn if any region pokes outside the UTCI grid footprint.
for (nm in expected) {
  ext <- do.call(rbind, regions[[nm]]$rings)
  if (min(ext[, 1]) < min(lon) || max(ext[, 1]) > max(lon) ||
      min(ext[, 2]) < min(lat) || max(ext[, 2]) > max(lat))
    warning("Region '", nm, "' extends beyond the UTCI grid footprint; ",
            "edge cells may be missing.")
}

# Cities: the single nearest cell to the region centroid.
city_idx <- tibble(GROUP = city_groups) %>%
  mutate(cell = map_int(GROUP, function(nm) {
    cen <- ring_centroid(regions[[nm]]$rings[[1]])
    which.min((cells$lon - cen["lon"])^2 + (cells$lat - cen["lat"])^2)
  }))

# Rural: every cell whose centre falls inside the region outline.
message("Building rural cell masks...")
rural_idx <- tibble(GROUP = rural_groups) %>%
  mutate(cell_ids = map(GROUP, function(nm)
    which(points_in_region(cells$lon, cells$lat, regions[[nm]]$rings))),
    n_cells  = map_int(cell_ids, length))

if (any(rural_idx$n_cells == 0))
  stop("Rural region(s) contain no grid cells: ",
       paste(rural_idx$GROUP[rural_idx$n_cells == 0], collapse = ", "))

message("Rural cells per region: ",
        paste(rural_idx$GROUP, rural_idx$n_cells, sep = "=", collapse = ", "))
message("City cells: ",
        paste(city_idx$GROUP,
              sprintf("(%.2f, %.2f)", cells$lon[city_idx$cell],
                      cells$lat[city_idx$cell]),
              sep = " ", collapse = ", "))

# --- Read one daily file -> flat min/max vectors --------------
read_utci_day <- function(path) {
  nc <- nc_open(path); on.exit(nc_close(nc), add = TRUE)
  out <- lapply(utci_vars, function(v) {
    vals <- ncvar_get(nc, v)
    if (length(dim(vals)) == 3) vals <- vals[, , 1]   # drop singleton time axis
    if (!identical(dim(vals), c(n_lon, n_lat)))
      stop("Unexpected grid in ", basename(path), ": ",
           paste(dim(vals), collapse = "x"))
    as.vector(vals) - kelvin_offset                   # Kelvin -> Celsius
  })
  names(out) <- names(utci_vars)
  out
}

# --- Extract one value per region per statistic for a day -----
extract_day <- function(path, date) {
  day <- read_utci_day(path)
  bind_rows(
    tibble(GROUP    = city_idx$GROUP,
           date     = date,
           utci_min = day$utci_min[city_idx$cell],
           utci_max = day$utci_max[city_idx$cell]),
    tibble(GROUP    = rural_idx$GROUP,
           date     = date,
           utci_min = map_dbl(rural_idx$cell_ids,
                              ~ mean(day$utci_min[.x], na.rm = TRUE)),
           utci_max = map_dbl(rural_idx$cell_ids,
                              ~ mean(day$utci_max[.x], na.rm = TRUE)))
  )
}

# --- Run over every day -> daily UTCI per region --------------
message("Processing ", nrow(files), " files...")
tick <- max(1L, nrow(files) %/% 20L)
results <- vector("list", nrow(files))
for (k in seq_len(nrow(files))) {
  results[[k]] <- extract_day(files$path[k], files$date[k])
  if (k %% tick == 0 || k == nrow(files))
    message("  ", k, "/", nrow(files), " (", round(100 * k / nrow(files)), "%)")
}
utci_daily <- bind_rows(results) %>% arrange(GROUP, date)

message("Daily UTCI rows: ", nrow(utci_daily))
print(utci_daily %>% group_by(GROUP) %>%
        summarise(n          = n(),
                  min_C      = round(min(utci_min,  na.rm = TRUE), 1),
                  max_C      = round(max(utci_max,  na.rm = TRUE), 1),
                  mean_max_C = round(mean(utci_max, na.rm = TRUE), 1),
                  .groups = "drop"))

write_csv(utci_daily, "utci_daily_by_region.csv")

# --- Build the DLNM analysis frame ----------------------------
# DLNM needs one row per GROUP per day, in date order, with NO GAPS:
# the crossbasis lag structure assumes consecutive rows are consecutive
# days.  CVD_data is long (one row per cvd_type) and only contains
# region-days on which a death was recorded, so it is pivoted wide and
# the date grid is completed.
cvd_dated <- CVD_data %>% mutate(date = mdy(DOD_STRING, quiet = TRUE))

n_baddate <- sum(is.na(cvd_dated$date))
if (n_baddate)
  warning(n_baddate, " CVD row(s) have an unparseable DOD_STRING; ",
          "check its format (expected M/D/YYYY).")

cvd_wide <- cvd_dated %>%
  filter(!is.na(date)) %>%
  select(GROUP, date, cvd_type, deaths) %>%
  pivot_wider(names_from = cvd_type, values_from = deaths,
              values_fn = sum, values_fill = 0)

death_cols <- setdiff(names(cvd_wide), c("GROUP", "date"))
message("Death columns: ", paste(death_cols, collapse = ", "))

# Complete the daily grid over the period covered by BOTH datasets.
start_date <- max(min(cvd_wide$date), min(utci_daily$date))
end_date   <- min(max(cvd_wide$date), max(utci_daily$date))
message("Analysis period: ", format(start_date), " to ", format(end_date))

full_grid <- expand_grid(
  GROUP = sort(unique(utci_daily$GROUP)),
  date  = seq(start_date, end_date, by = "day"))

dlnm_data <- full_grid %>%
  left_join(utci_daily, by = c("GROUP", "date")) %>%
  left_join(cvd_wide,   by = c("GROUP", "date")) %>%
  mutate(
    # A region-day absent from CVD_data means no deaths were recorded,
    # i.e. a true zero - not a missing observation.
    across(all_of(death_cols), ~ replace_na(.x, 0L)),
    year  = year(date),
    month = month(date),
    doy   = yday(date),
    dow   = wday(date, label = TRUE, week_start = 1),
    time  = as.integer(date - start_date) + 1L   # linear time index for the spline
  ) %>%
  arrange(GROUP, date) %>%
  select(GROUP, date, year, month, doy, dow, time,
         utci_min, utci_max, all_of(death_cols))

# --- Report data quality before modelling ---------------------
n_filled <- nrow(dlnm_data) - nrow(distinct(cvd_wide, GROUP, date))
message("Rows in analysis frame: ", nrow(dlnm_data),
        "  (", n_filled, " region-days filled as zero deaths)")

gaps <- dlnm_data %>%
  group_by(GROUP) %>%
  summarise(n_days       = n(),
            missing_utci = sum(is.na(utci_max)),
            .groups = "drop")
print(gaps)

if (any(gaps$missing_utci > 0))
  warning("Some region-days have no UTCI value. DLNM will drop these rows, ",
          "which BREAKS the lag structure. Inspect before modelling: ",
          "dlnm_data %>% filter(is.na(utci_max))")

# Each region must be a contiguous daily series for lags to be valid.
stopifnot(
  "Duplicate GROUP-date rows in dlnm_data." =
    !any(duplicated(dlnm_data[c("GROUP", "date")])),
  "Non-contiguous dates within a region." =
    dlnm_data %>% group_by(GROUP) %>%
    summarise(ok = all(diff(date) == 1), .groups = "drop") %>%
    pull(ok) %>% all()
)

write_csv(dlnm_data, "CVD-UTCI_dlnm_data.csv")
message("Wrote CVD-UTCI_dlnm_data.csv")

dlnm_data

# ============================================================
#  Next step: DLNM
#  Fit ONE REGION AT A TIME - a crossbasis computes lags down the
#  rows it is given, so a pooled frame would lag one region's
#  exposure into the next region's deaths.
#
#  library(dlnm); library(splines)
#  d  <- dlnm_data %>% filter(GROUP == "Omaha")
#  cb <- crossbasis(d$utci_max, lag = 21,
#                   argvar = list(fun = "ns", df = 4),
#                   arglag = list(fun = "ns", df = 4))
#  m  <- glm(total_cvd_deaths ~ cb + ns(time, df = 7 * length(unique(d$year))) +
#              dow, family = quasipoisson(), data = d)
#  pred <- crosspred(cb, m, cen = median(d$utci_max, na.rm = TRUE))
#  plot(pred, "overall")
#
#  Then repeat per region and pool with mvmeta if wanted.
# ============================================================