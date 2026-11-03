# =====================================================================
# Stage 3 -- attributable CVD deaths from UTCI
# Nebraska, 2005-2025
#
# Uses Gasparrini's attrdl() (attributable risk from DLNMs) with
# Monte Carlo simulation for 95% CIs (~1000 sims).
#
# Design (locked with collaborator):
#   - OVERALL attribution: single POOLED curve applied to every region's
#     days (pooled coef/vcov from stage 2), summed over all region-days.
#   - cold/heat split: ABSOLUTE UTCI cutpoints (Broede scale):
#       cold stress  = UTCI < 9
#       heat stress  = UTCI > 26
#       (no-stress 9-26 contributes ~baseline; excluded from cold/heat)
#   - by-category: full Broede 10-category scale, all bands except
#     no-stress.
#   - urban/rural: pooled-urban and pooled-rural curves (from stage 2),
#     each applied to its own group's region-days.
#   - per cause: same three breakdowns for heart disease, stroke,
#     hypertension, other vascular (NOT total).
#
# Requires stage 1 (per-region cb + coefs) AND stage 2 (pooled/urban/
# rural meta fits) to have been run for each OUTCOME at this EXPOSURE.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dlnm)
library(splines)
library(mixmeta)
library(dplyr)

if (!exists("EXPOSURE")) EXPOSURE <- "utci_max"
INFILE   <- "CVD-UTCI_dlnm_data.csv"
OUTDIR   <- "stage3_out"; dir.create(OUTDIR, showWarnings = FALSE)

REGIONS <- c("Omaha","Lincoln","Tri-Cities","West","Southeast","Northeast")
URBAN   <- c("Omaha","Lincoln","Tri-Cities")
RURAL   <- c("West","Southeast","Northeast")

OVERALL_CAUSE <- "total_cvd_deaths"
CAUSES <- c("heart_disease_deaths","stroke_deaths",
            "hypertension_deaths","other_vascular_deaths")
CAUSE_LAB <- c(total_cvd_deaths="Total CVD",
               heart_disease_deaths="Heart disease", stroke_deaths="Stroke",
               hypertension_deaths="Hypertension",
               other_vascular_deaths="Other vascular")

LAG <- 21
NSIM <- 1000
set.seed(1)

# UTCI Broede categories (lower bound, label). Ordered cold -> hot.
UTCI_CATS <- data.frame(
  lo  = c(-Inf, -40, -27, -13,   0,   9,  26,  32,  38,  46),
  hi  = c(-40, -27, -13,   0,   9,  26,  32,  38,  46, Inf),
  lab = c("Extreme cold","Very strong cold","Strong cold","Moderate cold",
          "Slight cold","No stress","Moderate heat","Strong heat",
          "Very strong heat","Extreme heat"),
  stringsAsFactors = FALSE)

## --- attrdl(): source Gasparrini's attributable-risk function -------
# Not exported by dlnm; ships with the 2014 BMC MRM paper. Try a local
# copy first, then GitHub. If neither is reachable, stop with guidance.
if (!exists("attrdl")) {
  local_attrdl <- "attrdl.R"
  if (file.exists(local_attrdl)) {
    source(local_attrdl)
  } else {
    ok <- tryCatch({
      source("https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/attrdl.R")
      TRUE
    }, error = function(e) FALSE)
    if (!ok || !exists("attrdl"))
      stop("attrdl() not found. Download attrdl.R from\n",
           "  github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata\n",
           "and place it in the working directory.")
  }
}

## --- load data + stage 1/2 objects ----------------------------------
dat <- read.csv(INFILE, stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date); dat$dow <- factor(dat$dow)

# basis spec (shared across causes; UTCI range is cause-independent)
bspec <- readRDS(file.path("stage1_out",
                           paste0("basis_", OVERALL_CAUSE, "_", EXPOSURE, ".rds")))

# rebuild the per-region cross-basis EXACTLY as stage 1 did, so attrdl
# can map coefficients onto each region's own exposure series.
build_cb <- function(x) {
  crossbasis(x, lag = LAG,
             argvar = list(fun="ns", knots=bspec$knots, Bound=bspec$bound),
             arglag = list(fun="ns", knots=logknots(LAG,3)))
}

# helper: return the FULL simulation vector of attributable numbers for
# a coef/vcov curve applied to a region's days, restricted to an EXPOSURE
# RANGE (attrdl has no `subset`; it restricts by exposure value via
# `range`). Retaining draws lets us sum across regions BEFORE taking
# percentiles -> exact joint CI.
#   rng = c(lo, hi) exposure interval; cases counted as days with x in rng
attr_sim <- function(coef, vcov, cb, x, cases, cen, rng, simseed) {
  set.seed(simseed)  # SAME seed across regions -> shared coef draws ->
  # correct joint (perfectly-dependent) pooled CI
  an  <- attrdl(x, cb, cases, coef = coef, vcov = vcov, type = "an",
                dir = "forw", cen = cen, range = rng, sim = FALSE)
  sims <- attrdl(x, cb, cases, coef = coef, vcov = vcov, type = "an",
                 dir = "forw", cen = cen, range = rng, sim = TRUE, nsim = NSIM)
  inrng <- x >= rng[1] & x <= rng[2]
  list(an = as.numeric(an), sims = as.numeric(sims),
       cases = sum(cases[inrng], na.rm = TRUE))
}

# Summarise a list of per-region attr_sim results into one row:
# sum point ANs, sum sim vectors elementwise, percentile once.
combine_regions <- function(parts) {
  an   <- sum(vapply(parts, `[[`, numeric(1), "an"))
  sims <- Reduce(`+`, lapply(parts, `[[`, "sims"))   # elementwise sum
  cases<- sum(vapply(parts, `[[`, numeric(1), "cases"))
  c(AN    = an,
    AN_lo = as.numeric(quantile(sims, .025, na.rm = TRUE)),
    AN_hi = as.numeric(quantile(sims, .975, na.rm = TRUE)),
    AF    = 100 * an / cases,
    AF_lo = 100 * as.numeric(quantile(sims, .025, na.rm=TRUE)) / cases,
    AF_hi = 100 * as.numeric(quantile(sims, .975, na.rm=TRUE)) / cases)
}

## --- core: attribute one cause with a given curve set ---------------
# curve_mode: "pooled" (one shared curve for all regions) or "urbanrural"
#   (urban curve for urban regions, rural for rural). Because a shared
#   curve is used within each pool, summing per-region sim draws (same
#   seed) yields the exact joint CI for that pool.
# Returns a data.frame: one row per breakdown (split/category), already
# combined across the relevant regions, with exact CIs.
attribute_cause <- function(cause, curve_mode = "pooled") {
  
  tag <- paste0(cause, "_", EXPOSURE)
  meta <- readRDS(file.path("stage2_out", paste0("meta_fits_", tag, ".rds")))
  pooled_coef <- coef(meta$mv);  pooled_vcov <- vcov(meta$mv)
  urb_coef <- coef(meta$mv_urb); urb_vcov <- vcov(meta$mv_urb)
  rur_coef <- coef(meta$mv_rur); rur_vcov <- vcov(meta$mv_rur)
  cen <- meta$pool_mmt
  
  # Breakdowns as OPEN specifications (lo, hi) on the UTCI scale; the
  # open ends (-Inf / Inf) get clamped below to each POOL's own observed
  # UTCI range, so attribution never runs past the data that curve saw.
  # Cold = below no-stress (9); Heat = above no-stress (26); plus each
  # non-no-stress Broede band.
  breakdowns <- c(
    list(`Cold (<9)`  = c(-Inf, 9),
         `Heat (>26)` = c(26, Inf)),
    setNames(
      lapply(which(UTCI_CATS$lab != "No stress"),
             function(k) c(UTCI_CATS$lo[k], UTCI_CATS$hi[k])),
      UTCI_CATS$lab[UTCI_CATS$lab != "No stress"])
  )
  
  # for pooled mode: combine over all regions (one pool).
  # for urbanrural mode: combine within urban and within rural separately.
  if (curve_mode == "pooled") groups <- list(All = REGIONS)
  else groups <- list(Urban = URBAN, Rural = RURAL)
  
  out <- list()
  for (gname in names(groups)) {
    regs <- groups[[gname]]
    
    # OBSERVED UTCI range for THIS pool (full range, all days in these
    # regions). Open-ended breakdown bounds are clamped to this so no
    # extrapolated exposure contributes to the attribution.
    pool_x   <- dat[[EXPOSURE]][dat$GROUP %in% regs]
    pool_lo  <- min(pool_x, na.rm = TRUE)
    pool_hi  <- max(pool_x, na.rm = TRUE)
    
    for (bn in names(breakdowns)) {
      raw <- breakdowns[[bn]]
      # clamp to the pool's observed range
      rng <- c(max(raw[1], pool_lo), min(raw[2], pool_hi))
      # skip degenerate ranges (category lies entirely outside the data)
      if (rng[1] >= rng[2]) next
      parts <- list()
      for (g in regs) {
        d <- dat[dat$GROUP == g, ]; d <- d[order(d$date), ]
        x <- d[[EXPOSURE]]; cases <- d[[cause]]
        cb <- build_cb(x)
        if (curve_mode == "pooled") { cf <- pooled_coef; vc <- pooled_vcov }
        else if (g %in% URBAN)      { cf <- urb_coef;   vc <- urb_vcov }
        else                        { cf <- rur_coef;   vc <- rur_vcov }
        # skip if no days fall in this range for this region
        if (sum(x >= rng[1] & x <= rng[2], na.rm = TRUE) == 0) next
        # seed keyed to cause+breakdown (NOT region) so all regions in a
        # pool share the same coefficient draws -> exact joint CI.
        sd <- sum(utf8ToInt(paste(cause, bn))) + 1000
        parts[[g]] <- attr_sim(cf, vc, cb, x, cases, cen, rng, sd)
      }
      if (length(parts) == 0) next
      cr <- combine_regions(parts)
      out[[paste(gname, bn)]] <- data.frame(
        cause = CAUSE_LAB[cause], curve = curve_mode, group = gname,
        split = bn, t(cr), check.names = FALSE)
    }
  }
  do.call(rbind, out)
}

## --- run: overall (total CVD) and each cause ------------------------
message("Attributing OVERALL (", OVERALL_CAUSE, ") ...")
res_overall_pooled <- attribute_cause(OVERALL_CAUSE, "pooled")
res_overall_ur     <- attribute_cause(OVERALL_CAUSE, "urbanrural")

res_causes_pooled <- do.call(rbind, lapply(CAUSES, function(c1) {
  message("Attributing ", c1, " (pooled) ..."); attribute_cause(c1, "pooled")
}))
res_causes_ur <- do.call(rbind, lapply(CAUSES, function(c1) {
  message("Attributing ", c1, " (urban/rural) ..."); attribute_cause(c1, "urbanrural")
}))

## --- results already carry EXACT CIs (draws summed before percentile)
# attribute_cause() returns pre-combined rows with columns:
#   cause, curve, group, split, AN, AN_lo, AN_hi, AF, AF_lo, AF_hi
# pooled mode -> group == "All"; urbanrural mode -> group in Urban/Rural.

# save raw combined results
write.csv(res_overall_pooled, file.path(OUTDIR, paste0("overall_pooled_", EXPOSURE, ".csv")), row.names=FALSE)
write.csv(res_overall_ur,     file.path(OUTDIR, paste0("overall_urbanrural_", EXPOSURE, ".csv")), row.names=FALSE)
write.csv(res_causes_pooled,  file.path(OUTDIR, paste0("causes_pooled_", EXPOSURE, ".csv")), row.names=FALSE)
write.csv(res_causes_ur,      file.path(OUTDIR, paste0("causes_urbanrural_", EXPOSURE, ".csv")), row.names=FALSE)

cat("\nDone. CSVs in ", normalizePath(OUTDIR), "\n", sep="")
cat("  overall_pooled_", EXPOSURE, ".csv       (Total CVD, pooled curve)\n", sep="")
cat("  overall_urbanrural_", EXPOSURE, ".csv   (Total CVD, urban/rural curves)\n", sep="")
cat("  causes_pooled_", EXPOSURE, ".csv        (each cause, pooled curve)\n", sep="")
cat("  causes_urbanrural_", EXPOSURE, ".csv    (each cause, urban/rural curves)\n", sep="")
cat("\nCIs are exact: within each pool all regions share one curve, so the\n")
cat("simulation draws are summed across regions BEFORE taking the 2.5/97.5\n")
cat("percentiles. Intervals may be asymmetric and may cross zero -- that is\n")
cat("the honest result for a stratum whose burden is not distinguishable\n")
cat("from zero.\n")