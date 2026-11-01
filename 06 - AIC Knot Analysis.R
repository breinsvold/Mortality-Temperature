# =====================================================================
# Stage 1b -- QAIC selection of cross-basis knots / df
# Nebraska CVD vs UTCI
#
# Purpose: justify the argvar (exposure-response) and arglag (lag-
# response) specifications with QAIC instead of asserting them.
#
# QAIC = -2 * logLik_Poisson / phi_hat + 2*k
#   phi_hat = overdispersion, estimated ONCE from a reference (most
#   complex) model and held FIXED across all candidates so the criterion
#   is comparable. (Standard Gasparrini/Peng approach.)
#
# IMPORTANT for stage 2: the six regions must ultimately share a COMMON
# basis to be poolable. So this script reports:
#   (A) per-region QAIC-optimal specs   -> diagnostic / sensitivity
#   (B) QAIC SUMMED across regions for each spec, using a COMMON basis
#       -> this is what you use to pick the single shared specification.
# Pick the common spec by (B); use (A) only to see how much regions
# disagree about their ideal flexibility.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dlnm)
library(splines)

INFILE   <- "CVD-UTCI_dlnm_data.csv"
OUTDIR   <- "stage1_qaic"
dir.create(OUTDIR, showWarnings = FALSE)

EXPOSURE <- "utci_max"
OUTCOME  <- "total_cvd_deaths"
LAG      <- 21
TAG      <- paste0(OUTCOME, "_", EXPOSURE)

REGIONS  <- c("Omaha", "Lincoln", "Tri-Cities",
              "West", "Southeast", "Northeast")

## --- load ------------------------------------------------------------
dat <- read.csv(INFILE, stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date)
dat$dow  <- factor(dat$dow)

## --- candidate grid --------------------------------------------------
# argvar (exposure-response): natural spline knot placements to test.
# Expressed as percentile sets of the POOLED distribution so the tested
# bases are common across regions.
argvar_cands <- list(
  "ns_50"            = c(.50),
  "ns_10_90"         = c(.10, .90),
  "ns_25_50_75"      = c(.25, .50, .75),
  "ns_10_75_90"      = c(.10, .75, .90),   # the current default
  "ns_5_35_65_95"    = c(.05, .35, .65, .95),
  "ns_10_50_75_90"   = c(.10, .50, .75, .90)
)

# arglag (lag-response): number of internal log-spaced knots to test.
arglag_cands <- list(
  "lag_2" = 2,
  "lag_3" = 3,   # current default
  "lag_4" = 4
)

# common knots/bounds from pooled distribution (shared across regions)
common_bound <- range(dat[[EXPOSURE]], na.rm = TRUE)
pooled_x     <- dat[[EXPOSURE]]

## --- helper: fit one region with a given spec, return Poisson dev ----
# returns components needed for a FIXED-phi QAIC:
#   dev_pois = Poisson deviance (for logLik proxy), k = # params, n, df.res
fit_one <- function(d, var_knots_pct, lag_nknots) {
  
  x <- d[[EXPOSURE]]; y <- d[[OUTCOME]]
  
  argvar <- list(fun = "ns",
                 knots = as.numeric(quantile(pooled_x, var_knots_pct,
                                             na.rm = TRUE)),
                 Bound = common_bound)
  arglag <- list(fun = "ns", knots = logknots(LAG, lag_nknots))
  cb <- crossbasis(x, lag = LAG, argvar = argvar, arglag = arglag)
  
  n_years  <- length(unique(d$year))
  trend_df <- max(2, round(n_years / 10))
  seas  <- ns(d$doy,  df = 4)
  trend <- ns(d$time, df = trend_df)
  yearf <- factor(d$year)
  dow   <- d$dow
  
  m <- tryCatch(
    glm(y ~ cb + seas * yearf + trend + dow, family = quasipoisson()),
    error = function(e)
      glm(y ~ cb + seas + trend + dow, family = quasipoisson())
  )
  
  # Poisson deviance from the SAME fitted means (quasi and Poisson share
  # the point estimates; only the dispersion differs). Use the model
  # deviance, which for the (quasi)poisson family is the Poisson deviance.
  list(dev = m$deviance,
       k   = length(coef(m)),          # params incl. cross-basis
       n   = nrow(d),
       dfres = m$df.residual,
       pearson = sum(residuals(m, type = "pearson")^2, na.rm = TRUE),
       fitted_ok = m$converged)
}

## --- reference phi (overdispersion), estimated ONCE ------------------
# Use the most complex spec (finest argvar + most lag knots) on the
# POOLED data as the reference for phi_hat, then hold phi fixed for all
# QAIC calcs. This is the comparability requirement.
ref_var <- argvar_cands[["ns_10_50_75_90"]]
ref_lag <- arglag_cands[["lag_4"]]

ref_fits <- lapply(REGIONS, function(g)
  fit_one(dat[dat$GROUP == g, ], ref_var, ref_lag))
phi_hat <- sum(sapply(ref_fits, `[[`, "pearson")) /
  sum(sapply(ref_fits, `[[`, "dfres"))
cat(sprintf("Fixed overdispersion phi_hat = %.3f  (from reference spec)\n\n",
            phi_hat))

## --- QAIC for a fit, given fixed phi --------------------------------
# QAIC = dev/phi + 2*k   (dev is the Poisson deviance; the -2logLik and
# constants cancel in model COMPARISON, so deviance/phi + 2k ranks
# identically to the full QAIC. We report this comparable form.)
qaic <- function(fit, phi) fit$dev / phi + 2 * fit$k

## --- (A) per-region optimal spec ------------------------------------
per_region <- list()
for (g in REGIONS) {
  d <- dat[dat$GROUP == g, ]
  rows <- list()
  for (av in names(argvar_cands)) for (al in names(arglag_cands)) {
    f <- fit_one(d, argvar_cands[[av]], arglag_cands[[al]])
    rows[[paste(av, al)]] <- data.frame(
      region = g, argvar = av, arglag = al,
      k = f$k, qaic = round(qaic(f, phi_hat), 1),
      converged = f$fitted_ok)
  }
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$qaic), ]
  tab$dQAIC <- round(tab$qaic - min(tab$qaic), 1)
  per_region[[g]] <- tab
}
per_region_df <- do.call(rbind, per_region)
write.csv(per_region_df,
          file.path(OUTDIR, paste0("qaic_per_region_", TAG, ".csv")),
          row.names = FALSE)

## --- (B) COMMON spec: QAIC summed across regions --------------------
# For each (argvar x arglag) combination, sum QAIC over the six regions.
# The minimum-sum spec is the single shared basis best supported overall.
combo_rows <- list()
for (av in names(argvar_cands)) for (al in names(arglag_cands)) {
  qsum <- 0; kk <- NA; ok <- TRUE
  for (g in REGIONS) {
    f <- fit_one(dat[dat$GROUP == g, ], argvar_cands[[av]],
                 arglag_cands[[al]])
    qsum <- qsum + qaic(f, phi_hat)
    kk <- f$k; ok <- ok && f$fitted_ok
  }
  combo_rows[[paste(av, al)]] <- data.frame(
    argvar = av, arglag = al, k_per_region = kk,
    qaic_sum = round(qsum, 1), all_converged = ok)
}
combo <- do.call(rbind, combo_rows)
combo <- combo[order(combo$qaic_sum), ]
combo$dQAIC_sum <- round(combo$qaic_sum - min(combo$qaic_sum), 1)
write.csv(combo, file.path(OUTDIR, paste0("qaic_common_", TAG, ".csv")),
          row.names = FALSE)

## --- report ----------------------------------------------------------
cat("===== (B) COMMON-basis QAIC, summed across regions =====\n")
cat("    (pick the shared stage-1 spec from the top row)\n\n")
print(combo, row.names = FALSE)

cat("\n===== (A) per-region best spec (top of each) =====\n")
for (g in REGIONS) {
  best <- per_region[[g]][1, ]
  cat(sprintf("  %-11s -> %s / %s   (QAIC %.1f)\n",
              g, best$argvar, best$arglag, best$qaic))
}

cat("\nInterpretation:\n")
cat(" - Use the TOP row of (B) as your common stage-1 basis.\n")
cat(" - dQAIC_sum < ~2 vs the default => specs are effectively tied;\n")
cat("   prefer the conventional 10/75/90 + lag_3 for comparability and\n")
cat("   better tail extrapolation (per France-study reasoning).\n")
cat(" - (A) shows whether regions disagree about ideal flexibility;\n")
cat("   wide disagreement is itself worth noting, but stage 2 still\n")
cat("   requires the single common spec from (B).\n")
cat("\nSaved to ", normalizePath(OUTDIR), "\n", sep = "")