# =====================================================================
# Stage 1 -- region-specific DLNM of CVD mortality vs UTCI heat stress
# Nebraska, 6 regions, 2005-2025
#
# Fits an independent DLNM per region, locates the minimum-mortality
# UTCI (MMT), re-centres there, and saves the reduced overall
# exposure-response (coefs + vcov) for the stage-2 meta-analysis.
#
# Design choices (locked with collaborator):
#   exposure   : utci_max primary, utci_min secondary  (EXPOSURE var)
#   lag        : 21 days
#   argvar     : natural spline, knots at 10/75/90 pctiles of POOLED
#                UTCI (common basis across regions -> poolable in stage 2)
#   arglag     : natural spline w/ intercept, 3 log-spaced knots
#   season/trend: ns(doy) x factor(year) interaction + ns(time) trend
#   family     : quasipoisson
#   centering  : region-specific MMT
# =====================================================================

## --- 0. setup --------------------------------------------------------
# custom ncdf4 lib path not needed here, but harmless to keep consistent
.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dlnm)
library(splines)
library(dplyr)

set.seed(1)

INFILE   <- "CVD-UTCI_dlnm_data.csv"
OUTDIR   <- "stage1_out"
dir.create(OUTDIR, showWarnings = FALSE)

# ---- exposure: honor a value set by a driver script, else default ---
# A purrr::walk driver can set EXPOSURE before source()-ing this file;
# if nothing is set, fall back to the primary analysis (utci_max).
if (!exists("EXPOSURE")) EXPOSURE <- "utci_max"
if (!exists("OUTCOME"))  OUTCOME  <- "total_cvd_deaths"
LAG      <- 21
TAG      <- paste0(OUTCOME, "_", EXPOSURE)   # names output files

REGIONS  <- c("Omaha", "Lincoln", "Tri-Cities",
              "West", "Southeast", "Northeast")

## --- 1. load ---------------------------------------------------------
dat <- read.csv(INFILE, stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date)
# dow shipped as an ordered factor in the CSV build; re-assert
dat$dow  <- factor(dat$dow)

# sanity: contiguous dates per region are guaranteed by the build script,
# but the cross-basis depends on it absolutely, so cheap re-check here.
chk <- dat %>% group_by(GROUP) %>%
  summarise(gaps = any(diff(sort(date)) != 1), .groups = "drop")
if (any(chk$gaps)) stop("Non-contiguous dates in: ",
                        paste(chk$GROUP[chk$gaps], collapse = ", "))

## --- 1b. COMMON exposure basis (computed once, from all regions) -----
# Knots and boundary are set from the POOLED UTCI distribution and reused
# for every region. This is what makes the region-specific reduced
# coefficients commensurable, so stage 2 can pool them. Do NOT switch
# these back to per-region quantiles -- that breaks the meta-analysis.
COMMON_KNOTS <- as.numeric(
  quantile(dat[[EXPOSURE]], c(.10, .75, .90), na.rm = TRUE))
COMMON_BOUND <- range(dat[[EXPOSURE]], na.rm = TRUE)

cat("Common argvar knots (", EXPOSURE, "): ",
    paste(round(COMMON_KNOTS, 2), collapse = ", "),
    "  | bounds: ", paste(round(COMMON_BOUND, 2), collapse = " to "),
    "\n", sep = "")

## --- 2. containers ---------------------------------------------------
coef_list <- vector("list", length(REGIONS)); names(coef_list) <- REGIONS
vcov_list <- vector("list", length(REGIONS)); names(vcov_list) <- REGIONS
mmt_vec   <- setNames(rep(NA_real_, length(REGIONS)), REGIONS)
diag_rows <- list()

## --- 3. loop ---------------------------------------------------------
for (g in REGIONS) {
  
  d <- dat[dat$GROUP == g, ]
  d <- d[order(d$date), ]
  x <- d[[EXPOSURE]]
  y <- d[[OUTCOME]]
  
  # ---- 3a. diagnostics (decide later how far to trust each region) ---
  diag_rows[[g]] <- data.frame(
    region      = g,
    n_days      = nrow(d),
    mean_deaths = round(mean(y), 3),
    pct_zero    = round(100 * mean(y == 0), 1),
    utci_sd     = round(sd(x), 2),
    utci_p01    = round(quantile(x, .01), 1),
    utci_p99    = round(quantile(x, .99), 1)
  )
  
  # ---- 3b. cross-basis ----------------------------------------------
  # COMMON basis across regions: same knots + boundary for every region
  # (computed once above from pooled UTCI). This is what lets stage 2
  # pool the reduced coefficients. argvar identical across regions;
  # only the data (x) differs.
  argvar <- list(fun   = "ns",
                 knots = COMMON_KNOTS,
                 Bound = COMMON_BOUND)
  arglag <- list(fun = "ns", knots = logknots(LAG, 3))  # 3 log-spaced
  
  cb <- crossbasis(x, lag = LAG, argvar = argvar, arglag = arglag)
  
  # ---- 3c. confounder terms -----------------------------------------
  # Seasonality: cyclic-ish flexible shape allowed to differ by year.
  #   ns(doy, 4 df) x factor(year)  -> relaxes constant-season assumption
  # Long-term trend: slow ns(time), ~1 df / decade.
  n_years <- length(unique(d$year))
  trend_df <- max(2, round(n_years / 10))   # ~1 per decade, floor 2
  
  # Confounder terms as matrices in THIS environment. We do NOT pack cb
  # or these splines into a data.frame: the crossbasis column names must
  # survive intact so crosspred() can match them to model coefficients.
  # glm() resolves bare object names (cb, seas, trend) from the calling
  # frame, which is the idiom used in the dlnm vignette.
  seas  <- ns(d$doy,  df = 4)      # seasonal shape
  trend <- ns(d$time, df = trend_df)  # slow long-term trend, ~1 df/decade
  yearf <- factor(d$year)
  dow   <- d$dow
  
  # formula: cb + (seasonal spline * year) + slow trend + dow
  form  <- y ~ cb + seas * yearf + trend + dow
  form2 <- y ~ cb + seas + trend + dow   # fallback: drop season x year
  
  m <- tryCatch(
    glm(form, family = quasipoisson()),
    error = function(e) {
      message("  [", g, "] full model failed (", conditionMessage(e),
              "); retrying without season x year interaction")
      glm(form2, family = quasipoisson())
    }
  )
  
  # ---- 3d. find MMT --------------------------------------------------
  # predict across the observed UTCI range, restricted to a sensible
  # central band so a noisy tail dip can't masquerade as the optimum.
  pctl   <- quantile(x, c(.01, .99), na.rm = TRUE)
  atseq  <- seq(pctl[1], pctl[2], length.out = 100)
  p0     <- crosspred(cb, m, at = atseq)
  # search MMT only within 1st-99th pctile (already the range), take the
  # UTCI minimising the overall cumulative RR
  mmt    <- atseq[which.min(p0$allfit)]
  mmt_vec[g] <- mmt
  
  # ---- 3e. reduce to overall cumulative ER, centred at MMT ----------
  red <- crossreduce(cb, m, type = "overall", cen = mmt)
  coef_list[[g]] <- coef(red)
  vcov_list[[g]] <- vcov(red)
  
  # ---- 3f. per-region prediction object for plotting/checking -------
  pred <- crosspred(cb, m, cen = mmt, by = 0.5)
  saveRDS(pred, file.path(OUTDIR, paste0("pred_", TAG, "_", g, ".rds")))
  
  message(sprintf("[%s] done | mean deaths=%.2f  %%zero=%.1f  MMT=%.1f",
                  g, mean(y), 100*mean(y==0), mmt))
}

## --- 4. save stage-2 inputs -----------------------------------------
# coefs as a matrix (rows = regions), vcov as a list, aligned by name.
coef_mat <- do.call(rbind, coef_list)
rownames(coef_mat) <- REGIONS

saveRDS(coef_mat,  file.path(OUTDIR, paste0("coef_", TAG, ".rds")))
saveRDS(vcov_list, file.path(OUTDIR, paste0("vcov_", TAG, ".rds")))
saveRDS(mmt_vec,   file.path(OUTDIR, paste0("mmt_",  TAG, ".rds")))

# save the common basis spec so stage 2 rebuilds an IDENTICAL basis
# (knots + bounds + lag) without hand-copying any numbers.
basis_spec <- list(knots = COMMON_KNOTS, bound = COMMON_BOUND,
                   lag = LAG, exposure = EXPOSURE, outcome = OUTCOME)
saveRDS(basis_spec, file.path(OUTDIR, paste0("basis_", TAG, ".rds")))

diag <- do.call(rbind, diag_rows)
write.csv(diag, file.path(OUTDIR, paste0("diagnostics_", TAG, ".csv")),
          row.names = FALSE)

cat("\n==== stage 1 diagnostics (", TAG, ") ====\n", sep = "")
print(diag, row.names = FALSE)
cat("\nMMT per region:\n"); print(round(mmt_vec, 1))
cat("\nSaved to ", normalizePath(OUTDIR), "\n", sep = "")