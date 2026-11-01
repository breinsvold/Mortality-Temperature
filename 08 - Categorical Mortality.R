# =====================================================================
# Cause analysis -- one pooled 6-region curve PER CVD cause, compared
# Nebraska CVD vs UTCI
#
# Design (Interpretation A): each cause is run through the SAME stage-1
# -> stage-2 pipeline as the main analysis (6 independent region DLNMs,
# meta-analysed into one pooled state-wide curve). Regions are
# independent, so each cause's pooling is clean. We then compare the
# cause-specific POOLED curves descriptively -- we do NOT meta-analyse
# the causes against each other, because the causes are computed from
# the same region-days (they partition total_cvd_deaths) and are not
# independent, which would bias a between-cause I2.
#
# Produces, per cause: the usual stage1_out/ + stage2_out/ files, tagged
# by cause. Plus a cross-cause comparison figure + table in cause_out/.
#
# Fixed exposure here (set EXPOSURE below). Run once per exposure if you
# want both utci_max and utci_min cause breakdowns.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(purrr)
library(dlnm)

STAGE1 <- "04 - Regional DLNM.R"
STAGE2 <- "05 - BLUP Analysis.R"

# exposure: honor a value set by a driver (EXPOSURE), else default here.
if (exists("EXPOSURE")) EXPOSURE_FIXED <- EXPOSURE else EXPOSURE_FIXED <- "utci_max"
OUTDIR <- "cause_out"; dir.create(OUTDIR, showWarnings = FALSE)

# causes to analyse. total_cvd_deaths included as the anchor/reference.
CAUSES <- c("heart_disease_deaths", "stroke_deaths",
            "hypertension_deaths", "other_vascular_deaths",
            "total_cvd_deaths")

# pretty labels for plotting
CAUSE_LAB <- c(
  heart_disease_deaths  = "Heart disease",
  stroke_deaths         = "Stroke",
  hypertension_deaths   = "Hypertension",
  other_vascular_deaths = "Other vascular",
  total_cvd_deaths      = "Total CVD")

## --- run the pipeline for each cause --------------------------------
run_script <- function(path, exposure, outcome) {
  env <- new.env(parent = globalenv())
  assign("EXPOSURE", exposure, envir = env)
  assign("OUTCOME",  outcome,  envir = env)
  message(sprintf(">> %s   [%s | %s]", path, outcome, exposure))
  sys.source(path, envir = env)
  invisible(NULL)
}
run_cause <- function(cause) {
  run_script(STAGE1, EXPOSURE_FIXED, cause)
  run_script(STAGE2, EXPOSURE_FIXED, cause)
}

walk(CAUSES, run_cause)

## --- gather the cause-specific pooled curves ------------------------
# stage 2 saved pooled_pred_<OUTCOME>_<EXPOSURE>.rds in stage2_out/.
# Load each, align on the common UTCI grid, and overlay.
load_pool <- function(cause) {
  tag <- paste0(cause, "_", EXPOSURE_FIXED)
  f   <- file.path("stage2_out", paste0("pooled_pred_", tag, ".rds"))
  if (!file.exists(f)) { warning("missing pooled pred for ", cause); return(NULL) }
  readRDS(f)
}
pools <- setNames(lapply(CAUSES, load_pool), CAUSES)
pools <- pools[!vapply(pools, is.null, logical(1))]

# also pull each cause's daily-count diagnostics (mean, %zero) so the
# rare subtypes are clearly flagged as noisy.
load_diag <- function(cause) {
  tag <- paste0(cause, "_", EXPOSURE_FIXED)
  f   <- file.path("stage1_out", paste0("diagnostics_", tag, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f)
  data.frame(cause = cause,
             mean_deaths_per_region_day = round(mean(d$mean_deaths), 3),
             max_pct_zero = round(max(d$pct_zero), 1))
}
diag_tab <- do.call(rbind, lapply(CAUSES, load_diag))

## --- plotting window: trim to observed data range -------------------
# CIs explode beyond the data (few days at UTCI extremes); restrict all
# cause plots to the widest p01-p99 across regions (from any cause's
# stage-1 diagnostics -- the UTCI range is the same regardless of cause).
plt <- tryCatch({
  # use total_cvd diagnostics if present, else the first available
  dgf <- file.path("stage1_out",
                   paste0("diagnostics_total_cvd_deaths_", EXPOSURE_FIXED, ".csv"))
  if (!file.exists(dgf)) {
    cand <- list.files("stage1_out",
                       pattern = paste0("diagnostics_.*_", EXPOSURE_FIXED, ".csv"),
                       full.names = TRUE)
    dgf <- cand[1]
  }
  dg <- read.csv(dgf)
  c(min(dg$utci_p01), max(dg$utci_p99))
}, error = function(e) range(lapply(pools, function(p) range(p$predvar))))
cat(sprintf("Plot x-window (data p01-p99): %.1f to %.1f C\n", plt[1], plt[2]))

# helper: y-range of a set of preds restricted to the data window
yr_in <- function(preds, lo, hi, cap_hi = 3.0, cap_lo = 0.4) {
  vals <- unlist(lapply(preds, function(p) {
    m <- p$predvar >= lo & p$predvar <= hi
    c(p$allRRlow[m], p$allRRhigh[m])
  }))
  yr <- range(vals, finite = TRUE)
  c(max(yr[1], cap_lo), min(yr[2], cap_hi))
}

## --- cross-cause comparison figure ----------------------------------
# each pooled crosspred carries its own predvar grid ($predvar) and RRfit.
cols <- c(heart_disease_deaths="#1b4965", stroke_deaths="#bb3e03",
          hypertension_deaths="#8a5a44", other_vascular_deaths="#5f7470",
          total_cvd_deaths="#333333")

pdf(file.path(OUTDIR, paste0("cause_comparison_", EXPOSURE_FIXED, ".pdf")),
    width = 8, height = 5.5)
xr <- plt
yr <- yr_in(pools, plt[1], plt[2])

plot(NA, xlim = xr, ylim = yr,
     xlab = paste0(EXPOSURE_FIXED, " (UTCI, \u00B0C)"),
     ylab = "RR of death (vs cause-specific pooled MMT)",
     main = "UTCI-CVD pooled curves by cause")
abline(h = 1, col = "grey70")
for (cs in names(pools)) {
  p <- pools[[cs]]
  lines(p$predvar, p$allRRfit, col = cols[cs],
        lwd = if (cs == "total_cvd_deaths") 3 else 2,
        lty = if (cs == "total_cvd_deaths") 2 else 1)
}
legend("top", legend = CAUSE_LAB[names(pools)],
       col = cols[names(pools)],
       lwd = 2, lty = ifelse(names(pools) == "total_cvd_deaths", 2, 1),
       bty = "n", ncol = 2, cex = 0.9)
dev.off()

## --- faceted: one panel per cause, WITH CI band ---------------------
# each cause's pooled curve + 95% CI, with total CVD overlaid in grey
# behind it for reference. Shared axes so causes are comparable.
tot <- pools[["total_cvd_deaths"]]          # reference (may be NULL)

# shared limits from all causes (clamped so a sparse cause can't blow up)
facet_x <- plt
facet_y <- yr_in(pools, plt[1], plt[2])

ncause <- length(pools)
nc <- 3; nr <- ceiling(ncause / nc)

pdf(file.path(OUTDIR, paste0("cause_facets_", EXPOSURE_FIXED, ".pdf")),
    width = 3.1 * nc, height = 3.0 * nr)
op <- par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

for (cs in names(pools)) {
  p <- pools[[cs]]
  plot(NA, xlim = facet_x, ylim = facet_y,
       xlab = paste0(EXPOSURE_FIXED, " (\u00B0C)"),
       ylab = "RR of death",
       main = CAUSE_LAB[cs])
  
  # total CVD reference behind (skip on the total panel itself)
  if (!is.null(tot) && cs != "total_cvd_deaths") {
    lines(tot$predvar, tot$allRRfit, col = "grey55", lwd = 1.3, lty = 2)
  }
  
  # this cause's CI band + line
  polygon(c(p$predvar, rev(p$predvar)),
          c(p$allRRlow, rev(p$allRRhigh)),
          col = adjustcolor(cols[cs], 0.20), border = NA)
  lines(p$predvar, p$allRRfit, col = cols[cs], lwd = 2.4)
  
  abline(h = 1, col = "grey70")
}
mtext(sprintf("Pooled UTCI-CVD curve by cause, with 95%% CI  (grey dashed = total CVD)  |  %s",
              EXPOSURE_FIXED), outer = TRUE, cex = 0.9)
par(op)
dev.off()

## --- cross-cause RR summary table -----------------------------------
# read each cause's pooled curve at the data-window edges (not the full
# basis range, where CIs are meaningless extrapolation).
hot  <- plt[2]
cold <- plt[1]
at_idx <- function(p, v) which.min(abs(p$predvar - v))

rr_tab <- do.call(rbind, lapply(names(pools), function(cs) {
  p <- pools[[cs]]
  ih <- at_idx(p, hot); ic <- at_idx(p, cold)
  data.frame(
    cause      = CAUSE_LAB[cs],
    RR_hot     = round(p$allRRfit[ih],  3),
    RR_hot_lo  = round(p$allRRlow[ih],  3),
    RR_hot_hi  = round(p$allRRhigh[ih], 3),
    RR_cold    = round(p$allRRfit[ic],  3),
    RR_cold_lo = round(p$allRRlow[ic],  3),
    RR_cold_hi = round(p$allRRhigh[ic], 3))
}))

if (!is.null(diag_tab))
  rr_tab <- merge(rr_tab,
                  transform(diag_tab, cause = CAUSE_LAB[cause]),
                  by = "cause", all.x = TRUE)

write.csv(rr_tab,
          file.path(OUTDIR, paste0("cause_RR_table_", EXPOSURE_FIXED, ".csv")),
          row.names = FALSE)

cat("\n===== cause comparison (", EXPOSURE_FIXED, ") =====\n", sep = "")
cat(sprintf("Read at cold = %.1f C, hot = %.1f C (common window)\n\n",
            cold, hot))
print(rr_tab, row.names = FALSE)
cat("\nNOTE: interpret hypertension / other_vascular cautiously -- low\n")
cat("daily counts (see max_pct_zero) make those pooled curves noisy.\n")
cat("Causes are NOT independent (they partition total_cvd_deaths), so\n")
cat("this is a DESCRIPTIVE comparison of pooled curves, not a formal\n")
cat("between-cause test.\n")
cat("\nFigures in ", OUTDIR, "/:\n", sep = "")
cat("  cause_comparison_", EXPOSURE_FIXED, ".pdf  (all causes overlaid, no CI)\n", sep = "")
cat("  cause_facets_",     EXPOSURE_FIXED, ".pdf  (one panel per cause, WITH CI)\n", sep = "")
cat("Per-cause 6-region panels are in stage2_out/blup_facets_<cause>_",
    EXPOSURE_FIXED, ".pdf\n", sep = "")
cat("\nSaved to ", normalizePath(OUTDIR), "\n", sep = "")