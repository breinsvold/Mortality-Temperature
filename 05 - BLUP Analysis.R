# =====================================================================
# Stage 2 -- multivariate meta-analysis of the 6 region DLNM curves
# Nebraska CVD vs UTCI, 2005-2025
#
# Consumes stage-1 output (all in stage1_out/):
#   coef_<TAG>.rds    6 x k matrix of reduced overall-cumulative coefs
#   vcov_<TAG>.rds    list of 6 (k x k) covariance matrices
#   mmt_<TAG>.rds     named vector of region MMTs
#   basis_<TAG>.rds   common basis spec (knots, bounds, lag) <- imported
#                     so the basis here is IDENTICAL to stage 1
#
# Produces (in stage2_out/):
#   pooled_curve_<TAG>.pdf        state-wide random-effects curve
#   blup_curves_<TAG>.pdf         6 shrunk region curves overlaid  [Comp A]
#   urban_vs_rural_<TAG>.pdf      3-urban pool vs 3-rural pool      [Comp B']
#   region_RR_table_<TAG>.csv     RR at hot/cold marks per region
#   pooled_pred / blup_pred .rds  objects for stage 3 (attribution)
#
# Only 6 units => intercept-only random-effects pool + BLUPs. No
# meta-regression with covariates (would over-fit). Urban/rural is
# answered by pooling the 3+3 separately, not by a 6-point indicator.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dlnm)
library(splines)
library(mixmeta)

# ---- exposure: honor a value set by a driver script, else default ---
# Must match the stage-1 run being pooled. A purrr::walk driver sets this
# before source()-ing; otherwise defaults to the primary analysis.
if (!exists("EXPOSURE")) EXPOSURE <- "utci_max"
if (!exists("OUTCOME"))  OUTCOME  <- "total_cvd_deaths"
TAG      <- paste0(OUTCOME, "_", EXPOSURE)

INDIR    <- "stage1_out"
OUTDIR   <- "stage2_out"
dir.create(OUTDIR, showWarnings = FALSE)

REGIONS  <- c("Omaha", "Lincoln", "Tri-Cities",
              "West", "Southeast", "Northeast")
URBAN    <- c("Omaha", "Lincoln", "Tri-Cities")
RURAL    <- c("West", "Southeast", "Northeast")

## --- 1. load stage-1 outputs ----------------------------------------
coef_mat  <- readRDS(file.path(INDIR, paste0("coef_",  TAG, ".rds")))
vcov_list <- readRDS(file.path(INDIR, paste0("vcov_",  TAG, ".rds")))
mmt_vec   <- readRDS(file.path(INDIR, paste0("mmt_",   TAG, ".rds")))
bspec     <- readRDS(file.path(INDIR, paste0("basis_", TAG, ".rds")))

# enforce a known region order everywhere
coef_mat  <- coef_mat[REGIONS, , drop = FALSE]
vcov_list <- vcov_list[REGIONS]
mmt_vec   <- mmt_vec[REGIONS]

## --- 1b. rebuild the EXACT common basis from stage 1 ----------------
# Because stage 1 saved knots + bounds, the basis here is guaranteed to
# match. onebasis over a UTCI grid; the reduced coefficients live in
# exactly this space.
#
# The prediction grid spans the full basis bounds (needed so the basis
# columns evaluate correctly), but everything we INTERPRET or minimise
# over is restricted below to avoid the dataless extrapolation tails
# that produced the spurious -48.8 C MMT.
axis_lo   <- bspec$bound[1]
axis_hi   <- bspec$bound[2]
utci_grid <- seq(axis_lo, axis_hi, by = 0.5)

ob <- onebasis(utci_grid, fun = "ns",
               knots = bspec$knots, Bound = bspec$bound)

# guard: basis dimension must equal the coefficient dimension
if (ncol(ob) != ncol(coef_mat))
  stop("Basis dim (", ncol(ob), ") != coef dim (", ncol(coef_mat),
       "). Stage 1 and stage 2 bases are out of sync.")

# ---- MMT search window: UTCI "no thermal stress" band --------------
# Standard UTCI assessment scale (Broede et al. 2012): 9-26 C = no
# thermal stress. Constrain the minimum-mortality search to this
# physiologically-neutral band so the optimum can't land in an
# extrapolated tail. A binding constraint at 9 or 26 is itself
# informative -- it means the true optimum lies outside the no-stress
# band (plausible for cold-dominated CVD).
NO_STRESS <- c(9, 26)
mmt_ok    <- utci_grid >= NO_STRESS[1] & utci_grid <= NO_STRESS[2]

# ---- plotting window: trim sparse tails so curves are readable ------
# stage 1 wrote per-region 1st/99th UTCI percentiles; use the widest
# p01 and p99 across regions as the x-limits for all plots. Falls back
# to the basis bounds if the diagnostics file is absent.
plt <- tryCatch({
  dg <- read.csv(file.path(INDIR, paste0("diagnostics_", TAG, ".csv")))
  c(min(dg$utci_p01), max(dg$utci_p99))
}, error = function(e) bspec$bound)
cat(sprintf("Plot x-window (data p01-p99): %.1f to %.1f C\n", plt[1], plt[2]))

## --- 2. pooled random-effects meta-analysis -------------------------
# multivariate meta: pools the k-vector of reduced coefs across regions,
# estimating between-region (co)variance by REML.
mv <- mixmeta(coef_mat, vcov_list, method = "reml")

cat("\n===== pooled multivariate meta (", TAG, ") =====\n", sep = "")
print(summary(mv))

# heterogeneity: multivariate Cochran Q, plus a scalar I2 summary
q  <- qtest(mv)
cat("\n--- heterogeneity (Cochran Q) ---\n"); print(q)
i2 <- max(0, (q$Q[1] - q$df[1]) / q$Q[1])
cat(sprintf("Approx. I2 = %.1f%%  (share of variation that is between-region)\n",
            100 * i2))

## --- 3. pooled state-wide curve -------------------------------------
# centre at the pooled MMT, searched ONLY within the no-stress band.
pool_lp   <- as.numeric(ob %*% coef(mv))
pool_mmt  <- utci_grid[mmt_ok][ which.min(pool_lp[mmt_ok]) ]

# diagnostic: where would the min have fallen if searched over all data?
uncon_mmt <- utci_grid[which.min(pool_lp)]
cat(sprintf("\nPooled MMT (constrained to 9-26 C) = %.1f C\n", pool_mmt))
cat(sprintf("Unconstrained min over full grid  = %.1f C", uncon_mmt))
if (pool_mmt %in% NO_STRESS)
  cat("  <-- NOTE: MMT is AT a no-stress boundary; true optimum",
      "likely lies outside the band.")
cat("\n")

pool_pred <- crosspred(ob, coef = coef(mv), vcov = vcov(mv),
                       model.link = "log", at = utci_grid, cen = pool_mmt)
saveRDS(pool_pred, file.path(OUTDIR, paste0("pooled_pred_", TAG, ".rds")))

pdf(file.path(OUTDIR, paste0("pooled_curve_", TAG, ".pdf")),
    width = 7, height = 5)
plot(pool_pred, "overall", lwd = 2, col = "#1b4965", ci = "area",
     ci.arg = list(col = adjustcolor("#1b4965", 0.15)),
     xlim = plt,
     xlab = paste0(EXPOSURE, " (UTCI, \u00B0C)"),
     ylab = "RR of CVD death (vs pooled MMT)",
     main = sprintf("Pooled Nebraska UTCI-CVD curve\nMMT = %.1f \u00B0C   I2 = %.0f%%",
                    pool_mmt, 100 * i2))
abline(h = 1, col = "grey70"); abline(v = pool_mmt, lty = 3)
dev.off()

## --- 4. BLUPs: shrunk region curves  [Comparison A] -----------------
# each region's coefs pulled toward the pool by its own imprecision.
bl <- blup(mv, vcov = TRUE)                 # list aligned to coef_mat rows
names(bl) <- REGIONS

# turn each BLUP into a curve, all centred at the SAME (pooled) MMT so
# they are directly comparable on one axis.
blup_pred <- lapply(REGIONS, function(g) {
  crosspred(ob, coef = bl[[g]]$blup, vcov = bl[[g]]$vcov,
            model.link = "log", at = utci_grid, cen = pool_mmt)
})
names(blup_pred) <- REGIONS
saveRDS(blup_pred, file.path(OUTDIR, paste0("blup_pred_", TAG, ".rds")))

cols <- c(Omaha="#1b4965", Lincoln="#3a7ca5", `Tri-Cities`="#5fa8d3",
          West="#b5651d", Southeast="#d98c3f", Northeast="#e8a95f")

pdf(file.path(OUTDIR, paste0("blup_curves_", TAG, ".pdf")),
    width = 7.5, height = 5.5)
inwin <- utci_grid >= plt[1] & utci_grid <= plt[2]   # data-range mask
ylim <- range(vapply(blup_pred,
                     function(p) range(p$allRRfit[inwin]), numeric(2)))
plot(utci_grid, blup_pred[[1]]$allRRfit, type = "n", ylim = ylim,
     xlim = plt,
     xlab = paste0(EXPOSURE, " (UTCI, \u00B0C)"),
     ylab = "RR of CVD death (vs pooled MMT)",
     main = "Region-specific BLUP curves (urban solid, rural dashed)")
abline(h = 1, col = "grey70"); abline(v = pool_mmt, lty = 3, col = "grey70")
for (g in REGIONS) {
  lines(utci_grid, blup_pred[[g]]$allRRfit, col = cols[g],
        lwd = if (g %in% URBAN) 2.4 else 1.8,
        lty = if (g %in% URBAN) 1   else 2)
}
legend("top", legend = REGIONS, col = cols[REGIONS], lwd = 2,
       lty = ifelse(REGIONS %in% URBAN, 1, 2), ncol = 2, bty = "n", cex = 0.9)
dev.off()

## --- 4b. per-region small multiples, WITH CIs  [Comparison A, detail]
# One panel per region: its BLUP curve + 95% CI band, with the pooled
# state-wide curve overlaid in grey behind it. A region whose band hugs
# the pooled line is indistinguishable from the state average -- this is
# the visual analogue of the low I2, now shown per region with the
# uncertainty the overlay plot could not display.
inwin <- utci_grid >= plt[1] & utci_grid <= plt[2]

# shared y-limits across panels = widest CI band in-window across regions
yl <- range(
  vapply(blup_pred, function(p)
    range(c(p$allRRlow[inwin], p$allRRhigh[inwin])), numeric(2)),
  pool_pred$allRRlow[inwin], pool_pred$allRRhigh[inwin]
)

pdf(file.path(OUTDIR, paste0("blup_facets_", TAG, ".pdf")),
    width = 9, height = 6)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

for (g in REGIONS) {
  p <- blup_pred[[g]]
  plot(utci_grid, p$allRRfit, type = "n", xlim = plt, ylim = yl,
       xlab = paste0(EXPOSURE, " (\u00B0C)"),
       ylab = "RR CVD death",
       main = sprintf("%s (%s)", g,
                      ifelse(g %in% URBAN, "urban", "rural")))
  
  # pooled curve behind, grey, with a light band
  polygon(c(utci_grid, rev(utci_grid)),
          c(pool_pred$allRRlow, rev(pool_pred$allRRhigh)),
          col = adjustcolor("grey60", 0.18), border = NA)
  lines(utci_grid, pool_pred$allRRfit, col = "grey45", lwd = 1.4, lty = 2)
  
  # this region's BLUP band + line
  polygon(c(utci_grid, rev(utci_grid)),
          c(p$allRRlow, rev(p$allRRhigh)),
          col = adjustcolor(cols[g], 0.20), border = NA)
  lines(utci_grid, p$allRRfit, col = cols[g], lwd = 2.4)
  
  abline(h = 1, col = "grey70")
  abline(v = pool_mmt, lty = 3, col = "grey70")
}
mtext(sprintf("Region BLUP curves (colour) vs pooled state curve (grey dashed)  |  I2 = %.0f%%",
              100 * i2), outer = TRUE, cex = 0.95)
par(op)
dev.off()

## --- 5. urban vs rural, pooled separately  [Comparison B'] ----------
# cleaner than a 6-point urban indicator: two 3-unit random-effects pools.
mv_urb <- mixmeta(coef_mat[URBAN, , drop = FALSE], vcov_list[URBAN], method = "reml")
mv_rur <- mixmeta(coef_mat[RURAL, , drop = FALSE], vcov_list[RURAL], method = "reml")

pred_urb <- crosspred(ob, coef = coef(mv_urb), vcov = vcov(mv_urb),
                      model.link = "log", at = utci_grid, cen = pool_mmt)
pred_rur <- crosspred(ob, coef = coef(mv_rur), vcov = vcov(mv_rur),
                      model.link = "log", at = utci_grid, cen = pool_mmt)

pdf(file.path(OUTDIR, paste0("urban_vs_rural_", TAG, ".pdf")),
    width = 7, height = 5)
plot(pred_urb, "overall", col = "#1b4965", lwd = 2.4, ci = "area",
     ci.arg = list(col = adjustcolor("#1b4965", 0.15)),
     xlim = plt,
     xlab = paste0(EXPOSURE, " (UTCI, \u00B0C)"),
     ylab = "RR of CVD death (vs pooled MMT)",
     main = "Urban vs rural pooled UTCI-CVD curves")
lines(pred_rur, "overall", col = "#b5651d", lwd = 2.4, ci = "area",
      ci.arg = list(col = adjustcolor("#b5651d", 0.15)))
abline(h = 1, col = "grey70")
legend("top", c("Urban pool (3)", "Rural pool (3)"),
       col = c("#1b4965", "#b5651d"), lwd = 2.4, bty = "n")
dev.off()

## --- 6. RR summary table per region ---------------------------------
# read each BLUP curve at fixed hot/cold marks (edges of the data window).
hot  <- plt[2]
cold <- plt[1]
idx  <- function(v) which.min(abs(utci_grid - v))

tab <- do.call(rbind, lapply(REGIONS, function(g) {
  p <- blup_pred[[g]]
  data.frame(
    region     = g,
    type       = ifelse(g %in% URBAN, "urban", "rural"),
    mmt        = round(mmt_vec[g], 1),
    RR_hot     = round(p$allRRfit[idx(hot)],   3),
    RR_hot_lo  = round(p$allRRlow[idx(hot)],   3),
    RR_hot_hi  = round(p$allRRhigh[idx(hot)],  3),
    RR_cold    = round(p$allRRfit[idx(cold)],  3),
    RR_cold_lo = round(p$allRRlow[idx(cold)],  3),
    RR_cold_hi = round(p$allRRhigh[idx(cold)], 3)
  )
}))
write.csv(tab, file.path(OUTDIR, paste0("region_RR_table_", TAG, ".csv")),
          row.names = FALSE)

# save fitted meta objects for the record / stage 3
saveRDS(list(mv = mv, mv_urb = mv_urb, mv_rur = mv_rur,
             pool_mmt = pool_mmt, i2 = i2),
        file.path(OUTDIR, paste0("meta_fits_", TAG, ".rds")))

cat("\n===== region RR summary (", TAG, ") =====\n", sep = "")
print(tab, row.names = FALSE)
cat(sprintf("\nPooled MMT = %.1f  |  I2 = %.0f%%\n", pool_mmt, 100 * i2))
cat("\nPlots written:\n")
cat("  pooled_curve_",   TAG, ".pdf   (state-wide curve)\n", sep = "")
cat("  blup_curves_",    TAG, ".pdf   (6 regions overlaid, no CI)\n", sep = "")
cat("  blup_facets_",    TAG, ".pdf   (6 regions, 1 panel each, WITH CI)\n", sep = "")
cat("  urban_vs_rural_", TAG, ".pdf   (3+3 split pools)\n", sep = "")
cat("Saved plots + tables + objects to ", normalizePath(OUTDIR), "\n", sep = "")