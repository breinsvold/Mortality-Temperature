# =====================================================================
# 11 - I2 (within-pool heterogeneity) across outcomes, pools, exposures
# Nebraska CVD vs UTCI
#
# For every outcome (Total CVD + 4 causes) and every pool (Overall 6-
# region, Urban 3-region, Rural 3-region), extracts I2 from the stage-2
# meta fits and renders a grouped bar chart, for BOTH utci_max and
# utci_min.
#
# I2 = share of the variation BETWEEN the individual region time series
# WITHIN a pool that is real heterogeneity rather than sampling noise.
# It is a WITHIN-pool measure: it does NOT compare urban vs rural.
#
# Reads stage2_out/meta_fits_<outcome>_<exposure>.rds (written by stage 2
# / the cause analysis). Run those for both exposures first.
#
# CAUTION printed on the chart: urban/rural I2 come from only 3 regions
# and are imprecise; rare-cause I2 can be unstable.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(mixmeta)

INDIR  <- "stage2_out"
OUTDIR <- "stage2_out"

OUTCOMES <- c("total_cvd_deaths","heart_disease_deaths","stroke_deaths",
              "hypertension_deaths","other_vascular_deaths")
OUT_LAB  <- c(total_cvd_deaths="Total CVD", heart_disease_deaths="Heart disease",
              stroke_deaths="Stroke", hypertension_deaths="Hypertension",
              other_vascular_deaths="Other vascular")
EXPOSURES <- c("utci_max","utci_min")
POOLS     <- c("Overall","Urban","Rural")

## --- I2 from a fitted mixmeta object --------------------------------
# I2 = max(0, (Q - df) / Q), with df = (#units - 1) * (#coefs). qtest()
# gives the multivariate Q and its df.
get_i2 <- function(fit) {
  q <- tryCatch(qtest(fit), error = function(e) NULL)
  if (is.null(q)) return(NA_real_)
  Q  <- as.numeric(q$Q[1]); df <- as.numeric(q$df[1])
  if (!is.finite(Q) || Q <= 0) return(0)
  max(0, 100 * (Q - df) / Q)
}

## --- gather I2 for every outcome x pool x exposure ------------------
rows <- list()
for (exp in EXPOSURES) for (oc in OUTCOMES) {
  f <- file.path(INDIR, paste0("meta_fits_", oc, "_", exp, ".rds"))
  if (!file.exists(f)) { message("missing: ", f, " (skipped)"); next }
  m <- readRDS(f)
  rows[[paste(exp, oc)]] <- data.frame(
    exposure = exp, outcome = OUT_LAB[oc],
    Overall = get_i2(m$mv), Urban = get_i2(m$mv_urb), Rural = get_i2(m$mv_rur),
    stringsAsFactors = FALSE)
}
i2tab <- do.call(rbind, rows)
if (is.null(i2tab) || nrow(i2tab) == 0)
  stop("No meta_fits_*.rds found in ", INDIR,
       ". Run stage 2 / cause analysis for both exposures first.")

write.csv(i2tab, file.path(OUTDIR, "i2_all_pools.csv"), row.names = FALSE)
cat("\nI2 table (%):\n"); print(i2tab, row.names = FALSE, digits = 1)

## --- chart: one panel per exposure, grouped bars ---------------------
# x = outcomes, 3 bars per outcome (Overall / Urban / Rural).
outcome_levels <- OUT_LAB[OUTCOMES]
cols <- c(Overall = "#1b4965", Urban = "#5fa8d3", Rural = "#b5651d")

pdf(file.path(OUTDIR, "i2_all_pools.pdf"), width = 10, height = 5.5)
op <- par(mfrow = c(1, length(EXPOSURES)), mar = c(8, 4, 3, 1), oma = c(0,0,2,0))

for (exp in EXPOSURES) {
  sub <- i2tab[i2tab$exposure == exp, ]
  sub <- sub[match(outcome_levels, sub$outcome), ]   # order
  M <- t(as.matrix(sub[, c("Overall","Urban","Rural")]))
  colnames(M) <- sub$outcome
  M[is.na(M)] <- 0
  bp <- barplot(M, beside = TRUE, col = cols, ylim = c(0, 100),
                las = 2, ylab = "I² (%)",
                main = exp, border = NA, cex.names = 0.8)
  # reference bands for low/moderate/high
  abline(h = c(25, 50, 75), col = "grey80", lty = 3)
  legend("topright", legend = names(cols), fill = cols, border = NA,
         bty = "n", cex = 0.85)
}
mtext("Within-pool heterogeneity (I²) by outcome, pool, and exposure",
      outer = TRUE, cex = 1.05, font = 2)
par(op)
dev.off()

cat("\nSaved:\n  i2_all_pools.csv\n  i2_all_pools.pdf\n")
cat("\nI2 = within-pool disagreement among the pool's region time series\n")
cat("(NOT urban-vs-rural). Urban/Rural use only 3 regions -> imprecise;\n")
cat("rare-cause I2 can be unstable. Read as low/moderate/high, not exact.\n")