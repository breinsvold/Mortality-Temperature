# =====================================================================
# 09 - Combined max/min attribution tables (presentation only)
# Nebraska CVD vs UTCI
#
# Rebuilds the two stage-3 tables so each CELL shows BOTH exposures,
# stacked: first line = utci_max, second line = utci_min (each with its
# 95% CI). This is a side-by-side presentation of two DIFFERENT
# exposures -- it is NOT an average (the two analyses use the same
# deaths related to different daily exposures, so averaging would double
# count; stacking lets the reader see both honestly).
#
# Reads the CSVs that stage 3 (script 08) already wrote for each
# exposure -- no re-computation. Run 08 for BOTH utci_max and utci_min
# first (e.g. via the 06 driver).
#
# Note: utci_min rarely reaches heat stress (>26), so its HEAT cells are
# typically empty; those show "\u2014" on the min line.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

INDIR  <- "stage3_out"
OUTDIR <- "stage3_out"

CAUSE_ORDER <- c("Heart disease","Stroke","Hypertension","Other vascular")

## --- load both exposures' CSVs --------------------------------------
rd <- function(stub, exp) {
  f <- file.path(INDIR, paste0(stub, "_", exp, ".csv"))
  if (!file.exists(f)) stop("missing: ", f,
                            "\nRun script 08 for exposure '", exp, "' first.")
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

overall_pooled_max <- rd("overall_pooled",     "utci_max")
overall_ur_max     <- rd("overall_urbanrural", "utci_max")
causes_pooled_max  <- rd("causes_pooled",      "utci_max")
causes_ur_max      <- rd("causes_urbanrural",  "utci_max")

overall_pooled_min <- rd("overall_pooled",     "utci_min")
overall_ur_min     <- rd("overall_urbanrural", "utci_min")
causes_pooled_min  <- rd("causes_pooled",      "utci_min")
causes_ur_min      <- rd("causes_urbanrural",  "utci_min")

## --- helpers ---------------------------------------------------------
# significant if CI excludes zero; bold the string if so.
is_sig <- function(lo, hi) (lo > 0 & hi > 0) | (lo < 0 & hi < 0)
fmt1 <- function(an, lo, hi) {
  if (is.na(an)) return("\u2014")
  s <- sprintf("%.0f (%.0f, %.0f)", an, lo, hi)
  if (is_sig(lo, hi)) paste0("<strong>", s, "</strong>") else s
}
# stack max over min in one cell: "max<br/>min", each labelled.
stack_cell <- function(an_mx, lo_mx, hi_mx, an_mn, lo_mn, hi_mn) {
  mx <- fmt1(an_mx, lo_mx, hi_mx)
  mn <- fmt1(an_mn, lo_mn, hi_mn)
  paste0("max: ", mx, "<br/>min: ", mn)
}

# pull one value (AN, lo, hi) for a given filter from a df; NA if absent.
pick <- function(df, ...) {
  r <- df %>% filter(...)
  if (nrow(r) == 0) return(c(AN=NA, lo=NA, hi=NA))
  c(AN = r$AN[1], lo = r$AN_lo[1], hi = r$AN_hi[1])
}

## --- TABLE 1: OVERALL (Total CVD), max/min stacked ------------------
# rows: cold total, heat total, each cold/heat UTCI category, urban/rural
# cold + heat. Cold cells get both exposures; heat cells too (min heat
# usually blank).

# category labels present in the overall pooled files (cold + heat bands)
cats_present <- unique(c(
  overall_pooled_max$split, overall_pooled_min$split))
# keep only the non-total splits, ordered cold->hot
cat_order <- c("Extreme cold","Very strong cold","Strong cold",
               "Moderate cold","Slight cold",
               "Moderate heat","Strong heat","Very strong heat","Extreme heat")
cat_rows <- cat_order[cat_order %in% cats_present]

mkrow <- function(label, split_val, src = "pooled", grp = NULL) {
  if (src == "pooled") {
    mx <- pick(overall_pooled_max, split == split_val)
    mn <- pick(overall_pooled_min, split == split_val)
  } else {
    mx <- pick(overall_ur_max, split == split_val, group == grp)
    mn <- pick(overall_ur_min, split == split_val, group == grp)
  }
  data.frame(Category = label,
             `Attributable deaths (max / min, 95% CI)` =
               stack_cell(mx["AN"],mx["lo"],mx["hi"], mn["AN"],mn["lo"],mn["hi"]),
             check.names = FALSE)
}

t1 <- bind_rows(
  data.frame(Category="\u2014 Cold vs heat (whole curve) \u2014",
             `Attributable deaths (max / min, 95% CI)`="", check.names=FALSE),
  mkrow("Cold (<9)",  "Cold (<9)"),
  mkrow("Heat (>26)", "Heat (>26)"),
  data.frame(Category="\u2014 By UTCI category \u2014",
             `Attributable deaths (max / min, 95% CI)`="", check.names=FALSE),
  bind_rows(lapply(cat_rows, function(cc) mkrow(cc, cc))),
  data.frame(Category="\u2014 Urban / rural pools \u2014",
             `Attributable deaths (max / min, 95% CI)`="", check.names=FALSE),
  mkrow("Urban - Cold (<9)",  "Cold (<9)",  "ur", "Urban"),
  mkrow("Urban - Heat (>26)", "Heat (>26)", "ur", "Urban"),
  mkrow("Rural - Cold (<9)",  "Cold (<9)",  "ur", "Rural"),
  mkrow("Rural - Heat (>26)", "Heat (>26)", "ur", "Rural")
)

k1 <- kbl(t1, caption = paste0("Attributable Total-CVD deaths from UTCI, ",
                               "max and min stacked per cell (max = daytime, min = night-time ",
                               "felt-temp). Cold < 9\u00B0C, Heat > 26\u00B0C. Bold = 95% CI excludes ",
                               "zero. \u2014 = no days in range (min rarely reaches heat stress)."),
          booktabs = TRUE, format = "html", escape = FALSE) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover")) %>%
  row_spec(which(grepl("^\u2014", t1$Category)), bold = TRUE, background = "#f0f0f0")
save_kable(k1, file.path(OUTDIR, "combined_table1_overall_maxmin.html"))

## --- TABLE 2: per-cause, max/min stacked ----------------------------
# columns: Cold (whole), Heat (whole), Urban cold, Urban heat,
#          Rural cold, Rural heat. Each cell stacks max over min.
cell_pooled <- function(cause_lab, split_val) {
  mx <- pick(causes_pooled_max, cause == cause_lab, split == split_val)
  mn <- pick(causes_pooled_min, cause == cause_lab, split == split_val)
  stack_cell(mx["AN"],mx["lo"],mx["hi"], mn["AN"],mn["lo"],mn["hi"])
}
cell_ur <- function(cause_lab, split_val, grp) {
  mx <- pick(causes_ur_max, cause == cause_lab, split == split_val, group == grp)
  mn <- pick(causes_ur_min, cause == cause_lab, split == split_val, group == grp)
  stack_cell(mx["AN"],mx["lo"],mx["hi"], mn["AN"],mn["lo"],mn["hi"])
}

causes_in <- CAUSE_ORDER[CAUSE_ORDER %in%
                           unique(c(causes_pooled_max$cause, causes_pooled_min$cause))]

t2 <- do.call(rbind, lapply(causes_in, function(cl) {
  data.frame(
    Cause         = cl,
    `Cold (<9)`   = cell_pooled(cl, "Cold (<9)"),
    `Heat (>26)`  = cell_pooled(cl, "Heat (>26)"),
    `Urban cold`  = cell_ur(cl, "Cold (<9)",  "Urban"),
    `Urban heat`  = cell_ur(cl, "Heat (>26)", "Urban"),
    `Rural cold`  = cell_ur(cl, "Cold (<9)",  "Rural"),
    `Rural heat`  = cell_ur(cl, "Heat (>26)", "Rural"),
    check.names = FALSE)
}))

k2 <- kbl(t2, caption = paste0("Attributable deaths by CVD cause from UTCI, ",
                               "max and min stacked per cell. Cold < 9\u00B0C, Heat > 26\u00B0C. ",
                               "Values: AN (95% CI); bold = significant. \u2014 = no days in range."),
          booktabs = TRUE, format = "html", escape = FALSE) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover")) %>%
  add_header_above(c(" "=1, "Whole curve"=2, "Urban pool"=2, "Rural pool"=2))
save_kable(k2, file.path(OUTDIR, "combined_table2_by_cause_maxmin.html"))

cat("\nDone. Combined max/min tables in ", normalizePath(OUTDIR), "\n", sep="")
cat("  combined_table1_overall_maxmin.html\n")
cat("  combined_table2_by_cause_maxmin.html\n")
cat("\nEach cell: 'max:' line (daytime UTCI) over 'min:' line (night-time).\n")
cat("These are two DIFFERENT exposures shown side by side, NOT averaged.\n")