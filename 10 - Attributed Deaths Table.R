# =====================================================================
# 10 - Attribution tables, max/min stacked (presentation only)
# Nebraska CVD vs UTCI
#
# Builds TWO tables. Each cell stacks the two exposures: first line =
# utci_max (daytime), second line = utci_min (night-time), each with its
# 95% CI. This is a side-by-side presentation of two DIFFERENT exposures
# -- NOT an average (both analyses use the same deaths, so averaging
# would double count).
#
#   Table 1 (Total CVD): rows = total cold, total heat, then each UTCI
#     thermal-stress category (all except no-stress); columns =
#     Combined (whole pooled curve), Urban, Rural.
#   Table 2 (by cause): rows = CVD causes; columns = cold/heat for the
#     whole curve, urban pool and rural pool.
#
# Values are attributable deaths AN (95% CI); bold = CI excludes zero;
# "\u2014" = no days in range for that pool/exposure (e.g. utci_min
# rarely reaches heat stress).
#
# Reads stage-3 output (script 08) for BOTH exposures -- no
# recomputation. Run script 08 for utci_max and utci_min first (e.g. via
# the 06 driver).
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(dplyr)
library(knitr)
library(kableExtra)

INDIR  <- "stage3_out"
OUTDIR <- "stage3_out"

CAUSE_ORDER <- c("Heart disease","Stroke","Hypertension","Other vascular")

# category order (cold -> hot), excluding no-stress
CAT_ORDER <- c("Extreme cold","Very strong cold","Strong cold",
               "Moderate cold","Slight cold",
               "Moderate heat","Strong heat","Very strong heat","Extreme heat")

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
  if (length(an) == 0 || is.na(an)) return("\u2014")
  s <- sprintf("%.0f (%.0f, %.0f)", an, lo, hi)
  if (is_sig(lo, hi)) paste0("<strong>", s, "</strong>") else s
}

# pull one value (AN, lo, hi) for a given filter from a df; NA if absent.
pick <- function(df, ...) {
  r <- df %>% filter(...)
  if (nrow(r) == 0) return(c(AN=NA, lo=NA, hi=NA))
  c(AN = r$AN[1], lo = r$AN_lo[1], hi = r$AN_hi[1])
}

## --- TABLE 1: Total CVD by UTCI category, Combined / Urban / Rural --
# Cells stack max over min with NO per-line labels (explained in caption).
stack_unlab <- function(mx, mn) {
  paste0(fmt1(mx["AN"], mx["lo"], mx["hi"]),
         "<br/>",
         fmt1(mn["AN"], mn["lo"], mn["hi"]))
}
# Combined = whole pooled curve (group "All" in overall_pooled files)
cell_comb <- function(split_val) {
  stack_unlab(pick(overall_pooled_max, split == split_val),
              pick(overall_pooled_min, split == split_val))
}
cell_pool <- function(grp, split_val) {
  stack_unlab(pick(overall_ur_max, group == grp, split == split_val),
              pick(overall_ur_min, group == grp, split == split_val))
}

# keep only categories that appear for at least one pool/exposure
present <- unique(c(overall_pooled_max$split, overall_pooled_min$split,
                    overall_ur_max$split,     overall_ur_min$split))
cats <- CAT_ORDER[CAT_ORDER %in% present]

t1_splits <- c("Cold (<9)", "Heat (>26)", cats)
t1_labels <- c("Total cold (<9)", "Total heat (>26)", cats)

t1 <- data.frame(
  `UTCI category` = t1_labels,
  Combined = vapply(t1_splits, cell_comb, character(1)),
  Urban    = vapply(t1_splits, function(s) cell_pool("Urban", s), character(1)),
  Rural    = vapply(t1_splits, function(s) cell_pool("Rural", s), character(1)),
  check.names = FALSE, row.names = NULL)

note1 <- paste0(
  "Attributable Total-CVD deaths from UTCI by thermal-stress category. ",
  "Combined = whole pooled curve; Urban/Rural = pool-specific curves. ",
  "Each cell: first line = utci_max (daytime), second line = utci_min ",
  "(night-time). Cold < 9\u00B0C, Heat > 26\u00B0C. Values: AN (95% CI); ",
  "bold = 95% CI excludes zero; \u2014 = no days in range for that ",
  "pool/exposure (min rarely reaches heat stress).")

k1 <- kbl(t1, booktabs = TRUE, format = "html", escape = FALSE) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover")) %>%
  footnote(general = note1, general_title = "", footnote_as_chunk = TRUE) %>%
  pack_rows("Total", 1, 2) %>%
  pack_rows("By UTCI category", 3, nrow(t1))
save_kable(k1, file.path(OUTDIR, "categorical_mortality.html"))

## --- TABLE 2: per-cause, max/min stacked ----------------------------
# columns: Combined cold, Combined heat, Urban cold, Urban heat,
#          Rural cold, Rural heat. Each cell stacks max over min
#          (no per-line labels; explained in the footnote).
cell_pooled <- function(cause_lab, split_val) {
  stack_unlab(pick(causes_pooled_max, cause == cause_lab, split == split_val),
              pick(causes_pooled_min, cause == cause_lab, split == split_val))
}
cell_ur <- function(cause_lab, split_val, grp) {
  stack_unlab(pick(causes_ur_max, cause == cause_lab, split == split_val, group == grp),
              pick(causes_ur_min, cause == cause_lab, split == split_val, group == grp))
}

causes_in <- CAUSE_ORDER[CAUSE_ORDER %in%
                           unique(c(causes_pooled_max$cause, causes_pooled_min$cause))]

t2 <- do.call(rbind, lapply(causes_in, function(cl) {
  data.frame(
    Cause            = cl,
    `Combined cold`  = cell_pooled(cl, "Cold (<9)"),
    `Combined heat`  = cell_pooled(cl, "Heat (>26)"),
    `Urban cold`     = cell_ur(cl, "Cold (<9)",  "Urban"),
    `Urban heat`     = cell_ur(cl, "Heat (>26)", "Urban"),
    `Rural cold`     = cell_ur(cl, "Cold (<9)",  "Rural"),
    `Rural heat`     = cell_ur(cl, "Heat (>26)", "Rural"),
    check.names = FALSE)
}))

note2 <- paste0(
  "Attributable deaths by CVD cause from UTCI. ",
  "Combined = whole pooled curve; Urban/Rural = pool-specific curves. ",
  "Each cell: first line = utci_max (daytime), second line = utci_min ",
  "(night-time). Cold < 9\u00B0C, Heat > 26\u00B0C. Values: AN (95% CI); ",
  "bold = 95% CI excludes zero; \u2014 = no days in range for that ",
  "pool/exposure (min rarely reaches heat stress).")

k2 <- kbl(t2, booktabs = TRUE, format = "html", escape = FALSE) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover")) %>%
  footnote(general = note2, general_title = "", footnote_as_chunk = TRUE)
save_kable(k2, file.path(OUTDIR, "causational_mortality.html"))

cat("\nDone. Tables in ", normalizePath(OUTDIR), "\n", sep = "")
cat("  categorical_mortality.html  (Total CVD: totals + categories; Combined/Urban/Rural)\n")
cat("  causational_mortality.html     (each cause: cold/heat; Combined/Urban/Rural)\n")
cat("\nEach cell stacks utci_max (daytime) over utci_min (night-time).\n")
cat("These are two DIFFERENT exposures shown side by side, NOT averaged.\n")