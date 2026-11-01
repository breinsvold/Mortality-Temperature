# =====================================================================
# Driver -- full pipeline for BOTH exposures
# Nebraska CVD vs UTCI
#
# For each exposure in c("utci_max","utci_min"), runs:
#   07 - Cause Analysis      (internally runs stage 1 + stage 2 for every
#                             cause, including total_cvd_deaths, and the
#                             cross-cause comparison figures)
#   08 - Attributable Deaths (stage 3: attributable deaths, cold/heat, by
#                             UTCI category, urban/rural, per cause, with
#                             exact Monte Carlo 95% CIs and kableExtra
#                             tables)
#
# 07 already covers stage 1 (04) and stage 2 (05) for all causes, so we
# do NOT call 04/05 separately -- that would just redo the total_cvd run.
#
# Each script runs in a fresh child environment with EXPOSURE pre-set;
# the scripts honor a pre-set EXPOSURE (guarded assignments), so this
# wins over their internal defaults.
# =====================================================================

.libPaths(c("~/R/x86_64-module-4.5/library", .libPaths()))

library(purrr)

CAUSE_ANALYSIS <- "08 - Categorical Mortality.R"
ATTRIBUTION    <- "09 - Attributed Deaths.R"
EXPOSURES      <- c("utci_max", "utci_min")

# run one script in a fresh child env, with EXPOSURE bound so the
# script's `if (exists/!exists("EXPOSURE"))` guard picks it up.
run_script <- function(path, exposure) {
  env <- new.env(parent = globalenv())
  assign("EXPOSURE", exposure, envir = env)
  message("\n=====================================================")
  message(sprintf(">> %s   [EXPOSURE = %s]", path, exposure))
  message("=====================================================")
  sys.source(path, envir = env)
  invisible(NULL)
}

# per exposure: cause analysis (stage 1+2 for all causes) THEN attribution
run_pipeline <- function(exposure) {
  run_script(CAUSE_ANALYSIS, exposure)   # -> stage1_out/, stage2_out/, cause_out/
  run_script(ATTRIBUTION,    exposure)   # -> stage3_out/
}

walk(EXPOSURES, run_pipeline)

message("\nAll exposures complete.")
message("Stage 1/2 per cause : stage1_out/ , stage2_out/  (tagged by cause & exposure)")
message("Cause comparison    : cause_out/                 (tagged by exposure)")
message("Attributable deaths : stage3_out/                (tagged by exposure)")