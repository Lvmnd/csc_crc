# =============================================================================
# 00_set_seed.R — Global random seed initialisation for R scripts
#
# Purpose:
#   Every R script in this project sources this file FIRST to ensure all
#   stochastic analyses are fully reproducible.  The global seed (42) is
#   documented in logs/decisions.md (entry 2026-07-10-01).
#
# Usage:
#   source("scripts/00_set_seed.R")
#
# Effects:
#   - Calls set.seed(42) immediately
#   - Registers a .Random.seed in the global environment if not present
#   - Logs the seed to the console on source
# =============================================================================

GLOBAL_SEED <- 42L

set.seed(GLOBAL_SEED)

cat("[00_set_seed.R] Global random seed set to", GLOBAL_SEED, "\n")

# ── Helper: wrapper for functions that need an explicit seed argument ──────
with_global_seed <- function(expr) {
  set.seed(GLOBAL_SEED)
  expr
}

# If you want to enforce that every script *must* source this file,
# uncomment the line below.  see also scripts/00_set_seed.py
# invisible()
