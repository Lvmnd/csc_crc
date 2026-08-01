# =============================================================================
# 00_set_seed.py — Global random seed initialisation for Python scripts
#
# Purpose:
#   Every Python script in this project imports this module FIRST to ensure
#   all stochastic analyses are fully reproducible.  The global seed (42) is
#   documented in logs/decisions.md (entry 2026-07-10-01).
#
# Usage:
#   from scripts.00_set_seed import GLOBAL_SEED, set_global_seed
#   set_global_seed()
#
# Effects:
#   - Seeds random, numpy, and scikit-learn (if available) in one call.
# =============================================================================

import os
import random
import warnings

GLOBAL_SEED = 42


def set_global_seed(seed: int = GLOBAL_SEED) -> None:
    """Seed all stochastic libraries with a single value.

    Parameters
    ----------
    seed : int
        The seed value (default: ``GLOBAL_SEED = 42``).
    """
    # Python built-in
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)

    # NumPy
    try:
        import numpy as np

        np.random.seed(seed)
    except ModuleNotFoundError:
        pass

    # scikit-learn
    try:
        from sklearn.utils import check_random_state

        check_random_state(seed)
    except (ModuleNotFoundError, ImportError):
        pass

    # PyTorch (optional — uncomment if torch is used)
    # try:
    #     import torch
    #     torch.manual_seed(seed)
    #     torch.cuda.manual_seed_all(seed)
    #     torch.backends.cudnn.deterministic = True
    #     torch.backends.cudnn.benchmark = False
    # except ModuleNotFoundError:
    #     pass

    print(f"[00_set_seed.py] Global random seed set to {seed}")


# ── Auto-seed on import ────────────────────────────────────────────────────
set_global_seed(GLOBAL_SEED)
