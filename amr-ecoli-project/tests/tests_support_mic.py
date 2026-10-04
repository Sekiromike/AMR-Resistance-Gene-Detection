"""Shared synthetic MIC data for model tests (not a test module)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate(n: int = 3000, seed: int = 0):
    """log2 MIC = -3 + 3*gyrA + 2*qnr + 0*noise_gene + N(0, 0.7); panel 0.25..8 mg/L."""
    rng = np.random.default_rng(seed)
    X = rng.binomial(1, [0.3, 0.2, 0.5], size=(n, 3)).astype(float)
    latent = -3.0 + X @ np.array([3.0, 2.0, 0.0]) + rng.normal(0, 0.7, n)
    d = np.ceil(latent)
    lo = np.where(d <= -2, -np.inf, np.where(d > 3, 4.0, d))
    hi = np.where(d <= -2, -2.0, np.where(d > 3, np.inf, d))
    return X, lo, hi


def small_cohort(seed: int = 0):
    X, lo, hi = simulate(n=600, seed=seed)
    rows, features = [], []
    for i in range(600):
        sign, value = ("<=", 0.25) if np.isneginf(lo[i]) else ((">", 8.0) if np.isposinf(hi[i]) else ("=", 2.0 ** hi[i]))
        rows.append({"isolate_id": f"I{i}", "antibiotic": "ciprofloxacin", "evaluation_split": "development",
                     "lineage_group": f"SLV:ST{i % 60}", "intended_use_population": "human_clinical",
                     "measurement_sign": sign, "ast_value": value})
        for j, name in enumerate(("gyrA_S83L", "qnrS1", "noise")):
            if X[i, j]:
                features.append({"isolate_id": f"I{i}", "element_symbol": name})
    # A determinant seen only in outer fold 0 must never enter fold 0's model.
    features += [{"isolate_id": f"I{i}", "element_symbol": "fold0_only"} for i in range(600) if (i % 60) % 5 == 0]
    cohort = pd.DataFrame(rows)
    splits = pd.DataFrame({"isolate_id": cohort["isolate_id"], "evaluation_split": "development",
                           "cv_group": cohort["lineage_group"],
                           "cv_fold": [(i % 60) % 5 for i in range(600)]})
    return cohort, splits, pd.DataFrame(features)
