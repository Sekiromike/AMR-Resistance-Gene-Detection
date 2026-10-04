"""Known-mechanism model: ridge-penalised interval-censored regression of log2 MIC.

For each antibiotic, log2 MIC ~ Normal(b0 + X beta, sigma^2), where X holds
presence of AMRFinderPlus AMR determinants (scripts/build_amr_features.py).
Reference MICs enter as intervals (a dilution d means the latent value is in
(d-1, d]; comparators open the interval), so censored MICs are used as
information, never imputed. beta carries an L2 penalty; b0 and sigma do not.

Nested cross-validation (docs/EVALUATION_DESIGN_AMENDMENT.md):
- outer folds come from the frozen split table and only report;
- inside each outer training set, features are kept when present in at least
  MIN_FEATURE_COUNT training isolates, and the penalty is chosen by inner
  folds over whole CV units (cv_group) on inner essential agreement;
- the model is refitted on the whole outer training set and predicts the
  held-out fold once.

Development only. The external set is not read.

Run from the repository root:
    python scripts/run_censored_regression.py --cohort ... --splits ... \\
        --features ... --population human_clinical --output ... --manifest ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import log_ndtr

try:
    from scripts.evaluate_mic_predictions import predicted_dilution, reference_bounds
except ModuleNotFoundError:  # Direct execution
    from evaluate_mic_predictions import predicted_dilution, reference_bounds  # type: ignore[no-redef]

POPULATIONS = ("human_clinical", "all_sources")
PENALTY_GRID = (0.1, 1.0, 10.0, 100.0, 1000.0)
INNER_FOLDS = 4
MIN_FEATURE_COUNT = 5
LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


class ModelError(ValueError):
    """Raised when inputs violate the modelling contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _log_interval_probability(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """log(Phi(b) - Phi(a)) for standardized bounds, stable in both tails."""
    out = np.empty_like(a)
    upper_open = np.isposinf(b)
    lower_open = np.isneginf(a) & ~upper_open
    both = ~upper_open & ~lower_open
    out[upper_open] = log_ndtr(-a[upper_open])
    out[lower_open] = log_ndtr(b[lower_open])
    aa, bb = a[both], b[both]
    upper_tail = aa > 0
    hi = np.where(upper_tail, log_ndtr(-aa), log_ndtr(bb))
    lo = np.where(upper_tail, log_ndtr(-bb), log_ndtr(aa))
    out[both] = hi + np.log1p(-np.exp(np.minimum(lo - hi, -1e-300)))
    return out


def negative_log_likelihood(theta: np.ndarray, X: np.ndarray, lo: np.ndarray, hi: np.ndarray,
                            penalty: float) -> tuple[float, np.ndarray]:
    """Penalised interval-censored normal NLL and its gradient.

    theta = [b0, log sigma, beta...]; latent interval (lo - 1, hi].
    """
    b0, log_sigma, beta = theta[0], theta[1], theta[2:]
    sigma = np.exp(log_sigma)
    mu = b0 + X @ beta
    a = (lo - 1.0 - mu) / sigma
    b = (hi - mu) / sigma
    log_p = _log_interval_probability(a, b)

    def density_ratio(z: np.ndarray) -> np.ndarray:  # phi(z) / P, zero at +-inf
        finite = np.isfinite(z)
        ratio = np.zeros_like(z)
        ratio[finite] = np.exp(-0.5 * z[finite] ** 2 - LOG_SQRT_2PI - log_p[finite])
        return ratio

    ra, rb = density_ratio(a), density_ratio(b)
    d_mu = (ra - rb) / sigma
    za = np.where(np.isfinite(a), a, 0.0)
    zb = np.where(np.isfinite(b), b, 0.0)
    d_log_sigma = za * ra - zb * rb
    value = -float(log_p.sum()) + penalty * float(beta @ beta)
    gradient = np.concatenate([[-d_mu.sum(), -d_log_sigma.sum()], -(X.T @ d_mu) + 2.0 * penalty * beta])
    return value, gradient


# sigma is bounded to [0.05, 64] log2 units. Unbounded, the line search could
# try log sigma large enough to overflow exp() when features were badly scaled
# (Evo 2 PCA scores, model round 64983779), yielding NaN gradients.
LOG_SIGMA_BOUNDS = (float(np.log(0.05)), float(np.log(64.0)))


def fit_checked(X: np.ndarray, lo: np.ndarray, hi: np.ndarray, penalty: float) -> tuple[np.ndarray, bool, str]:
    """Fit and report whether L-BFGS-B converged, with its message."""
    finite = np.concatenate([lo[np.isfinite(lo)], hi[np.isfinite(hi)]])
    if finite.size == 0:
        raise ModelError("Every training interval is unbounded")
    theta0 = np.concatenate([[float(np.median(finite)), np.log(2.0)], np.zeros(X.shape[1])])
    bounds = [(None, None), LOG_SIGMA_BOUNDS] + [(None, None)] * X.shape[1]
    result = minimize(negative_log_likelihood, theta0, args=(X, lo, hi, penalty), jac=True,
                      method="L-BFGS-B", bounds=bounds, options={"maxiter": 5000, "gtol": 1e-6})
    if not np.all(np.isfinite(result.x)):
        raise ModelError(f"Censored regression failed: {result.message}")
    return result.x, bool(result.success), str(result.message)


def fit(X: np.ndarray, lo: np.ndarray, hi: np.ndarray, penalty: float) -> np.ndarray:
    return fit_checked(X, lo, hi, penalty)[0]


def predict(theta: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, float]:
    return theta[0] + X @ theta[2:], float(np.exp(theta[1]))


def essential_agreement(mu: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> float:
    p = predicted_dilution(mu)
    error = np.maximum(np.maximum(p - hi, lo - p), 0.0)
    return float((error <= 1).mean())


def inner_folds(groups: pd.Series, n_folds: int, seed: int) -> np.ndarray:
    """Deterministic group folds: largest groups first, each to the smallest fold."""
    sizes = groups.value_counts()
    order = sorted(sizes.index, key=lambda g: (-sizes[g], hashlib.sha256(f"{seed}:{g}".encode()).hexdigest()))
    load = np.zeros(n_folds)
    assignment = {}
    for group in order:
        fold = int(np.argmin(load))
        assignment[group] = fold
        load[fold] += sizes[group]
    return groups.map(assignment).to_numpy()


def design(isolates: pd.Series, features: pd.DataFrame, vocabulary: list[str]) -> np.ndarray:
    index = {symbol: j for j, symbol in enumerate(vocabulary)}
    row_of = {isolate: i for i, isolate in enumerate(isolates)}
    X = np.zeros((len(isolates), len(vocabulary)))
    subset = features[features["isolate_id"].isin(row_of) & features["element_symbol"].isin(index)]
    X[subset["isolate_id"].map(row_of).to_numpy(), subset["element_symbol"].map(index).to_numpy()] = 1.0
    return X


def run_model(cohort: pd.DataFrame, splits: pd.DataFrame, features: pd.DataFrame, population: str,
              seed: int) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if population not in POPULATIONS:
        raise ModelError(f"population must be one of {POPULATIONS}")
    development = cohort.loc[cohort["evaluation_split"].astype(str).eq("development")].copy()
    if population == "human_clinical":
        development = development.loc[development["intended_use_population"].astype(str).eq("human_clinical")]
    folds = splits.loc[splits["evaluation_split"].astype(str).eq("development"), ["isolate_id", "cv_fold", "cv_group"]]
    development = development.merge(folds, on="isolate_id", how="left", validate="many_to_one")
    if development["cv_fold"].isna().any() or development["cv_group"].isna().any():
        raise ModelError("Every development isolate needs a cv_fold and cv_group")
    development["cv_fold"] = development["cv_fold"].astype(int)
    bounds = [reference_bounds(str(s), float(v)) for s, v in zip(development["measurement_sign"], development["ast_value"])]
    development["_lo"] = [b[0] for b in bounds]
    development["_hi"] = [b[1] for b in bounds]

    predictions: list[pd.DataFrame] = []
    records: list[dict[str, Any]] = []
    for antibiotic, drug in development.groupby("antibiotic", sort=True):
        for fold in sorted(drug["cv_fold"].unique()):
            train = drug.loc[drug["cv_fold"] != fold].reset_index(drop=True)
            test = drug.loc[drug["cv_fold"] == fold].reset_index(drop=True)
            present = features[features["isolate_id"].isin(set(train["isolate_id"]))]
            counts = present.groupby("element_symbol")["isolate_id"].nunique()
            vocabulary = sorted(counts[counts >= MIN_FEATURE_COUNT].index)
            X_train = design(train["isolate_id"], features, vocabulary)
            lo, hi = train["_lo"].to_numpy(float), train["_hi"].to_numpy(float)

            inner = inner_folds(train["cv_group"].astype(str), INNER_FOLDS, seed)
            scores = {}
            for penalty in PENALTY_GRID:
                agree = []
                for k in range(INNER_FOLDS):
                    fit_rows, eval_rows = inner != k, inner == k
                    theta = fit(X_train[fit_rows], lo[fit_rows], hi[fit_rows], penalty)
                    mu, _ = predict(theta, X_train[eval_rows])
                    agree.append(essential_agreement(mu, lo[eval_rows], hi[eval_rows]) * eval_rows.sum())
                scores[penalty] = sum(agree) / len(train)
            best = max(PENALTY_GRID, key=lambda p: (scores[p], p))  # ties: stronger penalty

            theta = fit(X_train, lo, hi, best)
            mu, sigma = predict(theta, design(test["isolate_id"], features, vocabulary))
            predictions.append(test.assign(predicted_log2_mic=mu, predicted_log2_sd=sigma, selected_penalty=best))
            coefficients = dict(zip(vocabulary, theta[2:]))
            top = sorted(coefficients.items(), key=lambda item: -abs(item[1]))[:10]
            records.append({
                "antibiotic": str(antibiotic), "cv_fold": int(fold), "n_train": int(len(train)),
                "n_test": int(len(test)), "n_features": len(vocabulary), "inner_ea_by_penalty": scores,
                "selected_penalty": best, "intercept_log2": float(theta[0]), "sigma_log2": sigma,
                "top_coefficients": [{"feature": f, "log2_effect": float(c)} for f, c in top],
            })
    columns = ["isolate_id", "antibiotic", "evaluation_split", "lineage_group", "intended_use_population",
               "cv_fold", "measurement_sign", "ast_value", "predicted_log2_mic", "predicted_log2_sd",
               "selected_penalty"]
    output = pd.concat(predictions, ignore_index=True)[columns]
    return output.sort_values(["antibiotic", "isolate_id"]).reset_index(drop=True), records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--population", choices=POPULATIONS, required=True)
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    seed = int(json.loads(args.study_config.read_text(encoding="utf-8"))["splitting"]["seed"])
    cohort = pd.read_csv(args.cohort, dtype={"isolate_id": str, "lineage_group": str}, low_memory=False)
    splits = pd.read_csv(args.splits, dtype={"isolate_id": str, "cv_group": str})
    features = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
    try:
        output, records = run_model(cohort, splits, features, args.population, seed)
    except ModelError as exc:
        raise SystemExit(f"Censored regression failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    payload = {
        "schema_version": "1.0.0",
        "operation": "run_censored_regression",
        "model": "ridge-penalised interval-censored normal regression on AMRFinderPlus AMR determinants",
        "generated_at_utc": utc_now(),
        "population": args.population,
        "nested_cv": {"penalty_grid": list(PENALTY_GRID), "inner_folds": INNER_FOLDS,
                      "min_feature_count": MIN_FEATURE_COUNT, "selection_metric": "inner essential agreement",
                      "seed": seed},
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in
                   (("cohort", args.cohort), ("splits", args.splits), ("features", args.features),
                    ("study_config", args.study_config))},
        "folds": records,
        "outputs": {"predictions": {"path": str(args.output), "sha256": sha256_file(args.output),
                                    "rows": int(len(output))}},
        "scientific_boundary": {"external_rows_read": False},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
