"""Post-hoc read-out of docs/POSTHOC_FAMILY_ANALYSIS.md (after scoring only).

Among external isolates with ceftriaxone right-censored above the panel,
reports for each scored prediction file the share placed on the resistant
side (> 2 mg/L), separately for carriers of the resistance variants that were
absent or near-absent from development (fixed list below, from
docs/EXTERNAL_RESULTS.md) and for carriers of common CTX-M variants.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scripts.evaluate_mic_predictions import score_rows
except ModuleNotFoundError:  # Direct execution
    from evaluate_mic_predictions import score_rows  # type: ignore[no-redef]

UNSEEN = ("blaCTX-M-2", "blaCTX-M-8", "blaCTX-M-3", "blaCTX-M-24", "blaCTX-M-65", "blaIMP-1")
COMMON = ("blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-27", "blaCTX-M-55")
RESISTANT_SIDE_LOG2 = 1.0  # > 2 mg/L


def read_out(scored: pd.DataFrame, features: pd.DataFrame) -> dict[str, float | int]:
    unseen = set(features.loc[features["element_symbol"].isin(UNSEEN), "isolate_id"])
    common = set(features.loc[features["element_symbol"].isin(COMMON), "isolate_id"]) - unseen
    s = score_rows(scored)
    high = s[(s["antibiotic"] == "ceftriaxone") & np.isposinf(s["ref_hi"])]
    a, b = high[high["isolate_id"].isin(unseen)], high[high["isolate_id"].isin(common)]
    return {"unseen_n": int(len(a)), "unseen_resistant_side": float((a["predicted_dilution"] > RESISTANT_SIDE_LOG2).mean()),
            "common_n": int(len(b)), "common_resistant_side": float((b["predicted_dilution"] > RESISTANT_SIDE_LOG2).mean())}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scored", nargs="+", required=True, help="name=scored_predictions.csv")
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    features = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
    result = {}
    for spec in args.scored:
        name, _, path = spec.partition("=")
        result[name] = read_out(pd.read_csv(path, dtype={"isolate_id": str, "lineage_group": str}), features)
        print(name, json.dumps(result[name]))
    args.output.write_text(json.dumps({"post_hoc": True, "unseen_variants": list(UNSEEN), "common_variants": list(COMMON),
                                       "results": result}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
