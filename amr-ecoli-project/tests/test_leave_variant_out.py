from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from scripts import build_hierarchical_features as HIER
from scripts import leave_variant_out as MODULE


def family_world(n: int = 600, seed: int = 3):
    """Ceftriaxone driven by any blaX allele; blaX-3 is the held-out variant."""
    rng = np.random.default_rng(seed)
    rows, feats = [], []
    for i in range(n):
        allele = rng.choice(["none", "blaX-1", "blaX-2", "blaX-3"], p=[0.55, 0.2, 0.15, 0.1])
        high = allele != "none"
        rows.append({"isolate_id": f"I{i}", "antibiotic": "ceftriaxone", "evaluation_split": "development",
                     "lineage_group": f"SLV:ST{i % 40}", "intended_use_population": "human_clinical",
                     "measurement_sign": ">" if high else "<=", "ast_value": "64" if high else "1"})
        if allele != "none":
            feats.append({"isolate_id": f"I{i}", "cohort": "development", "element_symbol": allele, "subtype": "AMR",
                          "class": "BETA-LACTAM", "subclass": "CEPHALOSPORIN", "method": "EXACTX"})
    cohort = pd.DataFrame(rows)
    splits = pd.DataFrame({"isolate_id": cohort["isolate_id"], "evaluation_split": "development",
                           "cv_group": cohort["lineage_group"], "cv_fold": [(i % 40) % 5 for i in range(n)]})
    allele = pd.DataFrame(feats)
    return cohort, splits, allele, HIER.add_hierarchy(allele)


class LeaveVariantOutTests(unittest.TestCase):
    def test_held_out_allele_never_trains_and_family_grouping_recovers_it(self) -> None:
        cohort, splits, allele, family = family_world()
        train, test = MODULE.split_on_allele(cohort, allele, "blaX-3")
        self.assertFalse(set(train.isolate_id) & set(allele[allele.element_symbol == "blaX-3"].isolate_id))
        predictions, records = MODULE.run_target("blaX-3", cohort, splits, allele, family, None, seed=7,
                                                 feature_sets=("allele", "family"))
        self.assertEqual({r["n_test"] for r in records}, {len(test)})
        detected = {m: MODULE.resistant_side(p)["detected"].mean() for m, p in predictions.groupby("model")}
        self.assertLess(detected["ridge_censored:allele"], 0.2)
        self.assertGreater(detected["ridge_censored:family"], 0.9)
        self.assertGreater(detected["xgboost_aft:family"], 0.9)

    def test_readout_handles_isolates_carrying_two_targets(self) -> None:
        # Regression: read-out job 65051617 refused isolates counted under two targets.
        cohort, splits, allele, family = family_world()
        both = allele[allele.element_symbol == "blaX-3"].assign(element_symbol="blaY-1", subclass="BETA-LACTAM")
        allele2 = pd.concat([allele, both], ignore_index=True)
        family2 = HIER.add_hierarchy(allele2)
        frames = [MODULE.run_target(t, cohort, splits, allele2, family2, None, seed=7,
                                    feature_sets=("allele", "family"))[0] for t in ("blaX-3", "blaY-1")]
        predictions = pd.concat(frames, ignore_index=True)
        fam = predictions[predictions["model"].str.endswith(":family")]
        predictions = pd.concat([predictions, fam.assign(model=fam["model"].str.replace("family", "evo2"))],
                                ignore_index=True)
        result = MODULE.readout(predictions, iterations=100, seed=1)
        self.assertEqual(set(result["per_target"]), {"blaX-3", "blaY-1"})

    def test_paired_bootstrap_on_identical_rows(self) -> None:
        a = pd.DataFrame({"target": "t", "isolate_id": [f"I{i}" for i in range(20)],
                          "lineage_group": [f"L{i % 5}" for i in range(20)], "detected": 1.0})
        b = a.assign(detected=[1.0] * 10 + [0.0] * 10)
        result = MODULE.paired_bootstrap(a, b, iterations=200, seed=1)
        self.assertAlmostEqual(result["difference"], 0.5)
        self.assertLessEqual(result["ci_low"], 0.5)
        self.assertGreaterEqual(result["ci_high"], 0.5)
