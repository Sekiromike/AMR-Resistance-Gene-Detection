"""Add curated gene-family and drug-class indicators to AMRFinderPlus features.

docs/POSTHOC_FAMILY_ANALYSIS.md (post hoc). Input and output share the
build_amr_features.py long format; each isolate gains one `family:<name>` row
per gene family and one `subclass:<name>` row per AMRFinderPlus Subclass,
alongside its unchanged allele-level rows. Genotype only.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

BETA_LACTAMASE = re.compile(r"^(bla[A-Za-z]+(?:-[A-Za-z]+)?)(?:-\d+[A-Za-z]*)?$")
TRAILING_ALLELE = re.compile(r"^(.*?[^\d])(\d+)$")
ROMAN_SUBCLASS = re.compile(r"^(.*-[IVX]+)[a-z]$")


def gene_family(symbol: str) -> str | None:
    """Family name for a gene symbol; None for point mutations."""
    if "_" in symbol:
        return None
    match = BETA_LACTAMASE.match(symbol)
    if match:
        return match.group(1)
    family = symbol
    trailing = TRAILING_ALLELE.match(family)
    if trailing:
        family = trailing.group(1)
    roman = ROMAN_SUBCLASS.match(family)
    if roman:
        family = roman.group(1)
    return family


def add_hierarchy(features: pd.DataFrame) -> pd.DataFrame:
    rows = [features]
    families = features.assign(family=features["element_symbol"].map(gene_family)).dropna(subset=["family"])
    rows.append(families.assign(element_symbol="family:" + families["family"])
                [features.columns].drop_duplicates(["isolate_id", "element_symbol"]))
    # AMRFinderPlus joins several subclasses with "/" (e.g. AMIKACIN/KANAMYCIN).
    subclasses = features[features["subclass"].str.strip().ne("")]
    expanded = subclasses.assign(item=subclasses["subclass"].str.split("/")).explode("item")
    expanded = expanded[expanded["item"].str.strip().ne("")]
    rows.append(expanded.assign(element_symbol="subclass:" + expanded["item"].str.strip())
                [features.columns].drop_duplicates(["isolate_id", "element_symbol"]))
    return pd.concat(rows, ignore_index=True).sort_values(["cohort", "isolate_id", "element_symbol"]).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    features = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
    out = add_hierarchy(features)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, lineterminator="\n")
    added = out["element_symbol"].str.extract(r"^(family|subclass):")[0].value_counts().to_dict()
    print({"input_rows": len(features), "output_rows": len(out), "added": added})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
