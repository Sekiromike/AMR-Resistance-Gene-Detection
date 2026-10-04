"""Render complete comparison tables from the verified exports."""
from pathlib import Path
import csv

PAPER = Path(__file__).resolve().parents[1]


def build(tables):
    lines = ['---', 'title: "Supplementary results: curated resistance gene families and Evo 2"',
             'author: "Rushath Rajeev"', 'date: "3 October 2026"', 'documentclass: article',
             'fontsize: 10pt', 'geometry: margin=0.8in', '---', '',
             'All values below were checked against the original completed-run predictions. '
             'EA is essential agreement against an interval-censored MIC reference; EA exact '
             'uses exact references only. Brackets give 95% lineage-bootstrap intervals '
             '(2,000 draws, fixed predictions). Models with family features are post hoc. '
             'No model fitting was performed during manuscript verification.', '']
    names = {'constant': 'Constant', 'neighbor': 'Nearest relative',
             'ridge_censored': 'Ridge: allele', 'ridge_censored_amr': 'Ridge: allele',
             'ridge_censored_evo2': 'Ridge: Evo 2', 'ridge_censored_amr_evo2': 'Ridge: allele + Evo 2',
             'xgboost_aft': 'XGBoost: allele', 'xgboost_aft_amr': 'XGBoost: allele',
             'xgboost_aft_evo2': 'XGBoost: Evo 2', 'xgboost_aft_amr_evo2': 'XGBoost: allele + Evo 2',
             'ridge_censored_family': 'Ridge: allele + family + class*',
             'xgboost_aft_family': 'XGBoost: allele + family + class*'}
    for section, key in [('S1. Locked external evaluation and exploratory family models', 'full_external'),
                         ('S2. Final development out-of-fold evaluations', 'final_development')]:
        lines += ['# ' + section, '']
        for population in ('human_clinical', 'all_sources'):
            lines += ['## ' + population.replace('_', ' ').capitalize(), '']
            for drug in ('ceftriaxone', 'ciprofloxacin', 'gentamicin'):
                data = [r for r in tables[key] if r['population'] == population and r['drug'] == drug]
                lines += ['### ' + drug.capitalize(), '', f"All models: n = {data[0]['n']}; exact-reference n = {data[0]['n_exact']}.", '',
                          '| Model | EA [95% interval] | EA exact |', '|:---|:---|---:|']
                for r in sorted(data, key=lambda r: (r['evidence_tier'].endswith('posthoc') or r['evidence_tier']=='posthoc_external', r['model'])):
                    lines.append(f"| {names[r['model']]} | {r['ea']:.3f} [{r['ci_low']:.3f}, {r['ci_high']:.3f}] | {r['ea_exact']:.3f} |")
                lines += ['', '*Family models were developed after inspection of the original external results.*', '']
    lines += ['# S3. Variant withholding: all carriers and above-panel carriers', '',
              'For each target, all its development carriers were withheld. EA uses all held-out carriers; '
              'the >2 mg/L share uses right-censored references only. The four primary CTX-M targets '
              'contribute 568 target–isolate observations from 565 unique isolates and 45 lineages.', '']
    with (PAPER / 'verification/variant_results.csv').open() as handle:
        variants = list(csv.DictReader(handle))
    for target in sorted({r['target'] for r in variants}):
        lines += ['## ' + target, '', '| Learner / features | Held out | Above panel | EA | >2 mg/L share |',
                  '|:---|---:|---:|---:|---:|']
        for r in variants:
            if r['target'] != target:
                continue
            label = r['model'].replace('ridge_censored:', 'Ridge / ').replace('xgboost_aft:', 'XGBoost / ').replace('amr+evo2', 'allele + Evo 2').replace('family', 'allele + family + class')
            lines.append(f"| {label} | {r['n_test']} | {r['n_right_censored']} | {float(r['essential_agreement']):.3f} | {float(r['resistant_side_share']):.3f} |")
        lines += ['']
    lines += ['# S4. Numerical convergence and computational provenance', '',
              'The bounded ridge development rerun recorded the following nonconverged '
              'fits, counting inner selection and outer final fits together. These counts '
              'do not imply that every selected final fit failed to converge.', '',
              '| Population | Features | Nonconverged / total fits |', '|:---|:---|---:|']
    with (PAPER / 'verification/selected_settings.csv').open() as handle:
        settings = list(csv.DictReader(handle))
    for population in ('human_clinical', 'all_sources'):
        for feature, name in [('', 'Allele'), ('_evo2', 'Evo 2'), ('_amr_evo2', 'Allele + Evo 2')]:
            suffix = f'results/models/65006995/manifest_ridge_censored{feature}_{population}.json'
            selected = [r for r in settings if r['manifest'] == suffix]
            lines.append(f"| {population.replace('_', ' ')} | {name} | {sum(int(r['not_converged']) for r in selected)} / {sum(int(r['fits']) for r in selected)} |")
    lines += ['', 'For final external prediction (including development-fold selection), '
              'the corresponding counts were 7/318 for human-clinical allele + Evo 2, '
              '16/318 for all-source allele + Evo 2, and 3/318 for all-source Evo 2. '
              'All other external ridge configurations and both family ridge runs recorded '
              'zero nonconverged fits. Full per-fold settings and messages are in `selected_settings.csv`.', '',
              'Production embedding extraction used four completed H100 tasks, with '
              'recorded wall times totaling approximately 3.6 GPU-hours. Software versions '
              'and the 64 embedding-shard checksums accompany the release. This count '
              'excludes training the foundation model, CPU modelling, genome processing, '
              'engineering attempts, and transfers.', '',
              '# Machine-readable supplements', '',
              'The release includes `full_external.csv`, `final_development.csv`, '
              '`all_evaluations.csv` (including superseded development runs, identified by run path), '
              '`reference_and_lineage_strata.csv`, `paired_comparisons.csv`, '
              '`source_sensitivity.csv`, `variant_results.csv`, and `selected_settings.csv`. '
              'Source exclusions re-evaluate existing out-of-fold predictions without retraining. '
              'No clinical S/I/R, clinical very-major-error, independent replication, or '
              'laboratory-validation result is claimed.', '']
    (PAPER / 'supplement.md').write_text('\n'.join(lines))
