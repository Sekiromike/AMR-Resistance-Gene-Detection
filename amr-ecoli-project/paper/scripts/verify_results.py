"""Verify original completed runs and export full-precision supplementary tables.

Recomputes scores of frozen predictions only. Does not fit or select models.
Run from the project directory after extracting the release assets.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import scipy
from scripts.evaluate_mic_predictions import evaluate, score_rows
from scripts.compare_mic_models import paired_difference, source_robustness
from scripts.leave_variant_out import readout, resistant_side, PRIMARY_TARGETS
from scripts.unseen_variant_analysis import read_out as unseen_readout


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def compare(expected, actual, context=''):
    """Check every originally recorded value; allow newly added reporting keys."""
    if isinstance(expected, dict):
        for key, value in expected.items():
            assert key in actual, (context, 'missing', key)
            compare(value, actual[key], context + '/' + key)
    elif isinstance(expected, list):
        assert len(expected) == len(actual), context
        for i, (x, y) in enumerate(zip(expected, actual)):
            compare(x, y, context + '/' + str(i))
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        assert (math.isnan(expected) and math.isnan(actual)) or math.isclose(expected, actual, rel_tol=1e-10, abs_tol=1e-10), (context, expected, actual)
    else:
        assert expected == actual, (context, expected, actual)


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', type=Path, default=ROOT / 'paper/artifacts')
    parser.add_argument('--output', type=Path, default=ROOT / 'paper/verification')
    args = parser.parse_args()
    artifacts = args.artifacts.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'status': 'RUNNING', 'model_fits': 0, 'external_selection': False,
              'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
              'hashes_checked': [], 'unavailable_hash_targets': [], 'prediction_locks': [],
              'evaluations': [], 'comparisons': []}
    files = {str(p.relative_to(artifacts)): p for p in artifacts.rglob('*') if p.is_file()}
    release_path = artifacts / 'release_manifest.json'
    distributed = {}
    if release_path.exists():
        distributed = {r['file']: r for r in json.loads(release_path.read_text())['files']}
        for name, entry in distributed.items():
            assert (artifacts / name).is_file() and sha(artifacts / name) == entry['published_sha256'], ('distributed checksum mismatch', name)
    def original_hash_matches(path, expected):
        actual = sha(path)
        if actual == expected:
            return True
        relative = str(path.relative_to(artifacts)) if path.is_relative_to(artifacts) else ''
        entry = distributed.get(relative, {})
        return entry.get('redacted') is True and entry.get('original_sha256') == expected and entry.get('published_sha256') == actual


    def resolve(name, parent=None):
        name = str(name).replace('${PROJECT_ROOT}/', '')
        if name.startswith('./'):
            name = name[2:]
        for prefix in ('results/', 'provenance/', 'workflow/'):
            if prefix in name:
                name = prefix + name.split(prefix, 1)[1]
                break
        candidates = [artifacts / name]
        if parent:
            candidates.insert(0, parent / name)
        candidates.append(ROOT / name)
        return next((p for p in candidates if p.is_file()), None)

    # Validate recorded artifact digests, including all embedding shards.
    for name, manifest in sorted(files.items()):
        if '.sha256' not in manifest.name:
            continue
        for line in manifest.read_text().splitlines():
            fields = line.split(maxsplit=1)
            if len(fields) != 2 or len(fields[0]) != 64:
                continue
            digest, target = fields
            path = resolve(target.lstrip('*'), manifest.parent)
            if path is None:
                report['unavailable_hash_targets'].append({'manifest': name, 'target': target.split('amr-ecoli-project/')[-1]})
                continue
            actual = sha(path)
            assert original_hash_matches(path, digest), ('artifact hash mismatch', name, target, digest, actual)
            report['hashes_checked'].append({'manifest': name, 'file': str(path.relative_to(artifacts)) if path.is_relative_to(artifacts) else str(path.relative_to(ROOT)), 'sha256': actual})
    for lock in ('external_prediction_lock.tsv', 'external_posthoc_lock.tsv'):
        with (ROOT / 'config' / lock).open() as handle:
            entries = list(csv.DictReader(handle, delimiter='\t'))
        for entry in entries:
            lock_parent = artifacts / ('results/external/predictions/65019273' if lock == 'external_prediction_lock.tsv' else 'results/posthoc/family/65022960')
            path = resolve(entry['file'], lock_parent)
            assert path is not None and sha(path) == entry['sha256'], ('prediction lock mismatch', entry)
            report['prediction_locks'].append({'lock': lock, **entry})
    assert len(report['prediction_locks']) == 18
    print('PASS: artifact hashes and 18 committed prediction locks', flush=True)

    config = json.loads((ROOT / 'config/study.json').read_text())['evaluation']
    iterations, seed = int(config['cluster_bootstrap_iterations']), int(config['seed'])
    cohort = pd.read_csv(artifacts / 'results/cohort/64948707/research_cohort.csv', low_memory=False)
    assert len(cohort) == 29619 and cohort.isolate_id.nunique() == 11520
    sources = cohort.drop_duplicates('isolate_id').set_index('isolate_id')['bioproject_accession']
    cohorts = []
    for population, part in [('all_source_development', cohort[cohort.evaluation_split.eq('development')]),
                             ('human_clinical_development', cohort[cohort.evaluation_split.eq('development') & cohort.intended_use_population.eq('human_clinical')]),
                             ('external', cohort[cohort.evaluation_split.eq('external')])]:
        cohorts.append({'population': population, 'isolates': part.isolate_id.nunique(), **part.antibiotic.value_counts().to_dict()})
    report['cohort'] = cohorts
    write_csv(output / 'cohort.csv', cohorts)
    scores, strata_rows, comparison_rows, source_rows = [], [], [], []
    for name, path in sorted(files.items()):
        if path.suffix != '.json' or not (path.name.startswith(('evaluation', 'dev_evaluation', 'comparison')) or '/comparisons/' in name):
            continue
        data = json.loads(path.read_text())
        if data.get('operation') != 'evaluate_mic_predictions':
            continue
        prediction_path = resolve(data['inputs']['predictions']['path'])
        assert prediction_path is not None, name
        assert original_hash_matches(prediction_path, data['inputs']['predictions']['sha256']), name
        frame = pd.read_csv(prediction_path, dtype={'isolate_id': str, 'lineage_group': str})
        verified = evaluate(frame, iterations=iterations, seed=seed, locked_external=data['locked_external'], dominant_lineage=config['dominant_lineage'])['per_antibiotic']
        compare(data['per_antibiotic'], verified, name)
        # Independently check EA with scalar interval arithmetic.
        for drug, part in frame.groupby('antibiotic'):
            success = 0
            for row in part.itertuples():
                d = round(math.log2(row.ast_value))
                p = math.floor(row.predicted_log2_mic + .5)
                bounds = {'=': (d,d), '==': (d,d), '<=': (-math.inf,d), '<': (-math.inf,d-1), '>=': (d,math.inf), '>': (d+1,math.inf)}
                lo, hi = bounds[row.measurement_sign]
                success += max(lo-p, p-hi, 0) <= 1
            compare(success / len(part), verified[drug]['essential_agreement'], name + '/independent_EA')
        report['evaluations'].append({'file': name, 'predictions': str(prediction_path.relative_to(artifacts)), 'rows': len(frame), 'sha256': sha(path)})
        for drug, metrics in verified.items():
            row = {'evaluation': name, 'drug': drug, 'locked_external': data['locked_external']}
            row.update({k: v for k, v in metrics.items() if not isinstance(v, (dict, list))})
            row.update(ci_low=metrics['bootstrap_95ci']['essential_agreement'][0], ci_high=metrics['bootstrap_95ci']['essential_agreement'][1])
            scores.append(row)
            for kind in ('by_reference_class', 'strata'):
                for group, values in metrics[kind].items():
                    stratum = {'evaluation': name, 'drug': drug, 'group_type': kind, 'group': group, **{k: v for k,v in values.items() if not isinstance(v, (dict, list))}}
                    for k, v in values.items():
                        if isinstance(v, list):
                            stratum[k + '_low'], stratum[k + '_high'] = v
                    for k, v in values.get('bootstrap_95ci', {}).items():
                        if isinstance(v, list):
                            stratum[k + '_95ci_low'], stratum[k + '_95ci_high'] = v
                    strata_rows.append(stratum)
        print('Verified', name, flush=True)

    for name, path in sorted(files.items()):
        if path.suffix != '.json' or not (path.name.startswith('comparison') or '/comparisons/' in name):
            continue
        data = json.loads(path.read_text())
        if data.get('operation') != 'compare_mic_models':
            continue
        model_name = data['model']
        model_path = resolve(data['inputs'][model_name]['path'])
        assert model_path is not None, name
        model = pd.read_csv(model_path)
        for drug, detail in data['per_antibiotic'].items():
            m = model[model.antibiotic.eq(drug)]
            for baseline, expected in detail['comparisons'].items():
                bp = resolve(data['inputs'][baseline]['path'])
                assert bp is not None, (name, baseline)
                b = pd.read_csv(bp)
                b = b[b.antibiotic.eq(drug)]
                actual = paired_difference(m, b, iterations, seed)
                compare(expected, actual, name + '/' + drug + '/' + baseline)
                comparison_rows.append({'evaluation': name, 'drug': drug, 'model': model_name, 'baseline': baseline, **{k: v for k,v in actual.items() if not isinstance(v, list)}, 'ci_low': actual['difference_95ci'][0], 'ci_high': actual['difference_95ci'][1]})
                if baseline in detail.get('source_robustness', {}):
                    computed = source_robustness(m, b, sources, iterations, seed)
                    compare(detail['source_robustness'][baseline], computed, name + '/source/' + baseline)
                    for project, v in computed.items():
                        source_rows.append({'evaluation': name, 'drug': drug, 'model': model_name, 'baseline': baseline, 'source_excluded': project, **{k:x for k,x in v.items() if not isinstance(x,list)}, 'ci_low': v.get('difference_95ci',['',''])[0], 'ci_high': v.get('difference_95ci',['',''])[1]})
        report['comparisons'].append({'file': name, 'sha256': sha(path)})

    family_dir = artifacts / 'results/posthoc/family/65022960/scoring-65023163'
    unseen = json.loads((family_dir / 'unseen_variants.json').read_text())
    features = pd.read_csv(artifacts / 'results/models/64960056/amr_features.tsv', sep='\t')
    for name, expected in unseen['results'].items():
        learner, feature = name.split('_', 1)
        learner = 'ridge_censored' if learner == 'ridge' else 'xgboost_aft'
        if feature == 'family':
            prediction_path = family_dir / f'scored_predictions_{learner}_family_human_clinical.csv'
        else:
            feature = 'amr' if feature == 'allele' else feature
            prediction_path = artifacts / f'results/external/evaluation/65020832/scored_predictions_{learner}_{feature}_human_clinical.csv'
        compare(expected, unseen_readout(pd.read_csv(prediction_path), features), 'unseen/' + name)
    (output / 'external_variants.json').write_text(json.dumps(unseen, indent=2) + '\n')
    report['external_variant_models_verified'] = len(unseen['results'])

    variant_dir = artifacts / 'results/leave-variant-out/65051616'
    variants = pd.concat([pd.read_csv(p) for p in sorted(variant_dir.glob('predictions_*.csv'))], ignore_index=True)
    variant_results = readout(variants, iterations, seed)
    compare(json.loads((variant_dir / 'readout.json').read_text()), variant_results, 'variant_readout')
    high = resistant_side(variants[variants.model.eq('ridge_censored:family')])
    primary = high[high.target.isin(PRIMARY_TARGETS)]
    report['variant_primary'] = {'target_isolate_observations': len(primary), 'unique_isolates': primary.isolate_id.nunique(), 'lineages': primary.lineage_group.nunique()}
    write_csv(output / 'variant_results.csv', [{'target': target, 'model': model, **value} for target, entries in variant_results['per_target'].items() for model, value in entries.items()])
    (output / 'variant_readout.json').write_text(json.dumps(variant_results, indent=2) + '\n')
    write_csv(output / 'all_evaluations.csv', scores)
    write_csv(output / 'reference_and_lineage_strata.csv', strata_rows)
    write_csv(output / 'paired_comparisons.csv', comparison_rows)
    write_csv(output / 'source_sensitivity.csv', source_rows)
    # Original binary/content hashes are retained even if public manifests redact private paths.
    write_csv(output / 'original_artifacts.csv', [{'file': name, 'bytes': p.stat().st_size, 'sha256': sha(p)} for name, p in sorted(files.items())])
    selected_settings = []
    for path in sorted(artifacts.rglob('*manifest*.json')):
        data = json.loads(path.read_text())
        if not isinstance(data, dict):
            continue
        records = data.get('folds', data.get('per_antibiotic', []))
        if not isinstance(records, list):
            continue
        for row in records:
            optimiser = row.get('optimiser', {})
            selected_settings.append({'manifest': str(path.relative_to(artifacts)), 'drug': row.get('antibiotic', ''), 'fold': row.get('cv_fold', 'full_development'), 'selected': json.dumps(row.get('selected', {}), sort_keys=True), 'fits': optimiser.get('fits', ''), 'not_converged': optimiser.get('not_converged', ''), 'n_train': row.get('n_train', ''), 'n_test': row.get('n_test', row.get('n_external', '')), 'n_features': row.get('n_features', ''), 'messages': json.dumps(optimiser.get('not_converged_messages', []))})
    write_csv(output / 'selected_settings.csv', selected_settings)
    report['status'] = 'PASS'
    report['counts'] = {'evaluations': len(report['evaluations']), 'comparison_files': len(report['comparisons']), 'score_rows': len(scores), 'paired_comparison_rows': len(comparison_rows), 'source_exclusion_rows': len(source_rows), 'variant_rows': len(variants), 'artifact_files': len(files)}
    (output / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], **report['counts'], 'variant_primary': report['variant_primary']}, indent=2))


if __name__ == '__main__':
    main()
