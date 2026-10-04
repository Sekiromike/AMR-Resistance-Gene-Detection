"""Build paper tables from verified, full-precision machine-readable results."""
import csv
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
DRUGS = ('ceftriaxone', 'ciprofloxacin', 'gentamicin')
LEARNERS = {'ridge': 'ridge_censored', 'xgboost': 'xgboost_aft'}
EXT = 'results/external/evaluation/65020832/'
FAMILY = 'results/posthoc/family/65022960/scoring-65023163/'


def rows(name):
    with (PAPER / 'verification' / (name + '.csv')).open() as handle:
        return list(csv.DictReader(handle))


def one(data, **match):
    found = [r for r in data if all(r.get(k) == v for k, v in match.items())]
    assert len(found) == 1, (match, len(found))
    return found[0]


def extract():
    audit = json.loads((PAPER / 'verification/audit.json').read_text())
    assert audit['status'] == 'PASS'
    evaluations, strata = rows('all_evaluations'), rows('reference_and_lineage_strata')
    comparisons, sensitivity = rows('paired_comparisons'), rows('source_sensitivity')
    tables = {'locked_external': [], 'posthoc_family_external': [], 'posthoc_external_variants': [],
              'leave_variant_out': [], 'leave_variant_out_hypotheses': [], 'cohort': [], 'gentamicin_gain': [],
              'full_external': [], 'final_development': []}
    for drug in DRUGS:
        for learner, full in LEARNERS.items():
            path = EXT + f'evaluation_predictions_{full}_amr_human_clinical.json'
            value = one(evaluations, evaluation=path, drug=drug)
            right = one(strata, evaluation=path, drug=drug, group_type='by_reference_class', group='right_censored')
            tables['locked_external'].append(dict(drug=drug, learner=learner,
                ea=float(value['essential_agreement']), ci_low=float(value['ci_low']), ci_high=float(value['ci_high']),
                undercall=float(right['under_call_2plus']), undercall_ci_low=float(right['under_call_2plus_95ci_low']),
                undercall_ci_high=float(right['under_call_2plus_95ci_high']), evidence_tier='locked_external', source=path))
            path = FAMILY + f'evaluation_{full}_family_human_clinical.json'
            value = one(evaluations, evaluation=path, drug=drug)
            diff = one(comparisons, evaluation=FAMILY + f'comparison_{full}_family_human_clinical.json', drug=drug, baseline=full + '_allele')
            tables['posthoc_family_external'].append(dict(drug=drug, learner=learner,
                ea=float(value['essential_agreement']), ci_low=float(value['ci_low']), ci_high=float(value['ci_high']),
                difference=float(diff['difference']), difference_ci_low=float(diff['ci_low']), difference_ci_high=float(diff['ci_high']),
                evidence_tier='posthoc_external', source=path))
    external = json.loads((PAPER / 'verification/external_variants.json').read_text())
    for name, r in external['results'].items():
        learner, feature = name.split('_', 1)
        tables['posthoc_external_variants'].append(dict(features='amr+evo2' if feature == 'amr_evo2' else feature,
            learner=learner, n=r['unseen_n'], share=r['unseen_resistant_side'], common_n=r['common_n'],
            common_share=r['common_resistant_side'], evidence_tier='posthoc_external', source=FAMILY + 'unseen_variants.json'))
    variant = json.loads((PAPER / 'verification/variant_readout.json').read_text())
    for target, models in variant['per_target'].items():
        if target == 'blaOXA-1':
            continue
        for name, r in models.items():
            learner, feature = name.split(':')
            tables['leave_variant_out'].append(dict(target=target, features=feature,
                learner='ridge' if learner == 'ridge_censored' else 'xgboost', n_above=r['n_right_censored'],
                n_test=r['n_test'], share=r['resistant_side_share'], ea=r['essential_agreement'],
                evidence_tier='development_variant_withholding', source='results/leave-variant-out/65051616/readout.json'))
    for learner, hypotheses in variant['hypotheses'].items():
        for name, r in hypotheses.items():
            tables['leave_variant_out_hypotheses'].append(dict(learner='ridge' if learner == 'ridge_censored' else 'xgboost',
                comparison=name.split('_', 1)[1], difference=r['difference'], ci_low=r['ci_low'], ci_high=r['ci_high'],
                n_target_isolate_observations=r['n'], n_unique_isolates=audit['variant_primary']['unique_isolates'],
                evidence_tier='development_variant_withholding', source='results/leave-variant-out/65051616/readout.json'))
    for r in audit['cohort']:
        tables['cohort'].append({**r, 'source': 'results/cohort/64948707/research_cohort.csv'})
    for learner, full in LEARNERS.items():
        path = f'results/comparisons/65016028/{learner}_amr_evo2_human_clinical.json'
        baseline = learner + '_amr'
        for stage, r in [
            ('Development', one(comparisons, evaluation=path, drug='gentamicin', baseline=baseline)),
            ('Development excluding source', one(sensitivity, evaluation=path, drug='gentamicin', baseline=baseline, source_excluded='PRJNA1297298')),
            ('Locked external', one(comparisons, evaluation=EXT + f'comparison_{full}_amr_evo2_human_clinical.json', drug='gentamicin', baseline=baseline))]:
            tables['gentamicin_gain'].append(dict(stage=stage, learner=learner, difference=float(r['difference']),
                ci_low=float(r['ci_low']), ci_high=float(r['ci_high']), source=r['evaluation']))
    for r in evaluations:
        path = r['evaluation']
        population = 'human_clinical' if 'human_clinical' in path else 'all_sources'
        if path.startswith(EXT) or path.startswith(FAMILY):
            table, tier = 'full_external', 'posthoc_external' if path.startswith(FAMILY) else 'locked_external'
            name = Path(path).stem.replace('evaluation_predictions_', '').replace('evaluation_', '')
        elif ('/65006995/' in path or '/64983779/evaluation_xgboost' in path
              or '/64960056/evaluation_constant_' in path or '/64960056/evaluation_neighbor_' in path
              or '/dev_evaluation_' in path):
            table, tier = 'final_development', 'development_posthoc' if '/dev_evaluation_' in path else 'development'
            name = Path(path).stem.replace('dev_evaluation_', '').replace('evaluation_', '')
        else:
            continue
        name = name.removesuffix('_' + population)
        tables[table].append(dict(model=name, population=population, drug=r['drug'], n=int(r['n']),
            ea=float(r['essential_agreement']), ci_low=float(r['ci_low']), ci_high=float(r['ci_high']),
            n_exact=int(r['n_exact_reference']), ea_exact=float(r['essential_agreement_exact_only']),
            mean_log_likelihood=float(r['mean_interval_log_likelihood']) if r['mean_interval_log_likelihood'] else '',
            evidence_tier=tier, source=path))
    assert len(tables['full_external']) == 54 and len(tables['final_development']) == 54
    return tables
