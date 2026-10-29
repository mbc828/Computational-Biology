# Public source snapshot: replace path placeholders before running.
from pathlib import Path
import json
import math
import csv
import numpy as np
from scipy import stats

ROOT = Path(r'<PROJECT_ROOT>')
SRC = ROOT / 'outputs/elisa_redraw_20261026/elisa_recalculated_20261026.json'
OUT = ROOT / 'outputs/manuscript_v21_github_release/review_v85/analysis'
OUT.mkdir(parents=True, exist_ok=True)
data = json.loads(SRC.read_text(encoding='utf-8'))

def values(analyte, group):
    return np.asarray(data[analyte]['plate']['groups'][group]['concentration'], dtype=float)

def mean_sd(x):
    return float(np.mean(x)), float(np.std(x, ddof=1)), len(x)

def contrast(analyte):
    # Additional baicalein effect in the IKK-16 background minus its effect without IKK-16.
    groups = [values(analyte, g) for g in (
        'HKSA+BAI(24)+IKK-16', 'HKSA+IKK-16', 'HKSA+BAI(24)', 'HKSA')]
    signs = [1.0, -1.0, -1.0, 1.0]
    est = sum(s * np.mean(x) for s, x in zip(signs, groups))
    se2 = sum(np.var(x, ddof=1) / len(x) for x in groups)
    se = math.sqrt(se2)
    df = se2 ** 2 / sum((np.var(x, ddof=1) / len(x)) ** 2 / (len(x) - 1) for x in groups)
    crit = stats.t.ppf(0.975, df)
    return {
        'estimate': float(est), 'se': se, 'df': float(df),
        'ci_low': float(est - crit * se), 'ci_high': float(est + crit * se),
        'p_welch': float(2 * stats.t.sf(abs(est / se), df)),
        'extra_baicalein_without_ikk16': float(np.mean(groups[2]) - np.mean(groups[3])),
        'extra_baicalein_with_ikk16': float(np.mean(groups[0]) - np.mean(groups[1])),
    }

def trend(analyte):
    xs, ys = [], []
    for dose in (6, 12, 24, 48):
        vals = values(analyte, f'HKSA+BAI({dose})')
        xs.extend([float(dose)] * len(vals))
        ys.extend(vals.tolist())
    x, y = np.asarray(xs), np.asarray(ys)
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    slope = np.sum((x - xm) * (y - ym)) / sxx
    intercept = ym - slope * xm
    resid = y - (intercept + slope * x)
    df = len(x) - 2
    mse = np.sum(resid ** 2) / df
    se = math.sqrt(mse / sxx)
    crit = stats.t.ppf(0.975, df)
    return {
        'slope_pg_ml_per_um': float(slope), 'se': se, 'df': df,
        'ci_low': float(slope - crit * se), 'ci_high': float(slope + crit * se),
        'p_linear_trend_unpaired': float(2 * stats.t.sf(abs(slope / se), df)),
        'n_measurements': len(x),
        'group_means': {str(d): mean_sd(values(analyte, f'HKSA+BAI({d})'))[0] for d in (6, 12, 24, 48)},
    }

results = {
    'source': str(SRC),
    'scope_note': 'Exploratory calculations treating each of five cell experiments as independent. Culture-batch IDs and explicit cross-treatment pairing were not present; no batch-adjusted or paired confirmatory model was fitted.',
    'contrast_definition': '(HKSA+BAI24+IKK16 - HKSA+IKK16) - (HKSA+BAI24 - HKSA). Positive values indicate attenuation of the additional baicalein effect in the IKK-16 background.',
    'analytes': {},
}
for analyte in ('IL-1β', 'IL-6', 'TNF-α'):
    results['analytes'][analyte] = {
        'interaction_contrast': contrast(analyte),
        'linear_trend_6_48um': trend(analyte),
        'group_summary': {g: dict(zip(('mean', 'sd', 'n'), mean_sd(values(analyte, g))))
                          for g in data[analyte]['plate']['groups']},
    }

(OUT / 'elisa_exploratory_statistics.json').write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
with (OUT / 'elisa_exploratory_statistics.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(['analyte', 'analysis', 'estimate', 'ci_low', 'ci_high', 'df', 'p_exploratory', 'unit'])
    for analyte, d in results['analytes'].items():
        c = d['interaction_contrast']
        w.writerow([analyte, 'IKK-16 background attenuation contrast', c['estimate'], c['ci_low'], c['ci_high'], c['df'], c['p_welch'], 'pg/mL'])
        t = d['linear_trend_6_48um']
        w.writerow([analyte, 'unpaired linear dose trend', t['slope_pg_ml_per_um'], t['ci_low'], t['ci_high'], t['df'], t['p_linear_trend_unpaired'], 'pg/mL per uM'])

for analyte, d in results['analytes'].items():
    c, t = d['interaction_contrast'], d['linear_trend_6_48um']
    print(analyte)
    print('  attenuation contrast = %.3f pg/mL (95%% CI %.3f to %.3f; df=%.2f; exploratory P=%.3g)' %
          (c['estimate'], c['ci_low'], c['ci_high'], c['df'], c['p_welch']))
    print('  linear trend = %.4f pg/mL per uM (95%% CI %.4f to %.4f; exploratory P=%.3g)' %
          (t['slope_pg_ml_per_um'], t['ci_low'], t['ci_high'], t['p_linear_trend_unpaired']))
print(results['scope_note'])
