# Public source snapshot: replace path placeholders before running.
#!/usr/bin/env python3
"""只读整合结构结果；不运行轨迹计算或MMGBSA。所有新文件仅写入指定产物目录。"""
from __future__ import annotations
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import datetime

sys.dont_write_bytecode = True
BASE = Path('<PROJECT_ROOT>')
PARENT = BASE / 'outputs'
OLD = PARENT / 'MD分析汇总_20260920'
REPAIR = Path('<MD_DATA_ROOT>/postprocess_repair_20260921_001533_74632')
LEGACY = Path('<MD_DATA_ROOT>/analysis_all_completed')
OUT = PARENT / '综合分析汇总_20260921' / 'structural'
if not PARENT.is_dir():
    raise RuntimeError(f'父目录不存在：{PARENT}')
OUT.mkdir(parents=True, exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(OUT / 'matplotlib_cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Microsoft YaHei', 'axes.unicode_minus': False,
                     'font.size': 9, 'pdf.fonttype': 42})
MAPDIR = OUT / '结构残基映射'
MAPDIR.mkdir(exist_ok=True)
PARSER = BASE / 'MD_SHXD_PJI/05_analysis/parse_contact_occupancy.py'
FITS = ['protein_bb', 'TLR4_bb', 'MD2_bb', 'pocket_bb']
IDENT = ['chain', 'input_resnr', 'insertion_code', 'resname', 'moltype_name']
COLORS = ['#187b9d', '#db803a', '#7853a2']
issues, provenance, bindings, numeric, identities, qc_rows = [], {}, [], [], [], []
long_rows, residue_rmsf, residue_contact, blocks, drifts, trajectories = [], [], [], [], [], []
series = {}
fit_rows = []
reference_rows = []
pocket_selection = {}
contact_time_checks, contact_plot_checks = [], []


def issue(code, context, detail):
    issues.append({'code': code, 'context': str(context), 'detail': str(detail)})


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def register(p, role='结构输入'):
    p = Path(p)
    key = str(p)
    if key not in provenance:
        if not p.is_file():
            provenance[key] = {'path': key, 'role': role, 'status': 'MISSING'}
            issue('MISSING_INPUT', p, role)
        else:
            st = p.stat()
            skip = p.suffix.lower() in {'.xtc', '.trr', '.tng'} or st.st_size > 64 * 1024 * 1024
            provenance[key] = {'path': key, 'role': role, 'size': st.st_size,
                               'mtime_ns': st.st_mtime_ns, 'sha256': None if skip else sha(p),
                               'status': 'SKIPPED_LARGE_NOT_REHASHED' if skip else 'HASHED_NOW'}
    return provenance[key]


def read_json(p):
    register(p)
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def read_tsv(p):
    register(p)
    with Path(p).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(name, rows, keys=None):
    rows = list(rows)
    keys = keys or list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: 'NA' if r.get(k) is None else r.get(k, '') for k in keys})


def write_json(name, d):
    (OUT / name).write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def winpath(s):
    if s.startswith('/mnt/'):
        return Path(s[5].upper() + ':' + s[6:])
    return Path(s)


def num(x):
    if x is None or str(x).strip().lower() in {'', 'na', 'nan', 'none'}:
        return None
    v = float(x)
    if not math.isfinite(v):
        raise ValueError(f'非有限值 {x}')
    return v


def stats(v):
    v = [float(x) for x in v if x is not None]
    return {'mean': statistics.mean(v) if v else None,
            'sd_between_replicates': statistics.stdev(v) if len(v) > 1 else None,
            'n_replicates': len(v)}


def idkey(r):
    return tuple(str(r[k]).strip() for k in IDENT)


def check_binding(p, expected, evidence, size=None, mtime=None):
    p = Path(p)
    current = register(p, '历史绑定输入')
    status = ('MISSING_FILE' if current['status'] == 'MISSING' else
              'MISSING_EXPECTED_HASH' if not expected else
              'SKIPPED_LARGE_NOT_REHASHED' if current['sha256'] is None else
              'MATCH' if current['sha256'] == expected else 'HASH_MISMATCH')
    record = {'evidence': str(evidence), 'path': str(p), 'expected_sha256': expected,
              'actual_sha256': current.get('sha256'), 'status': status,
              'size_match': None if size is None else current.get('size') == size,
              'mtime_match': None if mtime is None else current.get('mtime_ns') == mtime}
    bindings.append(record)
    if status in {'MISSING_FILE', 'MISSING_EXPECTED_HASH', 'HASH_MISMATCH'}:
        issue(status, evidence, p)
    if record['size_match'] is False:
        issue('SIZE_MISMATCH', evidence, p)
    return record


def walk_bindings(obj, evidence):
    if isinstance(obj, dict):
        if 'path' in obj:
            check_binding(winpath(obj['path']), obj.get('sha256'), evidence,
                          obj.get('size'), obj.get('mtime_ns'))
        for k, v in obj.items():
            if k.endswith('_sha256') and isinstance(v, str) and k[:-7] in obj:
                check_binding(winpath(obj[k[:-7]]), v, evidence)
            elif isinstance(v, (dict, list)):
                walk_bindings(v, evidence)
    elif isinstance(obj, list):
        for v in obj:
            walk_bindings(v, evidence)


def xvg(p, default='ps'):
    register(p)
    text = Path(p).read_text(encoding='utf-8', errors='replace')
    a = np.array([[float(x) for x in l.split()] for l in text.splitlines()
                  if l.strip() and not l.lstrip().startswith(('#', '@'))])
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError(f'XVG格式或数值错误：{p}')
    unit = 'ns' if 'Time (ns)' in text else default
    return a, unit, text


def metadata(system, rep, source):
    return {'system': system, 'target': system.split('_')[1], 'ligand': system.split('_')[2],
            'rep': str(rep), 'trajectory_id': f'{system}/rep{rep}', 'source': source}


def metric_row(meta, metric, s, source_file, scope, unit='nm', fit='none', note='', qc=''):
    s = s or {}
    r = dict(meta, metric=metric, definition_scope=scope, fit=fit, unit=unit,
             window_start_ns=50, window_end_ns=100, mean=num(s.get('mean')),
             sd_within=num(s.get('sd', s.get('sd_population', s.get('sd_within')))),
             within_sd_definition='population_ddof0' if 'sd_population' in s else
                 ('source1_preserved' if 'sd_within' in s else 'sample_ddof1'),
             min=num(s.get('min')), max=num(s.get('max')), n_frames=s.get('n', s.get('n_frames', 0)),
             missing_reason='' if s.get('mean') is not None else ('apo不适用' if meta['ligand'] == 'apo' else '源缺失'),
             source_file=str(source_file), qc_layer=qc, note=note)
    long_rows.append(r)
    return r


def add_series(meta, metric, t, y, expected, p):
    t, y = np.asarray(t, dtype=float).copy(), np.asarray(y, dtype=float)
    if len(t) != len(y) or not np.isfinite(t).all() or not np.isfinite(y).all():
        raise ValueError(f'时间序列非有限或维数不一致：{p} {metric}')
    # 旧XVG ns输出在整十ns处有约1e-5的格式误差，只吸附这些已知边界。
    for edge in [50, 60, 70, 80, 90, 100]:
        t[np.abs(t - edge) <= 1e-4] = edge
    mask = (t >= 50) & (t <= 100)
    tw, yw = t[mask], y[mask]
    if len(tw) == 0 or np.any(np.diff(tw) <= 0):
        raise ValueError(f'窗口为空或时间不递增：{p} {metric}')
    mu, n = float(yw.mean()), len(yw)
    within_sample = float(yw.std(ddof=1)) if n > 1 else None
    within_population = float(yw.std(ddof=0))
    expected_mean = expected.get('mean') if expected else None
    expected_n = expected.get('n_frames') if expected else None
    tolerance = max(1e-5, abs(mu) * 2e-7)
    ok = (expected_mean is None or abs(mu - expected_mean) <= tolerance) and (not expected_n or n == int(expected_n))
    source_sd = expected.get('sd_within') if expected else None
    sd_definition = expected.get('within_sd_definition') if expected else None
    sd_recomputed = within_population if sd_definition == 'population_ddof0' else within_sample
    sd_ok = source_sd is None or abs(source_sd - sd_recomputed) <= tolerance
    if sd_definition == 'source1_preserved':
        sd_ok = source_sd is None or min(abs(source_sd - within_sample), abs(source_sd - within_population)) <= tolerance
    minmax_ok = all(expected.get(k) is None or abs(expected[k] - float(f(yw))) <= tolerance
                    for k, f in [('min', np.min), ('max', np.max)]) if expected else True
    ok = ok and sd_ok and minmax_ok
    numeric.append(dict(meta, metric=metric, source_file=str(p), n_frames=n,
                        first_ns=float(tw[0]), last_ns=float(tw[-1]), mean_recomputed=mu,
                        sd_sample_recomputed=within_sample, sd_population_recomputed=within_population,
                        sd_source=source_sd, sd_check=sd_ok, minmax_check=minmax_ok,
                        mean_source=expected_mean, absolute_difference=None if expected_mean is None else abs(mu - expected_mean),
                        mean_tolerance=tolerance, status='PASS' if ok else 'FAIL'))
    if not ok:
        issue('NUMERIC_MISMATCH', meta['trajectory_id'], metric)
    series[(meta['system'], meta['rep'], metric)] = (tw, yw)
    br = []
    for i, start in enumerate(range(50, 100, 10)):
        sel = (tw >= start) & ((tw < start + 10) if i < 4 else (tw <= 100))
        v, ts = yw[sel], tw[sel]
        r = dict(meta, metric=metric, block=i + 1, start_ns=start, end_ns=start + 10,
                 right_closed=i == 4, n_frames=len(v), first_observed_ns=float(ts[0]) if len(v) else None,
                 last_observed_ns=float(ts[-1]) if len(v) else None,
                 mean=float(v.mean()) if len(v) else None,
                 sd_within_block_sample=float(v.std(ddof=1)) if len(v) > 1 else None,
                 note='相关时间块，仅描述漂移；不是独立重复')
        blocks.append(r)
        br.append(r)
    assert sum(r['n_frames'] for r in br) == n
    vals = [r['mean'] for r in br]
    complete = all(v is not None for v in vals)
    drifts.append(dict(meta, metric=metric, last_minus_first=vals[-1] - vals[0] if complete else None,
                       range_of_block_means=max(vals) - min(vals) if complete else None,
                       slope_per_ns=float(np.polyfit([55, 65, 75, 85, 95], vals, 1)[0]) if complete else None,
                       last_observed_ns=float(tw[-1]), n_blocks_observed=sum(v is not None for v in vals)))


def ndx_groups(p):
    register(p)
    groups, current = {}, None
    for line in Path(p).read_text().splitlines():
        line = line.split(';', 1)[0].strip()
        if not line:
            continue
        if line.startswith('['):
            current = line.strip('[] ')
            groups[current] = []
        else:
            groups[current].extend(map(int, line.split()))
    return groups


def get_old_mapping(meta, per, run, parser):
    tpr = run / 'md.tpr'
    tpr_hash = register(tpr, '旧轨迹拓扑')['sha256']
    cache = MAPDIR / (meta['system'] + '_rep' + meta['rep'] + '_protein_atoms.tsv')
    stamp = cache.with_suffix('.json')
    parser_hash = register(PARSER, '只读拓扑解析器')['sha256']
    groups = ndx_groups(per / 'analysis.ndx')
    expected = {'tpr_sha256': tpr_hash, 'parser_sha256': parser_hash,
                'analysis_ndx_sha256': register(per / 'analysis.ndx')['sha256']}
    if cache.is_file() and stamp.is_file():
        prior = json.loads(stamp.read_text(encoding='utf-8'))
        if all(prior.get(k) == v for k, v in expected.items()) and prior.get('mapping_sha256') == sha(cache):
            register(stamp, '新建映射绑定')
            return read_tsv(cache), groups
    wsl_tpr = '/mnt/' + tpr.drive[0].lower() + tpr.as_posix()[2:]
    cmd = ['wsl.exe', '-d', 'Ubuntu-24.04', '--cd', '/tmp', '--', 'env', 'OMP_NUM_THREADS=1',
           '/opt/gromacs-2026.2/bin/gmx', 'dump', '-s', wsl_tpr]
    result = subprocess.run(cmd, capture_output=True, encoding='utf-8', errors='replace', check=True)
    mapping = parser.parse_topology_dump_text(result.stdout)
    rows = [asdict(mapping.atoms[a]) for a in groups['Protein']]
    write_tsv(cache.relative_to(OUT), rows)
    stamp.write_text(json.dumps(dict(expected, source_tpr=str(tpr), mapping_sha256=sha(cache),
                                    command=cmd, identity_method='逐rep原始TPR展开，与Protein及C-alpha索引对应'),
                               ensure_ascii=False, indent=2), encoding='utf-8')
    register(stamp, '新建映射绑定')
    register(cache, '新建蛋白原子映射')
    return rows, groups


def check_contact_times(times, backbone_times, n_frames):
    times = np.asarray(times, dtype=float)
    backbone_times = np.asarray(backbone_times, dtype=float)
    assert times.ndim == backbone_times.ndim == 1 and len(times) == len(backbone_times) == n_frames > 0, '接触与骨架帧数不一致'
    assert np.isfinite(times).all() and np.isfinite(backbone_times).all(), '接触或骨架时间非有限'
    assert np.all(np.diff(times) > 0) and np.all(np.diff(backbone_times) > 0), '时间不严格递增'
    assert times.min() >= 50 - 1e-4 and times.max() <= 100 + 1e-4, '接触时间超出窗口'
    # 允许旧XVG时间格式精度误差，不使用随时间放大的相对容差。
    assert np.allclose(times, backbone_times, rtol=0, atol=1e-4), '接触NDX与该重复骨架实际帧时间不匹配'
    return float(np.max(np.abs(times - backbone_times)))


def old_residues(meta, per, run, parser, rmsf_description):
    atoms, groups = get_old_mapping(meta, per, run, parser)
    amap = {int(r['atom_index1']): r for r in atoms}
    residues = {int(r['global_residue_index0']): r for r in atoms}
    keys = [idkey(r) for r in residues.values()]
    if len(keys) != len(set(keys)):
        raise ValueError(f'残基身份不唯一：{meta}')
    a, _, text = xvg(per / 'rmsf_ca_steady.xvg')
    ca = [amap[i] for i in groups['C-alpha']]
    if len(a) != len(ca) or any(int(x[0]) != int(r['input_resnr']) for x, r in zip(a, ca)):
        raise ValueError(f'RMSF残基轴与CA映射不一致：{meta}')
    for x, r in zip(a, ca):
        residue_rmsf.append(dict(meta, **{k: str(r[k]).strip() for k in IDENT},
                                chain_source=r['chain_source'], fit='legacy_CA_gmx_default', rmsf_nm=float(x[1]),
                                source_file=str(per / 'rmsf_ca_steady.xvg'), window_ns='50-100',
                                fit_semantics='GROMACS rmsf默认拟合/所选CA；与修复骨架拟合分开'))
    assert abs(float(a[:, 1].mean()) - float(rmsf_description['mean_across_positions'])) < 1e-7
    identities.append(dict(meta, kind='RMSF', n_residues=len(ca), status='PASS',
                           method='逐rep TPR展开+CA原子索引+原始残基轴逐项匹配'))
    if meta['ligand'] == 'apo':
        return
    p = per / 'contact_residues_dynamic.ndx'
    register(p)
    text = p.read_text(encoding='utf-8')
    frames = parser.parse_dynamic_ndx_text(text)
    times = [float(t) / 1000 for t in re.findall(r'_t([0-9.]+)\s*\]', text)]
    backbone_times, _ = series[(meta['system'], meta['rep'], 'protein_bb_self_RMSD_nm')]
    max_time_difference = check_contact_times(times, backbone_times, len(frames))
    contact_time_checks.append(dict(meta, status='PASS', n_frames=len(frames),
        first_ns=times[0], last_ns=times[-1], max_absolute_difference_ns=max_time_difference,
        absolute_tolerance_ns=1e-4, relative_tolerance=0, source_file=str(p),
        backbone_source_file=str(per / 'rmsd_backbone_full.xvg')))
    counts = Counter()
    for frame in frames:
        assert frame <= amap.keys(), f'接触NDX含非蛋白原子：{meta}'
        counts.update({int(amap[i]['global_residue_index0']) for i in frame})
    old_contacts = read_tsv(per / 'contact_occupancy_residue.tsv')
    assert len({int(r['resind']) for r in old_contacts}) == len(old_contacts)
    assert {int(r['resind']) for r in old_contacts} == set(counts)
    for r in old_contacts:
        idx = int(r['resind'])
        rr = residues[idx]
        assert (r['resname'], int(r['resnr']), r['insertion_code'].strip()) == (
            rr['resname'], int(rr['input_resnr']), str(rr['insertion_code']).strip())
        assert int(r['hit_frames']) == counts[idx] and int(r['total_frames']) == len(frames)
        assert abs(float(r['occupancy']) - counts[idx] / len(frames)) <= 5.1e-7
    # 只有完整动态NDX逐帧与TPR残基全集都核验后，才把未命中残基记为真零。
    for idx, r in residues.items():
        residue_contact.append(dict(meta, **{k: str(r[k]).strip() for k in IDENT},
            chain_source=r['chain_source'], atom_mode='all_atom', geometry='PBC',
            method='legacy_gmx_select', hit_frames=counts[idx], total_frames=len(frames),
            occupancy=counts[idx] / len(frames), cutoff_nm=0.45,
            zero_evidence='完整动态NDX逐帧未命中' if counts[idx] == 0 else '源表与动态NDX一致',
            source_file=str(p)))
    identities.append(dict(meta, kind='contacts', n_residues=len(residues), status='PASS',
                           method='TPR全集+动态NDX重计数；与骨架实际帧逐项匹配；旧稀疏表逐项核对', n_frames=len(frames)))


LEGACY_METRICS = {
    'rmsd_backbone_full': ('protein_bb_self_RMSD_nm', '蛋白骨架；legacy GROMACS质量加权拟合', 'nm', 'legacy_protein_bb'),
    'rmsd_ligand_full': ('ligand_all_atom_RMSD_nm', '全原子配体；蛋白骨架拟合，不与修复重原子RMSD混合', 'nm', 'legacy_protein_bb'),
    'rg_protein_full': ('protein_mass_Rg_nm', '蛋白质量加权Rg', 'nm', 'none'),
    'sasa_protein_full': ('sasa_protein_only_nm2', 'protein-only surface', 'nm^2', 'none'),
    'mindist_protein_ligand': ('protein_lig_min_all_PBC_nm', '蛋白-配体全原子MIC最小距离', 'nm', 'none'),
    'contact_count_045nm': ('legacy_mindist_contact_count', 'GROMACS mindist原生接触计数；不混作唯一残基数', 'count', 'none'),
    'hbond_protein_ligand': ('hbond_protein_ligand', 'GROMACS蛋白-配体氢键', 'count', 'none'),
}


def identity_check(meta, rows, atoms, extras, kind):
    atom_by_index = {int(r['solute_index0']): r for r in atoms}
    assert len(atom_by_index) == len(atoms), (meta, '主原子映射索引重复')
    residue_by_index = {}
    for atom in atoms:
        idx = int(atom['global_residue_index0'])
        if idx in residue_by_index:
            prior = residue_by_index[idx]
            assert idkey(atom) == idkey(prior) and int(atom['molecule_index0']) == int(prior['molecule_index0']), (meta, '主残基映射身份冲突', idx)
        residue_by_index[idx] = atom
    assert kind in {'RMSF', 'contacts'}, kind
    seen = set()
    for r in rows:
        key = idkey(r)
        if kind == 'RMSF':
            idx = int(r['solute_index0'])
            assert idx in atom_by_index, (meta, kind, '原子索引不存在', idx)
            mapped = atom_by_index[idx]
            assert mapped['atom_name'] == 'CA', (meta, kind, '非CA原子', idx)
        else:
            idx = int(r['global_residue_index0'])
            assert idx in residue_by_index, (meta, kind, '残基索引不存在', idx)
            mapped = residue_by_index[idx]
        assert mapped['moltype_name'].startswith('Protein') and key == idkey(mapped), (meta, kind, '索引与五字段身份不符', idx, key)
        for field in ['global_residue_index0', 'molecule_index0']:
            assert int(r[field]) == int(mapped[field]), (meta, kind, field, r[field], mapped[field])
        full = key + tuple(r[k] for k in extras)
        assert full not in seen, (meta, kind, '重复身份', full)
        seen.add(full)
    identities.append(dict(meta, kind=kind, n_rows=len(rows), n_residues=len({idkey(r) for r in rows}),
        status='PASS', method=('solute_index0回查CA原子、五字段及global_residue_index0/molecule_index0' if kind == 'RMSF'
                              else 'global_residue_index0回查五字段及molecule_index0')))


def repaired(meta, d):
    q = read_json(d / 'QC.json')
    m = read_json(d / 'metrics.json')
    reference = read_json(d / 'reference.json')
    assert q['status'] == 'PASS' and m['gromacs']['status'] == 'PASS'
    accepted = d / 'supplemental/SUPPLEMENTAL_ACCEPTED.json'
    review = d / 'boundary_review/MANUAL_BOUNDARY_REVIEW_PASS.json'
    ac, rv = None, None
    if accepted.is_file():
        ac = read_json(accepted)
        assert ac['status'] == 'PASS'
        sp = (d / 'supplemental' / ac['attempt'] / 'summary.json').resolve()
        walk_bindings(ac, accepted)
    elif review.is_file():
        rv = read_json(review)
        assert rv['status'] == 'PASS_WITH_CUTOFF_BOUNDARY_REVIEW'
        sp = winpath(rv['evidence']['source_supplemental_summary']).resolve()
        walk_bindings(rv, review)
        audit = read_json(winpath(rv['evidence']['audit_summary']))
        assert rv['audit_metrics']['residue_identity_difference_in_original_summary'] == {'all_atom': 0, 'heavy_atom': 0}
    else:
        raise ValueError(f'缺少明确补充summary选择证据：{d}')
    assert sp.is_relative_to((d / 'supplemental').resolve())
    s = read_json(sp)
    assert winpath(s['rep']).resolve() == d.resolve()
    assert s['status'] == ('PASS' if ac else 'FAIL')
    walk_bindings(read_json(sp.parent / s['bindings_file']), sp.parent / s['bindings_file'])
    evidence = dict(meta, main_qc=q['status'], automatic_supplemental=s['status'],
        boundary_review=rv['status'] if rv else 'NOT_REQUIRED',
        selected_by=str(accepted if ac else review), supplemental_summary=str(sp),
        main_metrics=str(d / 'metrics.json'), main_qc_file=str(d / 'QC.json'),
        interpretation='主QC通过；自动补充失败保留；边界复核限定接触阈值' if rv else '主QC与自动补充通过')
    qc_rows.append(evidence)
    qc_documents.append(dict(trajectory_id=meta['trajectory_id'], main_qc=q, automatic_supplemental={
        k: s[k] for k in ['status', 'trajectory', 'main_numeric_validation', 'sasa', 'contacts', 'acceptance_scope']},
        boundary_review=rv, accepted_marker=str(accepted) if ac else None, supplemental_summary=str(sp)))
    ts = read_tsv(d / 'timeseries.tsv')
    t = np.array([float(r['time_ns']) for r in ts])
    assert len(ts) == m['frames'] and len(ts) == s['trajectory']['frames_read']
    assert np.allclose(t * 1000, [float(r['time_ps']) for r in ts])
    trajectories.append(dict(meta, directory=str(d), full_first_ns=float(t[0]), full_last_ns=float(t[-1]),
        total_frames=len(ts), window_frames=m['window_50_100ns_frames'],
        replaced_old_TLR4=True, main_qc=q['status'], automatic_supplemental=s['status'],
        boundary_review=evidence['boundary_review']))
    layer = 'main_PASS;automatic_' + s['status'] + ';review_' + evidence['boundary_review']
    for metric, v in m['statistics_50_100ns'].items():
        fit = next((f for f in FITS if metric.startswith(f + '_')), 'none')
        scope = ('修复连续构象；0ps参考；无权重骨架拟合；配体不自拟合' if fit != 'none'
                 else '蛋白质量加权Rg' if metric == 'protein_mass_Rg_nm' else '修复连续构象；' + metric)
        r = metric_row(meta, metric, v, d / 'metrics.json', scope, fit=fit, qc=layer)
        if v is not None:
            add_series(meta, metric, t, [float(x[metric]) for x in ts], r, d / 'timeseries.tsv')
        if 'ligand_heavy_RMSD' in metric:
            fit_rows.append(r)
    for key, value in m['gromacs'].items():
        if not isinstance(value, dict) or 'columns' not in value:
            continue
        p = d / ('sasa_50_100ns.xvg' if key == 'sasa_complex_and_protein' else key + '_50_100ns.xvg')
        a, unit, _ = xvg(p)
        for col, v in enumerate(value['columns'], 1):
            if key == 'sasa_complex_and_protein':
                metric = ['sasa_complex_total_nm2', 'sasa_protein_component_in_complex_nm2'][col - 1]
                scope, u = metric, 'nm^2'
            else:
                metric = key if col == 1 else key + '_non_hbond_pairs'
                scope, u = ('GROMACS氢键；第一数值列' if col == 1 else 'GROMACS氢键输出第二列：距离内非氢键对'), 'count'
            r = metric_row(meta, metric, v, d / 'metrics.json', scope, unit=u, qc='main_PASS')
            add_series(meta, metric, a[:, 0] / (1000 if unit == 'ps' else 1), a[:, col], r, p)
    assert s['sasa']['status'] == 'PASS' and s['rmsf']['status'] == 'PASS'
    r = metric_row(meta, 'sasa_protein_only_nm2', s['sasa']['protein_only_surface'], sp,
                   'protein-only surface', unit='nm^2', qc=layer)
    p = sp.parent / 'sasa_protein_only_50_100ns.xvg'
    a, unit, _ = xvg(p)
    add_series(meta, r['metric'], a[:, 0] / (1000 if unit == 'ps' else 1), a[:, 1], r, p)
    metric_row(meta, 'sasa_ligand_shielding_delta_nm2', s['sasa']['ligand_shielding_delta'], sp,
               'protein-only减复合物蛋白分量；蛋白侧遮蔽面积，非总埋藏面积', unit='nm^2', qc=layer)
    atoms = read_tsv(d / 'solute_atom_mapping.tsv')
    atom_by_index = {int(r['solute_index0']): r for r in atoms}
    for fit, indices in reference['fit_solute_indices0'].items():
        fit_atoms = [atom_by_index[i] for i in indices]
        assert all(a['atom_name'] in {'N', 'CA', 'C'} for a in fit_atoms)
        residue_keys = sorted({idkey(a) for a in fit_atoms})
        if fit == 'pocket_bb':
            pocket_selection[(meta['system'], meta['rep'])] = set(residue_keys)
        reference_rows.append(dict(meta, fit=fit, n_fit_atoms=len(indices), n_fit_residues=len(residue_keys),
            fit_identity_sha256=hashlib.sha256(json.dumps(residue_keys).encode()).hexdigest(),
            reference=reference['reference'], weighting=reference['fit_weight'],
            pocket_definition=reference['pocket_definition'],
            residue_identities=json.dumps(residue_keys, ensure_ascii=False), source_file=str(d / 'reference.json')))
    rr = read_tsv(sp.parent / 'rmsf_ca_50_100ns.tsv')
    identity_check(meta, rr, atoms, ['fit'], 'RMSF')
    for r in rr:
        residue_rmsf.append(dict(meta, **{k: str(r[k]).strip() for k in IDENT},
            chain_source=r['chain_source'], fit=r['fit'], rmsf_nm=float(r['rmsf_nm']),
            source_file=str(sp.parent / 'rmsf_ca_50_100ns.tsv'), window_ns='50-100',
            fit_semantics='第0帧参考；对应骨架N/CA/C无权重拟合', qc_layer=layer))
    for fit, v in s['rmsf']['groups'].items():
        vals = [float(r['rmsf_nm']) for r in rr if r['fit'] == fit]
        if vals:
            assert v['fit_solute_indices0'] == reference['fit_solute_indices0'][fit]
            assert len(vals) == v['ca_atoms']
            assert abs(statistics.mean(vals) - v['rmsf_nm']['mean']) < 1e-10
    if meta['ligand'] == 'apo':
        metric_row(meta, 'hbond_protein_ligand', None, d / 'metrics.json',
                   '无配体不适用', unit='count', qc=layer)
        return
    cr = read_tsv(sp.parent / 'contact_residue_occupancy_50_100ns.tsv')
    identity_check(meta, cr, atoms, ['atom_mode', 'geometry'], 'contacts')
    universe = {idkey(r) for r in rr if r['fit'] == 'protein_bb'}
    for mode in ['all_atom', 'heavy_atom']:
        for geometry in ['direct', 'PBC']:
            selected = [r for r in cr if r['atom_mode'] == mode and r['geometry'] == geometry]
            assert {idkey(r) for r in selected} == universe
            for r in selected:
                assert int(r['total_frames']) == m['window_50_100ns_frames']
                assert 0 <= int(r['hit_frames']) <= int(r['total_frames'])
                assert abs(float(r['occupancy']) - int(r['hit_frames']) / int(r['total_frames'])) < 1e-12
                residue_contact.append(dict(meta, **{k: str(r[k]).strip() for k in IDENT},
                    chain_source=r['chain_source'], atom_mode=mode, geometry=geometry,
                    method='repaired_supplemental', hit_frames=int(r['hit_frames']),
                    total_frames=int(r['total_frames']), occupancy=float(r['occupancy']), cutoff_nm=0.45,
                    zero_evidence='源全残基表显式零' if int(r['hit_frames']) == 0 else '源表显式非零',
                    source_file=str(sp.parent / 'contact_residue_occupancy_50_100ns.tsv'), qc_layer=layer))
    ct = read_tsv(sp.parent / 'contacts_50_100ns.tsv')
    for key, v in s['contacts']['statistics'].items():
        metric = 'contact_' + key
        r = metric_row(meta, metric, v, sp, '0.45nm；' + key, unit='count', qc=layer)
        add_series(meta, metric, [float(x['time_ps']) / 1000 for x in ct], [float(x[key]) for x in ct],
                   r, sp.parent / 'contacts_50_100ns.tsv')
        if key.endswith('_protein_residues'):
            mode, geometry = key.rsplit('_protein_residues', 1)[0].rsplit('_', 1)
            occsum = sum(float(x['occupancy']) for x in cr if x['atom_mode'] == mode and x['geometry'] == geometry)
            assert abs(occsum - v['mean']) < 1e-9


def aggregate(rows, group_keys, value):
    groups = defaultdict(list)
    for r in rows:
        groups[tuple(r[k] for k in group_keys)].append(r)
    result = []
    for key, rr in sorted(groups.items()):
        assert len({r['rep'] for r in rr}) == len(rr), ('重复计入', key)
        vals = {r['rep']: r[value] for r in rr if r[value] is not None}
        system = rr[0]['system']
        expected = sorted(system_reps[system])
        not_selected = []
        if value == 'rmsf_nm' and rr[0]['fit'] == 'pocket_bb':
            not_selected = [rep for rep in expected if idkey(rr[0]) not in pocket_selection[(system, rep)]]
            assert not set(vals).intersection(not_selected), ('未选入口袋却有RMSF', key)
        missing = [rep for rep in expected if rep not in vals and rep not in not_selected]
        n = len(vals)
        result.append(dict(zip(group_keys, key), **stats(list(vals.values())),
            n_expected=len(expected), n_selected=len(expected) - len(not_selected),
            missing_reps=';'.join(missing), not_selected_reps=';'.join(not_selected),
            completeness=('MISSING_NOT_ZERO_FILLED' if missing else
                          'COMPLETE_FOR_SELECTED_REPLICATES' if not_selected else 'COMPLETE'),
            inference=('n=0无有效数据；不作推断' if n == 0 else 'n=1探索性' if n == 1 else
                       f'n={n}独立重复描述性；不以帧数作为n'),
            replicate_availability=';'.join(f'rep{rep}=' + ('OBSERVED' if rep in vals else
                'NOT_SELECTED_IN_POCKET' if rep in not_selected else 'DATA_MISSING') for rep in expected),
            replicate_values=';'.join(f'rep{rep}={vals[rep]:.10g}' if rep in vals else f'rep{rep}=NA' for rep in expected)))
    return result


def savefig(fig, name):
    fig.savefig(OUT / (name + '.png'), dpi=180, bbox_inches='tight')
    fig.savefig(OUT / (name + '.pdf'), bbox_inches='tight')
    plt.close(fig)


def contact_bars(ax, rows, color):
    means = np.asarray([r['mean'] for r in rows], dtype=float)
    sds = np.asarray([r['sd_between_replicates'] or 0 for r in rows], dtype=float)
    bars = ax.barh(np.arange(len(rows)), means, xerr=sds, color=color, alpha=.85)
    low = min(0., float(np.min(means - sds))) if len(rows) else 0.
    high = max(1., float(np.max(means + sds))) if len(rows) else 1.
    padding = .05 * (high - low)
    ax.set_xlim(low - padding, high + padding)
    for boundary in [0, 1]:
        ax.axvline(boundary, color='#475569', linestyle='--', linewidth=.8)
    return bars


def plots(summary, rmsf_summary, contact_summary):
    fig, axes = plt.subplots(3, 4, figsize=(16, 11), layout='constrained')
    choices = [('protein_bb_self_RMSD_nm', '蛋白骨架RMSD (nm)'),
               ('protein_mass_Rg_nm', '蛋白Rg (nm)'), ('sasa_protein_only_nm2', 'Protein-only SASA (nm²)'),
               ('hbond_protein_ligand', '蛋白–配体氢键数')]
    for row, target in enumerate(['IKKb', 'NLRP3', 'TLR4MD2']):
        systems = sorted([s for s in system_reps if s.split('_')[1] == target], key=lambda s: (not s.endswith('_apo'), s))
        for ax, (metric, label) in zip(axes[row], choices):
            for i, system in enumerate(systems):
                vals = [r['mean'] for r in long_rows if r['system'] == system and r['metric'] == metric and r['mean'] is not None]
                if vals:
                    ax.scatter(i + np.linspace(-.13, .13, len(vals)), vals, s=22, color=COLORS[row], zorder=3)
                    ax.errorbar(i, statistics.mean(vals), yerr=statistics.stdev(vals) if len(vals) > 1 else None,
                                fmt='_', color='#202b36', capsize=4, markersize=12)
            ax.set_xticks(range(len(systems)), [s.split('_')[2] + (' (n=1)' if len(system_reps[s]) == 1 else '') for s in systems], rotation=45, ha='right')
            ax.set_title(target + ' / ' + label)
            ax.grid(axis='y', alpha=.2)
    fig.suptitle('50–100 ns：独立重复均值与样本SD；apo配体指标不适用\n各靶点分别解释；TLR4使用修复结果，拟合定义见长表', fontsize=14)
    savefig(fig, '图1_三靶点结构指标')
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), layout='constrained')
    ligands = ['Baicalein', 'OroxylinA', 'Wogonin']
    for ax, fit in zip(axes, FITS):
        for i, ligand in enumerate(ligands):
            vals = [r['mean'] for r in fit_rows if r['ligand'] == ligand and r['fit'] == fit]
            ax.scatter(i + np.linspace(-.1, .1, 3), vals, color=COLORS[i], s=28)
            ax.errorbar(i, statistics.mean(vals), yerr=statistics.stdev(vals), fmt='_', color='black', capsize=4)
        ax.set_xticks(range(3), ligands, rotation=25, ha='right')
        ax.set_title(fit)
        ax.set_ylabel('配体重原子RMSD (nm)')
    fig.suptitle('TLR4–MD2：四种独立骨架拟合，统一0ps参考；点为独立重复')
    savefig(fig, '图2_TLR4配体重原子RMSD_四种拟合')
    drift_metrics = ['protein_bb_self_RMSD_nm', 'protein_mass_Rg_nm', 'sasa_protein_only_nm2', 'hbond_protein_ligand']
    for target in ['IKKb', 'NLRP3', 'TLR4MD2']:
        systems = sorted(s for s in system_reps if s.split('_')[1] == target)
        fig, axes = plt.subplots(len(systems), 4, figsize=(15, 2.3 * len(systems)), squeeze=False, layout='constrained')
        for i, system in enumerate(systems):
            for j, metric in enumerate(drift_metrics):
                ax = axes[i, j]
                for rep in sorted(system_reps[system]):
                    rr = [r for r in blocks if r['system'] == system and r['rep'] == rep and r['metric'] == metric]
                    if rr:
                        ax.plot([r['start_ns'] + 5 for r in rr], [r['mean'] for r in rr], 'o-', ms=3, label='rep' + rep)
                ax.set_title(system.split('_')[2] + ' / ' + choices[j][1], fontsize=9)
                ax.set_xticks([55, 65, 75, 85, 95])
                ax.grid(alpha=.2)
                if j == 0:
                    ax.legend(fontsize=7)
        fig.suptitle(target + '：五个10ns块均值；块相关，仅描述漂移，不作为独立重复')
        savefig(fig, '图3_' + target + '_五块漂移')
    fig, axes = plt.subplots(4, 2, figsize=(14, 12), layout='constrained')
    for i, ligand in enumerate(['apo'] + ligands):
        system = 'v2_TLR4MD2_' + ligand
        for ax, chain in zip(axes[i], ['A', 'C']):
            rr = sorted([r for r in rmsf_summary if r['system'] == system and r['fit'] == 'protein_bb' and r['chain'] == chain], key=lambda r: int(r['input_resnr']))
            x = np.array([int(r['input_resnr']) for r in rr])
            y = np.array([r['mean'] for r in rr])
            sd = np.array([r['sd_between_replicates'] if r['sd_between_replicates'] is not None else np.nan for r in rr])
            ax.plot(x, y, lw=1, color=COLORS[0])
            ax.fill_between(x, y - sd, y + sd, alpha=.25, color=COLORS[0])
            ax.set_title(ligand + ' / chain ' + chain + ' / protein_bb fit')
            ax.set_ylabel('CA RMSF (nm)')
            ax.set_xlabel('输入残基编号')
    fig.suptitle('TLR4–MD2逐残基RMSF：五字段身份匹配后均值±重复间样本SD')
    savefig(fig, '图4_TLR4_RMSF')
    fig, axes = plt.subplots(3, 4, figsize=(18, 12), layout='constrained')
    for i, ligand in enumerate(ligands):
        for j, (mode, geometry) in enumerate([(a, g) for a in ['all_atom', 'heavy_atom'] for g in ['direct', 'PBC']]):
            ax = axes[i, j]
            rr = sorted([r for r in contact_summary if r['system'] == 'v2_TLR4MD2_' + ligand and r['atom_mode'] == mode and r['geometry'] == geometry], key=lambda r: r['mean'] or 0, reverse=True)[:12]
            rr.reverse()
            y = np.arange(len(rr))
            bars = contact_bars(ax, rr, COLORS[i])
            ax.set_yticks(y, [r['chain'] + ':' + r['resname'] + r['input_resnr'] + r['insertion_code'] for r in rr], fontsize=7)
            left, right = ax.get_xlim()
            segments = bars.errorbar.lines[2][0].get_segments()
            assert all(left < segment[:, 0].min() <= segment[:, 0].max() < right for segment in segments)
            contact_plot_checks.append(dict(ligand=ligand, atom_mode=mode, geometry=geometry,
                n_rows=len(rr), xmin=float(left), xmax=float(right), status='PASS',
                formerly_clipped_upper=sum(r['mean'] + (r['sd_between_replicates'] or 0) > 1.15 for r in rr),
                formerly_clipped_lower=sum(r['mean'] - (r['sd_between_replicates'] or 0) < 0 for r in rr)))
            ax.set_title(ligand + ' / ' + mode + ' / ' + geometry, fontsize=9)
            ax.set_xlabel('接触占有率（0–1）')
    fig.suptitle('TLR4–MD2接触占有率：每种定义独立列示前12残基；误差线为三重复样本SD')
    savefig(fig, '图5_TLR4_接触占有率')


def self_tests():
    out = []
    def check(name, condition):
        assert condition, name
        out.append({'test': name, 'status': 'PASS'})
    check('样本SD使用ddof1', abs(stats([1, 2, 3])['sd_between_replicates'] - 1) < 1e-12)
    check('单重复SD不可设零', stats([1])['sd_between_replicates'] is None)
    check('缺失值不可设零', stats([None, 2])['mean'] == 2 and stats([None, 2])['n_replicates'] == 1)
    a = dict(chain='A', moltype_name='Protein_chain_A', input_resnr='1', insertion_code='', resname='ALA')
    for key in IDENT:
        b = dict(a)
        b[key] += 'X'
        check('身份匹配包含_' + key, idkey(a) != idkey(b))
    t = np.array([50, 59.99, 60, 70, 80, 90, 99.5, 100])
    masks = [(t >= start) & ((t < start + 10) if start < 90 else (t <= 100)) for start in range(50, 100, 10)]
    check('五块边界无重叠无漏帧', np.all(np.sum(masks, axis=0) == 1))
    check('NA识别', num('NA') is None)

    def rejects(name, action):
        try:
            action()
        except AssertionError:
            check(name, True)
        else:
            raise AssertionError(name + '：未拦截错误')

    backbone = np.array([50., 60., 70., 80., 90., 99.5])
    check('接触允许真实不足100ns的完整帧序列', check_contact_times(backbone, backbone, 6) == 0)
    rounded = backbone.copy()
    rounded[1] += 1e-5
    check('接触允许XVG格式精度误差', check_contact_times(rounded, backbone, 6) <= 1e-4)
    rejects('接触缺帧必须拒绝', lambda: check_contact_times(backbone[:-1], backbone, 5))
    for name, wrong in [('内部递增但错位', [50, 61, 70, 80, 90, 99.5]),
                        ('窗口内整体偏移', backbone + .001),
                        ('不得使用相对容差', backbone + .0005),
                        ('重复帧', [50, 60, 60, 80, 90, 99.5]),
                        ('伪补100ns', [50, 60, 70, 80, 90, 100])]:
        rejects('接触' + name + '必须拒绝', lambda wrong=wrong: check_contact_times(wrong, backbone, 6))

    meta = metadata('test_TLR4MD2_ligand', '1', 'synthetic')
    atom = dict(a, solute_index0='0', global_residue_index0='0', molecule_index0='0', atom_name='CA')
    other = dict(atom, solute_index0='2', global_residue_index0='1', molecule_index0='1',
                 chain='C', input_resnr='2', insertion_code='B', resname='GLY', moltype_name='Protein_chain_C')
    atoms = [atom, dict(atom, solute_index0='1', atom_name='N'), other]
    rmsf = dict(atom, fit='protein_bb', rmsf_nm=1.)
    contact = dict(atom, atom_mode='all_atom', geometry='PBC')
    before_identities = len(identities)
    try:
        identity_check(meta, [rmsf], atoms, ['fit'], 'RMSF')
        identity_check(meta, [contact], atoms, ['atom_mode', 'geometry'], 'contacts')
        check('RMSF和接触索引回查正常记录通过', len(identities) == before_identities + 2)
        for field, wrong in [('solute_index0', '1'), ('solute_index0', '2'), ('solute_index0', '99'),
                             ('global_residue_index0', '1'), ('molecule_index0', '1')]:
            rejects('RMSF拒绝错配_' + field + '_' + wrong,
                    lambda field=field, wrong=wrong: identity_check(meta, [dict(rmsf, **{field: wrong})], atoms, ['fit'], 'RMSF'))
        for kind, row, extras in [('RMSF', rmsf, ['fit']), ('contacts', contact, ['atom_mode', 'geometry'])]:
            swapped = dict(row, **{k: other[k] for k in IDENT})
            rejects(kind + '拒绝合法五字段整体配错索引',
                    lambda: identity_check(meta, [swapped], atoms, extras, kind))
            for field in IDENT:
                bad = dict(row, **{field: other[field]})
                rejects(kind + '回查包含_' + field,
                        lambda: identity_check(meta, [bad], atoms, extras, kind))
            rejects(kind + '拒绝重复身份', lambda: identity_check(meta, [row, row], atoms, extras, kind))
        for field, wrong in [('global_residue_index0', '1'), ('global_residue_index0', '99'), ('molecule_index0', '1')]:
            rejects('接触拒绝错配_' + field + '_' + wrong,
                    lambda field=field, wrong=wrong: identity_check(meta, [dict(contact, **{field: wrong})], atoms, ['atom_mode', 'geometry'], 'contacts'))
    finally:
        del identities[before_identities:]

    system = meta['system']
    assert system not in system_reps and not any(k[0] == system for k in pocket_selection)
    system_reps[system] = {'1', '2', '3'}
    group_keys = ['system', 'fit'] + IDENT
    try:
        for n in [0, 1, 2, 3]:
            rows = [{**meta, **rmsf, 'rep': str(i + 1), 'rmsf_nm': float(i + 1) if i < n else None} for i in range(3)]
            result = aggregate(rows, group_keys, 'rmsf_nm')[0]
            check(f'汇总按实际n={n}推断', result['n_replicates'] == n and result['inference'].startswith(f'n={n}'))
            check(f'实际n={n}保留样本SD定义', result['sd_between_replicates'] == stats(list(range(1, n + 1)))['sd_between_replicates'])
        pocket_row = {**meta, **rmsf, 'fit': 'pocket_bb'}
        for rep in ['1', '2', '3']:
            pocket_selection[(system, rep)] = {idkey(pocket_row)} if rep != '3' else set()
        result = aggregate([pocket_row], group_keys, 'rmsf_nm')[0]
        check('口袋选择差异与真正缺失分开', result['missing_reps'] == '2' and result['not_selected_reps'] == '3'
              and result['n_selected'] == 2 and result['completeness'] == 'MISSING_NOT_ZERO_FILLED'
              and result['replicate_availability'] == 'rep1=OBSERVED;rep2=DATA_MISSING;rep3=NOT_SELECTED_IN_POCKET')
        result = aggregate([pocket_row, dict(pocket_row, rep='2', rmsf_nm=3.)], group_keys, 'rmsf_nm')[0]
        check('口袋仅选入两个重复不报数据缺失', result['n_replicates'] == 2 and result['missing_reps'] == ''
              and result['not_selected_reps'] == '3' and result['completeness'] == 'COMPLETE_FOR_SELECTED_REPLICATES'
              and result['mean'] == 2. and result['sd_between_replicates'] == stats([1, 3])['sd_between_replicates'])
        rejects('未选入口袋却有数值必须拒绝', lambda: aggregate([dict(pocket_row, rep='3')], group_keys, 'rmsf_nm'))
    finally:
        del system_reps[system]
        for rep in ['1', '2', '3']:
            pocket_selection.pop((system, rep), None)

    fig, ax = plt.subplots()
    try:
        rows = [stats([0, 0, 1]), stats([0, 1, 1]), stats([1])]
        bars = contact_bars(ax, rows, COLORS[0])
        segments = bars.errorbar.lines[2][0].get_segments()
        expected = np.array([[r['mean'] - (r['sd_between_replicates'] or 0),
                              r['mean'] + (r['sd_between_replicates'] or 0)] for r in rows])
        actual = np.array([segment[:, 0] for segment in segments])
        check('图5绘制真实完整mean±SD', np.allclose(actual, expected, rtol=0, atol=1e-12))
        left, right = ax.get_xlim()
        check('图5上下误差线及0和1均不裁切', left < min(0, expected.min()) and right > max(1, expected.max()))
        check('图5确实绘制0和1参考线', {tuple(line.get_xdata()) for line in ax.lines} == {(0, 0), (1, 1)})
        check('图5单重复不虚构SD', rows[-1]['sd_between_replicates'] is None)
    finally:
        plt.close(fig)
    return out


def report(summary, rmsf_summary, contact_summary, validation):
    def fmt(v, digits=4):
        return 'NA' if v is None else f'{v:.{digits}f}'
    def cell(system, metric):
        rr = [r for r in summary if r['system'] == system and r['metric'] == metric]
        if not rr or rr[0]['mean'] is None:
            return '不适用'
        r = rr[0]
        return fmt(r['mean']) + (' ± ' + fmt(r['sd_between_replicates']) if r['n_replicates'] > 1 else '（n=1探索性）')
    lines = ['# 结构综合分析汇总（2026-09-21）', '',
        '分析窗口为50–100 ns。纳入41条不重复轨迹、15个体系：29条既有非TLR4记录与12条修复TLR4–MD2记录；旧TLR4指标全部替换。', '',
        '## 数据、统计与证据口径', '',
        '- 每个重复先计算窗口均值，再对独立重复等权求均值和样本SD（ddof=1）；不按轨迹帧数加权。13个体系各有3条轨迹，IKKb_Apigenin和IKKb_OroxylinA各有1条；每项推断按实际有效重复数n标注，n=1仅探索性且SD为NA，n=0不作推断。',
        '- pocket_bb RMSF逐重复按reference.json固定口袋集合核对：未选入列为NOT_SELECTED_IN_POCKET，已选入而无值列为DATA_MISSING，两者均不补零。n_expected保留体系轨迹数，n_selected为选入数，n_replicates为实际有效值数；选入集合可能跨重复不同。',
        '- apo的配体RMSD、蛋白–配体氢键、配体距离和接触占有率不适用；不会产生虚假配体零值。',
        '- TLR4主QC共12条PASS；自动补充4条PASS、8条FAIL；后8条单列PASS_WITH_CUTOFF_BOUNDARY_REVIEW，自动FAIL未被改写。边界复核仅限定接触阈值一致性问题。',
        '- TLR4主metrics.json的statistics_50_100ns和gromacs均已读取；补充summary严格由ACCEPTED.attempt或边界复核evidence.source_supplemental_summary指定。',
        '- TLR4重原子配体RMSD采用0ps参考，分别按protein_bb、TLR4_bb、MD2_bb、pocket_bb骨架无权重拟合；配体不自拟合。旧非TLR4配体RMSD是全原子/旧GROMACS定义，单列，不与重原子定义合并。',
        '- Rg为蛋白质量加权。protein-only SASA、复合物总SASA和复合物蛋白分量单列；蛋白侧遮蔽差值不等于总埋藏面积，也不代表结合能。apo的遮蔽差值为同构象表面定义对照，不是虚拟配体结果。',
        '- 氢键使用GROMACS输出第一数值列；第二数值列非氢键对单列。距离与接触的全原子/重原子、direct/PBC定义保持独立。',
        '- RMSF与接触按chain/moltype_name/input_resnr/insertion_code/resname匹配，拟合、原子选择及几何定义另作分组键。TLR4 RMSF逐条用solute_index0回查CA原子、五字段、global_residue_index0和molecule_index0；接触用global_residue_index0回查五字段及molecule_index0。旧残基身份由逐重复原始TPR只读导出恢复。',
        '- 缺失不盲目补零：修复接触使用全残基源表显式零；旧接触仅在动态NDX与该重复骨架源时序50–100 ns实际帧逐项匹配（绝对容差1e-4 ns、相对容差0）、TPR残基全集核验后，将逐帧未命中残基记为真零；不补足真实不足100 ns的轨迹。',
        '- 图5按实际绘制的mean±样本SD上下界扩展坐标，并画0与1参考线。SD保留原值；误差线可超出占有率物理范围，不代表负占有率或大于1的观测值。',
        '- 五块使用[50,60)、[60,70)、[70,80)、[80,90)、[90,100]，边界帧只计一次。块相关，只作漂移描述，不能当作5个独立重复或收敛证明。', '',
        '## 真实时间覆盖', '']
    for r in validation['short_trajectories']:
        lines.append(f"- {r['trajectory_id']}：实际末帧 {r['full_last_ns']:.8g} ns，保留观测窗口；没有补足到100 ns。")
    lines += ['', '## 各体系核心指标', '',
              '|体系|n|骨架RMSD (nm)|Rg (nm)|Protein-only SASA (nm²)|蛋白–配体氢键|',
              '|---|---:|---:|---:|---:|---:|']
    for system in sorted(system_reps):
        lines.append('|' + '|'.join([system, str(len(system_reps[system]))] + [cell(system, m) for m in
            ['protein_bb_self_RMSD_nm', 'protein_mass_Rg_nm', 'sasa_protein_only_nm2', 'hbond_protein_ligand']]) + '|')
    lines += ['', '## TLR4配体重原子RMSD', '', '|配体|protein_bb|TLR4_bb|MD2_bb|pocket_bb|', '|---|---:|---:|---:|---:|']
    for lig in ['Baicalein', 'OroxylinA', 'Wogonin']:
        lines.append('|' + '|'.join([lig] + [cell('v2_TLR4MD2_' + lig, f + '_ligand_heavy_RMSD_nm') for f in FITS]) + '|')
    lines += ['', '单位nm，均值±重复间样本SD（n=3）。不同拟合参考组反映不同相对运动；口袋由每条轨迹0ps配体邻域定义，可能跨重复不同，不能当作完全相同的拟合原子集合。', '',
              '## 五块漂移与残基接触', '', '|TLR4配体/对照|重复|protein_bb RMSD末块−首块 (nm)|Rg末块−首块 (nm)|', '|---|---|---:|---:|']
    for system in sorted(s for s in system_reps if 'TLR4MD2' in s):
        for rep in sorted(system_reps[system]):
            vals = {r['metric']: r['last_minus_first'] for r in drifts if r['system'] == system and r['rep'] == rep}
            lines.append(f"|{system.split('_')[2]}|{rep}|{fmt(vals['protein_bb_self_RMSD_nm'])}|{fmt(vals['protein_mass_Rg_nm'])}|")
    lines += ['', '以下仅列重原子PBC定义下平均占有率最高的5个残基；其他定义完整保留于接触表。', '']
    for lig in ['Baicalein', 'OroxylinA', 'Wogonin']:
        top = sorted([r for r in contact_summary if r['system'] == 'v2_TLR4MD2_' + lig and r['atom_mode'] == 'heavy_atom' and r['geometry'] == 'PBC'], key=lambda r: r['mean'], reverse=True)[:5]
        lines.append('- ' + lig + '：' + '；'.join(f"{r['chain']}:{r['resname']}{r['input_resnr']} {r['mean']:.3f}±{r['sd_between_replicates']:.3f}" for r in top) + '。')
    lines += ['', '较低RMSD或较高接触占有率只能支持相应几何描述，不能据此作亲和力、药效排序。末块−首块与块均值范围应结合逐重复曲线阅读，不能单凭窗口均值断言已经收敛。', '',
              '## 核验与文件', '',
              f"- 重算时序核验：{len(numeric)}项；残基身份/接触重计数核验：{len(identities)}项。",
              '- 完整结果见01_统一长表.tsv、02_体系均值与样本SD.tsv、03_轨迹清单.tsv、04_QC分层.tsv及05–19专题表。',
              '- validation.json给出测试、覆盖、绑定和只读核验；input_small_file_hashes.tsv记录实际重算的小文件SHA-256；historical_binding_checks.tsv保存与历史绑定的比较。',
              '- 大轨迹未重算哈希，也未重读轨迹；其历史绑定哈希仅原样保留，当前仅核对文件元数据，明确标记SKIPPED_LARGE_NOT_REHASHED。本次数值重算基于小文件时序。',
              '- 未运行MMGBSA，未修改原始分析文件或其他脚本。只读TPR导出不涉及分子动力学重算。', '']
    for p in sorted(OUT.glob('图*.png')):
        lines += [f'![{p.stem}]({p.name})', '']
    (OUT / '结构综合分析报告.md').write_text('\n'.join(lines), encoding='utf-8')


qc_documents = []
system_reps = defaultdict(set)


def main():
    tests = self_tests()
    register(Path(__file__), '本次汇总脚本')
    register(BASE / 'MD_SHXD_PJI/analyze_all_completed.sh', '旧指标定义只读证据')
    old = read_tsv(OLD / '02_逐重复指标.tsv')
    coverage = read_tsv(OLD / '03_时间覆盖与数值核验.tsv')
    inventory = read_tsv(OLD / '04_完成与文件清单.tsv')
    rmsf_desc = read_tsv(OLD / '05_RMSF逐重复描述统计.tsv')
    run_status = read_tsv(LEGACY / 'summary/input_status.tsv')
    key = lambda r: (r['system'], str(r['rep']))
    inv = {key(r): r for r in inventory}
    assert len(inv) == len(inventory) == 41
    assert set(inv) == {key(r) for r in old} == {key(r) for r in coverage}
    status = {key(r): r for r in run_status}
    desc = {key(r): r for r in rmsf_desc}
    old_groups = defaultdict(list)
    for r in old:
        old_groups[key(r)].append(r)
    repaired_keys = {(d.parent.name, d.name[3:]) for d in REPAIR.glob('v2_*/rep*') if d.is_dir()}
    old_tlr = {k for k in inv if 'TLR4MD2' in k[0]}
    assert repaired_keys == old_tlr and len(repaired_keys) == 12
    for system, rep in inv:
        system_reps[system].add(rep)
    assert len(system_reps) == 15 and Counter(map(len, system_reps.values())) == {3: 13, 1: 2}
    spec = importlib.util.spec_from_file_location('structural_mapping_parser', PARSER)
    parser = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = parser
    spec.loader.exec_module(parser)
    for idx, ((system, rep), item) in enumerate(sorted(inv.items()), 1):
        print(f'[{idx}/41] {system}/rep{rep}', flush=True)
        if 'TLR4MD2' in system:
            repaired(metadata(system, rep, 'TLR4_repair_20260921'), REPAIR / system / ('rep' + rep))
            continue
        meta = metadata(system, rep, item['source'])
        per = LEGACY / 'per_rep' / item['source'] / system / ('rep' + rep)
        run = winpath(status[(system, rep)]['path'])
        assert not item['missing']
        register(per / 'done.flag')
        register(per / 'metrics.tsv')
        qlist = [r for r in coverage if key(r) == (system, rep)]
        backbone_cov = next(r for r in qlist if r['metric'] == 'rmsd_backbone_full')
        trajectories.append(dict(meta, directory=str(per), full_first_ns=float(backbone_cov['first_ns']),
            full_last_ns=float(backbone_cov['last_ns']), total_frames=int(backbone_cov['n_total']),
            window_frames=int(backbone_cov['n_50_100']), replaced_old_TLR4=False,
            main_qc='LEGACY_NUMERIC_VERIFIED', automatic_supplemental='NOT_APPLICABLE', boundary_review='NOT_APPLICABLE'))
        qc_rows.append(dict(meta, main_qc='LEGACY_NUMERIC_VERIFIED', automatic_supplemental='NOT_APPLICABLE',
            boundary_review='NOT_APPLICABLE', selected_by=str(OLD / '04_完成与文件清单.tsv'),
            interpretation='既有覆盖表与时序重算核验；非TLR4修复QC'))
        for r in old_groups[(system, rep)]:
            metric, scope, unit, fit = LEGACY_METRICS[r['metric']]
            row = metric_row(meta, metric, r, OLD / '02_逐重复指标.tsv', scope, unit, fit, qc='legacy_numeric')
            p = per / (r['metric'] + '.xvg')
            a, time_unit, _ = xvg(p)
            times = a[:, 0] / (1000 if time_unit == 'ps' else 1)
            cv = next(x for x in qlist if x['metric'] == r['metric'])
            assert abs(times[-1] - float(cv['last_ns'])) < 1e-7
            assert len(times) == int(cv['n_total'])
            add_series(meta, metric, times, a[:, 1], row, p)
        if meta['ligand'] == 'apo':
            for metric in ['ligand_all_atom_RMSD_nm', 'protein_lig_min_all_PBC_nm', 'hbond_protein_ligand', 'legacy_mindist_contact_count']:
                metric_row(meta, metric, None, OLD / '02_逐重复指标.tsv', '无配体不适用',
                           unit='count' if 'hbond' in metric or 'count' in metric else 'nm', qc='legacy_numeric')
        old_residues(meta, per, run, parser, desc[(system, rep)])
    assert len(trajectories) == len({r['trajectory_id'] for r in trajectories}) == 41
    for r in long_rows:
        if r['mean'] is not None and (r['system'], r['rep'], r['metric']) in series:
            t, _ = series[(r['system'], r['rep'], r['metric'])]
            r['observed_first_ns'], r['observed_last_ns'] = float(t[0]), float(t[-1])
    groups = ['system', 'target', 'ligand', 'metric', 'definition_scope', 'fit', 'unit']
    summary = aggregate(long_rows, groups, 'mean')
    rmsf_summary = aggregate(residue_rmsf, ['system', 'target', 'ligand', 'fit'] + IDENT, 'rmsf_nm')
    contact_summary = aggregate(residue_contact, ['system', 'target', 'ligand', 'atom_mode', 'geometry', 'method'] + IDENT, 'occupancy')
    for r in long_rows:
        if r['ligand'] == 'apo' and (('lig' in r['metric'] and 'shielding' not in r['metric']) or r['metric'] == 'legacy_mindist_contact_count'):
            assert r['mean'] is None, ('apo虚假配体值', r)
    assert not any(r['ligand'] == 'apo' for r in residue_contact)
    assert all(r['source'] == 'TLR4_repair_20260921' for r in long_rows if r['target'] == 'TLR4MD2')
    for r in summary:
        if len(system_reps[r['system']]) == 1:
            assert r['sd_between_replicates'] is None and r['inference'].startswith(f"n={r['n_replicates']}")
    short = [r for r in trajectories if r['full_last_ns'] < 99.9999]
    assert len(short) == 2
    assert Counter(r['automatic_supplemental'] for r in qc_rows if r['target'] == 'TLR4MD2') == {'PASS': 4, 'FAIL': 8}
    write_tsv('01_统一长表.tsv', long_rows)
    write_tsv('02_体系均值与样本SD.tsv', summary)
    write_tsv('03_轨迹清单.tsv', trajectories)
    write_tsv('04_QC分层.tsv', qc_rows)
    write_json('QC分层原始证据.json', qc_documents)
    write_tsv('05_TLR4配体重原子RMSD_逐重复.tsv', fit_rows)
    write_tsv('06_TLR4配体重原子RMSD_体系汇总.tsv', [r for r in summary if 'ligand_heavy_RMSD' in r['metric']])
    write_tsv('07_Rg_SASA_氢键_距离_逐重复.tsv', [r for r in long_rows if 'RMSD' not in r['metric']])
    write_tsv('08_Rg_SASA_氢键_距离_体系汇总.tsv', [r for r in summary if 'RMSD' not in r['metric']])
    write_tsv('09_RMSF_残基逐重复.tsv', residue_rmsf)
    write_tsv('10_RMSF_残基均值与样本SD.tsv', rmsf_summary)
    write_tsv('11_接触占有率_残基逐重复.tsv', residue_contact)
    write_tsv('12_接触占有率_残基均值与样本SD.tsv', contact_summary)
    write_tsv('13_五个10ns块.tsv', blocks)
    write_tsv('14_逐重复漂移描述.tsv', drifts)
    write_tsv('15_时序重算核验.tsv', numeric)
    write_tsv('16_残基身份核验.tsv', identities)
    write_tsv('17_旧TLR4替换清单.tsv', [dict(system=k[0], rep=k[1], old_source=inv[k]['source'],
        replacement=str(REPAIR / k[0] / ('rep' + k[1])), counted_once=True) for k in sorted(old_tlr)])
    write_tsv('18_输入时间覆盖表_原样保留.tsv', coverage)
    write_tsv('19_源指标定义.tsv', [dict(old_metric=k, integrated_metric=v[0], definition=v[1], unit=v[2], fit=v[3]) for k, v in LEGACY_METRICS.items()])
    write_tsv('20_TLR4参考拟合原子集合.tsv', reference_rows)
    print('生成PNG/PDF图表与中文报告', flush=True)
    plots(summary, rmsf_summary, contact_summary)
    unchanged = []
    for p, record in provenance.items():
        if record['status'] == 'MISSING':
            continue
        st = Path(p).stat()
        ok = st.st_size == record['size'] and st.st_mtime_ns == record['mtime_ns']
        if record['sha256'] is not None:
            ok = ok and sha(p) == record['sha256']
        unchanged.append({'path': p, 'status': 'PASS' if ok else 'FAIL',
                          'check': 'size+mtime+SHA256' if record['sha256'] else 'size+mtime_only'})
        if not ok:
            issue('INPUT_CHANGED_DURING_RUN', p, '本次执行前后不一致')
    validation = dict(status='PASS' if not issues else 'FAIL', created_at=datetime.now().isoformat(),
        python=sys.executable, numpy=np.__version__, matplotlib=matplotlib.__version__,
        trajectory_count=len(trajectories), system_count=len(system_reps), replacement_count=12,
        retained_non_TLR4_count=29, system_replicates={k: sorted(v) for k, v in system_reps.items()},
        long_rows=len(long_rows), system_metric_rows=len(summary), rmsf_rows=len(residue_rmsf),
        rmsf_summary_rows=len(rmsf_summary), contact_rows=len(residue_contact), contact_summary_rows=len(contact_summary),
        block_rows=len(blocks), numeric_checks=len(numeric), identity_checks=len(identities),
        short_trajectories=short, synthetic_tests=tests,
        legacy_contact_time_checks=contact_time_checks, contact_plot_checks=contact_plot_checks,
        repaired_identity_source_rows=sum(r.get('n_rows', 0) for r in identities if r['target'] == 'TLR4MD2'),
        pocket_partial_replicate_rows=[r for r in rmsf_summary if r['fit'] == 'pocket_bb' and r['n_replicates'] < r['n_expected']],
        main_qc=dict(Counter(r['main_qc'] for r in qc_rows if r['target'] == 'TLR4MD2')),
        automatic_supplemental=dict(Counter(r['automatic_supplemental'] for r in qc_rows if r['target'] == 'TLR4MD2')),
        boundary_review=dict(Counter(r['boundary_review'] for r in qc_rows if r['target'] == 'TLR4MD2')),
        binding_status_counts=dict(Counter(r['status'] for r in bindings)),
        unique_small_files_hashed=sum(r['status'] == 'HASHED_NOW' for r in provenance.values()),
        large_files_not_rehashed=[r for r in provenance.values() if r['status'] == 'SKIPPED_LARGE_NOT_REHASHED'],
        input_unchanged_checks=unchanged, issues=issues,
        caveats=['大轨迹哈希未重验；历史哈希仅保留，当前元数据核对不等于内容校验',
                 '未读取大轨迹；时序数值从既有小文件重算', '边界复核不改变自动补充FAIL',
                 'n=1仅探索性；五时间块不构成独立重复', 'apo无配体指标为NA；缺失不盲目补零'],
        actions={'MMGBSA_run': False, 'other_scripts_modified': False, 'trajectory_read': False,
                 'legacy_tpr_dump_readonly': True})
    write_tsv('input_small_file_hashes.tsv', list(provenance.values()))
    write_tsv('historical_binding_checks.tsv', bindings)
    write_json('validation.json', validation)
    report(summary, rmsf_summary, contact_summary, validation)
    artifacts = [dict(path=str(p.relative_to(OUT)), size=p.stat().st_size, sha256=sha(p))
                 for p in sorted(OUT.rglob('*')) if p.is_file() and 'matplotlib_cache' not in p.parts and p.name != 'output_small_file_hashes.tsv']
    write_tsv('output_small_file_hashes.tsv', artifacts)
    print(json.dumps({k: validation[k] for k in ['status', 'trajectory_count', 'system_count', 'long_rows',
        'system_metric_rows', 'rmsf_rows', 'contact_rows', 'block_rows', 'numeric_checks', 'identity_checks',
        'binding_status_counts', 'unique_small_files_hashed', 'issues']}, ensure_ascii=False, indent=2), flush=True)
    if issues:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
