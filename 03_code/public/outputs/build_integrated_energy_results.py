# Public source snapshot: replace path placeholders before running.
#!/usr/bin/env python3
"""只读验收已有MM/GBSA结果；不导入runner，不启动任何计算。"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parent.parent
OUT = BASE / 'outputs/综合分析汇总_20260921/energy'
HISTORY = BASE / 'MD_SHXD_PJI/05_analysis/core_3rep_mmpbsa_pbc/summary/core3rep_mmpbsa_per_rep.tsv'
RUNNER = BASE / 'MD_SHXD_PJI/05_analysis/mmgbsa_decomp_20260920_runner.py'
SYSTEMS = ['v2_IKKb_Baicalein', 'v2_IKKb_Wogonin', 'v2_IKKb_Berberine',
           'v2_NLRP3_Baicalein', 'v2_NLRP3_Berberine']
HIST_SYSTEMS = SYSTEMS[:3] + ['v2_IKKb_Emodin']
JOBS = [(s, r) for s in SYSTEMS for r in (1, 2, 3)]
PARAMS = dict(igb=5, intdiel=1, extdiel=78.5, saltcon=0.15, temperature=298.15,
              idecomp=1, dec_verbose=1, print_res='within 6', PBRadii=3,
              surften=0.0072, surfoff=0, molsurf=0, entropy=False, frames=51,
              begin_ps=80000, end_ps=100000, interval_ps=400)
REQUIRED = ['FINAL_RESULTS_MMGBSA.dat', 'FINAL_DECOMP_MMGBSA.dat',
            'ENERGY_PER_FRAME.csv', 'DECOMP_PER_FRAME.csv',
            'ENERGY_3BLOCKS.csv', 'DECOMP_3BLOCKS.csv']
ENERGY_SECTIONS = ['Complex Energy Terms', 'Receptor Energy Terms',
                   'Ligand Energy Terms', 'Delta Energy Terms']
DECOMP_TERMS = ['Internal', 'van der Waals', 'Electrostatic', 'Polar Solvation',
                'Non-Polar Solv.', 'TOTAL']
TOKENS = ['Total Decomposition Contribution (TDC)',
          'Sidechain Decomposition Contribution (SDC)',
          'Backbone Decomposition Contribution (BDC)']
ENERGY_BLOCK_FIELDS = ['block', 'first_frame', 'last_frame', 'first_ns', 'last_ns',
                       'n', 'term', 'mean_kcal_mol', 'sample_sd']
DECOMP_BLOCK_FIELDS = ['block', 'decomp_type', 'residue', 'n', 'term',
                       'mean_kcal_mol', 'sample_sd']
FAILURES = (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError, csv.Error)
UNSPECIFIED_HASH = object()


def required_hash(value):
    require(isinstance(value, str) and re.fullmatch(r'[0-9a-fA-F]{64}', value), '缺少有效SHA256')
    return value


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    result = float(value)
    require(math.isfinite(result), f'非有限数值：{value}')
    return result


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def text(path):
    return Path(path).read_text(encoding='utf-8-sig')


def load_json(path):
    return json.loads(text(path))


def source_path(value, distro='Ubuntu-24.04'):
    value = str(value).replace('\\', '/')
    if os.name == 'nt':
        match = re.match(r'^/mnt/([a-z])/(.*)$', value)
        if match:
            return Path(match[1].upper() + ':/' + match[2])
        if value.startswith('/') and not value.startswith('//'):
            return Path('//wsl.localhost/' + distro + value)
    elif re.match(r'^[A-Za-z]:/', value):
        return Path('/mnt/' + value[0].lower() + value[2:])
    return Path(value)


def mean_sd(values):
    values = list(values)
    return (statistics.mean(values) if values else None,
            statistics.stdev(values) if len(values) > 1 else None)


def csv_rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.reader(stream):
            row = [cell.strip() for cell in row]
            while row and not row[-1]:
                row.pop()
            if row:
                yield row


def read_table(path, fields, delimiter=','):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        require(reader.fieldnames is not None and len(set(reader.fieldnames)) == len(reader.fieldnames),
                f'{path.name}：表头缺失或重复')
        require(set(fields).issubset(reader.fieldnames), f'{path.name}：缺少列 {set(fields) - set(reader.fieldnames)}')
        rows = list(reader)
        for row in rows:
            require(None not in row and all(row[k] is not None for k in fields), f'{path.name}：行宽不一致')
        return rows


def parse_dat(path):
    content = text(path)
    values = re.findall(r'^\s*ΔTOTAL\s+(\S+)', content, re.M)
    frames = re.findall(r'Calculations performed using\s+(\d+)\s+complex frames', content)
    require(len(values) == 1, f'{path.name}：必须恰有一个ΔTOTAL，实际{len(values)}')
    require(frames == ['51'], f'{path.name}：最终dat帧数不是51：{frames}')
    return finite(values[0])


def parse_energy(path):
    tables, section, header = {}, None, None
    for row in csv_rows(path):
        first = row[0]
        if 'Energy Terms' in first:
            require(first in ENERGY_SECTIONS and len(row) == 1, f'未知能量区段：{row}')
            require(first not in tables, f'重复能量区段：{first}')
            section, header = first, None
            tables[section] = []
        elif first == 'Frame #':
            require(section is not None and header is None, '能量表头位置错误或重复')
            require(len(set(row)) == len(row) and 'TOTAL' in row, '能量列重复或缺少TOTAL')
            header = row
        elif header:
            require(first.isdigit() and len(row) == len(header), f'能量数据行格式错误：{row[:3]}')
            tables[section].append(dict(zip(header, map(finite, row))))
    require(set(tables) == set(ENERGY_SECTIONS), '缺少四个能量区段')
    for section, rows in tables.items():
        require([r['Frame #'] for r in rows] == list(range(1, 52)), f'{section}：必须是完整有序的1–51帧')
    delta = tables['Delta Energy Terms']
    columns = list(delta[0])
    require(all(list(t[0]) == columns for t in tables.values()), '四区段能量列不一致')
    for i in range(51):
        for term in columns[1:]:
            difference = (tables[ENERGY_SECTIONS[0]][i][term] - tables[ENERGY_SECTIONS[1]][i][term]
                          - tables[ENERGY_SECTIONS[2]][i][term])
            require(abs(difference - delta[i][term]) <= 0.035, f'第{i + 1}帧{term}：Complex−Receptor−Ligand不一致')
    return delta


def parse_decomp(path):
    series = defaultdict(list)
    section, token, header = None, None, None
    sections, headers = set(), set()
    for row in csv_rows(path):
        first = row[0]
        if first in ['Complex:', 'Receptor:', 'Ligand:', 'DELTAS:']:
            require(first not in sections, f'重复分解区段：{first}')
            sections.add(first)
            section, token, header = first, None, None
        elif 'Decomposition Contribution' in first:
            require(section is not None and first in TOKENS, f'未知分解类型：{first}')
            token, header = first, None
        elif first == 'Frame #':
            require(section and token and (section, token) not in headers, '分解表头位置错误或重复')
            require(row[:2] == ['Frame #', 'Residue'] and set(row[2:]) == set(DECOMP_TERMS)
                    and len(row) == 8, f'分解CSV列不符：{row}')
            header = row
            headers.add((section, token))
        elif header:
            require(first.isdigit() and len(row) == len(header) and row[1], f'分解数据行格式错误：{row[:3]}')
            series[(section, token, row[1])].append((int(first), dict(zip(header[2:], map(finite, row[2:])))))
    delta = {}
    for (section, token, residue), rows in series.items():
        require([r[0] for r in rows] == list(range(1, 52)), f'分解{section}/{token}/{residue}：不是完整1–51帧')
        if section == 'DELTAS:':
            delta[(token, residue)] = [r[1] for r in rows]
    require(delta and {k[0] for k in delta} == set(TOKENS), 'DELTAS分解缺失TDC/SDC/BDC（dec_verbose=1）')
    sets = [{res for tok, res in delta if tok == token} for token in TOKENS]
    require(sets[0] == sets[1] == sets[2], 'TDC/SDC/BDC残基集合不一致')
    return delta


def energy_blocks(delta):
    output = []
    for block in range(1, 4):
        begin = (block - 1) * 17
        for term in list(delta[0])[1:]:
            avg, sd = mean_sd(row[term] for row in delta[begin:begin + 17])
            output.append(dict(block=block, first_frame=begin + 1, last_frame=begin + 17,
                               first_ns=80 + begin * .4, last_ns=80 + (begin + 16) * .4,
                               n=17, term=term, mean_kcal_mol=avg, sample_sd=sd))
    return output


def decomp_blocks(series):
    output = []
    for (token, residue), rows in series.items():
        for block in range(1, 4):
            for term in DECOMP_TERMS:
                avg, sd = mean_sd(r[term] for r in rows[(block - 1) * 17:block * 17])
                output.append(dict(block=block, decomp_type=token, residue=residue, n=17,
                                   term=term, mean_kcal_mol=avg, sample_sd=sd))
    return output


def validate_blocks(path, expected, fields, keys):
    rows = read_table(path, fields)
    expected_map = {tuple(str(row[k]) for k in keys): row for row in expected}
    actual = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        require(key not in actual, f'{path.name}：块键重复 {key}')
        actual[key] = row
    require(actual.keys() == expected_map.keys(), f'{path.name}：块/项/残基集合不一致，实际{len(actual)}，应为{len(expected_map)}')
    for key, row in actual.items():
        for field in set(fields) - set(keys):
            require(abs(finite(row[field]) - expected_map[key][field]) <= 1e-6,
                    f'{path.name}：{key}/{field}与51帧重算不符')
    return expected


class Audit:
    def __init__(self):
        self.checks = []
        self.warnings = []

    @property
    def ok(self):
        return all(item['ok'] for item in self.checks)

    def check(self, name, action):
        record = dict(check=name, ok=False)
        try:
            value = action()
            record['ok'] = True
            return value
        except FAILURES as exc:
            record['error'] = f'{type(exc).__name__}: {exc}'
        finally:
            self.checks.append(record)

    def hashed(self, path, expected=UNSPECIFIED_HASH, role='sha256'):
        path = Path(path)
        def operation():
            require(path.is_file(), f'缺少文件：{path}')
            before = path.stat()
            digest = sha(path)
            after = path.stat()
            require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), f'读取期间文件改变：{path}')
            if expected is not UNSPECIFIED_HASH:
                require(isinstance(expected, str) and re.fullmatch(r'[0-9a-fA-F]{64}', expected), f'无效SHA256：{path}')
                require(digest == expected.lower(), f'哈希不符：{path}；expected={expected} actual={digest}')
            return dict(path=str(path), size=after.st_size, mtime_ns=after.st_mtime_ns, sha256=digest)
        value = self.check(f'{role}: {path}', operation)
        if value:
            self.checks[-1].update(value)
        return value


def validate_params(params):
    for key, expected in PARAMS.items():
        require(key in params and params[key] == expected, f'参数{key}：应为{expected!r}，实际{params.get(key)!r}')
    require(params['entropy'] is False, 'entropy必须为JSON false')


def validate_manifest(folder, system, rep, audit, distro):
    manifest_file = folder / 'manifest.json'
    audit.hashed(manifest_file, role='manifest快照（ACCEPTED未给出该文件hash）')
    manifest = audit.check('读取源manifest', lambda: load_json(manifest_file))
    if not isinstance(manifest, dict):
        audit.check('manifest对象', lambda: require(False, 'manifest不是JSON对象'))
        return
    audit.check('manifest身份', lambda: require(manifest['system'] == system and manifest['rep'] == rep, 'manifest体系/重复不匹配'))
    audit.check('manifest参数', lambda: validate_params(manifest['params']))
    src = audit.check('manifest源目录', lambda: source_path(manifest['source'], distro))
    for filename, key in [('topol.top', 'source_topology_sha256'), ('md.tpr', 'source_tpr_sha256')]:
        digest = audit.check(key, lambda key=key: required_hash(manifest[key]))
        if src is not None and digest is not None:
            audit.hashed(src / filename, digest, '源文件')
    digest = audit.check('分析轨迹hash字段', lambda: required_hash(manifest['analyzed_trajectory_sha256']))
    if digest is not None:
        audit.hashed(folder / 'mmpbsa_80_100ns_dt400_pbc.xtc', digest, '分析轨迹')
    def trajectory_check():
        record = manifest['original_trajectory']
        path = source_path(record['path'], distro)
        st = path.stat()
        require(st.st_size == record['size'], f'源轨迹大小变化：{path}')
        require(abs(st.st_mtime_ns - record['mtime_ns']) <= 100, f'源轨迹mtime变化：{path}')
        return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                    verification='size/mtime；源manifest未提供原轨迹哈希')
    original = audit.check('原轨迹size/mtime（跨平台容差100ns）', trajectory_check)
    if original:
        audit.checks[-1].update(original)
    dependencies = audit.check('manifest依赖列表', lambda: manifest['dependencies'])
    if isinstance(dependencies, list) and dependencies:
        for i, dep in enumerate(dependencies):
            for role in ['source', 'staged']:
                pair = audit.check(f'依赖{i}/{role}字段', lambda dep=dep, role=role: (source_path(dep[role], distro), dep[role + '_sha256']))
                if pair:
                    audit.hashed(*pair, role='拓扑依赖' + role)
    else:
        audit.check('非空拓扑依赖', lambda: require(False, 'dependencies为空或格式错误'))
    def pbc_check():
        pbc = manifest['pbc']
        require(pbc['frames'] == 51, 'PBC帧数不是51')
        require(0 <= finite(pbc['max_direct_vs_minimage_nm']) <= .003, 'PBC direct/minimage偏差超过0.003nm')
        if pbc['max_raw_vs_final_minimage_nm'] is not None:
            require(0 <= finite(pbc['max_raw_vs_final_minimage_nm']) <= .003, 'PBC raw/final偏差超过0.003nm')
        if 'NLRP3' in system:
            require(pbc['max_raw_vs_final_minimage_nm'] is not None, 'NLRP3缺少raw/final校验')
        audit.warnings.extend(pbc.get('warnings', []))
        if 'NLRP3' in system and finite(pbc['direct_max_nm']) >= .5:
            audit.warnings.append('NLRP3存在远离接触区的帧（最小原子距离>=0.5nm）；能量不能单独证明原口袋稳定结合。')
    audit.check('manifest PBC校验', pbc_check)
    after = audit.check('manifest读取后稳定性', lambda: sha(manifest_file))
    initial = next((c.get('sha256') for c in audit.checks if c.get('path') == str(manifest_file)), None)
    audit.check('manifest快照一致', lambda: require(after == initial and initial is not None, 'manifest读取期间改变'))


def validate_rep(root, system, rep, distro):
    folder = root / 'per_rep' / system / f'rep{rep}'
    marker = folder / 'ACCEPTED.json'
    audit = Audit()
    result = dict(system=system, rep=rep, path=str(folder), status='PENDING', accepted_marker=marker.is_file(),
                  checks=audit.checks, warnings=audit.warnings, delta_total_kcal_mol=None, frames=None)
    if not marker.is_file():
        result['reason'] = '无ACCEPTED.json；不读取可能仍在写入的能量CSV，不纳入统计'
        result['calculation_started'] = (folder / 'CALCULATION_STARTED.json').is_file()
        result['calculation_finished'] = (folder / 'CALCULATION_FINISHED.json').is_file()
        if (folder / 'manifest.json').is_file():
            validate_manifest(folder, system, rep, audit, distro)
        result['preflight_ok'] = audit.ok
        return result, None
    initial = audit.hashed(marker, role='ACCEPTED快照')
    accepted = audit.check('解析ACCEPTED', lambda: load_json(marker))
    if not isinstance(accepted, dict):
        result.update(status='INVALID', reason='ACCEPTED无法解析为JSON对象')
        return result, None
    audit.check('ACCEPTED身份及帧块', lambda: require(
        accepted['system'] == system and accepted['rep'] == rep and accepted['frames'] == 51
        and accepted['blocks'] == [17, 17, 17], 'ACCEPTED体系/重复/51帧/3×17块不匹配'))
    outputs = accepted.get('outputs')
    audit.check('ACCEPTED输出清单', lambda: require(isinstance(outputs, dict) and set(REQUIRED).issubset(outputs), 'outputs未包含六个必需文件'))
    if isinstance(outputs, dict):
        for name, digest in outputs.items():
            def safe_path(name=name):
                path = (folder / name).resolve()
                require(path.is_relative_to(folder.resolve()), f'outputs路径越出重复目录：{name}')
                return path
            path = audit.check(f'产物路径：{name}', safe_path)
            if path:
                audit.hashed(path, digest, 'ACCEPTED全部输出hash')
    validate_manifest(folder, system, rep, audit, distro)
    total = audit.check('最终dat', lambda: parse_dat(folder / REQUIRED[0]))
    audit.check('最终分解dat非空', lambda: require((folder / REQUIRED[1]).stat().st_size >= 100, '最终分解dat不足100字节'))
    delta = audit.check('逐帧能量CSV', lambda: parse_energy(folder / 'ENERGY_PER_FRAME.csv'))
    decomp = audit.check('逐帧分解CSV', lambda: parse_decomp(folder / 'DECOMP_PER_FRAME.csv'))
    eb, db = [], []
    if delta:
        avg = statistics.mean(r['TOTAL'] for r in delta)
        audit.check('CSV/dat/ACCEPTED总能量一致', lambda: require(
            total is not None and abs(avg - total) <= .0200001
            and abs(finite(accepted['delta_total_kcal_mol']) - total) <= 1e-8
            and abs(finite(accepted['csv_mean_kcal_mol']) - avg) <= 1e-6,
            'CSV/dat/ACCEPTED总能量不一致；CSV与dat容差0.02 kcal/mol'))
        eb = audit.check('能量3×17块重算', lambda: validate_blocks(folder / 'ENERGY_3BLOCKS.csv',
                         energy_blocks(delta), ENERGY_BLOCK_FIELDS, ['block', 'term'])) or []
    if decomp:
        audit.check('ACCEPTED分解序列数', lambda: require(accepted['decomposition_series'] == len(decomp), '分解序列数不符'))
        db = audit.check('分解3×17块重算', lambda: validate_blocks(folder / 'DECOMP_3BLOCKS.csv',
                         decomp_blocks(decomp), DECOMP_BLOCK_FIELDS, ['block', 'decomp_type', 'residue', 'term'])) or []
    if isinstance(outputs, dict):
        for name, digest in outputs.items():
            path = (folder / name).resolve()
            if path.is_relative_to(folder.resolve()):
                audit.hashed(path, digest, '解析后输出hash复核')
    audit.hashed(marker, initial['sha256'] if initial else None, 'ACCEPTED读取后复核')
    result['status'] = 'VALID' if audit.ok else 'INVALID'
    result['reason'] = '全部验收通过' if audit.ok else '存在验收失败；不纳入数值汇总'
    if not audit.ok:
        return result, None
    avg, frame_sd = mean_sd(row['TOTAL'] for row in delta)
    result.update(delta_total_kcal_mol=total, csv_mean_kcal_mol=avg, frames=51,
                  within_trajectory_sample_sd=frame_sd, decomp_series=len(decomp))
    return result, dict(system=system, rep=rep, delta_total_kcal_mol=total,
                        csv_mean_kcal_mol=avg, within_trajectory_sample_sd=frame_sd,
                        frames=51, energy_blocks=eb, decomp_blocks=db, decomp=decomp)


def validate_package(root, distro):
    audit = Audit()
    audit.check('批次目录', lambda: require(root.is_dir(), f'批次目录不存在：{root}'))
    info = audit.hashed(root / 'package_manifest.json', role='批次manifest快照')
    package = audit.check('读取批次manifest', lambda: load_json(root / 'package_manifest.json'))
    if isinstance(package, dict):
        audit.check('批次参数', lambda: validate_params(package['params']))
        audit.check('批次15任务', lambda: require(Counter(map(tuple, package['jobs'])) == Counter(JOBS), '批次任务必须恰为指定5体系×3重复'))
        digest = audit.check('runner hash字段', lambda: required_hash(package['runner_sha256']))
        if digest is not None:
            audit.hashed(RUNNER, digest, '只读原runner hash')
    else:
        audit.check('批次manifest对象', lambda: require(False, '批次manifest不是对象'))
    if info:
        audit.hashed(root / 'package_manifest.json', info['sha256'], '批次manifest稳定性')
    return dict(ok=audit.ok, checks=audit.checks, warnings=audit.warnings)


def summarize_systems(records, systems):
    output = []
    for system in systems:
        selected = [r for r in records if r['system'] == system]
        avg, sd = mean_sd(r['delta_total_kcal_mol'] for r in selected)
        output.append(dict(system=system, status='COMPLETE' if len(selected) == 3 else 'PARTIAL',
                           n_reps=len(selected), expected_reps=3,
                           available_reps=','.join(str(r['rep']) for r in selected),
                           missing_reps=','.join(str(r) for r in (1, 2, 3) if r not in {v['rep'] for v in selected}),
                           mean_kcal_mol=avg, between_rep_sample_sd=sd))
    return output


def summarize_decomp(records):
    per_rep, union, intersection = [], [], []
    grouped = defaultdict(dict)
    valid_reps = defaultdict(set)
    for record in records:
        system, rep = record['system'], record['rep']
        valid_reps[system].add(rep)
        for (token, residue), rows in record['decomp'].items():
            for term in DECOMP_TERMS:
                avg, sd = mean_sd(row[term] for row in rows)
                per_rep.append(dict(system=system, rep=rep, decomp_type=token, residue=residue,
                                    term=term, frames=51, mean_kcal_mol=avg, within_trajectory_sample_sd=sd))
                grouped[(system, token, residue, term)][rep] = avg
    for (system, token, residue, term), values in sorted(grouped.items()):
        avg, sd = mean_sd(values.values())
        complete = len(valid_reps[system]) == 3
        row = dict(system=system, status='COMPLETE' if complete else 'PARTIAL', decomp_type=token,
                   residue=residue, term=term, n_reps_available=len(valid_reps[system]),
                   n_reps_observed=len(values), expected_reps=3, coverage_of_expected=len(values) / 3,
                   coverage_of_available=len(values) / len(valid_reps[system]),
                   observed_reps=','.join(map(str, sorted(values))),
                   not_selected_reps=','.join(map(str, sorted(valid_reps[system] - values.keys()))),
                   unavailable_reps=','.join(map(str, sorted({1, 2, 3} - valid_reps[system]))),
                   intersection_3rep=complete and len(values) == 3,
                   mean_observed_kcal_mol=avg, between_observed_rep_sample_sd=sd,
                   **{f'rep{r}_kcal_mol': values.get(r) for r in (1, 2, 3)})
        union.append(row)
        if row['intersection_3rep']:
            intersection.append(row.copy())
    return per_rep, union, intersection


def validate_history(path, distro, root):
    audit = Audit()
    audit.hashed(path, role='历史TSV快照')
    table = audit.check('读取历史TSV', lambda: read_table(path,
        ['system', 'rep', 'delta_total', 'n_frames', 'source'], delimiter='\t')) or []
    expected = [(s, r) for s in HIST_SYSTEMS for r in (1, 2, 3)]
    indexed = defaultdict(list)
    for i, row in enumerate(table):
        key = audit.check(f'历史TSV第{i + 1}行身份', lambda row=row: (row['system'], int(row['rep'])))
        if key:
            indexed[key].append(row)
    audit.check('历史12行唯一身份', lambda: require(Counter({k: len(v) for k, v in indexed.items()}) == Counter(expected), '历史表不是指定4体系×3重复，或存在重复行'))
    guard = {}
    guard_file = root / 'protected_old_results.json'
    if guard_file.is_file():
        audit.hashed(guard_file, role='runner历史保护清单快照')
        raw = audit.check('历史保护清单解析', lambda: load_json(guard_file))
        if isinstance(raw, dict):
            guard = {str(source_path(k, distro).resolve()): v for k, v in raw.items()}
    records, details = [], []
    for system, rep in expected:
        local = Audit()
        entry = dict(system=system, rep=rep, status='INVALID', checks=local.checks)
        rows = indexed.get((system, rep), [])
        if len(rows) != 1:
            local.check('历史行唯一性', lambda: require(False, '缺少或重复历史行'))
        else:
            row = rows[0]
            source = source_path(row['source'], distro)
            expected_path = path.parent.parent / 'per_rep' / system / f'rep{rep}' / 'FINAL_RESULTS_MMPBSA_PBC_51FRAMES.dat'
            local.check('历史来源路径', lambda: require(source.resolve() == expected_path.resolve(), f'历史source与体系/重复路径不匹配：{source}'))
            info = local.hashed(source, guard.get(str(source.resolve()), UNSPECIFIED_HASH), '历史最终dat')
            total = local.check('历史最终dat总能量', lambda: parse_dat(source))
            local.check('历史TSV与dat一致', lambda: require(total is not None and int(row['n_frames']) == 51
                        and abs(finite(row['delta_total']) - total) <= 1e-8, '历史TSV帧数或总能量与dat不一致'))
            inp = source.parent / 'mmpbsa_51frames.in'
            local.hashed(inp, role='历史参数输入')
            def historical_method():
                content = re.sub(r'!.*', '', text(inp))
                require(re.search(r'\bigb\s*=\s*5\b', content, re.I), '历史输入igb不是5')
                require(re.search(r'\bsaltcon\s*=\s*0\.150?\b', content, re.I), '历史saltcon不是0.15')
                require(not re.search(r'\b(?:qh_entropy|interaction_entropy|c2_entropy|nmode)\s*=\s*[1-9]', content, re.I), '历史输入启用熵计算')
                require(not re.search(r'&(?:decomp|nmode)\b', content, re.I), '历史输入含分解或nmode，需单独审查')
            local.check('历史方法（GB、无熵、无分解）', historical_method)
            if info:
                local.hashed(source, info['sha256'], '历史dat稳定性')
            if local.ok:
                entry.update(status='VALID', delta_total_kcal_mol=total, source=str(source))
                records.append(dict(system=system, rep=rep, delta_total_kcal_mol=total, frames=51,
                                    source=str(source), dat_sha256=info['sha256'],
                                    within_trajectory_sample_sd=finite(row['within_traj_sd']) if row.get('within_traj_sd') else None))
        details.append(entry)
    return records, dict(status='COMPLETE' if len(records) == 12 and audit.ok else 'PARTIAL',
                         valid_reps=len(records), expected_reps=12, checks=audit.checks, jobs=details)


def output_path(folder, name):
    path = (folder / name).resolve()
    require(path.is_relative_to(OUT.resolve()), f'拒绝在指定新输出目录外写入：{path}')
    return path


def save_text(folder, name, content):
    path = output_path(folder, name)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(content, encoding='utf-8')
    temp.replace(path)


def save_json(folder, name, value):
    save_text(folder, name, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def save_table(folder, name, rows, fields):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    for row in rows:
        writer.writerow({key: 'NA' if row.get(key) is None else row[key] for key in fields})
    save_text(folder, name, '\ufeff' + stream.getvalue())


def plot_results(folder, new, old, summaries, old_summaries, union, intersection, per_rep):
    os.environ['MPLCONFIGDIR'] = str(output_path(folder, 'matplotlib_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ['Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC'] if f in fonts), 'DejaVu Sans')
    plt.rcParams.update({'font.family': font, 'axes.unicode_minus': False, 'pdf.fonttype': 42,
                         'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    colors = ['#237d9b', '#d78532', '#8260a5']
    files = []
    def finish(fig, stem):
        for ext in ['png', 'pdf']:
            filename = stem + '.' + ext
            fig.savefig(output_path(folder, filename), dpi=180, bbox_inches='tight')
            files.append(filename)
        plt.close(fig)
    for records, stats, label, stem in [(new, summaries, '新批次', 'new_energy_points'),
                                       (old, old_summaries, '独立历史批次', 'historical_energy_points')]:
        fig, ax = plt.subplots(figsize=(10, 5.2))
        for i, summary in enumerate(stats):
            for rep in (1, 2, 3):
                matches = [r for r in records if r['system'] == summary['system'] and r['rep'] == rep]
                if matches:
                    ax.scatter(i + (rep - 2) * .13, matches[0]['delta_total_kcal_mol'], color=colors[rep - 1], s=45, zorder=3)
            if summary['mean_kcal_mol'] is not None:
                ax.errorbar(i + .30, summary['mean_kcal_mol'], yerr=summary['between_rep_sample_sd'],
                            fmt='D', color='#222222', capsize=4, markersize=4)
        for rep in (1, 2, 3):
            ax.scatter([], [], color=colors[rep - 1], label=f'rep{rep}')
        ax.plot([], [], 'D', color='#222222', label='mean ± 独立重复sample SD')
        ax.set_xticks(range(len(stats)), [r['system'].replace('v2_', '').replace('_', '\n') +
                       f"\nn={r['n_reps']}/3 {r['status']}" for r in stats])
        ax.set_ylabel('ΔTOTAL (kcal/mol)，无熵')
        ax.set_title(f"{label}能量点图 | 已核验 {len(records)}/{len(stats) * 3}")
        ax.axhline(0, color='#999999', lw=.7)
        ax.grid(axis='y', alpha=.2)
        if not records:
            ax.text(.5, .55, 'PARTIAL：0条通过验收，暂无可绘制能量', transform=ax.transAxes, ha='center', fontsize=13)
        ax.legend(loc='upper left', fontsize=8, ncols=2)
        fig.tight_layout()
        finish(fig, stem)
    fig, axes = plt.subplots(3, 2, figsize=(11, 11), squeeze=False)
    for ax, system in zip(axes.flat, SYSTEMS):
        selected = [r for r in new if r['system'] == system]
        for record in selected:
            blocks = [r for r in record['energy_blocks'] if r['term'] == 'TOTAL']
            ax.plot([r['block'] for r in blocks], [r['mean_kcal_mol'] for r in blocks], 'o-',
                    color=colors[record['rep'] - 1], label=f"rep{record['rep']}")
        ax.set_title(system.replace('v2_', '') + f' | n={len(selected)}/3')
        ax.set_xticks([1, 2, 3], ['80–86.4', '86.8–93.2', '93.6–100'])
        ax.set_xlabel('时间块 (ns)，每块17帧')
        ax.set_ylabel('ΔTOTAL (kcal/mol)')
        if selected:
            ax.legend()
        else:
            ax.text(.5, .5, 'PARTIAL：暂无已验收时间块', ha='center', transform=ax.transAxes)
    axes.flat[-1].axis('off')
    fig.suptitle('三个时间块仅描述稳定性，不作为独立重复', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .97))
    finish(fig, 'new_energy_3blocks')
    for rows, stem, label in [(union, 'new_decomp_top_union', '并集：仅已观测重复均值，标注覆盖n/3'),
                              (intersection, 'new_decomp_top_intersection', '三重复交集：仅覆盖3/3的残基')]:
        fig, axes = plt.subplots(3, 2, figsize=(13, 15), squeeze=False)
        for ax, system in zip(axes.flat, SYSTEMS):
            selected = [r for r in rows if r['system'] == system and r['decomp_type'] == TOKENS[0]
                        and r['term'] == 'TOTAL' and r['residue'].startswith('R:')]
            selected = sorted(selected, key=lambda r: abs(r['mean_observed_kcal_mol']), reverse=True)[:12]
            selected.sort(key=lambda r: r['mean_observed_kcal_mol'])
            if selected:
                vals = [r['mean_observed_kcal_mol'] for r in selected]
                ax.barh(range(len(vals)), vals, color=['#237d9b' if v < 0 else '#d78532' for v in vals])
                for i, row in enumerate(selected):
                    sd = row['between_observed_rep_sample_sd']
                    if sd is not None:
                        ax.errorbar(row['mean_observed_kcal_mol'], i, xerr=sd, color='#333333', capsize=2)
                ax.set_yticks(range(len(vals)), [f"{r['residue']} ({r['n_reps_observed']}/3)" for r in selected])
                ax.invert_yaxis()
            else:
                ax.text(.5, .5, '暂无满足覆盖条件的受体残基', ha='center', transform=ax.transAxes)
            n = sum(r['system'] == system for r in new)
            ax.set_title(system.replace('v2_', '') + f" | n={n}/3 {'COMPLETE' if n == 3 else 'PARTIAL'}")
            ax.axvline(0, color='#888888', lw=.7)
            ax.set_xlabel('TDC TOTAL (kcal/mol)')
        axes.flat[-1].axis('off')
        fig.suptitle('受体残基绝对贡献Top12 | ' + label + '\n误差线为独立重复sample SD；缺失不填0；负值有利', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .955))
        finish(fig, stem)
    fig, axes = plt.subplots(5, 3, figsize=(17, 21), squeeze=False)
    for i, system in enumerate(SYSTEMS):
        for rep in (1, 2, 3):
            ax = axes[i, rep - 1]
            selected = [r for r in per_rep if r['system'] == system and r['rep'] == rep and
                        r['decomp_type'] == TOKENS[0] and r['term'] == 'TOTAL' and r['residue'].startswith('R:')]
            selected = sorted(selected, key=lambda r: abs(r['mean_kcal_mol']), reverse=True)[:10]
            selected.sort(key=lambda r: r['mean_kcal_mol'])
            if selected:
                vals = [r['mean_kcal_mol'] for r in selected]
                ax.barh(range(len(vals)), vals, color=['#237d9b' if v < 0 else '#d78532' for v in vals])
                ax.set_yticks(range(len(vals)), [r['residue'] for r in selected])
                ax.invert_yaxis()
            else:
                ax.text(.5, .5, '无已验收分解数据', ha='center', transform=ax.transAxes)
            ax.set_title(system.replace('v2_', '') + f' rep{rep}', fontsize=9)
            ax.set_xlabel('TDC TOTAL (kcal/mol)')
            ax.axvline(0, color='#888888', lw=.6)
    fig.suptitle('逐重复受体残基绝对贡献Top10 | 51帧均值；不使用帧内SD冒充重复间SD', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .975))
    finish(fig, 'new_decomp_top_per_rep')
    return files


def report_text(validation, summaries, historical_summaries, union, intersection):
    def number(value):
        return 'NA' if value is None else f'{value:.3f}'
    lines = ['# 能量综合分析报告', '',
             f"生成时间：{validation['generated_at']}", '',
             f"新批次：`{validation['root']}`", '',
             f"**新批次状态：{validation['status']}；通过验收 {validation['valid_reps']}/15；"
             f"ACCEPTED标记 {validation['accepted_markers']}/15。**", '',
             '仅VALID结果进入统计；无ACCEPTED的pilot或计算中目录不读取能量CSV。0条时均值和SD均为NA，图中明确显示无数据。', '',
             '## 新批次总能量', '',
             '| 体系 | 状态 | 独立重复数 | 均值 kcal/mol | 重复间sample SD | 缺失重复 |',
             '|---|---|---:|---:|---:|---|']
    for row in summaries:
        lines.append(f"| {row['system']} | {row['status']} | {row['n_reps']}/3 | {number(row['mean_kcal_mol'])} | {number(row['between_rep_sample_sd'])} | {row['missing_reps'] or '无'} |")
    lines += ['', 'mean按已验收独立重复等权计算；sample SD使用ddof=1，n<2时为NA。n<3为PARTIAL，现有重复统计仅作阶段性描述。',
              '51帧来自80–100 ns、dt=0.4 ns；三个17帧块为80–86.4、86.8–93.2、93.6–100 ns。块内SD仅描述轨迹波动；三个块相关，不视为三个独立重复，不计算推断性p值。', '',
              '## 验收与解析', '',
              '逐条核对ACCEPTED中的全部outputs SHA256、体系/重复、51帧、[17,17,17]、最终dat ΔTOTAL，以及四区段ENERGY_PER_FRAME和DELTAS分解。',
              '从逐帧数据重算ENERGY_3BLOCKS与DECOMP_3BLOCKS全部均值及sample SD；检查重复键、缺帧、非有限值、行宽和表头。CSV与dat均值容差0.02 kcal/mol；已由CSV生成的块表容差1e-6。',
              '源manifest校验参数、源topol.top/TPR、分析轨迹、每个源/暂存拓扑依赖hash及PBC距离一致性。原轨迹manifest只提供size/mtime，按这两项核验（mtime跨平台容差100ns），不能宣称原轨迹独立hash已校验。',
              '原runner的ACCEPTED不含manifest哈希；本脚本记录manifest当前SHA256并核对内容/源文件及读取前后稳定性，不声称恢复了验收当时manifest签名。未验收目录中已有manifest只作阶段性预检。', '']
    errors = []
    for check in validation['package']['checks']:
        if not check['ok']:
            errors.append('批次：' + check['check'] + '；' + check['error'])
    for job in validation['jobs']:
        if job['status'] != 'VALID':
            lines.append(f"- {job['system']} rep{job['rep']}：{job['status']}；{job['reason']}")
        for check in job['checks']:
            if not check['ok']:
                errors.append(f"{job['system']} rep{job['rep']}：{check['check']}；{check['error']}")
        for warning in job['warnings']:
            lines.append(f"- {job['system']} rep{job['rep']} 提示：{warning}")
    lines += ['', '实际验收/解析错误：' + (str(len(errors)) + '项，详见下列及validation.json。' if errors else '未发现。未完成任务不等于解析错误；尚无真实已完成CSV时，仅合成CSV格式通过测试不能替代真实批次验收。')]
    lines += ['- ' + error for error in errors]
    lines += ['', '## 分解覆盖与Top贡献', '',
              'within 6是首帧界面6 Å选区，各重复选中的残基可不同。以软件完整残基标识（包含R/L、链、残基名及编号）对齐，保留TDC/SDC/BDC六项数据；不自动丢弃链或拼接不同编号。',
              '并集表逐项记录rep1/2/3值、观测重复数、相对预期3重复的覆盖率、相对可用重复的覆盖率、未入选重复和未验收重复。两类缺失均为NA，不填0。并集均值和SD仅针对观测重复，不能解释为完整三重复均值。',
              '三重复交集只在体系三个重复均通过验收、且该残基在三个重复均出现时输出；部分批次不生成伪三重复交集。',
              f'当前并集共有{len(union)}个残基/类型/能量项，三重复交集{len(intersection)}项。',
              'Top图按TDC TOTAL均值绝对值排序，同时保留有利负贡献和不利正贡献，仅绘制R:受体残基；配体贡献仍保留在表中。TDC、SDC、BDC不能相加为三份贡献。局部选区残基总和不要求等于全体系ΔTOTAL。', '',
              '## 独立历史结果及方法对照', '',
              f"历史状态：**{validation['history']['status']}，{validation['history']['valid_reps']}/12**。历史TSV逐条与12个最终dat的51帧及ΔTOTAL核对；存在runner保护清单时同时复核保护哈希。", '',
              '| 历史体系 | 状态 | n | mean kcal/mol | 独立重复sample SD |',
              '|---|---|---:|---:|---:|']
    for row in historical_summaries:
        lines.append(f"| {row['system']} | {row['status']} | {row['n_reps']} | {number(row['mean_kcal_mol'])} | {number(row['between_rep_sample_sd'])} |")
    lines += ['', '| 方法 | 历史批次 | 新分解批次 |', '|---|---|---|',
              '| 能量模型 | 最终dat为GB，输入igb=5、saltcon=0.15 | igb=5、intdiel=1、extdiel=78.5、saltcon=0.15、PBRadii=3 |',
              '| 参数证据 | 历史输入未显式给出intdiel/extdiel/PBRadii，不假定已验证与新批次相同 | 源manifest显式核验；temperature=298.15 K |',
              '| 熵 | 无熵，历史输入未启用熵或nmode | 无熵，entropy=false |',
              '| 分解 | 历史输入无decomp，不能补造逐残基结果 | idecomp=1、dec_verbose=1、within 6首帧选区 |',
              '| 抽帧 | 既有PBC-safe 51帧（文件名80–100ns dt400） | 80–100ns，dt=0.4ns，51帧及3×17块 |',
              '| 体系 | IKKβ四配体×3，含Emodin，无NLRP3 | IKKβ三配体及NLRP3两配体×3，不含Emodin |', '',
              '历史命名MMPBSA不代表其使用PB；最终dat标明Generalized Born。新旧分别输出表图；新批次IKKβ复用原轨迹，同体系同rep属于重分析，绝不合成6个独立重复。跨方法数值差异不能归因于亲和力改变。',
              'IKKβ沿用既有cluster轨迹并复核direct/min-image；NLRP3仅whole及protein-centered mol重映像，保留真实远离。', '',
              '## NLRP3远离口袋警告及解释边界', '',
              '**NLRP3配体远离原结合口袋时，即使ΔTOTAL仍为负，也不能据此证明原口袋稳定结合或更强亲和力。** 应结合口袋距离、接触占有率及结构证据解释；全蛋白最小距离不能替代特定口袋距离。',
              'manifest若出现最小原子距离>=0.5nm则逐条显示距离警告；未有验收结果时此处仅为解释限制，不能声称已验证本批具体远离程度。',
              '所有结果均为无熵MM/GBSA端点估计，不等同实验结合自由能，不据此计算KD或跨靶点直接排名。', '',
              '## 输出', '',
              'new_energy_per_rep.tsv / new_energy_system_summary.tsv：逐重复状态与独立重复统计。',
              'new_energy_3blocks.tsv / new_decomp_3blocks.tsv：51帧重算的三时间块。',
              'new_decomp_per_rep.tsv / new_decomp_union_coverage.tsv / new_decomp_intersection_3rep.tsv：逐重复、并集覆盖与严格三重复交集。',
              'historical_energy_per_rep.tsv / historical_energy_system_summary.tsv：独立历史结果。',
              'new_energy_points、historical_energy_points、new_energy_3blocks、new_decomp_top_per_rep、new_decomp_top_union、new_decomp_top_intersection均有PNG/PDF。',
              'validation.json记录全部校验及来源快照；output_manifest.json记录本次表图和报告的SHA256。', '']
    historical_errors = [c for c in validation['history']['checks'] if not c['ok']]
    historical_errors += [c for j in validation['history']['jobs'] for c in j['checks'] if not c['ok']]
    if historical_errors:
        lines += ['历史校验错误：'] + ['- ' + c['check'] + '；' + c['error'] for c in historical_errors]
    return '\n'.join(lines)


def build(root, folder=OUT, history=HISTORY, distro='Ubuntu-24.04', plots=True):
    root, folder, history = Path(root).resolve(), Path(folder).resolve(), Path(history).resolve()
    require(folder.is_relative_to(OUT.resolve()), '输出仅允许在指定energy目录内')
    require(not folder.is_relative_to(root) and not root.is_relative_to(folder), '输出与输入批次目录不得包含彼此')
    folder.mkdir(parents=True, exist_ok=True)
    package = validate_package(root, distro)
    jobs, records = [], []
    for system, rep in JOBS:
        result, record = validate_rep(root, system, rep, distro)
        if record and not package['ok']:
            result.update(status='INVALID', reason='单条验收通过但批次manifest校验失败，排除汇总', delta_total_kcal_mol=None)
            record = None
        jobs.append(result)
        if record:
            records.append(record)
    historical, hist_validation = validate_history(history, distro, root)
    summaries = summarize_systems(records, SYSTEMS)
    historical_summaries = summarize_systems(historical, HIST_SYSTEMS)
    per_rep, union, intersection = summarize_decomp(records)
    validation = dict(schema_version=1, generated_at=datetime.now().astimezone().isoformat(), root=str(root),
                      output=str(folder), status='COMPLETE' if len(records) == 15 else 'PARTIAL',
                      expected_reps=15, valid_reps=len(records), accepted_markers=sum(j['accepted_marker'] for j in jobs),
                      invalid_reps=sum(j['status'] == 'INVALID' for j in jobs), no_accepted_results=not records,
                      params=PARAMS, package=package, jobs=jobs, history=hist_validation,
                      statistics='独立重复等权mean，sample SD ddof=1；帧/时间块不是独立重复；缺失不填0',
                      safety=dict(source_read_only=True, mmgbsa_executed=False),
                      runtime=dict(python=sys.executable, version=sys.version), script_sha256=sha(__file__))
    summary_fields = ['system', 'status', 'n_reps', 'expected_reps', 'available_reps', 'missing_reps', 'mean_kcal_mol', 'between_rep_sample_sd']
    union_fields = ['system', 'status', 'decomp_type', 'residue', 'term', 'n_reps_available', 'n_reps_observed',
                    'expected_reps', 'coverage_of_expected', 'coverage_of_available', 'observed_reps',
                    'not_selected_reps', 'unavailable_reps', 'intersection_3rep', 'mean_observed_kcal_mol',
                    'between_observed_rep_sample_sd', 'rep1_kcal_mol', 'rep2_kcal_mol', 'rep3_kcal_mol']
    tables = [
        ('new_energy_per_rep.tsv', jobs, ['system', 'rep', 'status', 'accepted_marker', 'frames', 'delta_total_kcal_mol',
           'csv_mean_kcal_mol', 'within_trajectory_sample_sd', 'reason', 'path']),
        ('new_energy_system_summary.tsv', summaries, summary_fields),
        ('new_energy_3blocks.tsv', [dict(system=r['system'], rep=r['rep'], **b) for r in records for b in r['energy_blocks']],
           ['system', 'rep'] + ENERGY_BLOCK_FIELDS),
        ('new_decomp_3blocks.tsv', [dict(system=r['system'], rep=r['rep'], **b) for r in records for b in r['decomp_blocks']],
           ['system', 'rep'] + DECOMP_BLOCK_FIELDS),
        ('new_decomp_per_rep.tsv', per_rep, ['system', 'rep', 'decomp_type', 'residue', 'term', 'frames', 'mean_kcal_mol', 'within_trajectory_sample_sd']),
        ('new_decomp_union_coverage.tsv', union, union_fields),
        ('new_decomp_intersection_3rep.tsv', intersection, union_fields),
        ('historical_energy_per_rep.tsv', historical, ['system', 'rep', 'frames', 'delta_total_kcal_mol', 'within_trajectory_sample_sd', 'source', 'dat_sha256']),
        ('historical_energy_system_summary.tsv', historical_summaries, summary_fields)]
    for name, rows, fields in tables:
        save_table(folder, name, rows, fields)
    files = [item[0] for item in tables]
    if plots:
        files += plot_results(folder, records, historical, summaries, historical_summaries, union, intersection, per_rep)
    save_json(folder, 'validation.json', validation)
    save_text(folder, '中文能量报告.md', report_text(validation, summaries, historical_summaries, union, intersection))
    files += ['validation.json', '中文能量报告.md']
    save_json(folder, 'output_manifest.json', dict(generated_at=validation['generated_at'], root=str(root),
              outputs={name: sha(folder / name) for name in files}))
    return validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, help='新MMGBSA批次根目录（只读，必须显式指定）')
    parser.add_argument('--wsl-distro', default='Ubuntu-24.04', help='用于只读解析Linux绝对路径的WSL发行版')
    args = parser.parse_args()
    root = source_path(args.root, args.wsl_distro)
    result = build(root, distro=args.wsl_distro)
    print(json.dumps(dict(status=result['status'], valid_reps=result['valid_reps'], expected_reps=15,
                          accepted_markers=result['accepted_markers'], invalid_reps=result['invalid_reps'],
                          historical_status=result['history']['status'], historical_valid_reps=result['history']['valid_reps'],
                          output=str(OUT)), ensure_ascii=False))
    has_errors = (not result['package']['ok'] or result['invalid_reps'] or result['history']['status'] != 'COMPLETE'
                  or any(not c['ok'] for job in result['jobs'] for c in job['checks']))
    return 2 if has_errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
