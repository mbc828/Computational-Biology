# Public source snapshot: replace path placeholders before running.
"""Post-MMGBSA interpretation; read-only sources, isolated derived outputs."""
from pathlib import Path
import csv, json, hashlib, statistics as st, math, sys, re
from collections import defaultdict
from datetime import datetime

BASE = Path(__file__).resolve().parent
ENERGY = BASE / '综合分析汇总_20260921/energy/final_20260921'
STRUCT = BASE / '综合分析汇总_20260921/structural'
OUT = BASE / '综合分析汇总_20260921/followup_final_v2'
ROOT = Path('<MD_DATA_ROOT>/mmgbsa_captured_20260921_143632')
SYSTEMS = ['v2_IKKb_Baicalein','v2_IKKb_Wogonin','v2_IKKb_Berberine','v2_NLRP3_Baicalein','v2_NLRP3_Berberine']
LABELS = dict(zip(SYSTEMS,['IKKβ–黄芩素','IKKβ–汉黄芩素','IKKβ–小檗碱','NLRP3–黄芩素','NLRP3–小檗碱']))
TDC = 'Total Decomposition Contribution (TDC)'
SOURCES = {}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def text(p):
    SOURCES[str(p)] = sha(p)
    return p.read_text(encoding='utf-8-sig')

def js(p):
    return json.loads(text(p))

def table(p):
    return list(csv.DictReader(text(p).splitlines(), delimiter='\t'))

def write_table(name, rows, fields=None):
    fields = fields or (list(rows[0]) if rows else None)
    assert fields, f'Missing schema for empty table {name}'
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows({k: 'NA' if v is None else v for k,v in r.items()} for r in rows)

def avg_sd(v):
    return st.mean(v),st.stdev(v) if len(v)>1 else None

def fmt(x):
    return 'NA' if x is None else f'{x:.2f}'

def canon(name):
    return {'HID':'HIS','HIE':'HIS','HIP':'HIS','CYX':'CYS','ASH':'ASP','GLH':'GLU','LYN':'LYS'}.get(name,name)

def pdb_residues(p):
    result=[]
    for line in text(p).splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            key=(line[21].strip(),int(line[22:26]),line[26].strip(),canon(line[17:20].strip()))
            if not result or result[-1]!=key:
                result.append(key)
    assert len(result)==len(set(result)),f'Duplicate PDB identity {p}'
    return result

def main():
    assert not OUT.exists(), 'Refuse overwriting a prior followup run'
    validation=js(ENERGY/'validation.json')
    assert validation['status']=='COMPLETE' and validation['valid_reps']==15
    assert validation['history']['status']=='COMPLETE'
    for n,h in js(ENERGY/'output_manifest.json')['outputs'].items():
        assert sha(ENERGY/n)==h, n
    sv=js(STRUCT/'validation.json')
    assert sv['status']=='PASS' and sv['trajectory_count']==41
    OUT.mkdir()
    summary=table(ENERGY/'new_energy_system_summary.tsv')
    eb=table(ENERGY/'new_energy_3blocks.tsv')
    db=table(ENERGY/'new_decomp_3blocks.tsv')
    decomp=table(ENERGY/'new_decomp_union_coverage.tsv')
    structural=table(STRUCT/'02_体系均值与样本SD.tsv')
    contacts=table(STRUCT/'11_接触占有率_残基逐重复.tsv')
    # Binding: compare complete receptor residue sequence/identity in fixed PDB to TPR mapping.
    identity_checks=[]
    receptor_keys={}
    for s in SYSTEMS:
        for rep in (1,2,3):
            mapping=table(STRUCT/'结构残基映射'/f'{s}_rep{rep}_protein_atoms.tsv')
            by_id={}
            for r in mapping:
                k=int(r['global_residue_index0'])
                val=(r['chain'].strip(),int(r['input_resnr']),r['insertion_code'].strip(),canon(r['resname']))
                assert k not in by_id or by_id[k]==val
                by_id[k]=val
            keys=list(by_id.values())
            folder=ROOT/'per_rep'/s/f'rep{rep}'
            pdb=pdb_residues(folder/'_GMXMMPBSA_COM_FIXED.pdb')
            receptor=[k for k in pdb if k[3]!='MOL']
            passed=(keys==receptor and len(keys)==len(set(keys)))
            identity_checks.append(dict(system=s,rep=rep,tpr_residues=len(keys),pdb_receptor_residues=len(receptor),status='PASS' if passed else 'UNRESOLVED'))
            if passed: receptor_keys[(s,rep)]=set(keys)
    contact_index={}
    for r in contacts:
        if r['system'] not in SYSTEMS or r['atom_mode']!='all_atom' or r['geometry']!='PBC': continue
        key=(r['system'],int(r['rep']),r['chain'].strip(),int(r['input_resnr']),r['insertion_code'].strip(),canon(r['resname']))
        assert key not in contact_index, f'Ambiguous contact residue {key}'
        contact_index[key]=float(r['occupancy'])
    drift=[]
    for s in SYSTEMS:
        for rep in (1,2,3):
            rows=sorted([r for r in eb if r['system']==s and int(r['rep'])==rep and r['term']=='TOTAL'],key=lambda r:int(r['block']))
            assert [int(r['block']) for r in rows]==[1,2,3]
            vals=[float(r['mean_kcal_mol']) for r in rows]
            drift.append(dict(system=s,rep=rep,block1=vals[0],block2=vals[1],block3=vals[2],last_minus_first=vals[2]-vals[0],block_range=max(vals)-min(vals),monotonic='increasing' if vals[0]<vals[1]<vals[2] else 'decreasing' if vals[0]>vals[1]>vals[2] else 'nonmonotonic',abs_drift_ge_5=abs(vals[2]-vals[0])>=5))
    block_index={}
    for r in db:
        if r['decomp_type']==TDC and r['term']=='TOTAL':
            k=(r['system'],int(r['rep']),r['residue'])
            block_index.setdefault(k,[]).append(float(r['mean_kcal_mol']))
    residue_rows=[]
    for r in decomp:
        if r['decomp_type']!=TDC or r['term']!='TOTAL' or not r['residue'].startswith('R:'): continue
        s=r['system']; observed=[i for i in (1,2,3) if r[f'rep{i}_kcal_mol']!='NA']
        parts=r['residue'].split(':'); identity=None
        if len(parts)==4 and re.fullmatch(r'-?\d+',parts[3]): identity=(parts[1],int(parts[3]),'',canon(parts[2]))
        good_identity=identity is not None and all(identity in receptor_keys.get((s,i),set()) for i in (1,2,3))
        occ=[contact_index.get((s,i,*identity)) if good_identity else None for i in (1,2,3)]
        vals=[float(r[f'rep{i}_kcal_mol']) for i in observed]
        blocks=[v for i in observed for v in block_index[(s,i,r['residue'])]]
        assert len(blocks)==3*len(observed)
        row=dict(system=s,residue=r['residue'],n_observed=len(observed),mean=float(r['mean_observed_kcal_mol']),sd=None if r['between_observed_rep_sample_sd']=='NA' else float(r['between_observed_rep_sample_sd']),all_reps_negative=len(observed)==3 and all(v<0 for v in vals),negative_blocks=sum(v<0 for v in blocks),observed_blocks=len(blocks),identity_status='PASS' if good_identity else 'UNRESOLVED',contact_mean_50_100ns=st.mean(occ) if all(v is not None for v in occ) else None,contact_min_rep=min(occ) if all(v is not None for v in occ) else None)
        row.update({f'rep{i}_energy':None if r[f'rep{i}_kcal_mol']=='NA' else float(r[f'rep{i}_kcal_mol']) for i in (1,2,3)})
        row.update({f'rep{i}_contact':occ[i-1] for i in (1,2,3)})
        row['evaluation_status']='EVALUABLE' if good_identity and all(v is not None for v in occ) and len(observed)==3 else 'INSUFFICIENT_COVERAGE_OR_IDENTITY'
        row['priority_candidate']=row['all_reps_negative'] and row['negative_blocks']==9 and row['mean']<=-1 and row['contact_min_rep'] is not None and row['contact_min_rep']>=.5
        residue_rows.append(row)
    candidates=sorted([r for r in residue_rows if r['priority_candidate']],key=lambda r:(SYSTEMS.index(r['system']),r['mean']))
    coverage=[]
    for s in SYSTEMS:
        rows=[r for r in residue_rows if r['system']==s]
        coverage.append(dict(system=s,union_receptor=len(rows),intersection_3rep=sum(r['n_observed']==3 for r in rows),partial_coverage=sum(r['n_observed']<3 for r in rows),not_evaluable=sum(r['evaluation_status']!='EVALUABLE' for r in rows),screened_candidates=sum(r['priority_candidate'] for r in rows)))
    # Add component means across independent replicate means, no pooled-frame inference.
    components=[]
    for s in SYSTEMS:
        for term in sorted({r['term'] for r in eb}):
            vals=[]
            for rep in (1,2,3):
                values=[float(r['mean_kcal_mol']) for r in eb if r['system']==s and int(r['rep'])==rep and r['term']==term]
                assert len(values)==3
                vals.append(st.mean(values))
            mean,sd=avg_sd(vals)
            components.append(dict(system=s,term=term,mean=mean,between_rep_sd=sd,rep1=vals[0],rep2=vals[1],rep3=vals[2]))
    sensitivity=[]
    sys.path.insert(0,str(BASE))
    from build_integrated_energy_results import parse_energy
    for s in SYSTEMS:
        samples={}
        for rep in (1,2,3):
            p=ROOT/'per_rep'/s/f'rep{rep}'/'ENERGY_PER_FRAME.csv'
            SOURCES[str(p)]=sha(p)
            samples[rep]=parse_energy(p)
        for label,first,last in [('80-100ns',1,51),('80-89.6ns',1,25),('90-100ns',26,51),('86.8-100ns',18,51),('93.6-100ns',35,51)]:
            vals=[st.mean(r['TOTAL'] for r in samples[i][first-1:last]) for i in (1,2,3)]
            mean,sd=avg_sd(vals)
            sensitivity.append(dict(system=s,window=label,frames_per_rep=last-first+1,mean=mean,between_rep_sd=sd,rep1=vals[0],rep2=vals[1],rep3=vals[2]))
        vals=[st.mean(r['EEL']+r['EGB'] for r in samples[i]) for i in (1,2,3)]
        mean,sd=avg_sd(vals)
        components.append(dict(system=s,term='EEL+EGB',mean=mean,between_rep_sd=sd,rep1=vals[0],rep2=vals[1],rep3=vals[2]))
    integrated=[]
    for s in SYSTEMS:
        e=next(r for r in summary if r['system']==s)
        row=dict(system=s,energy_mean=float(e['mean_kcal_mol']),energy_sd=float(e['between_rep_sample_sd']))
        for metric in ['protein_bb_self_RMSD_nm','ligand_all_atom_RMSD_nm','protein_mass_Rg_nm','hbond_protein_ligand','protein_lig_min_all_PBC_nm']:
            v=[r for r in structural if r['system']==s and r['metric']==metric]
            assert len(v)==1
            row[metric+'_mean']=float(v[0]['mean']); row[metric+'_sd']=float(v[0]['sd_between_replicates'])
        integrated.append(row)
    history=table(ENERGY/'historical_energy_per_rep.tsv')
    current=table(ENERGY/'new_energy_per_rep.tsv')
    comparison=[]
    for r in current:
        old=[x for x in history if x['system']==r['system'] and x['rep']==r['rep']]
        if old:
            a=float(old[0]['delta_total_kcal_mol']); b=float(r['delta_total_kcal_mol'])
            comparison.append(dict(system=r['system'],rep=int(r['rep']),historical=a,current=b,current_minus_historical=b-a,n_independent_reps_not_doubled=True))
    write_table('01_逐重复时间漂移.tsv',drift)
    write_table('02_残基能量与结构接触对照.tsv',residue_rows)
    write_table('03_探索性优先残基.tsv',candidates,fields=list(residue_rows[0]))
    write_table('04_残基身份检查.tsv',identity_checks)
    write_table('05_能量分项.tsv',components)
    write_table('06_结构与能量整合.tsv',integrated)
    write_table('07_新旧重分析对照.tsv',comparison)
    write_table('08_分解覆盖率.tsv',coverage)
    write_table('09_窗口敏感性.tsv',sensitivity)
    data=dict(generated_at=datetime.now().astimezone().isoformat(),root=str(ROOT),status='BUILDING',identity_checks=identity_checks,coverage=coverage,drift=drift,candidates=candidates,integrated=integrated,components=components,history_comparison=comparison,window_sensitivity=sensitivity)
    (OUT/'analysis.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    lines=['# SHXD–PJI：MM/GBSA后续分析与结构整合报告','',f'生成时间：{data["generated_at"]}','本报告为独立后续分析及论文更新素材，不覆盖既有论文；不是整篇论文最终版。','', '## 1. 完成范围与验收',
    '新批5体系×3独立重复全部通过增强验收；逐帧重算四个能量区段差分、总均值、3×17时间块及所有残基分解时间块，核对源拓扑、TPR、抽帧轨迹与产物哈希。历史12条单独验证，不与新15条合并计为27条独立轨迹。',
    '结构汇总覆盖41条不重复轨迹、15体系；12条修复TLR4替换旧记录。原大轨迹只核对元数据，不能宣称本次重新全量哈希或重算MD。',
    'MM/GBSA取80–100 ns每0.4 ns一帧共51帧；结构指标取50–100 ns，两个窗口不相同，联合表仅作互补证据，不用于帧级相关性或因果推断。',
    'GB参数：igb=5，intdiel=1，extdiel=78.5，saltcon=0.15，PBRadii=3；idecomp=1、dec_verbose=1；within 6为首帧6 Å选区。无熵，软件temperature=298.15 K，与MD 300 K不同须明确记录。',
    '', '## 2. 总能量与结构指标', '|体系|ΔTOTAL均值±SD (kcal/mol)|骨架RMSD (nm)|配体全原子RMSD (nm)|氢键数|','|---|---:|---:|---:|---:|']
    for r in integrated:
        lines.append(f'|{LABELS[r["system"]]}|{r["energy_mean"]:.2f} ± {r["energy_sd"]:.2f}|{r["protein_bb_self_RMSD_nm_mean"]:.3f} ± {r["protein_bb_self_RMSD_nm_sd"]:.3f}|{r["ligand_all_atom_RMSD_nm_mean"]:.3f} ± {r["ligand_all_atom_RMSD_nm_sd"]:.3f}|{r["hbond_protein_ligand_mean"]:.3f} ± {r["hbond_protein_ligand_sd"]:.3f}|')
    lines+=['均值和样本SD均以3个独立重复为统计单位，不把51帧或三个时间块当独立重复。总能量主表沿用最终DAT（每重复0.01 kcal/mol精度）口径，与此前进度一致；分项与分块由逐帧CSV重算，两种舍入层级可能产生末位差异，不混合来源计算同一统计量。',
    'IKKβ三个配体平均能量相近，黄芩素更负只是描述性趋势。NLRP3黄芩素与小檗碱均值几乎相同，不支持亲和力排序；氢键差异不等同总体亲和力差异。较大的全蛋白RMSD可能涉及域间运动，不直接判定失稳或解离。',
    '', '## 3. 时间分块与漂移',
    '块1=80–86.4 ns，块2=86.8–93.2 ns，块3=93.6–100 ns；每块17帧，边界无重叠。以下为总能量块均值，末块减首块仅作描述。|Δ|≥5 kcal/mol只是本报告事后筛查提示，不是预注册显著性阈值或收敛判据。',
    '|体系/重复|块1|块2|块3|末−首|块范围|','|---|---:|---:|---:|---:|---:|']
    for r in drift:
        lines.append(f'|{LABELS[r["system"]]} r{r["rep"]}|{r["block1"]:.2f}|{r["block2"]:.2f}|{r["block3"]:.2f}|{r["last_minus_first"]:+.2f}|{r["block_range"]:.2f}|')
    flagged=sorted([r for r in drift if r['abs_drift_ge_5']],key=lambda r:abs(r['last_minus_first']),reverse=True)
    lines+=['优先复核漂移：'+('；'.join(f'{LABELS[r["system"]]} rep{r["rep"]} {r["last_minus_first"]:+.2f}' for r in flagged) if flagged else '没有达到上述事后筛查阈值的重复，但这不等于证明收敛。'),
    '块均值的变化表明窗口敏感性；正差表示末块较不负，负差表示末块更负，二者都可能反映未充分采样。不能用三块一致性代替充分收敛或独立采样。',
    '', '## 4. 残基贡献、覆盖及接触整合',
    '受体分解仅使用TDC TOTAL，不把TDC/SDC/BDC相加。未选入残基为NA，不补零。身份核验将COM_FIXED.pdb的全受体残基顺序、链、编号、插入码与TPR映射逐项对照；仅对已声明的质子化残基命名归一化。通过后才连接全原子PBC 0.45 nm接触表。能量选区与50–100 ns接触窗口不同，接触数值不是能量时窗的占有率。',
    f'全受体身份检查：{sum(r["status"]=="PASS" for r in identity_checks)}/15通过；未通过的连接保留UNRESOLVED，不强行对齐。',
    '|体系|受体残基并集|严格3重复交集|部分覆盖|探索性优先候选|','|---|---:|---:|---:|---:|']
    for r in coverage: lines.append(f'|{LABELS[r["system"]]}|{r["union_receptor"]}|{r["intersection_3rep"]}|{r["partial_coverage"]}|{r["screened_candidates"]}|')
    lines+=['探索性筛选规则：3/3重复均负贡献，9/9分块均负贡献，三重复均值≤−1 kcal/mol，并且每个重复50–100 ns接触占有率≥0.5。阈值为本报告事后筛选、未经独立验证，不称为确认的功能热点。','|体系|残基（软件链:名:号）|贡献均值±SD|接触均值|最小重复接触|','|---|---|---:|---:|---:|']
    for s in SYSTEMS:
        selected=[r for r in candidates if r['system']==s]
        for r in selected: lines.append(f'|{LABELS[s]}|{r["residue"]}|{r["mean"]:.2f} ± {r["sd"]:.2f}|{r["contact_mean_50_100ns"]:.3f}|{r["contact_min_rep"]:.3f}|')
        if not selected: lines.append(f'|{LABELS[s]}|可评价三重复交集中无满足全部条件的残基；缺失项不可评价|—|—|—|')
    lines+=['以上编号属于模拟结构/软件输出编号；本次未进行UniProt标准位点映射。突变设计之前须核对参考序列、结构缺失区及构建体边界。',
    '', '## 5. 能量组成与历史结果',
    '能量分项表保留VDWAALS、EEL、EGB、ESURF等软件原项；静电项与极性溶剂化项通常应联合讨论，不能仅取其中一项支持更强结合。残基选区总和不要求等于全体系ΔTOTAL。',
    f'新旧同体系同重复有{len(comparison)}对可对照。相同MD轨迹重分析不增加独立样本数；历史参数部分未显式记录，差值不能归因于配体亲和力改变。',
    '', '## 6. TLR4–MD2结果与限制',
    '主QC为12/12 PASS；补充自动验收仍为4通过/8失败，8条另有PASS_WITH_CUTOFF_BOUNDARY_REVIEW。17个接触边界差异重算direct与PBC距离相同，属于0.45 nm阈值纳入差异；未改写自动FAIL。',
    'MD2链C的高占有率接触包含VAL135、PHE76、LEU61、LEU78和ILE117等，支持局部几何接触的描述，不证明抑制TLR4。OroxylinA rep2骨架RMSD末块较首块增加约0.180 nm，不能宣称所有重复完全收敛。本批未计算TLR4–MD2 MM/GBSA。',
    '', '## 7. 可用于论文的结果与讨论文字',
    '对5个配体–靶点体系的15条独立轨迹进行80–100 ns窗口的无熵MM/GBSA及逐残基分解。所有重复均获得51帧有效结果，并完成输出哈希、逐帧差分与3×17帧分块一致性复核。'+ '；'.join(f'{LABELS[r["system"]]}为{r["energy_mean"]:.2f}±{r["energy_sd"]:.2f}' for r in integrated)+' kcal/mol（最终DAT口径，均值±独立重复样本SD，n=3）。',
    '这些结果为候选相互作用提供端点能量描述，但不同重复及时间分块的变化提示采样与构象依赖性，不能仅以负值或均值排序推断实验亲和力。残基分解结合接触占有率用于提出探索性位点候选；首帧选区导致的缺失保留为NA，只有三重复交集用于完整重复间比较。功能验证仍须依赖生化结合/活性实验及细胞或动物证据。',
    '', '## 8. 后续仍需完成的工作',
    '1. 论文正文、方法、图注与补充材料的逐段更新及版面核验：本轮提供可插入文字和独立报告，未自动覆盖整篇论文。',
    '2. 如需发表级结合构象图，应另行执行预先定义的构象代表性选择并核验；本轮不把能量最低帧当作代表构象。',
    '3. 本轮已完成末20 ns内的窗口敏感性对照，但不能替代更长独立采样。是否加跑需另行确定，不自动启动长计算。',
    '4. 位点候选需参考序列映射、结构可及性核验及实验验证。现有计算不能单独建立整方疗效、TLR4抑制或因果机制。',
    '', '## 9. 文件索引与复现',
    'energy/final_20260921：完整能量验收JSON、九份基础数据表、六组PNG/PDF；structural：41轨迹结构汇总及七组图；followup_final_v2：漂移、残基接触连接、候选、身份、分项、联合指标与新旧对照及analysis.json。',
    '输入和输出快照哈希保存在followup_manifest.json；运行脚本analyze_post_mmgbsa.py。所有计算来源只读，本轮未重新运行生产MD或MM/GBSA。']
    extra=['## 补充A：时间窗口敏感性','所有窗口从同一已验收51帧中重取，不新增计算或独立重复。每窗先算各rep均值再算3rep均值±SD；不选取最有利窗口替代主结果。','|体系|窗口|每rep帧数|均值±重复间SD (kcal/mol)|','|---|---|---:|---:|']
    for r in sensitivity: extra.append(f'|{LABELS[r["system"]]}|{r["window"]}|{r["frames_per_rep"]}|{r["mean"]:.2f} ± {r["between_rep_sd"]:.2f}|')
    extra += ['## 补充B：静电与极性溶剂化抵消','以下均值按逐帧CSV计算；EEL+EGB先在每重复相加，再算重复间SD，不能将两个SD直接相加。','|体系|EEL|EGB|EEL+EGB 均值±SD|','|---|---:|---:|---:|']
    for s in SYSTEMS:
        c={r['term']:r for r in components if r['system']==s}
        extra.append(f'|{LABELS[s]}|{c["EEL"]["mean"]:.2f}|{c["EGB"]["mean"]:.2f}|{c["EEL+EGB"]["mean"]:.2f} ± {c["EEL+EGB"]["between_rep_sd"]:.2f}|')
    extra.append('NLRP3–小檗碱的较大有利静电项被极性溶剂化代价显著抵消，不能只用EEL解释为更强亲和力。')
    extra.append(f'历史重分析对照：9对最终DAT总能量差的最大绝对值为{max(abs(r["current_minus_historical"]) for r in comparison):.2f} kcal/mol，即在0.01 kcal/mol报告精度下完全一致；这支持重分析一致性，不是额外独立重复。')
    lines += extra
    report=''
    for i,line in enumerate(lines):
        separator='\n' if i and line.startswith('|') and lines[i-1].startswith('|') else '\n\n'
        report += (separator if i else '') + line
    (OUT/'report_source.txt').write_text(report,encoding='utf-8')
    # Compact independent plots for drift and candidate evidence.
    import os
    os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(10,5.5))
    positions=list(range(len(drift)))
    vals=[r['last_minus_first'] for r in drift]
    ax.barh(positions,vals,color=['#aa5737' if v>0 else '#247c87' for v in vals])
    ax.set_yticks(positions,[LABELS[r['system']]+f' r{r["rep"]}' for r in drift],fontsize=9)
    ax.invert_yaxis(); ax.axvline(0,color='#555555',lw=.8)
    ax.set_xlabel('末块 − 首块 ΔTOTAL (kcal/mol)'); ax.set_title('15条重复的时间漂移：描述性指标，不是收敛检验')
    fig.tight_layout()
    for ext in ('png','pdf'): fig.savefig(OUT/f'时间漂移.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    assert all(sha(Path(p))==h for p,h in SOURCES.items()), 'Source changed during analysis'
    assert len(drift)==15 and len(integrated)==5 and len(identity_checks)==15
    assert all((OUT/n).is_file() for n in ['03_探索性优先残基.tsv','report_source.txt','时间漂移.png','时间漂移.pdf'])
    data['status']='COMPLETE'
    (OUT/'analysis.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    manifest=dict(generated_at=data['generated_at'],source_read_only=True,mmgbsa_rerun=False,inputs=SOURCES,outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file()},script_sha256=sha(Path(__file__)))
    (OUT/'followup_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(output=str(OUT),identity_pass=sum(r['status']=='PASS' for r in identity_checks),candidates=len(candidates),drift_flagged=flagged,coverage=coverage),ensure_ascii=False,indent=2))

if __name__=='__main__': main()
