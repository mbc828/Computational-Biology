# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""逐残基接触频率：解析 pocket_mol.pdb（PBC已校正），输出每残基接触频率 CSV
接触判据：该残基任一原子与 MOL 任一原子距离 < 4.5 Å
"""
import sys, csv
import numpy as np

def analyze(pdb_path, out_csv, cutoff=4.5):
    frames = []  # 每帧: (prot_reskeys, prot_coords list per res), mol_coords
    cur = {}
    mol = []
    def flush():
        if cur or mol:
            frames.append(({k: np.array(v) for k, v in cur.items()},
                           np.array(mol) if mol else None))
    with open(pdb_path, encoding='ascii', errors='ignore') as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'MODEL':
                cur, mol = {}, []
            elif rec in ('ATOM', 'HETATM'):
                resn = line[17:20].strip()
                resi = line[22:26].strip()
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                if resn == 'MOL':
                    mol.append((x, y, z))
                else:
                    cur.setdefault(f'{resn}{resi}', []).append((x, y, z))
            elif rec == 'ENDMDL':
                flush()
                cur, mol = {}, []
    flush()
    nf = len(frames)
    res_freq = {}
    for pres, mc in frames:
        if mc is None or len(mc) == 0:
            continue
        for rk, coords in pres.items():
            d = np.sqrt(((coords[:, None, :] - mc[None, :, :]) ** 2).sum(-1)).min()
            if d < cutoff:
                res_freq[rk] = res_freq.get(rk, 0) + 1
    res_freq = {k: v / nf for k, v in res_freq.items()}
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['residue', 'contact_freq'])
        for k in sorted(res_freq, key=lambda s: int(''.join(ch for ch in s if ch.isdigit()))):
            w.writerow([k, f'{res_freq[k]:.3f}'])
    return nf, len(res_freq)

if __name__ == '__main__':
    nf, nr = analyze(sys.argv[1], sys.argv[2])
    print(f'{sys.argv[1]}: frames={nf} pocket_residues={nr}')
