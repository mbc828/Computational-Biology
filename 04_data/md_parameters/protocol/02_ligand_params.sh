#!/usr/bin/env bash
# ============================================================
# 02_ligand_params.sh — 配体GAFF2参数化（AM1-BCC电荷）+ GROMACS拓扑
# 输入: ligands_raw/{Name}_{CID}_3d.sdf
# 输出: params/{Name}_GMX.itp, params/{Name}_GMX.top, params/{Name}_ini.gro
# 依赖: ambertools(antechamber+parmchk2), acpype, openbabel
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p params

declare -A LIGS=( [Baicalein]=5281605 [Wogonin]=5281703 [OroxylinA]=5320315 )

for NAME in "${!LIGS[@]}"; do
  CID=${LIGS[$NAME]}
  SDF="ligands_raw/${NAME}_${CID}_3d.sdf"
  echo "== ${NAME} (CID ${CID}) =="

  # 1) SDF -> mol2（三角锥/芳香性检查以mol2为准）
  obabel "${SDF}" -O "params/${NAME}.mol2" --title "${NAME}"

  # 2) antechamber: GAFF2原子类型 + AM1-BCC电荷
  antechamber -i "params/${NAME}.mol2" -fi mol2 \
              -o "params/${NAME}_bcc.mol2" -fo mol2 \
              -c bcc -s 2 -at gaff2 -nc 0 -m 1 -rn "${NAME^^:0:3}" \
              -pf yes -dr no

  # 3) parmchk2: 缺失力场参数 -> frcmod
  parmchk2 -i "params/${NAME}_bcc.mol2" -f mol2 \
           -o "params/${NAME}.frcmod" -s gaff2

  # 4) acpype: AMBER -> GROMACS 拓扑（含原子电荷）
  mkdir -p "params/${NAME}.acpype"
  acpype -i "params/${NAME}_bcc.mol2" -n 0 -a gaff2 -b "${NAME}" \
         -o gmx -c user 2>/dev/null || \
  acpype -i "params/${NAME}_bcc.mol2" -n 0 -b "${NAME}" -o gmx
  cp "${NAME}.acpype/${NAME}_GMX.itp" params/ 2>/dev/null || \
  cp "params/${NAME}.acpype/${NAME}_GMX.itp" params/
  cp "${NAME}.acpype/${NAME}_GMX.gro" "params/${NAME}_ini.gro" 2>/dev/null || true
  rm -rf "${NAME}.acpype"
done
echo "[OK] 配体参数化完成 -> params/"
echo "QA: 检查 params/*.frcmod 中 'ATTN, need revision' 行数，>0 需人工确认二面角罚分"
