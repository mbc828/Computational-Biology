#!/usr/bin/env bash
# ============================================================
# 04_build_and_md.sh — 构建模拟体系并跑 平衡+100ns生产（单体系单重复）
# 用法: bash scripts/04_build_and_md.sh <SYSTEM> <REPLICA>
#   例: bash scripts/04_build_and_md.sh S1_NLRP3_OroxylinA 1
#       bash scripts/04_build_and_md.sh S6_MD2_apo 1
# 力场: amber99sb-ildn(蛋白) + GAFF2(配体) + TIP3P水 + 0.15M NaCl
# ============================================================
set -euo pipefail
cd "$(dirname "$0")/.."
SYS=$1; REP=$2
source systems/systems.env   # 由 scripts/run_all.sh 生成；或手动定义下方变量
: "${PROT:?}" "${LIG:=none}" # PROT=NLRP3|MD2|IKKb; LIG=Baicalein|Wogonin|OroxylinA|none

RUN="runs/${SYS}/rep${REP}"
mkdir -p "${RUN}" && cd "${RUN}"
FF=amber99sb-ildn
WATER=tip3p

echo "== [${SYS} rep${REP}] 蛋白拓扑 =="
gmx pdb2gmx -f "../../../prepared/${PROT}_protein.pdb" \
    -o protein.gro -p topol_protein.top -i posre_protein.itp \
    -ff ${FF} -water ${WATER} -ignh -merge no

if [ "${LIG}" != "none" ]; then
  echo "== 合并配体（对接姿势） =="
  # 对接第1模型(配体坐标) 与 参数化gro 对齐后拼接
  python ../../../scripts/merge_complex.py \
      --protein_gro protein.gro \
      --ligand_pose "../../../docking/${SYS}/docked.pdb" \
      --ligand_gro "../../../params/${LIG}_ini.gro" \
      --out complex.gro
  cat topol_protein.top > topol.top
  sed -i "s/; Include forcefield parameters/; Include forcefield parameters\n#include \"..\/..\/..\/params\/${LIG}_GMX.itp\"/" topol.top
  echo "${LIG}   1" >> topol.top
  # 配体位置限制
  gmx genrestr -f "../../../params/${LIG}_ini.gro" -o posre_lig.itp -fc 1000 1000 1000 <<< 0
  sed -i "/#include \"..\/..\/..\/params\/${LIG}_GMX.itp\"/a #ifdef POSRES_LIG\n#include \"posre_lig.itp\"\n#endif" topol.top
else
  cp protein.gro complex.gro
  cp topol_protein.top topol.top
fi

echo "== 溶剂盒/离子（0.15 M NaCl） =="
BOXD=1.0   # 蛋白表面到盒壁最小距离 nm
gmx editconf -f complex.gro -o boxed.gro -c -d ${BOXD} -bt dodecahedron
gmx solvate -cp boxed.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f "../../../mdp/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15 <<< SOL

echo "== 能量最小化 =="
gmx grompp -f "../../../mdp/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -deffnm em -ntmpi 1

echo "== NVT 100 ps =="
gmx grompp -f "../../../mdp/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt -ntmpi 1

echo "== NPT 1 ns =="
gmx grompp -f "../../../mdp/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt -ntmpi 1

echo "== 生产 100 ns =="
gmx grompp -f "../../../mdp/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2
gmx mdrun -deffnm md -ntmpi 1 -cpi md.cpt

echo "[OK] ${SYS} rep${REP} 完成 -> ${RUN}/md.xtc"
