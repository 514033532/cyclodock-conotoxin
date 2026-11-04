#!/bin/bash
# V4 补对接 mode A 批量驱动 v2 (2026-11-02)
export LD_PRELOAD=/home/zhh2/miniconda3/envs/RAPiDock/lib/libstdc++.so.6
HADDOCK=/home/zhh2/miniconda3/envs/pep_prot_calc/bin/haddock3
TEMPLATE=/home/zhh2/software/DrugCLIP_screen/haddock3-main/sample/work_dir/docking_template-ss.cfg
PEPDIR=/media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/benchmark_structures/pep_clean
POCK=/media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/pocket_pdbs_from_data
NETDIR=/media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/DrugCLIP_screen_pdb/NET
OUT=/media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/benchmark_docking/modeA
mkdir -p "$OUT"
P_CaV2=$(ls $POCK/7MIY-cav2.2.pdb_chain_B/*pocket_1_chainB.pdb 2>/dev/null | head -1)
P_NaV17=$(ls $POCK/6N4I-nav1.7.pdb_chain_B/*pocket_1_chainB.pdb 2>/dev/null | head -1)
P_NaV14=$(ls $POCK/6AGF-NaV1.4.pdb_chain_B/*pocket_1_chainB.pdb 2>/dev/null | head -1)
P_NaV18=$(ls $POCK/9DBK-nav1.8.pdb_chain_B/*pocket_1_chainB.pdb 2>/dev/null | head -1)
P_a3b4=$(ls $POCK/6PV7-CHRNA3-B4.pdb_chain_B/*pocket_1_chainB.pdb 2>/dev/null | head -1)
P_a7=$POCK/7EKO_CHRNA7.pdb_chain_AB/7EKO_CHRNA7_pocket_ECD-Cloop_chainAB.pdb
P_a4b2=$POCK/6CNJ_CHRNA4_B2.pdb_chain_AB/6CNJ_CHRNA4_B2_pocket_ECD-Cloop_chainAB.pdb
P_5ht3=$POCK/8AXD_HTR3A.pdb_chain_AB/8AXD_HTR3A_pocket_ECD-p7_chainAB.pdb
PAIRS="
MVIIA|MVIIA_1DW4|P_CaV2
MVIIA|MVIIA_1DW4|P_NaV17
MVIIA|MVIIA_1DW4|P_NaV14
GVIA|GVIA_2CCO|P_CaV2
GVIA|GVIA_2CCO|P_NaV17
CVID|CVID_1TTK|P_CaV2
KIIIA|KIIIA_2LXG|P_NaV17
KIIIA|KIIIA_2LXG|P_NaV14
GIIIA|GIIIA_1TCG|P_NaV14
GIIIA|GIIIA_1TCG|P_NaV17
CnIIIC|CnIIIC_2YEN|P_NaV14
CnIIIC|CnIIIC_2YEN|P_NaV17
MrVIB|MrVIB_1RMK|P_NaV18
MrVIB|MrVIB_1RMK|P_NaV17
ImI|ImI_1CNL|P_a7
ImI|ImI_1CNL|P_a4b2
ImI|ImI_1CNL|P_a3b4
GID|GID_1MTQ|P_a4b2
GID|GID_1MTQ|P_a3b4
AuIB|AuIB_1MXN|P_a3b4
AuIB|AuIB_1MXN|P_a7
TxIA|TxIA_2UZ6|P_a7
TxIA|TxIA_2UZ6|P_a3b4
Vc1.1|Vc1.1_2H8S|P_a7
Vc1.1|Vc1.1_2H8S|P_a4b2
MrIA|MrIA_2EW4|NET:NET1_DSP
MrIA|MrIA_2EW4|NET:NET2_BPP
MrIA|MrIA_2EW4|NET:NET3_68P
MrIA|MrIA_2EW4|NET:NET4_LNR
MrIA|MrIA_2EW4|NET:NET5_LDP
MrIA|MrIA_2EW4|NET:NET7_HJM
MrIA|MrIA_2EW4|NET:NET8_8PR
MrIA|MrIA_2EW4|NET:NET9_68P
MrIA|MrIA_2EW4|P_NaV17
GVIIIA|GVIIIA_9EBE|P_5ht3
"
JOBLIST=$(mktemp)
echo "$PAIRS" | while IFS='|' read -r tox pep tgt; do
  [ -z "$tox" ] && continue
  if [[ "$tgt" == NET:* ]]; then
    pdb="$NETDIR/${tgt#NET:}.pdb"; tname="NET_${tgt#NET:}"
  else
    pdb="${!tgt}"; tname="${tgt#P_}"
  fi
  if [ -z "$pdb" ] || [ ! -f "$pdb" ]; then echo "[跳过] $tox x $tgt (无结构: $pdb)"; continue; fi
  work="$OUT/${tox}__${tname}"
  if [ -d "$work/run1" ] && ls "$work/run1" | grep -q "clustfcc"; then echo "[已完成] $work"; continue; fi
  rm -rf "$work"; mkdir -p "$work"
  sed -e "s|__RUN_DIR__|$work/run1|g" -e "s|__PROTEIN_PDB__|$pdb|g" -e "s|__PEPTIDE_PDB__|$PEPDIR/${pep}.pdb|g" -e "s|^ncores = .*|ncores = 8|" -e "s|^tolerance = 10|tolerance = 30|" "$TEMPLATE" > "$work/run.cfg"
  echo "$work" >> "$JOBLIST"
done
NJOBS=$(wc -l < "$JOBLIST")
echo "共生成 $NJOBS 个作业"
FIRST=$(head -1 "$JOBLIST")
echo "冒烟测试: $FIRST"
( cd "$FIRST" && timeout 1800 $HADDOCK run.cfg > run.log 2>&1 )
if ! ls "$FIRST/run1" 2>/dev/null | grep -q clustfcc; then
  echo "❌ 冒烟测试失败/超时, 查看 $FIRST/run.log"; tail -25 "$FIRST/run.log"; exit 1
fi
echo "✅ 冒烟测试真完成, 启动剩余 $(($NJOBS-1)) 个作业 (3 并行)..."
tail -n +2 "$JOBLIST" | xargs -P 3 -I {} bash -c 'cd {} && { timeout 7200 /home/zhh2/miniconda3/envs/pep_prot_calc/bin/haddock3 run.cfg > run.log 2>&1; ls run1 2>/dev/null | grep -q clustfcc && echo "DONE {}" || echo "FAIL {}"; }'
echo "=== 全部完成 ==="
