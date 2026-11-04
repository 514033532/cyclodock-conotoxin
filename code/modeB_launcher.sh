#!/bin/bash
export LD_PRELOAD=/home/zhh2/miniconda3/envs/RAPiDock/lib/libstdc++.so.6
export HADDOCK=/home/zhh2/miniconda3/envs/pep_prot_calc/bin/haddock3
ls -d /media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/benchmark_docking/modeB/*__* > /tmp/modeB_list.txt
N=$(wc -l < /tmp/modeB_list.txt)
echo "模式B 共 $N 个作业, 3 并行启动"
cat /tmp/modeB_list.txt | xargs -P 3 -I {} bash -c 'cd {} && { rm -rf run1; timeout 7200 $HADDOCK run.cfg > run.log 2>&1; ls run1 2>/dev/null | grep -q clustfcc && echo "DONE {}" || echo "FAIL {}"; }'
echo "=== 模式B 全部完成 ==="
