#!/bin/bash
export LD_PRELOAD=/home/zhh2/miniconda3/envs/RAPiDock/lib/libstdc++.so.6
export HADDOCK=/home/zhh2/miniconda3/envs/pep_prot_calc/bin/haddock3
cat /media/zhh2/zhongxin4T01/CONO-PEPTIDE2026/benchmark_docking/ptm_sens/jobs.txt | xargs -P 2 -I {} bash -c 'cd {} && { [ -d run1 ] && ls run1 | grep -q clustfcc && echo "SKIP {}" && exit 0; rm -rf run1; timeout 7200 $HADDOCK run.cfg > run.log 2>&1; ls run1 2>/dev/null | grep -q clustfcc && echo "DONE {}" || echo "FAIL {}"; }'
echo "=== crossover 完成 ==="
