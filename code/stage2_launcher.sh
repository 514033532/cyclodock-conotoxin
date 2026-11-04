#!/bin/bash
export LD_PRELOAD=/home/zhh2/miniconda3/envs/RAPiDock/lib/libstdc++.so.6
export HADDOCK=/home/zhh2/miniconda3/envs/pep_prot_calc/bin/haddock3
cat /tmp/chrna3_list.txt | sort -u > /tmp/stage2_list.txt
N=$(wc -l < /tmp/stage2_list.txt)
echo "Stage2 共 $N 个作业, 4 并行"
cat /tmp/stage2_list.txt | xargs -P 4 -I {} bash -c 'cd {} && { rm -rf run1; timeout 7200 $HADDOCK run.cfg > run.log 2>&1; ls run1 2>/dev/null | grep -q clustfcc && echo "DONE {}" || echo "FAIL {}"; }'
echo "=== Stage2 完成 ==="
