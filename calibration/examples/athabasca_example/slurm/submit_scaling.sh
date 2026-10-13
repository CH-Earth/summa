#!/bin/bash
# Submit the scaling series; each run gets its own folder under my_work_athabasca/scaling_tests/<tag>/
# Usage:  bash submit_scaling.sh smoke | pilot | nodes
#   smoke : 1 node, 127 workers, 127 evaluations          -> scaling_tests/smoke_n1_w127/
#   pilot : 1 node, 16/32/64/127 workers, 256 evaluations   -> scaling_tests/pilot_n1_w016/ ... pilot_n1_w127/
#   nodes : 1/2/4/8 whole nodes (128n-1 workers), 5000 evals -> scaling_tests/nodes_n1_w127/ ... nodes_n8_w1023/
BASIN=<WORK_ROOT>/summa-spotpy/my_work_athabasca
cd $BASIN || exit 1
S=$BASIN/run_scaling.slurm
submit() {   # submit <nodes> <workers> <evals> <hours> <tag>
  local n=$1 w=$2 e=$3 h=$4 tag=$5 dir=$BASIN/scaling_tests/$5
  if [ -e "$dir/summary.json" ]; then echo "skip $tag: $dir already has results"; return; fi
  mkdir -p "$dir"
  sbatch -N $n --ntasks=$((w+1)) -t $(printf "%02d:00:00" $h) -J $tag \
         -o "$dir/slurm_%j.out" -e "$dir/slurm_%j.err" $S $e $tag
}
case "$1" in
  smoke) submit 1 127 127 2 smoke_n1_w127 ;;
  pilot) for w in 16 32 64 127; do
           case $w in 16) h=8;; 32) h=5;; 64) h=3;; 127) h=3;; esac   # limits from the develop smoke test
           submit 1 $w 256 $h $(printf "pilot_n1_w%03d" $w)
         done ;;
  nodes) for n in 1 2 4 8; do
           w=$((128*n-1))
           submit $n $w 5000 $(( 5000 / w * 40 / 60 + 4 )) nodes_n${n}_w${w}
         done ;;
  *) echo "usage: bash submit_scaling.sh smoke|pilot|nodes"; exit 1 ;;
esac
