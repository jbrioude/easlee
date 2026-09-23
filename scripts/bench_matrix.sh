#!/bin/sh
# BENCHMARK MATRIX for the EASLEE paper (armed 2026-09-22, Jerome: "once the bugs are
# fixed, start doing the benchmarks so the numbers in the paper are updated").
#
# Waits for a QUIET machine (no other julia), then runs the ETEX 1 M /
# 2-day protocol (docs/BENCHMARKS.md) sequentially — never two timing runs at once —
# and appends wall clock + cost centres per run to a summary. Each run records the
# commit, the backend actually selected (from its own log), PREFETCH and the governor.
J=${JULIA:-julia}; W=$(cd "$(dirname "$0")/.." && pwd); D=${EASLEE_RESULTS:-$W/results}   # override with JULIA= and EASLEE_RESULTS=
STAMP=$(date '+%Y%m%d'); OUT=$D/bench_matrix_$STAMP; mkdir -p $OUT; S=$OUT/summary.md; LOG=/tmp/bench_matrix.log
say() { echo "$(date '+%F %T') $*" >> $LOG; }
w() { while pgrep -f "julia .*configs/$1\.toml" > /dev/null; do sleep 30; done; }
# (the developer tree gated this on a running chemistry job and a marker file; the distribution starts at once)
# quiet machine: no other julia at all
while pgrep -x julia > /dev/null; do say "waiting: another julia is running"; sleep 300; done
cd $W; COMMIT=$(git rev-parse --short HEAD); GOV=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)
{ echo "# EASLEE benchmark matrix — $STAMP"; echo; echo "commit $COMMIT · governor $GOV · zion (Threadripper 3970X, RTX 4080) · ETEX 1 M particles, 2 days, 15 min steps"; echo;
  echo "| run | backend (log) | threads | PREFETCH | TOTAL(run) s | unit wall s | metread | advection | output | release |"; echo "|---|---|---|---|---|---|---|---|---|---|"; } > $S
say "starting matrix at commit $COMMIT, governor $GOV"
run_one() {  # tag base_cfg advance threads prefetch
  tag=$1; base=$2; adv=$3; thr=$4; pf=$5
  cfg=configs/bench_matrix/$tag.toml
  sed "s|^path *= *\"results/bench_[A-Za-z_0-9]*/\"|path = \"results/bm_$tag/\"|" configs/$base.toml > $cfg
  mkdir -p $D/bm_$tag; rm -f $D/bm_$tag/*; ln -sfn $D/bm_$tag results/bm_$tag
  t0=$(date +%s)
  systemd-run --user --unit=easlee-bm-$tag -p MemoryHigh=60G -p MemoryMax=80G \
    -E EASLEE_ADVANCE=$adv -E PREFETCH=$pf -E HDF5_USE_FILE_LOCKING=FALSE --working-directory=$W \
    sh -c "$J -t $thr --project=easlee scripts/run_config.jl $cfg > /tmp/bm_$tag.log 2>&1"
  sleep 60; w bench_matrix/$tag; t1=$(date +%s)
  L=/tmp/bm_$tag.log; mode=$(grep -oE 'mode = :[a-z_]+' $L | head -1 | sed 's/mode = ://')
  tot=$(grep -E '^\s*TOTAL\(run\)' $L | awk '{print $2}')
  g() { grep -E "^\s*$1 " $L | tail -1 | awk '{print $2}'; }
  ok=$(grep -c 'done →' $L)
  echo "| $tag | ${mode:-?} | $thr | $pf | ${tot:-FAILED} | $((t1-t0)) | $(g metread) | $(g advection) | $(g output) | $(g release) |" >> $S
  cp $L $OUT/$tag.log
  say "$tag: mode=$mode total=${tot:-FAILED} wall=$((t1-t0))s done=$ok"
}
# short runs first, the 1-thread run last, so most of the table exists early
run_one cpu_t64      bench_cpu_etex_1M threads 64   0
run_one cpu_t32      bench_cpu_etex_1M threads 32   0
run_one cpu_t16      bench_cpu_etex_1M threads 16   0
run_one cpu_t8       bench_cpu_etex_1M threads 8    0
run_one cpu_t4       bench_cpu_etex_1M threads 4    0
run_one gpu_t4       bench_gpu_etex_1M ka_gpu  4    0
run_one gpu_t16i     bench_gpu_etex_1M ka_gpu  16,1 0
run_one cpu_t16_pf   bench_cpu_etex_1M threads 16   1
run_one cpu_t64_pf   bench_cpu_etex_1M threads 64   1
run_one gpu_t16i_pf  bench_gpu_etex_1M ka_gpu  16,1 1
run_one gpu_t4_pf    bench_gpu_etex_1M ka_gpu  4    1
run_one cpu_t2       bench_cpu_etex_1M threads 2    0
run_one cpu_t1       bench_cpu_etex_1M threads 1    0
say "matrix finished -> $S"
