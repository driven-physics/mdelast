#!/bin/bash
# One-core Slurm array running benchmarks S, Z1 and L1 and verifying each.
# Submit from the repository root:   mkdir -p benchmark_output && sbatch benchmarks/slurm_example.sh
# Adapt the partition name and the environment setup to your cluster.
#SBATCH --job-name=mdelast_bench
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=benchmark_output/slurm_%a.out

export OMP_NUM_THREADS=1
BUILD=build/Release
OUT=benchmark_output
mkdir -p "$OUT"

case "$SLURM_ARRAY_TASK_ID" in
  0) python3 benchmarks/smoke/run_smoke.py --build-dir "$BUILD" --output-dir "$OUT/smoke" \
       && python3 benchmarks/verify.py smoke "$OUT/smoke" ;;
  1) python3 sample_case/run.py --build-dir "$BUILD" --output-dir "$OUT/zn06_300K" \
       && python3 benchmarks/verify.py zn06_300K "$OUT/zn06_300K" ;;
  2) python3 benchmarks/l12_cuzn3_5K/run.py --build-dir "$BUILD" --output-dir "$OUT/l12_cuzn3_5K" \
       && python3 benchmarks/verify.py l12_cuzn3_5K "$OUT/l12_cuzn3_5K" ;;
esac
