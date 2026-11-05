#!/bin/bash
# Theta sweep on Hamilton: one job-array task per (aggregation, start angle), 70 tasks.
#   tasks  1-35  mode-separated ('separated'),  task k = 5*k degrees
#   tasks 36-70  mode-collapsed ('smoothmax'),  task k = 5*(k-35) degrees
# Submit FROM the 'Unified Hashin' folder. logs/ must exist first, or Slurm cannot open the output file and the task fails without a log:
#     cd "Unified Hashin"; mkdir -p logs
#     sbatch --array=18,53 hamilton_thetaSweep.sh       # test: the 90 deg start for both aggregations
#     sbatch hamilton_thetaSweep.sh                     # the whole sweep
# Partition and module are copied from comp_min_FEA_stress_const/myjob1.sh. Memory and walltime are estimates: after the test read
#     sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State
# and set them to about 1.3 x the slowest task (override without editing: sbatch --mem=6G -t 0-06:00:00 ...).

#SBATCH -J thetaSweep
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-12:00:00
#SBATCH --array=1-70%10
#SBATCH -o logs/%x_%A_%a.out

module load matlab/R2023a
cd "$SLURM_SUBMIT_DIR"           # thetaSweep finds Functions/ and Results/ relative to the working directory
echo "task $SLURM_ARRAY_TASK_ID on $(hostname), $SLURM_CPUS_PER_TASK cpu"
matlab -nodisplay -singleCompThread -batch "thetaSweep"   # one thread (we asked for one core); -batch exits when done and exits non-zero on error