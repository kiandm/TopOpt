#!/bin/bash
# Separated vs smooth-max on the four-mode cantilever on Hamilton: one job-array task per (aggregation, load); 10 tasks with the defaults in scenarioPairs.m
#   tasks 1-5  = 'separated' at loads(1..5),   tasks 6-10 = 'smoothmax' at loads(1..5)     (edit loads / aggs / ndiv in scenarioPairs.m and keep --array in step)
# Submit FROM the folder that holds scenarioPairs.m and Functions/ ("Hamilton_Copy_of_Unified Hashin"). logs/ must exist first or the task fails without a log:
#     mkdir -p logs
#     sbatch --array=1,6 hamilton_scenarioPairs.sh        # test: the first load for both aggregations
#     sbatch hamilton_scenarioPairs.sh                    # all ten
# When they have finished, with Results/ in place:     matlab -nodisplay -batch "summarisePairs(50, {'separated','smoothmax'}, false)"
# Partition, module and the one-thread flags are copied from hamilton_thetaSweep.sh. Memory and walltime are estimates: after the test read
#     sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State
# and set them to about 1.3 x the slowest task (override without editing: sbatch --mem=4G -t 0-01:00:00 ...).

#SBATCH -J scenarioPairs
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 0-02:00:00
#SBATCH --array=1-10%10
#SBATCH -o logs/%x_%A_%a.out

module load matlab/R2023a
cd "$SLURM_SUBMIT_DIR"           # scenarioPairs finds Functions/ and Results/ relative to the working directory
echo "task $SLURM_ARRAY_TASK_ID on $(hostname), $SLURM_CPUS_PER_TASK cpu"
matlab -nodisplay -singleCompThread -batch "scenarioPairs"   # scenarioPairs.m sits in the top folder, so it is found without touching the path
