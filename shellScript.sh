#!/bin/bash
#SBATCH -J firstExample
#SBATCH -e /cluster/home/vlaferla/Logs/error_log.err
#SBATCH -o /cluster/home/vlaferla/Logs/logs.log
#SBATCH --partition=bigmem.24h
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

# Define the timing output file
TIMING_FILE="/cluster/home/vlaferla/fdsEulerJobs/64/job_timing.log"

# Record the start time
echo "Job started at $(date)" > "$TIMING_FILE"
start_time=$(date +%s)

# Change to the submission directory (if needed)
cd "$SLURM_SUBMIT_DIR"

# Ensure the executable has proper permissions
chmod +x /cluster/home/vlaferla/GSoC_FFD/bin/Release/net9.0/linux-x64/publish/GSoC_FFD

# Run the self-contained .NET application
/cluster/home/vlaferla/GSoC_FFD/bin/Release/net9.0/linux-x64/publish/GSoC_FFD

# Record the end time and calculate elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Job ended at $(date)" >> "$TIMING_FILE"
echo "Elapsed time: $(date -ud "@$elapsed_time" +'%H hours %M minutes %S seconds')" >> "$TIMING_FILE"
