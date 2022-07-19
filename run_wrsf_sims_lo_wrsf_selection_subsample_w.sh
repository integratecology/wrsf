#!/bin/bash
#SBATCH --job-name=wrsf_sims_sd    # name of the job
#SBATCH --partition=defq           # partition to be used (defq, gpu or intel)
#SBATCH --time=96:00:00            # walltime (up to 96 hours)
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks-per-node=1        # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=1          # number of cpus required to run the script
#SBATCH --mem-per-cpu=8G	   # memory required for process


module load gcc

export LD_LIBRARY_PATH="/home/alston92/software/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/gdal-3.3.0/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/proj-8.0.1/lib:$LD_LIBRARY_PATH"

ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/terra/libs/terra.so
ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/rgdal/libs/rgdal.so

module load R/3.6.3

cd /home/alston92/proj/wrsf   # where executable and data is located

date
echo "Initiating script"


if [ -f results/wrsf_sim_results_lo_wrsf_selection_subsample_w.csv ]; then
	echo "Results file already exists! continuing..."
else
	echo "creating results file"
	echo "sim_no,est,lcl,ucl" > results/wrsf_sim_results_lo_wrsf_selection_subsample_w.csv
fi

Rscript wrsf_sims_lo_wrsf_selection_subsample_w.R # name of script
echo "Script complete"
date
