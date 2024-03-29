#!/bin/bash
#SBATCH --job-name=caracal         # name of the job
#SBATCH --partition=defq           # partition to be used (defq, gpu or intel)
#SBATCH --time=96:00:00            # walltime (up to 96 hours)
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks-per-node=1        # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=1          # number of cpus required to run the script
#SBATCH --mem-per-cpu=32G         # memory required for process


module load gcc

export LD_LIBRARY_PATH="/home/alston92/software/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/gdal-3.3.0/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/alston92/software/proj-8.0.1/lib:$LD_LIBRARY_PATH"

ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/terra/libs/terra.so
ldd /home/alston92/R/x86_64-pc-linux-gnu-library/3.6/rgdal/libs/rgdal.so

module load R

cd /home/alston92/proj/wrsf   # where executable and data is located

date
echo "Initiating script"


if [ -f results/empirical_results.csv ]; then
        echo "Results file already exists! continuing..."
else
        echo "creating results file empirical_results.csv"
        echo "species,iid_pest,iid_plcl,iid_pucl,iid_uest,iid_ulcl,iid_uucl,w_pest,w_plcl,w_pucl,w_uest,w_ulcl,w_uucl" > results/empirical_results.csv
fi

Rscript caracal_example.R  # name of script
echo "Script complete"
date

