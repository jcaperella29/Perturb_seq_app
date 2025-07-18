#!/bin/bash
#SBATCH --job-name=perturb-shiny
#SBATCH --output=perturb-shiny_%j.log
#SBATCH --error=perturb-shiny_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=08:00:00
 #SBATCH --partition=debug
or compute, depending on your cluster

module load singularity   # If your cluster uses modules

IMAGE=/mnt/c/Users/jcape/Perturb_seq_app/CRISPR_APP/mixscape_all_deps.sif
APPDIR=/mnt/c/Users/jcape/Perturb_seq_app/CRISPR_APP

echo "Running on $(hostname)"
echo "Shiny app will be served at port 3838 on this node."

echo "If connecting from outside, SSH tunnel like:"
echo "ssh -N -L 3838:$(hostname):3838 your_user@cluster"

singularity run --bind $APPDIR:/app $IMAGE

