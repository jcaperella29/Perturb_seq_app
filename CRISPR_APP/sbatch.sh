#!/bin/bash
#SBATCH --job-name=shiny-mixscape
#SBATCH --output=shiny-mixscape_%j.log
#SBATCH --error=shiny-mixscape_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=interactive   # or compute, or whatever your cluster has

# Load singularity module if needed
module load singularity

# Set paths
IMAGE=/path/to/mixscape_all_deps.sif
APPDIR=/path/to/your/app

# If you want to bind a custom port, you can change 3838
PORT=3838

# Optional: print the node and port
echo "Running on $(hostname) at port $PORT"
echo "If using Shiny, tunnel with:"
echo "ssh -N -L 3838:$(hostname):3838 your_user@cluster"

# Run the container (foreground, don't detach)
singularity run --bind $APPDIR:/srv/shiny-server/app $IMAGE

