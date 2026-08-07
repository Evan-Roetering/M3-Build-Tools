#!/bin/bash
# Script to calculate area of undulating membrane from .gro topology snapshots
# Last updated by:
#   Evan Roetering
#   Klauda Lab
#   University of Maryland
#   8/7/2026

sbatch << UNDULATION-AREA
#!/bin/bash
#SBATCH --job-name=undulation_area
#SBATCH --open-mode=append
#SBATCH -o ./slurmlogs/undulation_area.out
#SBATCH -e ./slurmlogs/undulation_area.err
#SBATCH -A energybio-eng
#SBATCH -p scavenger
#SBATCH -t 2-00:00:00
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --mail-user=your-email@umd.edu
#SBATCH --mail-type=all

conda activate martiniglass
python3 ./undulation_area.py
UNDULATION-AREA