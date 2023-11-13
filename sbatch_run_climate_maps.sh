#!/bin/bash
#SBATCH --job-name=climate_maps
#SBATCH --output=climate_maps_%j.out
#SBATCH --error=climate_maps_%j.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=80
#SBATCH --partition=compute


python climate_maps.py merge-setups=[1] &
python climate_maps.py merge-setups=[2] &
python climate_maps.py merge-setups=[3] &

wait