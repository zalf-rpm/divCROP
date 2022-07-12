#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --partition=compute
#SBATCH --job-name=go_ascii
#SBATCH --time=10:00:00

SUB=$1

IMAGE_DIR_PYTHON=~/singularity/python
SINGULARITY_PYTHON_IMAGE=python3.7_2.0.sif
IMAGE_PYTHON_PATH=${IMAGE_DIR_PYTHON}/${SINGULARITY_PYTHON_IMAGE}
mkdir -p $IMAGE_DIR_PYTHON
if [ ! -e ${IMAGE_PYTHON_PATH} ] ; then
echo "File '${IMAGE_PYTHON_PATH}' not found"
cd $IMAGE_DIR_PYTHON
singularity pull docker://zalfrpm/python3.7:2.0
cd ~
fi

FOLDER=../agg_grids_clim/$SUB
mkdir -p $FOLDER/img
IMG=~/singularity/python/python3.7_2.0.sif
singularity run -B $FOLDER:/source,$FOLDER/img:/out $IMG python create_image_from_ascii.py path=cluster projection=25832

wait