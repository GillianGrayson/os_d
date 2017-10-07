#!/bin/sh

#SBATCH --time=200

#SBATCH --partition=gpu

export KMP_AFFINITY="granularity=fine,compact,1,0"
export MKL_NUM_THREADS=16

srun $1/../bin/cdiss_fb 16

