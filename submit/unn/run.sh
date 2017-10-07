#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --partition=gpu

export KMP_AFFINITY="granularity=fine,compact,1,0"
export MKL_NUM_THREADS=16

module load mkl2017

code_base=/home/yusipov_i/Work/os_d/source/cpp/CQdiss/bin
scratch=$1

cd $scratch

printf "We are in $(pwd) \n\n"
printf "$1 \n\n"

cat config.txt

srun $code_base/cdiss_fb 16

