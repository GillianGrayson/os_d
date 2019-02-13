#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --partition=gpu

module load mkl2017
export OMP_NUM_THREADS=32

code_base=/common/home/yusipov_i/Work/os_d/source/cpp/QJX/QJX
scratch=$1

cd $scratch

printf "We are in $(pwd) \n\n"
printf "$1 \n\n"

cat config.txt
cat params.txt

srun $code_base/qjx.out 32

