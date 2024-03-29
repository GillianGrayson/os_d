#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --partition=gpu

export KMP_AFFINITY="granularity=fine,compact,1,0"
export MKL_NUM_THREADS=16

code_base=/common/home/yusipov_i/Work/os_d/source/cpp/CQdiss_os_jcs_opt/bin
scratch=$1

cd $scratch

printf "We are in $(pwd) \n\n"
printf "$1 \n\n"

cat config.txt

srun $code_base/cdiss_os_jcs 16 

