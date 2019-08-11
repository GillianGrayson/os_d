#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=10000

export OMP_NUM_THREADS=1

scratch=/scratch/denysov/yusipov/os_d/$1
code_base=/home/denysov/yusipov/os_d/source/cpp/CQdiss_os_dimer/bin
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $1/config.txt .

cat config.txt

srun $code_base/cdiss_os_jcs 1 

cp -r $scratch/* $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*