#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

export OMP_NUM_THREADS=1

scratch=/scratch/ivanchen/yusipov/os_d/$1
code_base=/home/ivanchen/yusipov/os_d/source/cpp/QJX/QJX
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $1/config.txt .
cp $1/params.txt .

cat config.txt
cat params.txt

srun $code_base/qjx.out 1

cp -r $scratch/* $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*

