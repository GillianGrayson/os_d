# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=1500M,h_fsize=1000M,h_cpu=12:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

module load intel/2015.2

scratch=/scratch/yusipov/os_anderson/$1
code_base=$HOME/Work/os_anderson/source/cpp/CQdiss/bin

mkdir -p $scratch
mkdir -p $1

cd $scratch
cp $1/config.txt .

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt

$code_base/cdiss 16

cp -r $scratch/* $1
rm -r $scratch/*

