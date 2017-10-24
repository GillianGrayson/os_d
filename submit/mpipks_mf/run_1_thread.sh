# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=1500M,h_fsize=1000M,h_cpu=10:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

# --- Job Execution
# For faster disk access copy files to / scratch first .
module load mkl

data_path=/data/biophys/yusipov/os_d/mf
scratch=/scratch/yusipov/os_d/$1
code_base=$HOME/Work/os_d/source/cpp/MF/MF
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $1/config.txt .

# Execution - running the actual program .
# [ Remember : Don � t read or write to / home from here .]

echo " Running on $(hostname)"
echo " We are in $(pwd) "

cat config.txt
$code_base/MD.out

# Finish - Copy files back to your home directory , clean up .
cp -r $scratch/* $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*
#cd
#rm - rf $scratch
