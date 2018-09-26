# !/bin/bash -i
#$ -S /bin/bash

# --- Mandatory qsub arguments
# Hardware requirements.
#$ -l h_rss=4000M,h_fsize=4000M,h_cpu=80:00:00,hw=x86_64

# --- Optional qsub arguments
# Change working directory - your job will be run from the directory
# that you call qsub in. So stdout and stderr will end up there .
# $ -cwd

scratch=/scratch/yusipov/matlab_plot
code_base=/home/yusipov/Work/os_d/plot/mpipks/qjx/jcs
mkdir -p $scratch
cd $scratch
cp $code_base/* .

echo " Running on $(hostname)"
echo " We are in $(pwd) "

matlab -nodesktop -nosplash -r tbj_fit_color

rm -r $scratch/*


