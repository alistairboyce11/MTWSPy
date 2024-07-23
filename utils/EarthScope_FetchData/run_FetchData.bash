#!/bin/bash
#SBATCH --job-name=FetchD
#SBATCH --partition=Cascade-flix
#SBATCH --time=168:00:00
###########################################################
# USER PARAMETERS
# Cascade # LF
#SBATCH --nodes=1
# nombres de  processus MPI par noeuds Cascade
#SBATCH --ntasks-per-node=96

# nombres de processus MPI par processeurs (max=12) sachant qu'il y a 2 processeurs par noeuds
###########################################################
#precise le nombre de coeurs par processus MPI (pour etre sur que 1)
#SBATCH --cpus-per-task=1
#
# precise la memoire par coeur
#SBATCH --mem-per-cpu=4000mb
#
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt
# 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alistair.boyce@ens-lyon.fr
##########################################
source /usr/share/lmod/lmod/init/bash
module use /applis/PSMN/debian11/Cascade/modules/all
# Intel 2021.4.0
module purge
module load impi/2021.9.0-intel-compilers-2023.1.0
#

HOME=`echo ~`
CODE_HOME=`pwd`
# source ${HOME}/.bash_profile
conda activate MTWSPy  
export XDG_CACHE_HOME=${HOME}/tmp
today=`date +%Y-%m-%d`

year="2016"
months="1" #" 2 3 4 5 6" # "7 8 9 10 11 12" #
s_chan="LH"

###################################################################
# Check path for download: Specified in dl script.

###################################################################

cd $SLURM_SUBMIT_DIR

# Check we have the right number of cores requested in params file.

for month in $months; do
    python ./fetch_data.py ${year} ${month} ${s_chan}
done

####################################################################
exit
