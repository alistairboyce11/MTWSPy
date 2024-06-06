#!/bin/bash
#SBATCH --job-name=MTWSPy
#SBATCH --partition=Cascade-flix
#SBATCH --time=48:00:00
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
source ${HOME}/.bash_profile
conda activate MTWSPy
export PYTHONPATH=./MTWSPy:$PYTHONPATH
export XDG_CACHE_HOME=${HOME}/tmp
today=`date +%Y-%m-%d`

###################################################################

cd $SLURM_SUBMIT_DIR

# Check we have the right number of cores requested in params file.

p_cores=`grep "cores: " params_in.yaml | awk '{print $2}'`

if [ ${SLURM_NTASKS_PER_NODE} != ${p_cores} ]; then
    echo "Make sure number of cores is matching" 
    echo "in params_in.yaml and run_MTWSPy.bash" 
    echo "Currently: "${SLURM_NTASKS_PER_NODE}" & "${p_cores}
    exit
else
    echo "Number of cores requested are equal" 
    echo "Batch file: "${SLURM_NTASKS_PER_NODE}" & Params file: "${p_cores}
    echo "Continue..." 
fi

# ${HOME}/anaconda3/envs/MTWSPy/bin/python3.12 MTWSPy_main.py
MTWSPy_main

####################################################################
exit
