#!/bin/bash
#SBATCH --job-name=psf_nbd_r_nbf_r
#SBATCH --partition=Cascade-flix,Cascade
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
source ${HOME}/.bash_profile
DIR=$(pwd)
conda activate env3.12
CMTDIR=`echo ${HOME}/d_data_and_docs/gcmt`
export XDG_CACHE_HOME=${HOME}/tmp

cd $SLURM_SUBMIT_DIR


nbd=nbd_r
nbf=nbf_r
pac=pac_r
MODEL=model_r
YEAR=year_r
KCHAN=kchan_r


OUTDIR=`echo ${DIR}/CASE_${nbd}_${nbf}`

####################################################################
for run_num in $(seq ${nbd} ${nbf}); do
    # Get nb in 3 digit code:
    nb=`echo $run_num | awk '{ printf("%03d\n", $1) }'`
    
    echo "Start treatment of Earthquake N° ${nb}, $(date)"  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log

    num_sac_present=`echo ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/*sac | wc -w`

    if [ -d ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 ]; then 
        # Simulation has finished - ready to post process.
        echo "Directory present : OUTPUT_FILES_3D_seis"${nb}"_hd0"  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
        echo "Simulation produced ${num_sac_present} files..."  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log

        # Copy the HDX CMT solution into CMTSOLUTION
        cp ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION_run${nb} ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION

        # Do the Specfem postprocessing
        echo -e "starting SpecFEM postprocessing on $SLURM_NTASKS_PER_NODE processors in: "${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
        ${HOME}/anaconda3/envs/env3.12/bin/python ${HOME}/bin/parallel_proc_specfem_seis.py ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 ${SLURM_NTASKS_PER_NODE} ${KCHAN} >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
        echo -e "Finished SpecFEM postprocessing in: "${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log

        # Move all results to Results DIR
        MOVE_TO_DIR=`echo ${DIR}/py_formatted`

        if [ ! -d $MOVE_TO_DIR ]; then
            echo "Making move-to directory." >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
            mkdir -p $MOVE_TO_DIR
        fi

        # Move results into place.
        cp -r ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/py_formatted/${YEAR}* ${MOVE_TO_DIR}/

    else
        # Something wrong since not enough sac files present.
        echo  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
        echo "No SPECFEM OUTPUT Earthquake N° ${nb}, $(date)" >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log

    fi

    echo  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
    echo "End treatment of Earthquake N° ${nb}, $(date)" >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log
    echo  >> pyproc_sf_${nbd}_${nbf}_${SLURM_JOB_ID}.log


done
####################################################################
exit
