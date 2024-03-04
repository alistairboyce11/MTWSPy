#!/bin/bash
#SBATCH --job-name=psf_year_r
#SBATCH --partition=Cascade,Cascade-flix
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
source ${HOME}/.bash_profile
CMTDIR=`echo ${HOME}/d_data_and_docs/gcmt`
export XDG_CACHE_HOME=${HOME}/tmp

cd $SLURM_SUBMIT_DIR

N1=n1_r
N2=n2_r
PAC=pac_r
MODEL=model_r
YEAR=year_r
KCHAN=kchan_r

echo "Starting..." > date_${SLURM_JOB_ID}.txt

# Loop through Earthquakes and make packets of submissions.
for i in $(seq ${N1} ${PAC} ${N2}); do
    nbd=$i
    nbf=$(($((${nbd}+${PAC}))-1))
    if [ ${nbf} -gt ${N2} ]; then
        nbf=${N2}
    fi
    run_num=$(($(($((${nbd}+${PAC}))-1))/${PAC}))
    OUTDIR=`echo CASE_${nbd}_${nbf}`

    echo >> date_${SLURM_JOB_ID}.txt
    echo "Checking for the directory: "${OUTDIR} >> date_${SLURM_JOB_ID}.txt
    echo -e "For runs: nbd=${nbd}, nbf=${nbf}, run=${run_num}"  >> date_${SLURM_JOB_ID}.txt
    echo >> date_${SLURM_JOB_ID}.txt

    # Check MODEL/eYEAR directory:
    if [ ! -d ${OUTDIR} ]; then
        echo "Case directory: ./"${OUTDIR}" does not exist..."  >> date_${SLURM_JOB_ID}.txt
        echo "Exiting...." >> date_${SLURM_JOB_ID}.txt
        exit
    else
        echo "Case directory exists : ./"${OUTDIR}  >> date_${SLURM_JOB_ID}.txt
    fi

    ####################################################################
    for run_num in $(seq ${nbd} ${nbf}); do
        # Get nb in 3 digit code:
        nb=`echo $run_num | awk '{ printf("%03d\n", $1) }'`

        echo "Start treatment of Earthquake N° ${nb}, $(date)"  >> date_${SLURM_JOB_ID}.txt

        # Copy the HDX CMT solution into CMTSOLUTION
        cp ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION_run${nb} ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION

        # Do the Specfem postprocessing
        echo -e "starting SpecFEM postprocessing on $SLURM_NTASKS_PER_NODE processors in: "${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 >> date_${SLURM_JOB_ID}.txt
        ${HOME}/bin/parallel_proc_specfem_seis.sh ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 ${SLURM_NTASKS_PER_NODE} ${KCHAN} >> date_${SLURM_JOB_ID}.txt
        echo -e "Finished SpecFEM postprocessing in: "${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0 >> date_${SLURM_JOB_ID}.txt

        # Move all results to Results DIR
        MOVE_TO_DIR=`echo ./formatted`

        if [ ! -d $MOVE_TO_DIR ]; then
            echo "Making move-to directory." >> date_${SLURM_JOB_ID}.txt
            mkdir -p $MOVE_TO_DIR
        fi

        # Move results into place.
        cp -r ${OUTDIR}/OUTPUT_FILES_3D_seis${nb}_hd0/formatted/${YEAR}* ${MOVE_TO_DIR}/

        echo  >> date_${SLURM_JOB_ID}.txt
        echo "End treatment of Earthquake N° ${nb}, $(date)" >> date_${SLURM_JOB_ID}.txt
        echo  >> date_${SLURM_JOB_ID}.txt
    done
done
####################################################################
exit
