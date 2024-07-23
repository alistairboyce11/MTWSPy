#!/bin/bash
#SBATCH --job-name=sf_nbd_nbf
#SBATCH --partition=Cascade
#SBATCH --time=168:00:00
###########################################################
# USER PARAMETERS
# Cascade # LF
#SBATCH --nodes=4
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
# source ${HOME}/.bash_profile
CMTDIR=`echo ${HOME}/d_data_obs/gcmt`
export XDG_CACHE_HOME=${HOME}/tmp

cd $SLURM_SUBMIT_DIR

####################################################################
for run_num in $(seq nbd nbf); do
    # Get nb in 3 digit code:
    nb=`echo $run_num | awk '{ printf("%03d\n", $1) }'`
    echo "Start treatment of Earthquake N° ${nb}, $(date)"  >> date_${SLURM_JOB_ID}.txt

    # Only run the simulation if nonexistent.
    # Calculate number of sac files expected (three times the number of stations)
    num_sac=`cat DATA/STATIONS | wc -l | awk '{print $1*3}'`
    num_sac_present=`echo OUTPUT_FILES_3D_seis${nb}_hd0/*sac | wc -w`

    if [ -d OUTPUT_FILES_3D_seis${nb}_hd0 -a $num_sac -eq $num_sac_present ]; then 
        # It would appear that the simulation has already finished.
        # So we skip it here
        echo "Directory present : OUTPUT_FILES_3D_seis"${nb}"_hd0"  >> date_${SLURM_JOB_ID}.txt
        echo "Simulation appears complete..."  >> date_${SLURM_JOB_ID}.txt
        echo "Skipping..." >> date_${SLURM_JOB_ID}.txt
        continue
    else
        # If directory already exists there is a fault - so delete.
        if [ -d OUTPUT_FILES_3D_seis${nb}_hd0 ]; then
            # Probably a node crashed...
            echo "Faulty directory present : OUTPUT_FILES_3D_seis"${nb}"_hd0"  >> date_${SLURM_JOB_ID}.txt
            echo "Deleting..."                                          >> date_${SLURM_JOB_ID}.txt
            rm -rf OUTPUT_FILES_3D_seis${nb}_hd0
        fi

        # Copy Half-duration =0 file into CMTSOLUTION for SPECFEM, keep copies of full files.
        cp $CMTDIR/year_CMT/HD0/CMTSOLUTION_run${nb}_hd0 DATA/CMTSOLUTION
        cp $CMTDIR/year_CMT/HD0/CMTSOLUTION_run${nb}_hd0 DATA/CMTSOLUTION_run${nb}_hd0
        cp $CMTDIR/year_CMT/CMTSOLUTION_run${nb} DATA/CMTSOLUTION_run${nb} # ERROR 1
        ####################################################################

        BASEMPIDIR=`grep LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
        BASEMPIDIR_TMP=`grep LOCAL_TMP_PATH DATA/Par_file | cut -d = -f 2 `

        echo $BASEMPIDIR

        if [ ! -d $BASEMPIDIR ]; then
            echo "No Mesh DIR found at : "${BASEMPIDIR}  >> date_${SLURM_JOB_ID}.txt
            echo "Run will crash... so exiting...      " >> date_${SLURM_JOB_ID}.txt
            exit
        fi
        
        if [ ! -d $BASEMPIDIR_TMP ]; then
            echo "No TMP Mesh DIR found at : "${BASEMPIDIR_TMP}  >> date_${SLURM_JOB_ID}.txt
            echo "Making TMP directory...      " >> date_${SLURM_JOB_ID}.txt
            mkdir -p ${BASEMPIDIR_TMP}
        fi

        ####################################################################

        # script to run the mesher and the solver
        # read DATA/Par_file to get information about the run
        # compute total number of nodes needed
        NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
        NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
        NCHUNKS=`grep NCHUNKS DATA/Par_file | awk 'NR==1' | cut -d = -f 2`

        # total number of nodes is the product of the values read
        numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

        rm -r -f OUTPUT_FILES
        mkdir -p OUTPUT_FILES

        ####################################################################
        d=`date`
        echo -e "Starting compilation $(date)" >> date_${SLURM_JOB_ID}.txt
        ./configure FC=ifort CC=icc MPIFC=mpiifort
        make realclean
        make clean
        make xmeshfem3D
        make xcreate_header_file
        ./bin/xcreate_header_file
        make xspecfem3D
        d=`date`
        echo -e "Finished compilation $(date)" >> date_${SLURM_JOB_ID}.txt
        ####################################################################
        # backup files used for this simulation
        cp DATA/Par_file OUTPUT_FILES/
        cp DATA/STATIONS OUTPUT_FILES/
        cp DATA/CMTSOLUTION OUTPUT_FILES/
        cp DATA/CMTSOLUTION_run${nb}_hd0 OUTPUT_FILES/
        cp DATA/CMTSOLUTION_run${nb} OUTPUT_FILES/
        ##
        ## mesh generation
        ##
        echo ${SLURM_JOB_NODELIST} > OUTPUT_FILES/compute_nodes
        echo "${SLURM_JOB_ID}" > OUTPUT_FILES/jobid

        sleep 2

        # Attempt to read Mesh directly from scratch instead to save space....

        # echo -e "starting MPI mesher on $numnodes processors $(date)" >> date_${SLURM_JOB_ID}.txt
        # mpirun -np ${SLURM_NTASKS} -bootstrap ssh $PWD/bin/xmeshfem3D >> date_${SLURM_JOB_ID}.txt
        # echo -e "mesher done $(date)" >> date_${SLURM_JOB_ID}.txt

        ###################################################################
        ##
        ## forward simulation
        ##

        sync
        sleep 2

        echo -e "starting MPI solver on $numnodes processors in current directory $PWD $(date)" >> date_${SLURM_JOB_ID}.txt
        mpirun -np ${SLURM_NTASKS} -bootstrap ssh $PWD/bin/xspecfem3D >> date_${SLURM_JOB_ID}.txt
        echo -e "MPI solver finished successfully $(date)" >> date_${SLURM_JOB_ID}.txt
        #

        ####################################################################

        mv OUTPUT_FILES OUTPUT_FILES_3D_seis${nb}_hd0

        # Copy the HDX CMT solution into CMTSOLUTION
        cp OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION_run${nb} OUTPUT_FILES_3D_seis${nb}_hd0/CMTSOLUTION

    fi



    # HANDLE POSTPROCESSING IN ANOTHER SLURM SUBMISSION///

    echo  >> date_${SLURM_JOB_ID}.txt
    echo "End treatment of Earthquake N° ${nb}, $(date)" >> date_${SLURM_JOB_ID}.txt
    echo  >> date_${SLURM_JOB_ID}.txt
done
####################################################################
exit
