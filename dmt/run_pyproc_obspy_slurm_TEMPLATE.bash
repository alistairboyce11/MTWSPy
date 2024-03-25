#!/bin/bash
#SBATCH --job-name=pobs_year_r_month_r
#SBATCH --partition=Cascade,Cascade-flix
#SBATCH --time=48:00:00
#SBATCH --nodes=1
### SBATCH --ntasks=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#
## #SBATCH -o output_%j.txt
## #SBATCH -e stderr_%j.txt
# 
## #SBATCH --mail-type=BEGIN,END,FAIL
## #SBATCH --mail-user=alistair.boyce@ens-lyon.fr
#

source /usr/share/lmod/lmod/init/bash
module use /applis/PSMN/debian11/Cascade/modules/all

HOME=`echo ~`
source ${HOME}/.bash_profile
subdir=`echo ${HOME}/d_data_and_docs/dmt`
LOG_LOC=`echo ${subdir}"/LOG_FILES"`
export XDG_CACHE_HOME=${HOME}/tmp

# ############ For anaconda installation #################
export PATH=".:${HOME}:${HOME}/anaconda3/bin:$PATH"
source activate env3.12

cd $SLURM_SUBMIT_DIR

YEAR=year_r
MONTH=month_r
KCHAN=kchan_r

LOGFILE=`echo pyproc_obspy_${MONTH}_${YEAR}_${KCHAN}_${SLURM_JOB_ID}.log`

echo " " > $LOGFILE
echo "Start treatment of year ${YEAR}, month ${MONTH}, $(date)" >> $LOGFILE
echo " " >> $LOGFILE

evt_list=`echo ${SLURM_SUBMIT_DIR}/${YEAR}${MONTH}??_*.a`


num_evts=`echo $evt_list | wc -w`

if [ $num_evts -eq 0 ]; then
    echo " " >> $LOGFILE
    echo "No EQs to process for ${YEAR} ${MONTH}" >> $LOGFILE
    echo "Exiting..." >> $LOGFILE
    echo " " >> $LOGFILE
    exit
else
    echo " " >> $LOGFILE
    echo " We have ${num_evts} to process.... continue" >> $LOGFILE
    echo " " >> $LOGFILE    
fi

# Check for EVENTS-INFO
if [ ! -f ${SLURM_SUBMIT_DIR}/EVENTS-INFO/catalog.txt ]; then
    # No station file to read
    echo "No catalog.txt file to read for ${YEAR}..." >> $LOGFILE
    echo "Skip YEAR..." >> $LOGFILE
    echo "Exit..."  >> $LOGFILE
    exit
fi



# Move all results to Results DIR
MOVE_TO_DIR=`echo ${SLURM_SUBMIT_DIR}/py_formatted`

if [ ! -d $MOVE_TO_DIR ]; then
    echo "Making move-to directory." >> $LOGFILE
    mkdir -p $MOVE_TO_DIR
fi

for evt in $evt_list; do

    if [ -d ${evt} ]; then

        if [ ! -f ${evt}/info/station_event ]; then
            # No station file to read
            echo "No station_event file to read for ${evt}..." >> $LOGFILE
            echo "Skip event..." >> $LOGFILE
            echo "Continue..."  >> $LOGFILE
            continue
        fi

        if [ ! -d ${evt}/py_formatted ]; then
            echo "Making outfile directory." >> $LOGFILE
            mkdir -p ${evt}/py_formatted
        fi


        # Do the ObspyDMT postprocessing
        echo -e "starting SpecFEM postprocessing on $SLURM_NTASKS_PER_NODE processors for: "${evt} >> $LOGFILE
        ${HOME}/anaconda3/envs/env3.12/bin/python ${HOME}/bin/parallel_proc_obspyDMT_seis.py ${evt} ${SLURM_NTASKS_PER_NODE} ${KCHAN} >> $LOGFILE
        echo -e "Finished SpecFEM postprocessing for: "${evt} >> $LOGFILE

        # Move results into place.
        cp -r ${evt}/py_formatted/${YEAR}* ${MOVE_TO_DIR}/

    else
        # Something wrong
        echo  >> $LOGFILE
        echo "No obspy data for ${evt}, $(date)" >> $LOGFILE
    fi


done

echo  >> $LOGFILE
echo "End treatment of year ${YEAR}, month ${MONTH}, $(date)" >> $LOGFILE
echo  >> $LOGFILE


exit



