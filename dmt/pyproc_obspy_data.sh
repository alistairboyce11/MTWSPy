#!/bin/bash

if [ "$#" -lt 2 ];
then
    echo " Please specify YEAR and CHANNEL of obspyDMT Data to process"
    echo " USEAGE:   pyproc_obspy_data.sh <YEAR> <CHANNEL>"
    echo " EXAMPLE:  pyproc_obspy_data.sh 2008 LH"
    echo " "
    exit
else
    ######## ++++ SET YEAR AND CHANNEL +++++++++++++++#########
    YEAR=$1
    CHANNEL=$2
    echo " "
    echo "Processing obspyDMT data for: "$YEAR", assuming data channel: "$CHANNEL
    echo " "
fi

HOME=`echo ~`
source ${HOME}/.bash_profile
DIR=$(pwd)
export XDG_CACHE_HOME=${HOME}/tmp
today=`date +%Y-%m-%d`


if [ ! -d ${DIR}"/LOG_FILES" ]; then
    # Make log file directory
    mkdir -p ${DIR}"/LOG_FILES"
fi

LOG_FILE=`echo ${DIR}"/LOG_FILES/ProcessPy_obspyDMT_YR_"${YEAR}"_CH_"${CHANNEL}"_"${today}".log"`
echo Logfile name : $LOG_FILE
echo
echo Logfile name : $LOG_FILE  > $LOG_FILE
echo  >> $LOG_FILE

###################################################################

# First run a check to see if all modules/packages/codes are present
echo  >> $LOG_FILE
echo "Starting modules/packages/code check" >> $LOG_FILE

${HOME}/bin/check_obspyDMT_SPECFEM_install.sh # OUTFILE = ${HOME}/tmp.check-sh
mod_check=`tail -1 ${HOME}/tmp/tmp.check-sh | awk '{print $1}'`

if [ $mod_check == 'YES' ]; then
    echo "Finished modules/packages/code check" >> $LOG_FILE
    echo  >> $LOG_FILE
else
    echo "FAILED modules/packages/code check"  >> $LOG_FILE
    echo "Exiting..."  >> $LOG_FILE
    exit
fi

###################################################################

# Check Data directory
if [ ! -d e$YEAR ]; then
    echo "Data directory: e"$YEAR" not present..." >> $LOG_FILE
    echo "Exiting..." >> $LOG_FILE
    exit
else
    echo "Data directory e"$YEAR" exists" >> $LOG_FILE
fi

###################################################################

cd $DIR
OUTDIR=`echo ${DIR}/e${YEAR}`


echo >> $LOG_FILE
# Months
n1=1
n2=12
nbs=$(($((${n2}-${n1}))+1))
echo -e "Number of Simultaneous runs: ${nbs}"  >> $LOG_FILE

# Loop through Months and make packets of submissions.

for i in $(seq ${n1} ${n2}); do

    MONTH=`echo $i | awk '{ printf("%02d\n", $1) }'`

    echo "Making batch file for Processing obspyDMT submit..." >> $LOG_FILE

    # Replace options in batch file template before submission.
    sed -e "s/year_r/${YEAR}/g" run_pyproc_obspy_slurm_TEMPLATE.bash  > ${OUTDIR}/run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash 
    sed -i -e "s/month_r/${MONTH}/g" ${OUTDIR}/run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash 
    sed -i -e "s/kchan_r/${CHANNEL}/g" ${OUTDIR}/run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash 

    echo "Made batch file for Processing obspyDMT submit..." >> $LOG_FILE
    echo >> $LOG_FILE

    echo "Submitting batch file: ${OUTDIR}/run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash" >> $LOG_FILE

    cd ${OUTDIR}/
    echo sbatch run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash   >> $LOG_FILE
    sbatch run_pyproc_obspy_slurm_${MONTH}_${YEAR}_${CHANNEL}.bash   >> $LOG_FILE
    sleep 3
    cd ${DIR}

    echo "Submitted batch file" >> $LOG_FILE
    echo "Should find the submission in the following squeue call....:" >> $LOG_FILE

    echo >> $LOG_FILE
    squeue | grep "pobs_" >> $LOG_FILE
    echo >> $LOG_FILE

done

echo "Finished submission for: "$YEAR", assuming data channel: "$CHANNEL >> $LOG_FILE
echo >> $LOG_FILE



exit








