#!/bin/bash

if [ "$#" -lt 1 ];
then
    echo " Please specify YEAR of obspyDMT Data to fetch"
    echo " USEAGE:   get_obspy_data.sh <YEAR>"
    echo " EXAMPLE:  get_obspy_data.sh 2008"
    echo " "
    exit
else
    ######## ++++ SET YEAR AND CHANNEL +++++++++++++++#########
    # DATA_SOURCE=$1
    YEAR=$1
    # CHANNEL=$3
    echo " "
    echo "Requesting obspyDMT data for: "$YEAR # " from: "$DATA_SOURCE
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

LOG_FILE=`echo ${DIR}"/LOG_FILES/Request_obspyDMT_YR_"${YEAR}"_"${today}".log"`
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

# Make year directory:
if [ ! -d "./e"${YEAR} ]; then
    mkdir -p "./e"${YEAR}
    echo "Year directory created : ./e"${YEAR}  >> $LOG_FILE
else
    echo "Year directory exists : ./e"${YEAR}  >> $LOG_FILE
fi

###################################################################

echo "Making batch file for obspyDMT submit..." >> $LOG_FILE

# Replace options in batch file template before submission.
sed -e "s/YXXR/${YEAR}/g" run_obspy_slurm_TEMPLATE.bash  > run_obspy_slurm_${YEAR}.bash

echo "Made batch file for obspyDMT submit..." >> $LOG_FILE
echo >> $LOG_FILE

echo "Submitting batch file: run_obspy_slurm_${YEAR}.bash" >> $LOG_FILE

echo sbatch run_obspy_slurm_${YEAR}.bash    >> $LOG_FILE
sbatch run_obspy_slurm_${YEAR}.bash    >> $LOG_FILE
sleep 3

echo "Submitted batch file" >> $LOG_FILE
echo "Should find the submission in the following squeue call....:" >> $LOG_FILE

echo >> $LOG_FILE
squeue | grep "obs_" >> $LOG_FILE
echo >> $LOG_FILE

echo "Finished submission for: "$YEAR >> $LOG_FILE
echo >> $LOG_FILE

exit
