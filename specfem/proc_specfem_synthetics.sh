#!/bin/bash

if [ "$#" -lt 3 ];
then
    echo " Please specify MODEL, YEAR and CHANNEL of synthetics to process"
    echo " USEAGE:   proc_specfem_synthetics.sh <MODEL> <YEAR> <CHANNEL>"
    echo " EXAMPLE:  proc_specfem_synthetics.sh IASP91_CRUST1 2008 LH"
    echo " "
    exit
else
    ######## ++++ SET YEAR AND CHANNEL +++++++++++++++#########
    MODEL=$1
    YEAR=$2
    CHANNEL=$3
    echo " "
    echo "Processing SPECFEM SYNTHETICS for model: "$MODEL", for: "$YEAR", assuming data channel: "$CHANNEL
    echo " "
fi

HOME=`echo ~`
source ${HOME}/.bash_profile
DIR=$(pwd)
CMTDIR=`echo $HOME/d_data_and_docs/gcmt`
STADIR=`echo $HOME/d_data_and_docs/dmt`
export XDG_CACHE_HOME=${HOME}/tmp

today=`date +%Y-%m-%d`

if [ ! -d ${DIR}"/LOG_FILES" ]; then
    # Make log file directory
    mkdir -p ${DIR}"/LOG_FILES"
fi

LOG_FILE=`echo ${DIR}"/LOG_FILES/Process_SPECFEM_MOD_"${MODEL}"_YR_"${YEAR}"_CH_"${CHANNEL}"_"${today}".log"`
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

# ###################################################################

num_CMT_files=`echo ${CMTDIR}/${YEAR}_CMT/HD0/CMTSOLUTION_run* | wc -w`
echo "We have "${num_CMT_files}" CMT files to process..."  >> $LOG_FILE
echo  >> $LOG_FILE

# ###################################################################

cd $DIR

set_mod=`grep "MODEL" ./DATA/Par_file | awk 'NR==1' | awk '{print $3}'`
set_length=`grep "RECORD_LENGTH_IN_MINUTES" ./DATA/Par_file  | awk 'NR==1' | awk '{print $3}'`
set_stations=`wc -l ./DATA/STATIONS | awk '{print $1}'`

echo  >> $LOG_FILE
echo "Producing results for model: "${set_mod}", record length: "${set_length}" minutes, "${set_stations}" stations" >> $LOG_FILE
echo  >> $LOG_FILE

###################################################################

OUTDIR=`echo ${MODEL}/e${YEAR}`

# Check model/year directory:
if [ ! -d ${OUTDIR} ]; then
    echo "Year directory: ./"${OUTDIR}" does not exist..."  >> $LOG_FILE
    echo "Exiting...."
    exit
else
    echo "Year directory exists : ./"${OUTDIR}  >> $LOG_FILE
fi

echo >> $LOG_FILE
# Number of Earthquakes
n1=1
n2=$num_CMT_files
nbs=$(($((${n2}-${n1}))+1))
#nbs=66
echo -e "Number of Earthquakes: ${nbs}"  >> $LOG_FILE
# Packet size
pac=20
echo -e "Size of Packets: ${pac}"  >> $LOG_FILE
#Number of Runs
nbc=$(($((${nbs}/${pac}))+1))
echo -e "Number of Runs: ${nbc}"  >> $LOG_FILE
#
echo >> $LOG_FILE


echo "Making batch file for Processing SPECFEM submit..." >> $LOG_FILE

# Replace options in batch file template before submission.
sed -e "s/n1_r/${n1}/g" run_proc_specfem_slurm_TEMPLATE.bash  > ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash
sed -i -e "s/n2_r/${n2}/g" ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash
sed -i -e "s/pac_r/${pac}/g" ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash
sed -i -e "s/model_r/${MODEL}/g" ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash
sed -i -e "s/year_r/${YEAR}/g" ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash
sed -i -e "s/kchan_r/${CHANNEL}/g" ${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash

echo "Made batch file for Processing SPECFEM submit..." >> $LOG_FILE
echo >> $LOG_FILE

###################################################################
echo "Submitting batch file: "${OUTDIR}/run_proc_specfem_slurm_${n1}_${n2}.bash >> $LOG_FILE

cd ${OUTDIR}/
sbatch run_proc_specfem_slurm_${n1}_${n2}.bash >> $LOG_FILE
cd ${DIR}
sleep 10

echo "Submitted batch file" >> $LOG_FILE
echo "Should find the submission in the following squeue call....:" >> $LOG_FILE

echo >> $LOG_FILE
squeue | grep "psf_" >> $LOG_FILE
echo >> $LOG_FILE

echo "Finished the directory: "${OUTDIR} >> $LOG_FILE
echo -e "For runs: nbd=${n1}, nbf=${n2}, runs=${nbc}"  >> $LOG_FILE
echo >> $LOG_FILE
###################################################################


exit
