#!/bin/bash

if [ "$#" -lt 3 ];
then
    echo " Please specify MODEL, YEAR and CHANNEL of synthetics to process"
    echo " USEAGE:   pyproc_specfem_synthetics.sh <MODEL> <YEAR> <CHANNEL>"
    echo " EXAMPLE:  pyproc_specfem_synthetics.sh IASP91_CRUST1 2008 LH"
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
SCRATCH_DIR=`echo $DIR | sed 's|home|scratch\/Cascade|g'`
CMTDIR=`echo $HOME/d_data_obs/gcmt`
STADIR=`echo $HOME/d_data_obs/dmt`
export XDG_CACHE_HOME=${HOME}/tmp

today=`date +%Y-%m-%d`

if [ ! -d ${DIR}"/LOG_FILES" ]; then
    # Make log file directory
    mkdir -p ${DIR}"/LOG_FILES"
fi

LOG_FILE=`echo ${DIR}"/LOG_FILES/ProcessPy_SPECFEM_MOD_"${MODEL}"_YR_"${YEAR}"_"${today}".log"`
echo Logfile name : $LOG_FILE
echo
echo Logfile name : $LOG_FILE  >> $LOG_FILE
echo  >> $LOG_FILE

###################################################################

cd $DIR
OUTDIR=`echo ${DIR}/${MODEL}/e${YEAR}`

num_specfem_sims=`echo ${OUTDIR}/CASE_*/OUTPUT_FILES_3D_seis???_hd0 | wc -w`

###################################################################

echo >> $LOG_FILE
# Number of Earthquakes
n1=1
n2=$num_specfem_sims
nbs=$(($((${n2}-${n1}))+1))
echo -e "Number of Earthquakes: ${nbs}"  >> $LOG_FILE

# Packet size
pac=20
echo -e "Size of Packets: ${pac}"  >> $LOG_FILE
#Number of Runs
nbc=$(($((${nbs}/${pac}))+1))
echo -e "Number of Runs: ${nbc}"  >> $LOG_FILE

echo >> $LOG_FILE

###################################################################
# Loop through Earthquakes and make packets of submissions.
for i in $(seq ${n1} ${pac} ${n2}); do
    nbd=$i
    nbf=$(($((${nbd}+${pac}))-1))
    if [ ${nbf} -gt ${n2} ]; then
        nbf=${n2}
    fi
    run_num=$(($(($((${nbd}+${pac}))-1))/${pac}))

    echo "Making batch file for SPECFEM submit..." >> $LOG_FILE
    # Replace options in batch file template before submission.
    sed -e "s/nbd_r/${nbd}/g" run_pyproc_specfem_TEMPLATE.bash > ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash
    sed -i -e "s/nbf_r/${nbf}/g" ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash
    sed -i -e "s/pac_r/${pac}/g" ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash
    sed -i -e "s/model_r/${MODEL}/g" ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash
    sed -i -e "s/year_r/${YEAR}/g" ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash
    sed -i -e "s/kchan_r/${CHANNEL}/g" ${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash

    echo "Made batch file for SPECFEM submit..." >> $LOG_FILE
    echo >> $LOG_FILE

    ###################################################################
    echo "Submitting batch file: "${OUTDIR}/run_pyproc_specfem_${nbd}_${nbf}.bash >> $LOG_FILE

    cd ${OUTDIR}/
    sbatch run_pyproc_specfem_${nbd}_${nbf}.bash >> $LOG_FILE
    cd ${DIR}
    sleep 3

    echo "Submitted batch file" >> $LOG_FILE
    echo "Should find the submission in the following squeue call....:" >> $LOG_FILE

    echo >> $LOG_FILE
    squeue | grep "psf_" >> $LOG_FILE
    echo >> $LOG_FILE
    
    echo "Finished the directory: "${OUTDIR} >> $LOG_FILE
    echo -e "For runs: nbd=${nbd}, nbf=${nbf}, run =${run_num}"  >> $LOG_FILE
    echo >> $LOG_FILE
    ###################################################################
done

exit
