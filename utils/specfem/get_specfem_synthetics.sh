#!/bin/bash

if [ "$#" -lt 2 ];
then
    echo " Please specify MODEL and YEAR of synthetics to produce"
    echo " USEAGE:   get_specfem_synthetics.sh <MODEL> <YEAR>"
    echo " EXAMPLE:  get_specfem_synthetics.sh IASP91_CRUST1 2008"
    echo " "
    exit
else
    ######## ++++ SET MODEL and YEAR +++++++++++++++#########
    MODEL=$1
    YEAR=$2
    echo " "
    echo "Requesting SPECFEM SYNTHETICS for model: "$MODEL", for: "$YEAR
    echo " "
fi

HOME=`echo ~`
source ${HOME}/.bash_profile
DIR=$(pwd)
SCRATCH_DIR=`echo $DIR | sed 's|home|scratch\/Cascade|g'`
CMTDIR=`echo $HOME/d_data_and_docs/gcmt`
STADIR=`echo $HOME/d_data_and_docs/dmt`
export XDG_CACHE_HOME=${HOME}/tmp

today=`date +%Y-%m-%d`

if [ ! -d ${DIR}"/LOG_FILES" ]; then
    # Make log file directory
    mkdir -p ${DIR}"/LOG_FILES"
fi

LOG_FILE=`echo ${DIR}"/LOG_FILES/Request_SPECFEM_MOD_"${MODEL}"_YR_"${YEAR}"_"${today}".log"`
echo Logfile name : $LOG_FILE
echo
echo Logfile name : $LOG_FILE  >> $LOG_FILE
echo  >> $LOG_FILE

###################################################################
# Check all CMT files are present.

if [ ! -d $CMTDIR ]; then
    echo "CMT directory not found at: "$CMTDIR  >> $LOG_FILE
    echo "Exiting..."  >> $LOG_FILE
    exit
else
    echo "Moving to: "$CMTDIR" for CMT file check" >> $LOG_FILE
fi

cd $CMTDIR

if [ ! -f ${YEAR}_CMT/${YEAR}_data.txt ]; then
    echo "Need to download CMT files..."  >> $LOG_FILE
    if [ ! -f ./download_cmtsolutions.py ]; then
        echo "No Python script to make CMT files..."  >> $LOG_FILE
        echo "Exiting..."  >> $LOG_FILE
        exit
    else
        python3 download_cmtsolutions.py $YEAR  >> $LOG_FILE
        echo "CMT files ready..."  >> $LOG_FILE

    fi
else
    echo "CMT files already present..."  >> $LOG_FILE
    echo "e.g. ..."  >> $LOG_FILE
    echo  >> $LOG_FILE
    tail -20 ${YEAR}_CMT/${YEAR}_data.txt  >> $LOG_FILE
    echo  >> $LOG_FILE
    num_CMT_files=`echo ${YEAR}_CMT/HD0/CMTSOLUTION_run* | wc -w`
    echo "We have "${num_CMT_files}" CMT files to process..."  >> $LOG_FILE
    echo  >> $LOG_FILE
fi
###################################################################

cd $DIR
###################################################################
# Get the stations file:
if [ ! -d $STADIR ]; then
    echo "STATION directory not found at: "$STADIR  >> $LOG_FILE
    echo "Exiting..."  >> $LOG_FILE
    exit
else
    echo "Moving to: "$STADIR" for Station file check" >> $LOG_FILE
fi

cd $STADIR

if [ ! -f e${YEAR}/STATIONS_${YEAR}_obspy ]; then
    echo "Need to make STATION file..."  >> $LOG_FILE
    if [ ! -f ./mk_SPECFEM_STATIONS_obspy.py ]; then
        echo "No Python script to make STATION files..."  >> $LOG_FILE
        echo "Exiting..."  >> $LOG_FILE
        exit
    else
        python3 mk_SPECFEM_STATIONS_obspy.py LH?,BH?,HH? $YEAR  >> $LOG_FILE
        echo "STATION file ready..."  >> $LOG_FILE

    fi
else
    echo "STATION file already present..."  >> $LOG_FILE
    echo "e.g. ..."  >> $LOG_FILE
    echo  >> $LOG_FILE
    tail -20 e${YEAR}/STATIONS_${YEAR}_obspy  >> $LOG_FILE
    echo  >> $LOG_FILE
fi

# Copy stations file to SLURM DIRECTORY
echo  "Copying STATIONS file..." >> $LOG_FILE
echo "cp e${YEAR}/STATIONS_${YEAR}_obspy ${DIR}/DATA/STATIONS"  >> $LOG_FILE
cp e${YEAR}/STATIONS_${YEAR}_obspy ${DIR}/DATA/STATIONS
echo  >> $LOG_FILE

###################################################################

cd $DIR

set_mod=`grep "MODEL" ./DATA/Par_file | awk 'NR==1' | awk '{print $3}'`
set_length=`grep "RECORD_LENGTH_IN_MINUTES" ./DATA/Par_file  | awk 'NR==1' | awk '{print $3}'`
set_stations=`wc -l ./DATA/STATIONS | awk '{print $1}'`

echo  >> $LOG_FILE
echo "Producing results for model: "${set_mod}", record length: "${set_length}" minutes, "${set_stations}" stations" >> $LOG_FILE
echo  >> $LOG_FILE

###################################################################

# Make model/year directory:
if [ ! -d ${MODEL}"/e"${YEAR} ]; then
    mkdir -p ${MODEL}"/e"${YEAR}
    echo "Year directory created : ./"${MODEL}"/e"${YEAR}  >> $LOG_FILE
else
    echo "Year directory exists : ./"${MODEL}"/e"${YEAR}  >> $LOG_FILE
fi

# # Make directory on Scratch for temporary files
# if [ ! -d ${SCRATCH_DIR}"/"${MODEL}"/e"${YEAR} ]; then
#     mkdir -p ${SCRATCH_DIR}"/"${MODEL}"/e"${YEAR}
#     echo "Year directory created on scratch : "${SCRATCH_DIR}"/"${MODEL}"/e"${YEAR}  >> $LOG_FILE
# else
#     echo "Year directory exists on scratch : "${SCRATCH_DIR}"/"${MODEL}"/e"${YEAR}  >> $LOG_FILE
# fi


echo >> $LOG_FILE
# Number of Earthquakes
n1=201
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

###################################################################
# Loop through Earthquakes and make packets of submissions.
for i in $(seq ${n1} ${pac} ${n2}); do
    nbd=$i
    nbf=$(($((${nbd}+${pac}))-1))
    if [ ${nbf} -gt ${n2} ]; then
        nbf=${n2}
    fi
    run_num=$(($(($((${nbd}+${pac}))-1))/${pac}))
    OUTDIR=`echo ${MODEL}/e${YEAR}/CASE_${nbd}_${nbf}`

    echo >> $LOG_FILE
    echo "Creating the directory: "${OUTDIR} >> $LOG_FILE
    echo -e "For runs: nbd=${nbd}, nbf=${nbf}, run =${run_num}"  >> $LOG_FILE
    echo >> $LOG_FILE

    # Do not overwrite
    # rm -rf ${OUTDIR}
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi 

    ###################################################################
    echo "Starting rsync for SPECFEM directory files..." >> $LOG_FILE
   
    # All files checked for SpecFEM V8.0.0
    rsync -az ${DIR}/bin  ${OUTDIR}/
    rsync -az ${DIR}/DATA ${OUTDIR}/
    rsync -az ${DIR}/change_simulation_type.pl ${OUTDIR}/
    # rsync -az ${DIR}/compil.sh ${OUTDIR}/
    rsync -az ${DIR}/config.guess ${OUTDIR}/
    rsync -az ${DIR}/config.log ${OUTDIR}/
    rsync -az ${DIR}/config.status ${OUTDIR}/
    rsync -az ${DIR}/config.sub ${OUTDIR}/
    rsync -az ${DIR}/configure ${OUTDIR}/
    rsync -az ${DIR}/configure.ac ${OUTDIR}/
    # rsync -az ${DIR}/download_the_whole_topography_database_if_you_want_other_topographic_models.bash ${OUTDIR}/
    rsync -az ${DIR}/doc ${OUTDIR}/
    rsync -az ${DIR}/flags.guess ${OUTDIR}/
    rsync -az ${DIR}/install-sh ${OUTDIR}/
    rsync -az ${DIR}/LICENSE ${OUTDIR}/
    rsync -az ${DIR}/m4 ${OUTDIR}/
    rsync -az ${DIR}/Makefile ${OUTDIR}/
    rsync -az ${DIR}/Makefile.in ${OUTDIR}/
    rsync -az ${DIR}/obj ${OUTDIR}/
    # rsync -az ${DIR}/original ${OUTDIR}/
    rsync -az ${DIR}/SEM ${OUTDIR}/
    rsync -az ${DIR}/setup ${OUTDIR}/
    rsync -az ${DIR}/src ${OUTDIR}/
    rsync -az ${DIR}/utils ${OUTDIR}/

    echo "Finished rsync for SPECFEM directory files" >> $LOG_FILE
    echo >> $LOG_FILE

    echo "Making batch file for SPECFEM submit..." >> $LOG_FILE
    # Replace options in batch file template before submission.
    sed -e "s/nbd/${nbd}/g" run_specfem_slurm_TEMPLATE.bash > ${OUTDIR}/run_specfem_slurm_${nbd}_${nbf}.bash
    sed -i -e "s/nbf/${nbf}/g" ${OUTDIR}/run_specfem_slurm_${nbd}_${nbf}.bash
    sed -i -e "s/year/${YEAR}/g" ${OUTDIR}/run_specfem_slurm_${nbd}_${nbf}.bash

    echo "Made batch file for SPECFEM submit..." >> $LOG_FILE
    echo >> $LOG_FILE

    ###################################################################
    echo "Submitting batch file: "${OUTDIR}/run_specfem_slurm_${nbd}_${nbf}.bash >> $LOG_FILE

    cd ${OUTDIR}/
    sbatch run_specfem_slurm_${nbd}_${nbf}.bash >> $LOG_FILE
    cd ${DIR}
    sleep 15

    echo "Submitted batch file" >> $LOG_FILE
    echo "Should find the submission in the following squeue call....:" >> $LOG_FILE

    echo >> $LOG_FILE
    squeue | grep "sf_" >> $LOG_FILE
    echo >> $LOG_FILE
    
    echo "Finished the directory: "${OUTDIR} >> $LOG_FILE
    echo -e "For runs: nbd=${nbd}, nbf=${nbf}, run =${run_num}"  >> $LOG_FILE
    echo >> $LOG_FILE
    ###################################################################
done

exit
