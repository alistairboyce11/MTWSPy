#!/bin/bash

if [ "$#" -lt 3 ];
then
    echo " Please specifiy Start date and years  "
    echo " USEAGE:   parallel_proc_specfem_seis.sh <DATA_LOC> <PARALLEL_PROCESSES> <CHANNEL>"
    echo " EXAMPLE:  parallel_proc_specfem_seis.sh ./OUTPUT_FILES_3D_seis_hd0 96 LH"
    echo " "
    exit
else
    DATA_LOC=$1
    if [ ! -d $DATA_LOC ]; then
        echo
        echo "Directory given is not present.... exiting...."
        echo
        exit
    fi
     
    JOBS=$2
    CHANNEL=$3
    echo
    echo "Processing "$DATA_LOC", using "$JOBS" processes, assuming data channel: "$CHANNEL
fi

############# Log file ###########################

HOME=`echo ~`
source ${HOME}/.bash_profile
USER=`whoami`
export XDG_CACHE_HOME=${HOME}/tmp
subdir=`echo ${HOME}/s_seismology/specfem/sge.globe/C_source_codes/i_intel`
# LOG_LOC=`echo ${HOME}"/LOG_FILES"`
today=`date +%Y-%m-%d`

LOG_FILE=`echo $DATA_LOC"/process_seis_NP_"${JOBS}"_CH_"${CHANNEL}"_"${today}".log"`
echo Logfile name : $LOG_FILE
echo
echo Logfile name : $LOG_FILE  > $LOG_FILE
echo  >> $LOG_FILE

########### Set whether to use parallel or sequential processing.

DO_PAR=`which parallel | wc -l`

if [ $DO_PAR -eq 1 ]; then
    DO_FOR=0
    echo "Attempting run in parallel with "${JOBS}" cores..." >> $LOG_FILE
else
    DO_FOR=1
    echo "Attempting run in series..." >> $LOG_FILE
fi

######## ++++ SET event list to loop through ++#########
file_list=`echo $DATA_LOC"/*sac"`
num_file=`echo $file_list | wc -w`

# This will catch events with no or less than 1 processed files which are useless anyway.
if [ $num_file -le 1 ]; then
    echo  >> $LOG_FILE
    echo "No sac files to processing in :"$DATA_LOC"/"  >> $LOG_FILE
    echo "Exiting postprocessing..."
    echo  >> $LOG_FILE
    exit
else
    echo  >> $LOG_FILE
    echo "We have "$num_file" to process..."  >> $LOG_FILE
    echo  >> $LOG_FILE
fi


# Make directory for formatted data if it does not exist.
if [ ! -d $DATA_LOC"/formatted" ]; then
    echo "Making formatted data directory."  >> $LOG_FILE
	mkdir -p $DATA_LOC"/formatted"
fi

# Start by getting CMT params
CMT_FILE=$DATA_LOC"/CMTSOLUTION"

if [ ! -f $CMT_FILE ]; then
    echo "Some issue... "$DATA_LOC"/CMTSOLUTION not present..."
    echo "Exiting postprocessing..."
    exit    
fi

# Then STATION params
STAT_FILE=$DATA_LOC"/STATIONS"
if [ ! -f $STAT_FILE ]; then
    echo "Some issue... "$DATA_LOC"/STATIONS not present..."
    echo "Exiting postprocessing..."
    exit    
fi

######################################## GET CMT PARAMS #######################
echo  >> $LOG_FILE
echo "Collecting CMT params..."  >> $LOG_FILE
echo  >> $LOG_FILE

evt_params=`awk 'NR==1' $CMT_FILE` 
oyear=`echo $evt_params | awk '{print $2}'`
omonth=`echo $evt_params | awk '{ printf("%02d\n", $3) }'`
oday=`echo $evt_params | awk '{ printf("%02d\n", $4) }'`
ohr=`echo $evt_params | awk '{ printf("%02d\n", $5) }'`
omin=`echo $evt_params | awk '{ printf("%02d\n", $6) }'`
osec=`echo $evt_params | awk '{printf "%.3f\n", $7}' | awk -F. '{printf("%02d\n", $1)}'`
omsec=`echo $evt_params | awk '{printf "%.3f\n", $7}' | awk -F.  '{print $2}'`
evlat=`grep "latitude" $CMT_FILE | awk '{print $2}'`
evlon=`grep "longitude" $CMT_FILE | awk '{print $2}'`
evmag=`echo $evt_params | awk '{print $12}'`
evdep=`grep "depth" $CMT_FILE | awk '{print $2}'`
evname=`grep "event name" $CMT_FILE | awk '{print $3}'`
tshift=`grep "shift" $CMT_FILE | awk '{print $3}'`
# Try to make sure we dont get hdur = 0
hdur=`grep "duration" $CMT_FILE* | awk '{print $3}' | sort -n -r | awk 'NR==1'`
hdur_positive=`echo "$hdur > 0" | bc`

if [ $hdur_positive -eq 0 ]; then
    echo  >> $LOG_FILE
    echo "Half duration zero found..."  >> $LOG_FILE
    echo "Exiting postprocessing..."  >> $LOG_FILE
    echo  >> $LOG_FILE
    exit
else
    echo  >> $LOG_FILE
    echo "Postive half duration found..."  >> $LOG_FILE
    echo  >> $LOG_FILE
fi

echo  >> $LOG_FILE
echo oyear, omonth, oday, ohr, omin, osec.omsec, tshift, hdur, evname, evlat, evlon, evdep  >> $LOG_FILE
echo $oyear, $omonth, $oday, $ohr, $omin, $osec"."$omsec, $tshift, $hdur, $evname, $evlat, $evlon, $evdep  >> $LOG_FILE
echo  >> $LOG_FILE

# Get time shifted time
# echo python3 ${HOME}/bin/get_CMTSOLUTION_tshift.py $oyear"-"$omonth"-"$oday"T"$ohr":"$omin":"$osec"."$omsec" "$tshift

time_shifted=`${HOME}/anaconda3/bin/python ${HOME}/bin/get_CMTSOLUTION_tshift.py $oyear"-"$omonth"-"$oday"T"$ohr":"$omin":"$osec"."$omsec $tshift`

oyear1=`echo $time_shifted | awk '{print $1}'`
omonth1=`echo $time_shifted | awk '{ printf("%02d\n", $2) }'`
oday1=`echo $time_shifted | awk '{ printf("%02d\n", $3) }'`
ojday1=`echo $time_shifted | awk '{ printf("%03d\n", $4) }'`
ohr1=`echo $time_shifted | awk '{ printf("%02d\n", $5) }'`
omin1=`echo $time_shifted | awk '{ printf("%02d\n", $6) }'`
osec1=`echo $time_shifted | awk '{ printf("%02d\n", $7) }'`
omsec1=`echo $time_shifted | awk '{print substr($8,1,3)}'`

echo  >> $LOG_FILE
echo oyear1, omonth1, oday1, ojday1, ohr1, omin1, osec1.omsec1  >> $LOG_FILE
echo $oyear1, $omonth1, $oday1, $ojday1, $ohr1, $omin1, $osec1"."$omsec1  >> $LOG_FILE
echo  >> $LOG_FILE

echo $oyear > $CMT_FILE"_summary.temp"
echo $omonth >> $CMT_FILE"_summary.temp"
echo $oday >> $CMT_FILE"_summary.temp"
echo $ohr >> $CMT_FILE"_summary.temp"
echo $omin >> $CMT_FILE"_summary.temp"
echo $osec >> $CMT_FILE"_summary.temp"
echo $omsec >> $CMT_FILE"_summary.temp"
echo $evname >> $CMT_FILE"_summary.temp"
echo $evlat >> $CMT_FILE"_summary.temp"
echo $evlon >> $CMT_FILE"_summary.temp"
echo $evdep >> $CMT_FILE"_summary.temp"
echo $tshift >> $CMT_FILE"_summary.temp"
echo $hdur >> $CMT_FILE"_summary.temp"
echo $oyear1 >> $CMT_FILE"_summary.temp"
echo $omonth1 >> $CMT_FILE"_summary.temp"
echo $oday1 >> $CMT_FILE"_summary.temp"
echo $ojday1 >> $CMT_FILE"_summary.temp"
echo $ohr1 >> $CMT_FILE"_summary.temp"
echo $omin1 >> $CMT_FILE"_summary.temp"
echo $osec1 >> $CMT_FILE"_summary.temp"
echo $omsec1  >> $CMT_FILE"_summary.temp"

if [ ! -d $DATA_LOC"/formatted/"$oyear1$omonth1$oday1$ohr1$omin1$osec1 ]; then
	mkdir -p $DATA_LOC"/formatted/"$oyear1$omonth1$oday1$ohr1$omin1$osec1
fi

##############################################################

if [ $DO_FOR -eq 1 ]; then
    # For loop version of the code, sequential processing:
    echo "Processing files in :"$DATA_LOC"/ USING FOR-LOOP..."  >> $LOG_FILE
    echo  >> $LOG_FILE
    for file_dir in $file_list; do
        # echo $file_dir  >> $LOG_FILE
        # echo " "  >> $LOG_FILE
        file=`echo $file_dir | awk -F/ '{print $NF}'`
        # echo "${HOME}/bin/proc_specfem_seis_file.sh "$DATA_LOC $file $CHANNEL
        ${HOME}/bin/proc_specfem_seis_file.sh $DATA_LOC $file $CHANNEL  >> $LOG_FILE
    done
fi

if [ $DO_PAR -eq 1 ]; then
    # Parallel version of the code, Parallel processing:
    echo "Processing files in :"$DATA_LOC"/ USING PARALLEL..."  >> $LOG_FILE
    for file_dir in $file_list; do
        # echo $file_dir  >> $LOG_FILE
        # echo " "  >> $LOG_FILE
        file=`echo $file_dir| awk -F/ '{print $NF}'`
        # echo "${HOME}/bin/parallel --semaphore --id ${USER} --jobs "${JOBS}" ${HOME}/bin/proc_specfem_seis_file.sh "$DATA_LOC $file $CHANNEL
        ${HOME}/bin/parallel --semaphore --id ${USER} --jobs ${JOBS} --delay 2 --timeout 5000% ${HOME}/bin/proc_specfem_seis_file.sh $DATA_LOC $file $CHANNEL  >> $LOG_FILE
    done
    ${HOME}/bin/parallel --semaphore --wait --id ${USER}  >> $LOG_FILE
fi

#### Cleanup
rm $CMT_FILE"_summary.temp" 
rm -rf $XDG_CACHE_HOME"/parallel/semaphores/id-"${USER}
rm -rf $XDG_CACHE_HOME"/parallel/tmp"
rm -rf $HOME"/.parallel/semaphores/id-"${USER}
rm -rf tmp.log ttimes.lst

echo " "  >> $LOG_FILE
echo "Specfem seismogram processing complete for "$DATA_LOC" "  >> $LOG_FILE
echo " "  >> $LOG_FILE





