#!/bin/bash

if [ "$#" -lt 3 ];
then
    echo " Please specifiy Start date and years  "
    echo " USEAGE:   obspyDMT_data_par_proc_all.sh <DATA_LOC> <YEAR> <PARALLEL_PROCESSES> <CHANNEL>"
    echo " EXAMPLE:  obspyDMT_data_par_proc_all.sh /home/aboyce/d_data_and_docs/dmt/e2008 2008 96 LH"
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
     
    YEAR=$2
    JOBS=$3
    CHANNEL=$4
    echo
    echo "Processing "$DATA_LOC", year: "$YEAR", using "$JOBS" processes, assuming data channel: "$CHANNEL
fi

############# Log file ###########################

HOME=`echo ~`
source ${HOME}/.bash_profile
USER=`whoami`
export XDG_CACHE_HOME=${HOME}/tmp
subdir=`echo $HOME/d_data_and_docs/dmt`
LOG_LOC=`echo $subdir"/LOG_FILES"`
today=`date +%Y-%m-%d`

LOG_FILE=`echo $LOG_LOC"/process_e"${YEAR}"_NP_"${JOBS}"_CH_"${CHANNEL}"_"${today}".log"`
echo Logfile name : $LOG_FILE
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
event_list=`echo $DATA_LOC"/"$YEAR"????_??????.a"`
num_events=`echo $event_list | wc -w`

# This will catch events with no or less than 1 processed files which are useless anyway.
if [ $num_events -le 1 ]; then
    echo  >> $LOG_FILE
    echo "No events to processing in :"$DATA_LOC"/"$YEAR"*"  >> $LOG_FILE
    echo  >> $LOG_FILE
    exit
else
    echo  >> $LOG_FILE
    echo "We have "$num_events" to process..."  >> $LOG_FILE
    echo  >> $LOG_FILE
fi

# Make directory for formatted data if it does not exist.
if [ ! -d $DATA_LOC"/formatted" ]; then
    echo "Making formatted data directory."  >> $LOG_FILE
    mkdir $DATA_LOC"/formatted"		
fi


if [ $DO_FOR -eq 1 ]; then
    # For loop version of the code, sequential processing:
    echo "Processing events in :"$DATA_LOC"/"$YEAR" USING FOR-LOOP..."  >> $LOG_FILE
    echo  >> $LOG_FILE
    for event_dir in $event_list; do
        # echo $event_dir  >> $LOG_FILE
        # echo " "  >> $LOG_FILE
        event=`echo $event_dir | awk -F/ '{print $NF}'`
        # echo "${HOME}/bin/obspyDMT_data_proc_event.sh "$DATA_LOC $event $CHANNEL
        ${HOME}/bin/obspyDMT_data_proc_event.sh $DATA_LOC $event $CHANNEL  >> $LOG_FILE
    done
fi


if [ $DO_PAR -eq 1 ]; then
    # Parallel version of the code, Parallel processing:
    echo "Processing events in :"$DATA_LOC"/"$YEAR" USING PARALLEL..."  >> $LOG_FILE
    for event_dir in $event_list; do
        # echo $event_dir  >> $LOG_FILE
        # echo " "  >> $LOG_FILE
        event=`echo $event_dir| awk -F/ '{print $NF}'`
        # echo "${HOME}/bin/parallel --semaphore --id ${USER} --jobs "$JOBS" ${HOME}/bin/obspyDMT_data_proc_event.sh "$DATA_LOC $event $CHANNEL
        ${HOME}/bin/parallel --semaphore --id ${USER} --jobs ${JOBS} --delay 2 --timeout 5000% ${HOME}/bin/obspyDMT_data_proc_event.sh $DATA_LOC $event $CHANNEL  >> $LOG_FILE
    done
    ${HOME}/bin/parallel --semaphore --wait --id ${USER}  >> $LOG_FILE
fi

#### Cleanup
rm -rf $XDG_CACHE_HOME"/parallel/semaphores/id-"${USER}
rm -rf $XDG_CACHE_HOME"/parallel/tmp"
rm -rf $HOME"/.parallel/semaphores/id-"${USER}

echo " "  >> $LOG_FILE
echo "Instrument Correction and Rotation complete for "$DATA_LOC", year: "$YEAR  >> $LOG_FILE
echo " "  >> $LOG_FILE



