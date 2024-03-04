#!/bin/bash

if [ "$#" -lt 3 ];
then
	echo " Please specifiy Data location and event to process  "
	echo " USEAGE:   obspyDMT_data_proc_event.sh <DATA_LOC> <EVENT> <CHANNEL>"
	echo " EXAMPLE:  obspyDMT_data_proc_event.sh /home/aboyce/d_data_and_docs/dmt/e2008 20081230_194956.a LH"
	echo " "
	exit
else
    ######## ++++ SET data location and event +++++++++++++++#########
	DATA_LOC=$1
	EVENT=$2
	CHANNEL=$3
	echo "Processing "$EVENT", within: "$DATA_LOC" assuming data channel: "$CHANNEL
fi

HOME=`echo ~`
source ${HOME}/.bash_profile
subdir=`echo ${HOME}/d_data_and_docs/dmt`

catalog_table=`echo $DATA_LOC/EVENTS-INFO/catalog.txt`
# catalog_jday=`echo $DATA_LOC/EVENTS-INFO/catalog_jday.txt`
proc_event_dir=`echo $DATA_LOC/$EVENT/processed` # OR: processed

# Check if proc_event_directory exists
# if not there is no raw data usually... So exit...

if [ ! -d $proc_event_dir ]; then
	echo "No directory: "$proc_event_dir
	echo "Likely no raw data... Exiting..."
	echo
	exit
fi


############# Finding event details: ##########################

event_details=`grep $EVENT $catalog_table | uniq | awk 'NR==1'` # Assuming the top one is correct.
#number,event_id,datetime,latitude,longitude,depth,magnitude,magnitude_type,author,flynn_region,mrr,mtt,mpp,mrt,mrp,mtp,stf_func,stf_duration,t1,t2
KEVNM=`echo $event_details | awk -F, '{print $2}'`
evlat=`echo $event_details | awk -F, '{print $4}'`
evlon=`echo $event_details | awk -F, '{print $5}'`
evdep=`echo $event_details | awk -F, '{print $6}'`
mw=`echo $event_details | awk -F, '{print $7}'`
mwt=`echo $event_details | awk -F, '{print $8}'`

ev_timing_details=`echo $event_details | awk -F, '{print $3}' | sed 's/\:/ /g' | sed 's/\-/ /g' | sed 's/\T/ /g'`
Year=`echo $ev_timing_details | awk '{print $1}'`
Month=`echo $ev_timing_details | awk '{print $2}'`
Day=`echo $ev_timing_details | awk '{print $3}'`
Hour=`echo $ev_timing_details | awk '{print $4}'`
Minute=`echo $ev_timing_details | awk '{print $5}'`
Second=`echo $ev_timing_details | awk '{print $6}' | sed 's/\./ /g' | awk '{print $1}'`
Msec=`echo $ev_timing_details | awk '{print $6}' | sed 's/\./ /g' | awk '{print substr($2,1,3)}'`
# Use a new python script to grab the julian day. Havent compiled the C or fortran one.

Jday=`python3 ${HOME}/bin/get_julday.py $Year"-"$Month"-"$Day"T"$Hour":"$Minute":"$Second`

#echo "Working on event: " $Year $Jday $Hour $Minute $Second $Msec $evlat $evlon $evdep $mw

saclst stla stlo stel f $proc_event_dir/*.?H? > $proc_event_dir/station_locs.out

sac_file_list=`echo $proc_event_dir/*.?H?`
num_files=`echo $sac_file_list | wc -w`

# THis will catch events with no or less than 4 processed files which are useless anyway.
if [ $num_files -lt 1 ]; then
    echo "Not enough files to bother processing for :"$proc_event_dir
    exit
fi

# Make a directory for the formatted data using the usual naming convention.
form_event_dir=`echo $DATA_LOC/formatted/$Year$Month$Day$Hour$Minute$Second`
# Make directory for formatted data
if [ ! -d $form_event_dir ]; then
	mkdir -p $form_event_dir
fi

######### Loop over individual sac files #############

for sacfile in $sac_file_list; do
	
	#echo "Processing: "$sacfile
	
	# locate resp file
	respfile=`echo $sacfile | sed 's/processed\//resp\/STXML./g'` # OR: processed

	NTWRK=`echo $sacfile | awk -F/ '{print $NF}' | awk -F. '{print $1}'`
	KSTNM=`echo $sacfile | awk -F/ '{print $NF}' | awk -F. '{print $2}'`
	KSTLOC=`echo $sacfile | awk -F/ '{print $NF}' | awk -F. '{printf "%02d\n", $3}'`
	# KCHAN=$simresp_KCHAN
	KCHAN=`echo $sacfile | awk -F/ '{print $NF}' | awk -F. '{print $4}'`
	
	# Want stel in km. obspyDMT outputs in m so divide by 1000.
	stat_details=`grep $sacfile $proc_event_dir/station_locs.out`
	stlat=`echo $stat_details | awk '{print $2}'`
	stlon=`echo $stat_details | awk '{print $3}'`
	stel=`echo $stat_details | awk '{print $4}'` # Do not convert STEL to meters
	
	outfile=$form_event_dir/$Year$Month$Day$Hour$Minute$Second"_"$NTWRK"_"$KSTNM"."$KSTLOC"."$KCHAN
	
	# echo $NTWRK $KSTNM $KSTLOC $KCHAN $stlat $stlon $stel
	# echo " "
	# echo $sacfile $respfile $outfile
	
	# Now use obspy to deal with the instrument response using trace.simulate - OUTPUT TO DISPLACEMENT
	python3 ${HOME}/bin/obspyDMT_rem_inst_resp.py $sacfile $respfile $outfile $CHANNEL
	
	# Now make sure headers are correct.
	# keep stel in meters.

	SACMAC=$form_event_dir"/mac.sm"
	echo r $outfile > $SACMAC
	echo synchronize >> $SACMAC
	echo rmean >> $SACMAC
	echo rtrend >> $SACMAC
	echo taper type hanning width 0.015 >> $SACMAC
	if [ $CHANNEL = "LH" ]; then
		echo interpolate delta 1 >> $SACMAC
	elif [ $CHANNEL = "BH" ]; then
		echo interpolate delta 0.02 >> $SACMAC
	elif [ $CHANNEL = "HH" ]; then
		echo interpolate delta 0.005 >> $SACMAC
	fi
	# echo taper width 0.1 >> $SACMAC
	echo ch idep idisp >> $SACMAC		# remove resp - outputs to displacement
	echo ch o >> $SACMAC
	echo ch evla $evlat evlo $evlon evdp $evdep >> $SACMAC
	echo ch stla $stlat stlo $stlon stel $stel >> $SACMAC
	echo ch o gmt $Year $Jday $Hour $Minute $Second $Msec >> $SACMAC
	echo ch kevnm $KEVNM >> $SACMAC
	echo ch khole $KSTLOC >> $SACMAC
	echo ch lovrok true >> $SACMAC
	echo ch lcalda true >> $SACMAC
	echo w over >> $SACMAC
	echo quit >> $SACMAC
	echo '' >> $SACMAC
	
			
	echo "m "$form_event_dir"/mac.sm" | sac
		
	rm $form_event_dir/mac.sm	
done

rm $proc_event_dir/station_locs.out



echo " "
echo " Instrument correction complete for "$EVENT", within: "$DATA_LOC
echo " "

####################################################################
# 					Rotate components using obspy
####################################################################

python3 ${HOME}/bin/obspy_rotate_comps.py $form_event_dir $CHANNEL


echo " "
echo " Component rotation complete for "$EVENT", within: "$DATA_LOC
echo " "

####################################################################
# 					Rename components to match request
####################################################################

sac_file_list=`echo $form_event_dir"/*.?H?"`

for sacfile in $sac_file_list; do
	newfile=`echo $sacfile  | awk -v channel=$CHANNEL '{start=length($0)-2; end=length($0)-1; $0 = substr($0, 1, start-1) channel substr($0, end+1); print $0}'`
	mv $sacfile $newfile
done

echo " "
echo " Component renaming complete for "$EVENT", within: "$DATA_LOC
echo " "



