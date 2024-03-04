#!/bin/bash

if [ "$#" -lt 3 ];
then
	echo " Please specifiy Data location and file to process  "
	echo " USEAGE:   proc_specfem_seis_file.sh <DATA_LOC> <FILE> <CHANNEL>"
	echo " EXAMPLE:  proc_specfem_seis_file.sh  ./OUTPUT_FILES_3D_seis_hd0 II.AAK.BXZ.sem.sac LH"
	echo " "
	exit
else
    ######## ++++ SET data location and file +++++++++++++++#########
	DATA_LOC=$1
	FILE=$2
	CHANNEL=$3
    echo " "
	echo "Processing "$FILE", within: "$DATA_LOC" assuming data channel: "$CHANNEL
    echo " "
fi

HOME=`echo ~`
source ${HOME}/.bash_profile
# Start by getting CMT params
CMT_FILE=$DATA_LOC"/CMTSOLUTION"

if [ ! -f $CMT_FILE ]; then
    echo "Some issue... "$DATA_LOC"/CMTSOLUTION not present..."
    echo "Exiting...."
    exit    
fi

# Then STATION params
STAT_FILE=$DATA_LOC"/STATIONS"
if [ ! -f $STAT_FILE ]; then
    echo "Some issue... "$DATA_LOC"/STATIONS not present..."
    echo "Exiting...."
    exit    
fi


NTWRK=`echo $FILE | awk -F. '{print $1}'`
KSTNM=`echo $FILE | awk -F. '{print $2}'`
CHAN_IN=`echo $FILE | awk -F. '{print $3}' | awk '{print substr($0,1,2)}'`
COMP_IN=`echo $FILE | awk -F. '{print $3}' | awk '{print substr($0,3,1)}'`
KSTLOC="00"
KCMPNM=$CHANNEL$COMP_IN

# echo $NTWRK $KSTNM $CHAN_IN $COMP_IN $KCMPNM


############## READ CMT PARAMS ########################


if [ -f $CMT_FILE"_summary.temp" ]; then

    oyear=`awk 'NR==1' $CMT_FILE"_summary.temp"`
    omonth=`awk 'NR==2' $CMT_FILE"_summary.temp"`
    oday=`awk 'NR==3' $CMT_FILE"_summary.temp"`
    ohr=`awk 'NR==4' $CMT_FILE"_summary.temp"`
    omin=`awk 'NR==5' $CMT_FILE"_summary.temp"`
    osec=`awk 'NR==6' $CMT_FILE"_summary.temp"`
    omsec=`awk 'NR==7' $CMT_FILE"_summary.temp"`
    evname=`awk 'NR==8' $CMT_FILE"_summary.temp"`
    evlat=`awk 'NR==9' $CMT_FILE"_summary.temp"`
    evlon=`awk 'NR==10' $CMT_FILE"_summary.temp"`
    evdep=`awk 'NR==11' $CMT_FILE"_summary.temp"`
    tshift=`awk 'NR==12' $CMT_FILE"_summary.temp"`
    hdur=`awk 'NR==13' $CMT_FILE"_summary.temp"`
    oyear1=`awk 'NR==14' $CMT_FILE"_summary.temp"`
    omonth1=`awk 'NR==15' $CMT_FILE"_summary.temp"`
    oday1=`awk 'NR==16' $CMT_FILE"_summary.temp"`
    ojday1=`awk 'NR==17' $CMT_FILE"_summary.temp"`
    ohr1=`awk 'NR==18' $CMT_FILE"_summary.temp"`
    omin1=`awk 'NR==19' $CMT_FILE"_summary.temp"`
    osec1=`awk 'NR==20' $CMT_FILE"_summary.temp"`
    omsec1=`awk 'NR==21' $CMT_FILE"_summary.temp"`

else
    # No CMT summary file present...
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

    time_shifted=`${HOME}/anaconda3/bin/python ${HOME}/bin/get_CMTSOLUTION_tshift.py $oyear"-"$omonth"-"$oday"T"$ohr":"$omin":"$osec"."$omsec $tshift`

    oyear1=`echo $time_shifted | awk '{print $1}'`
    omonth1=`echo $time_shifted | awk '{ printf("%02d\n", $2) }'`
    oday1=`echo $time_shifted | awk '{ printf("%02d\n", $3) }'`
    ojday1=`echo $time_shifted | awk '{ printf("%03d\n", $4) }'`
    ohr1=`echo $time_shifted | awk '{ printf("%02d\n", $5) }'`
    omin1=`echo $time_shifted | awk '{ printf("%02d\n", $6) }'`
    osec1=`echo $time_shifted | awk '{ printf("%02d\n", $7) }'`
    omsec1=`echo $time_shifted | awk '{print substr($8,1,3)}'`


fi


# echo $oyear, $omonth, $oday, $ohr, $omin, $osec"."$omsec, $tshift, $hdur, $evname, $evlat, $evlon, $evdep
# echo $oyear1, $omonth1, $oday1, $ojday1, $ohr1, $omin1, $osec1"."$omsec1

# Make directory for formatted data if it does not exist.
if [ ! -d $DATA_LOC"/formatted" ]; then
	mkdir -p $DATA_LOC"/formatted"
fi

if [ ! -d $DATA_LOC"/formatted/"$oyear1$omonth1$oday1$ohr1$omin1$osec1 ]; then
	mkdir -p $DATA_LOC"/formatted/"$oyear1$omonth1$oday1$ohr1$omin1$osec1
fi


# # Get info frm sac headers
# headers=`saclst b delta dist gcarc f $DATA_LOC"/"$FILE`
# begin_time=`echo $headers | awk '{print $2}'`
# delta=`echo $headers | awk '{print $3}'`
# dist=`echo $headers | awk '{print $4}'`
# gcarc=`echo $headers | awk '{print $5}'`

# # echo $begin_time $delta $dist $gcarc


################### Get predicted travel times ##################

# p_pred_times=`phtimes.csh $evdep $gcarc P`
# s_pred_times=`phtimes.csh $evdep $gcarc S`

# Pph=`echo $p_pred_times | awk '{print $1}'`
# Ptime=`echo $p_pred_times | awk '{print $2}'`
# Sph=`echo $s_pred_times | awk '{print $1}'`
# Stime=`echo $s_pred_times | awk '{print $2}'`

# echo $Pph,$Ptime,$Sph,$Stime


############### Get station parameters ###########################

sta_params=`cat $STAT_FILE | awk  -v var1=$KSTNM -v var2=$NTWRK '$1 == var1 && $2 == var2 {print}'`

stlat=`echo $sta_params | awk '{print $3}'`
stlon=`echo $sta_params | awk '{print $4}'`
stel=`echo $sta_params | awk '{print $5}'` # in Meters
stdp=`echo $sta_params | awk '{print $6}'` # in Meters

# echo $stlat $stlon $stel $sta_dep

############### Set files names #####################
echo "Processing file:             "$DATA_LOC"/"$FILE
# outfile=$DATA_LOC"/formatted/"$FILE
# Not sure which of these it is yet............
# outfile=$DATA_LOC"/formatted/"$oyear$omonth$oday$ohr$omin$osec"/"$oyear$omonth$oday$ohr$omin$osec"_"$NTWRK"_"$KSTNM"."$KSTLOC"."$KCMPNM
outfile=$DATA_LOC"/formatted/"$oyear1$omonth1$oday1$ohr1$omin1$osec1"/"$oyear1$omonth1$oday1$ohr1$omin1$osec1"_"$NTWRK"_"$KSTNM"."$KSTLOC"."$KCMPNM
echo "To outfile:                  "$outfile

######### Now write sac macro to sort out the postprocessing ###########


SACMAC=$DATA_LOC"/formatted/mac_"$NTWRK"_"$KSTNM"_"$KCMPNM".sm"
#### Reading
echo r $DATA_LOC"/"$FILE > $SACMAC
echo synchronize >> $SACMAC
echo rtrend >> $SACMAC
echo rmean >> $SACMAC
echo taper type hanning width 0.015 >> $SACMAC
echo ch idep idisp >> $SACMAC # Specfem outputs displacement seismograms.
#### Make tshift
echo ch allt "-"$tshift >> $SACMAC	
echo ch nzyear $oyear1 nzjday $ojday1 nzhour $ohr1 nzmin $omin1 nzsec $osec1 nzmsec $omsec1 >> $SACMAC	
echo ch o 0 >> $SACMAC	
#### Event details
echo ch evla $evlat evlo $evlon evdp $evdep >> $SACMAC
echo ch kevnm $evname >> $SACMAC
echo ch lovrok true >> $SACMAC
echo ch lcalda true >> $SACMAC
#### Convolve STF 
echo convolve tri $hdur centered on >> $SACMAC
#### Station details
echo ch stla $stlat stlo $stlon stel $stel stdp $stdp >> $SACMAC
echo ch kstnm $KSTNM knetwk $NTWRK kcmpnm $KCMPNM >> $SACMAC
echo ch khole $KSTLOC >> $SACMAC

#### Filter
echo rtrend >> $SACMAC
echo rmean >> $SACMAC
echo taper type hanning width 0.015 >> $SACMAC
if [ $CHANNEL = "LH" ]; then
    echo bp bu co 0.01 0.4 n 2 p 2 >> $SACMAC
else
    echo bp bu co 0.01 6 n 2 p 2 >> $SACMAC
fi
echo rtrend >> $SACMAC
echo rmean >> $SACMAC
echo taper type hanning width 0.015 >> $SACMAC
#### Interpolate
if [ $CHANNEL = "LH" ]; then
    echo interpolate delta 1 >> $SACMAC
else
    echo interpolate delta 0.05 >> $SACMAC
fi
#### Set arrival time headers
# echo evaluate to tmp1 $Ptime - $tshift >> $SACMAC
# echo evaluate to tmp2 $Stime - $tshift >> $SACMAC
# echo ch t1 %tmp1% t2 %tmp2% >> $SACMAC
# echo ch kt1 $Pph kt2 $Sph >> $SACMAC
#### SAVE
echo w $outfile >> $SACMAC
echo quit >> $SACMAC
echo '' >> $SACMAC

#### Execute macro
echo "m "$DATA_LOC"/formatted/mac_"$NTWRK"_"$KSTNM"_"$KCMPNM".sm" | sac
# cat $DATA_LOC"/formatted/mac_"$NTWRK"_"$KSTNM"_"$KCMPNM".sm"

#### Cleanup
rm $DATA_LOC"/formatted/mac_"$NTWRK"_"$KSTNM"_"$KCMPNM".sm"	

echo " "
echo "Seismogram processing complete for "$FILE", within: "$DATA_LOC" assuming data channel: "$CHANNEL
echo " "








