#!/bin/bash
#SBATCH --job-name=obs_DMT
#SBATCH --partition=Cascade-flix,Cascade
#SBATCH --time=168:00:00
#SBATCH --nodes=1
### SBATCH --ntasks=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=4
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
# export PATH=".:${HOME}:${HOME}/anaconda3/bin:$PATH"
# source activate env3.12


year="YXXR"
#sources="AUSPASS ETH GEOFON GEONET ICGC INGV IPGP IRIS KNMI KOERI LMU NCEDC NIEP NOA ORFEUS RESIF SCEDC UIB-NORSAR USP"
sources="IRIS KNMI KOERI LMU NCEDC NIEP NOA ORFEUS RESIF SCEDC UIB-NORSAR USP"
evt_cat="GCMT_COMBO"
nwk=`echo \?\?`
cha=`echo LH\?,BH\?,HH\?`
loc=`echo \*`
sta=`echo \*`
p_set="60" # preset (seconds before epicentral time)
o_set="7200" # offset (seconds after epicentral time)
min_mg="5.5" # min magnitude
max_mg="7.5" # min magnitude
min_dp="0" # min depth
max_dp="700" # max depth
npr="8" # parallel requests
numcpus=$(( ${SLURM_NTASKS_PER_NODE} * ${SLURM_CPUS_PER_TASK} ))
npp=${numcpus} # parallel processes

echo " "   >> date_${SLURM_JOB_ID}.txt
echo "obspyDMT request with the following parameters:" >> date_${SLURM_JOB_ID}.txt  
echo " "   >> date_${SLURM_JOB_ID}.txt

echo "--datapath e"$year   >> date_${SLURM_JOB_ID}.txt
echo "--data_source "$sources  >> date_${SLURM_JOB_ID}.txt
echo "--event_catalog "$evt_cat  >> date_${SLURM_JOB_ID}.txt
echo "--net "$nwk  >> date_${SLURM_JOB_ID}.txt
echo "--cha "$cha  >> date_${SLURM_JOB_ID}.txt
echo "--loc "$loc  >> date_${SLURM_JOB_ID}.txt
echo "--sta "$sta  >> date_${SLURM_JOB_ID}.txt
echo "--preset "$p_set  >> date_${SLURM_JOB_ID}.txt
echo "--offset "$o_set  >> date_${SLURM_JOB_ID}.txt
echo "--min_mag "$min_mg  >> date_${SLURM_JOB_ID}.txt
echo "--max_mag "$max_mg  >> date_${SLURM_JOB_ID}.txt
echo "--min_depth "$min_dp  >> date_${SLURM_JOB_ID}.txt
echo "--max_depth "$max_dp  >> date_${SLURM_JOB_ID}.txt
echo "--min_date=$year-01-01T00:00:00"  >> date_${SLURM_JOB_ID}.txt
echo "--max_date=$year-12-31T00:00:00"  >> date_${SLURM_JOB_ID}.txt
echo "--waveform_format sac"  >> date_${SLURM_JOB_ID}.txt
echo "--response True"  >> date_${SLURM_JOB_ID}.txt
echo "--pre_process process_unit_default" >> date_${SLURM_JOB_ID}.txt
echo "--req_parallel"  >> date_${SLURM_JOB_ID}.txt
echo "--req_np "$npr >> date_${SLURM_JOB_ID}.txt
echo "--parallel_process"  >> date_${SLURM_JOB_ID}.txt
echo "--process_np "$npp >> date_${SLURM_JOB_ID}.txt
echo " " >> date_${SLURM_JOB_ID}.txt

# Make Data directory
if [ ! -d e$year ]; then
    echo "Making Data directory: e"$year >> date_${SLURM_JOB_ID}.txt
    mkdir "e"$year
else
    echo "Data directory e"$year" exists" >> date_${SLURM_JOB_ID}.txt
fi

echo " " >> date_${SLURM_JOB_ID}.txt
echo "Starting download..." >> date_${SLURM_JOB_ID}.txt
echo " " >> date_${SLURM_JOB_ID}.txt
# Run the request...
# echo obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-12-30T00:00:00" --max_date="$year-12-31T00:00:00" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> date_${SLURM_JOB_ID}.txt

for source in $sources; do
    echo obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-01-01T00:00:00" --max_date="$year-12-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> date_${SLURM_JOB_ID}.txt
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-01-01T00:00:00" --max_date="$year-01-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M01_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-02-01T00:00:00" --max_date="$year-02-28T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M02_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-03-01T00:00:00" --max_date="$year-03-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M03_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-04-01T00:00:00" --max_date="$year-04-30T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M04_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-05-01T00:00:00" --max_date="$year-05-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M05_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-06-01T00:00:00" --max_date="$year-06-30T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M06_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-07-01T00:00:00" --max_date="$year-07-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M07_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-08-01T00:00:00" --max_date="$year-08-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M08_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-09-01T00:00:00" --max_date="$year-09-30T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M09_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-10-01T00:00:00" --max_date="$year-10-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M10_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-11-01T00:00:00" --max_date="$year-11-30T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M11_${SLURM_JOB_ID}.log
    obspyDMT --datapath e$year --data_source $source --event_catalog $evt_cat --net $nwk --cha $cha --sta "*" --loc "*" --preset $p_set --offset $o_set --min_mag $min_mg --max_mag $max_mg --min_depth $min_dp --max_depth $max_dp --min_date="$year-12-01T00:00:00" --max_date="$year-12-31T23:59:59" --waveform_format "sac" --response "True" --pre_process "process_unit_default" --req_parallel --req_np $npr --parallel_process --process_np $npp >> $LOG_LOC/download_${source}_e${year}_M12_${SLURM_JOB_ID}.log

done
echo " " >> date_${SLURM_JOB_ID}.txt
echo "Finished download..." >> date_${SLURM_JOB_ID}.txt
echo " " >> date_${SLURM_JOB_ID}.txt





