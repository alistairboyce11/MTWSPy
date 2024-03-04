#!/bin/bash
#SBATCH --job-name=pobs_DMT
#SBATCH --partition=Cascade,Cascade-flix
#SBATCH --time=168:00:00
#SBATCH --nodes=1
### SBATCH --ntasks=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
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
sources="AUSPASS ETH GEOFON GEONET ICGC INGV IPGP IRIS KNMI KOERI LMU NCEDC NIEP NOA ORFEUS RESIF SCEDC UIB-NORSAR USP"
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
npr="4" # parallel requests
numcpus=$(( ${SLURM_NTASKS_PER_NODE} * ${SLURM_CPUS_PER_TASK} ))
npp=${numcpus} # parallel processes

echo " "   >> date_${SLURM_JOB_ID}.txt
echo "obspyDMT processing with the following parameters:" >> date_${SLURM_JOB_ID}.txt  
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

# Check Data directory
if [ ! -d e$year ]; then
    echo "Data directory: e"$year" not present..." >> date_${SLURM_JOB_ID}.txt
    echo "Exiting..."
    exit
else
    echo "Data directory e"$year" exists" >> date_${SLURM_JOB_ID}.txt
fi

echo " " >> date_${SLURM_JOB_ID}.txt
echo "Start treatment of year ${year}, $(date)" >> date_${SLURM_JOB_ID}.txt
echo " " >> date_${SLURM_JOB_ID}.txt

echo -e "Starting obspyDMT postprocessing on "${SLURM_NTASKS_PER_NODE}" processors in datapath: ./e"${year} >> date_${SLURM_JOB_ID}.txt
echo ${HOME}/bin/obspyDMT_data_par_proc_all.sh ${HOME}/d_data_and_docs/dmt/e${year} ${year} ${SLURM_NTASKS_PER_NODE} ${cha} >> date_${SLURM_JOB_ID}.txt
${HOME}/bin/obspyDMT_data_par_proc_all.sh ${HOME}/d_data_and_docs/dmt/e${year} ${year} ${SLURM_NTASKS_PER_NODE} ${cha} >> date_${SLURM_JOB_ID}.txt
echo -e "Finished obspyDMT postprocessing in datapath: ./e"${year} >> date_${SLURM_JOB_ID}.txt

echo  >> date_${SLURM_JOB_ID}.txt
echo "End treatment of year ${year}, $(date)" >> date_${SLURM_JOB_ID}.txt
echo  >> date_${SLURM_JOB_ID}.txt


exit



