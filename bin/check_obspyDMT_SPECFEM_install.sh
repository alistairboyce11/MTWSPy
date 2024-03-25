#!/bin/bash
#----------- REQUIREMENTS --------------#
# python3 ['obspy','numpy', 'scipy', 'matplotlib', 'sys', 'glob', 'shutil', 'os', 'time', 'warnings']
# obspyDMT : https://kasra-hosseini.github.io/obspyDMT/
# parallel : https://www.gnu.org/software/parallel/
# sac : http://www.iris.edu/dms/nodes/dmc/software/downloads/sac
# saclst : http://geophysics.eas.gatech.edu/people/zpeng/Software/sac_msc.tar.gz
# ttimes : https://ds.iris.edu/pub/programs/iaspei-tau/

#------ ObspyDMT Data processing. ------#
# parallel_proc_obspyDMT_seis.py
#------ Specfem Data processing. ------#
# parallel_proc_specfem_seis.py

# >> 


# Script to check for all scripts required to run the obspyDMT data processing
# conda activate env3.12

HOME=`echo ~`
source ${HOME}/.bash_profile
subdir=`echo ${HOME}/bin/`
SPECFEM_LOC=`echo ${HOME}/s_seismology/specfem/sge.globe/C_source_codes/i_intel`
if [ ! -d ${HOME}/tmp ]; then
    mkdir -p ${HOME}/tmp
    chmod 777 ${HOME}/tmp
fi
export XDG_CACHE_HOME=${HOME}/tmp
LOG_FILE_SH=${HOME}/tmp/tmp.check-sh
LOG_FILE_PY=${HOME}/tmp/tmp.check-py
LOG_FILE_PYDMT=${HOME}/tmp/tmp.check-pydmt
rm -rf $LOG_FILE_SH $LOG_FILE_PY $LOG_FILE_PYDMT

######################################################################################

# Locate python install:
py_loc=`which python`
py_ver=`$py_loc --version | awk -F. '{print $1"."$2}' | awk '{print $2}'`
py_ver_3=`echo "$py_ver > 2" |  bc`

echo  "------------------" >> $LOG_FILE_SH

if [ $py_ver_3 -eq 0 ]; then
    echo "Python3 not found..."  >> $LOG_FILE_SH
    exit
else
    echo "Python version "${py_ver}" found..."  >> $LOG_FILE_SH
fi
echo  "------------------" >> $LOG_FILE_SH

######################################################################################

### Check python packages necessary...
$py_loc ${subdir}/check_obspyDMT_SPECFEM_install.py > $LOG_FILE_PY
py_package_error=`tail -1 $LOG_FILE_PY`

if [ $py_package_error -eq 1 ]; then
    echo "Python package errors...."  >> $LOG_FILE_SH
    exit
else
    echo "Python packages found...."  >> $LOG_FILE_SH
fi
echo  "------------------" >> $LOG_FILE_SH

######################################################################################
### Check obspyDMT install:; 

obspyDMT_loc=`which obspyDMT`
obspyDMT_ver=`$obspyDMT_loc --version | wc -l`

if [ $obspyDMT_ver -eq 0 ]; then
    echo "obspyDMT not found...."  >> $LOG_FILE_SH
    exit
else
    echo "obspyDMT found at: "${obspyDMT_loc}  >> $LOG_FILE_SH
fi

$obspyDMT_loc --check  > $LOG_FILE_PYDMT

obspy_val=`grep 'obspy: not installed' $LOG_FILE_PYDMT | wc -l`
numpy_val=`grep 'numpy: not installed' $LOG_FILE_PYDMT | wc -l`
scipy_val=`grep 'scipy: not installed' $LOG_FILE_PYDMT | wc -l`
matpl_val=`grep 'matplotlib: not installed' $LOG_FILE_PYDMT | wc -l`

obDMT_val=`echo $(( $obspy_val + $numpy_val + $scipy_val + $matpl_val ))`

if [ $obDMT_val -ne 0 ]; then
    echo "Python package errors in obspyDMT...."  >> $LOG_FILE_SH
    exit
else
    echo "Python packages found in obspyDMT...."  >> $LOG_FILE_SH
fi

######################################################################################

# echo  "------------------" >> $LOG_FILE_SH
# ### Obspy Data processing:
# parallel_yes=`which parallel | wc -l`
# sac_yes=`which sac | wc -l`
# saclst_yes=`which saclst | wc -l`
# ttime_yes=`which ttimes | wc -l`

# if [ $parallel_yes -eq 0 ]; then
#     echo "Parallel not found....."  >> $LOG_FILE_SH
#     exit
# else
#     echo "Parallel found..."  >> $LOG_FILE_SH
# fi

# if [ $sac_yes -eq 0 ]; then
#     echo "SAC not found....."  >> $LOG_FILE_SH
#     exit
# else
#     echo "SAC found..."  >> $LOG_FILE_SH
# fi

# if [ $saclst_yes -eq 0 ]; then
#     echo "Saclst not found....."  >> $LOG_FILE_SH
#     exit
# else
#     echo "Saclst found..."  >> $LOG_FILE_SH
# fi

# if [ $ttime_yes -eq 0 ]; then
#     echo "TTimes not found....."  >> $LOG_FILE_SH
#     exit
# else
#     echo "TTimes found..."  >> $LOG_FILE_SH
# fi
# echo  "------------------" >> $LOG_FILE_SH

######################################################################################

### ObspyDMT Data processing.

if [ ! -f ${subdir}/parallel_proc_obspyDMT_seis.py ]; then
    echo "parallel_proc_obspyDMT_seis.py not found....."  >> $LOG_FILE_SH
    exit
else
    echo "obspyDMT processing scripts found....."  >> $LOG_FILE_SH
fi
echo  "------------------" >> $LOG_FILE_SH

######################################################################################

### Specfem Data processing.

if [ ! -f ${subdir}/parallel_proc_specfem_seis.py ]; then
    echo "parallel_proc_specfem_seis.py not found....."  >> $LOG_FILE_SH
    exit
else
    echo "SPECFEM processing scripts found....."  >> $LOG_FILE_SH
fi
echo  "------------------" >> $LOG_FILE_SH

echo "YES - Codes/Dependencies/Scripts checked... continue :)"   >> $LOG_FILE_SH



