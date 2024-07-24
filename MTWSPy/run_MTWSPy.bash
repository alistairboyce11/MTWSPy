#!/bin/bash

HOME=`echo ~`
CODE_HOME=`pwd`
# source ${HOME}/.bash_profile
conda activate MTWSPy
export PYTHONPATH=./MTWSPy:$PYTHONPATH

###################################################################

# Main code
python ./MTWSPy/MTWSPy_main.py

# Post Processing on *.tdl files
python ./MTWSPy/post_processing/process_tdl_files.py 

# Create inversion ready file on Linux (Case sensitive paths issue on OSX)

# TauP - installation required:
taup_exists=`which taup_path | wc -w`

if [[ `uname -s` !=  "Darwin" ]] && [[ $taup_exists -eq 1 ]]; then 
    python ./MTWSPy/post_processing/create_inv_files.py
fi 

# GMT Plotting - installation required:
gmt_exists=`which gmt | wc -w`

if [[ $gmt_exists -eq 1 ]]; then
    echo "Plotting results using GMT..."
    ./MTWSPy/post_processing/plot_phase.gmt S NO
fi

####################################################################
exit
