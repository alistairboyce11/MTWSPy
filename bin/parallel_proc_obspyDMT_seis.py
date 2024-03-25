# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:25:02 2024
@author: alistairboyce
if this fails: >> conda activate env3.12
"""

import time, sys, os, glob, shutil
import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace, UTCDateTime, read, Stream, read_inventory
import obspy.signal
import obspy.signal.rotate
import obspy.geodetics.base
import obspy.geodetics
from scipy.signal import convolve
import pandas as pd
import concurrent.futures

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def read_CMT_CATALOG(file):
    '''
    Read catalog.txt file into dict

    Parameters
    ----------
    file: str
        catalog.txt file name

    Returns
    ----------
    cmt_df : df 

    '''
    cmt_table = []

    f=open(file,'r')
    lines=f.readlines()
    f.close()


    for line in lines:
        if "----" in line:
            continue
        elif "Command" in line:
            continue
        elif "obspyDMT" in line:
            continue
        elif "#number" in line: 
            continue
        elif len(line.split(',')) == 1:
            continue
        else:
            if len(line.split('\n')[0].split(',')[1:]) == 19:
                cmt_table.append(line.split('\n')[0].split(',')[1:])
    
    column_headers = ['event_id','datetime','latitude','longitude','depth','magnitude','magnitude_type','author','flynn_region','mrr','mtt','mpp','mrt','mrp','mtp','stf_func','stf_duration','t1','t2']
    try:
        df = pd.DataFrame(cmt_table, columns = column_headers)
        df = df.drop_duplicates()

        ev_df = pd.DataFrame(columns = ['ev_id'], index = range(0,len(df)))

        for index, row in df.iterrows(): 
            eq_time = UTCDateTime(row['datetime'])
            ev_df['ev_id'][index] = f'{eq_time.year:4d}{eq_time.month:02d}{eq_time.day:02d}{eq_time.hour:02d}{eq_time.minute:02d}{eq_time.second:02d}'

        df = pd.merge(df, ev_df, left_index = True, right_index = True)    
    except:
        #Return empy dataframe if not possible.
        df = pd.DataFrame()

    return df


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def read_STATION_file(file):
    '''
    Read STATION file into pandas dataframe

    Parameters
    ----------
    file: str
        STATION file name
       
    Returns
    ----------
    stations_df : df 
    '''


    f=open(file,'r')
    lines=f.readlines()
    f.close()

    stat_table = []
    for line in lines:
        if len(line.split('\n')[0].split(',')) == 18: 
            stat_table.append(line.split('\n')[0].split(',')[0:17])

    column_headers = ['net_name', 'sta_name', 'sta_loc', 'sta_chan', 'sta_lat', 'sta_lon', 'sta_elev', 'sta_depth', 'data_source', 'evt_dir', 'evlat', 'evlon', 'evdep', 'evmag', 'orientation', 'dip', 'unknown']
    
    try:
        df = pd.DataFrame(stat_table, columns = column_headers)
        df = df.drop_duplicates()
    except:
        #Return empy dataframe if not possible.
        df = pd.DataFrame()

    return df




# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_sacfile_dict(input_sacfile, cmt_df, stations_df, channel):
    '''
    Read input_sacfilename get parameter for processing.

    Parameters
    ----------
    input_sacfile: str
        filename of input sac file
    cmt_df : dict
        CMT catalog dataframe
    stations_df : df 
        Stations dataframe
    channel : str
        desired data channel after processing
       
    Returns
    ----------
    sacfile_dict : dict 
    input_sacfile : input sacfile name
    ntwrk : network name
    kstnm : station name
    chan_in : input channel
    comp_in : input component
    kstloc : station location
     : station component
    output_sacfile : output sacfile name
    stlat : station lat
    stlon : station lon
    stel : station elevation
    stdp : station depth
    resp_file : input response file name
    '''

    sacfile_dict = {}
    sacfile_dict['input_sacfile'] = input_sacfile

    network_name = input_sacfile.split('/')[-1].split('.')[0]
    station_name = input_sacfile.split('/')[-1].split('.')[1]

    sacfile_dict['ntwrk'] = network_name 
    sacfile_dict['kstnm'] = station_name 
    sacfile_dict['chan_in'] = input_sacfile.split('/')[-1].split('.')[3][0:2]
    sacfile_dict['comp_in'] = input_sacfile.split('/')[-1].split('.')[3][2]
    loc_in = input_sacfile.split('/')[-1].split('.')[2]
    sacfile_dict['loc_in'] = loc_in

    sacfile_dict['kstloc'] = str('00')
    sacfile_dict['kcmpnm'] = str(channel) + sacfile_dict['comp_in']

    # $oyear1$omonth1$oday1$ohr1$omin1$osec1"_"$ntwrk"_"$kstnm"."$kstloc"."$kcmpnm
    sacfile_dict['output_sacfile'] = cmt_df['ev_id'].to_list()[0] + '_' + sacfile_dict['ntwrk'] + '_' + sacfile_dict['kstnm'] + '.' + sacfile_dict['kstloc'] + '.' + sacfile_dict['kcmpnm']

    net_df = stations_df.query("net_name == @network_name")
    sta_df = net_df.query("sta_name == @station_name")
    loc_df = sta_df.query("sta_loc == @loc_in")
    channel_name_in = f'{sacfile_dict['chan_in']}{sacfile_dict['comp_in']}'
    chan_df = loc_df.query("sta_chan == @channel_name_in")


    if len(chan_df) != 1: 
        print(f'Failed searching for : {network_name}, {station_name}, {channel_name_in} in STATIONS file')
        print('Return empty dict....')
        return {}
    else:
        sacfile_dict['stlat'] = np.round(float(chan_df['sta_lat'].to_numpy()[0]),3)
        sacfile_dict['stlon'] = np.round(float(chan_df['sta_lon'].to_numpy()[0]),3)
        sacfile_dict['stel'] = np.round(float(chan_df['sta_elev'].to_numpy()[0]),3)
        sacfile_dict['stdp'] = np.round(float(chan_df['sta_depth'].to_numpy()[0]),3)

    # Add resp file:
    resp_file = input_sacfile.split('processed/')[0]+'resp/STXML.'+input_sacfile.split('processed/')[1]
    if os.path.isfile(resp_file):
        sacfile_dict['resp_file'] = resp_file
    else:
        print(f'Failed searching for resp file for: {network_name}, {station_name}, {channel_name_in}')
        print('Return empty dict....')
        return {}

    return sacfile_dict


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_input_dicts(input_directory, functions, do_parallel, jobs, channel):
    '''
    Parallel processing needs a list of inputs for each independent iteration
    So put all input into a dict and create a list of dicts

    Parameters
    ----------
    input_directory: str
        input_directory name where all the data is.
    functions : list
        list of functions to apply to inputs
    do_parallel : bool
    jobs :  int
        number of paralle workers
    channel : str
        required output channel of seismograms
        
    Returns
    ----------
    input_dicts : list (of dictionaries as below)
        input_dict : dict : specific to each file
            input dictionary containing:
                file : file name to process.
                sacfile_dict : dict of parameters describing file
                output_directory : formatted output dir
                cmt_df : pd df of CMT catalog
                stations_df : STATION file from Specfem execution as pd df
                functions : list of functions to execute on sacfile
                do_parallel : boolean
                jobs : number of jobs
                channel : desired output data channel
    '''

    # Get CMT solutions
    cmt_catalog = input_directory + '/../EVENTS-INFO/catalog.txt'
    if  os.path.isfile(cmt_catalog):
        cmt_df = read_CMT_CATALOG(cmt_catalog)
        # print(f'Got CMT dict from {cmt_catalog}')

        event_id = input_directory.split('/')[-1]
        cmt_df = cmt_df.query("event_id == @event_id")
        if len(cmt_df) != 1:
            # Some issue:
            print(f'**WARNING: cannot find unique event at: {str(cmt_catalog):s}\n')
            print('Exiting...')
            sys.exit()

    else:
        print(f'**WARNING: cannot find CMT CATALOG at: {str(cmt_catalog):s}\n')
        print('Exiting...')
        sys.exit()

    # Get station file dataframe
    station_filepath = input_directory + '/info/station_event'
    if  os.path.isfile(station_filepath):
        stations_df = read_STATION_file(station_filepath)
        # print(f'Got STATIONS df from {station_filepath}')
    else:
        print(f'**WARNING: cannot find STATIONS at: {str(station_filepath):s}\n')
        print('Exiting...')
        sys.exit()

    # Make the output formatted directory
    output_directory =  input_directory + '/py_formatted/' + cmt_df['ev_id'].to_list()[0]
    if not os.path.exists(output_directory):
        # print(f'Making formatted dir: {output_directory}')
        os.makedirs(output_directory)

    
    # Finds files, makes input dictionary of files, functions, inputs

    # files = sorted(glob.glob(input_directory+'/processed/1M.D36.*.?H?'))

    if channel == 'LH':
        files = sorted(glob.glob(input_directory+'/processed/??.*.?H?'))
    elif channel == 'BH':
        files = sorted(glob.glob(input_directory+'/processed/??.*.HH?') + glob.glob(input_directory+'/processed/??.*.BH?'))
    elif channel == 'HH':
        files = sorted(glob.glob(input_directory+'/processed/??.*.HH?'))
        
    num_files = len(files)
    input_dicts = []

    if num_files == 0:
        # No files found...
        print('----------')
        print('----------////   NO-FILES-IN: '+str(input_directory)+'    ////----------')
    else:

        print(f'----------')
        print(f'Processing list of {num_files} files...')
        print(f'----------')

        for file in files:

            sacfile_dict = get_sacfile_dict(file, cmt_df, stations_df, channel)
            
            if sacfile_dict: # I.e. is not empty
                input_dict = {}
                input_dict['file'] = file
                input_dict['sacfile_dict'] = sacfile_dict
                input_dict['output_directory'] = output_directory
                input_dict['cmt_df'] = cmt_df
                input_dict['stations_df'] =  stations_df
                input_dict['functions'] = functions
                input_dict['do_parallel'] = do_parallel
                input_dict['jobs'] = jobs
                input_dict['channel'] = channel

                input_dicts.append(input_dict)

        print('----------\n')

    return input_dicts


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_input_dicts_rotate(input_directory, functions, do_parallel, jobs, channel):
    '''
    Parallel processing needs a list of inputs for each independent iteration
    So put all input into a dict and create a list of dicts

    Parameters
    ----------
    input_directory: str
        input_directory name where all the data is.
    functions : list
        list of functions to apply to inputs
    do_parallel : bool
    jobs :  int
        number of paralle workers
    channel : str
        required output channel of seismograms
        
    Returns
    ----------
    input_dicts : list (of dictionaries as below)
        input_dict : dict : specific to each file
            input dictionary containing:
                file_v : file name of vertical component to process.
                horizontals : list of horizontal file names to process
                output_directory : formatted output dir
                functions : list of functions to execute on sacfile
                do_parallel : boolean
                jobs : number of jobs
                channel : desired output data channel
    '''

    # Get CMT solutions
    cmt_catalog = input_directory + '/../EVENTS-INFO/catalog.txt'
    if  os.path.isfile(cmt_catalog):
        cmt_df = read_CMT_CATALOG(cmt_catalog)

        event_id = input_directory.split('/')[-1]
        cmt_df = cmt_df.query("event_id == @event_id")
        if len(cmt_df) != 1:
            # Some issue:
            print(f'**WARNING: cannot find unique event at: {str(cmt_catalog):s}\n')
            print('Exiting...')
            sys.exit()

    # Make the output formatted directory
    output_directory =  input_directory + '/py_formatted/' + cmt_df['ev_id'].to_list()[0]
    if not os.path.exists(output_directory):
        # Big issue - this should exit
        print(f'There is no output directory with formatted data... Exit..')
        sys.exit()

    
    # Finds files, makes input dictionary of files, functions, inputs
        
    vert_files = sorted(glob.glob(output_directory+'/??????????????_??_*.' + str(channel) + 'Z'))

    num_files = len(vert_files)
    input_dicts = []

    if num_files == 0:
        # No files found...
        print('----------')
        print('----------////   NO-VERT-FILES-IN: '+str(output_directory)+'    ////----------')
    else:

        print(f'----------')
        print(f'Processing list of {num_files} files...')
        print(f'----------')

        for file in vert_files:

            horizontals = sorted(glob.glob(file[:-1] + '*'))
            horizontals.remove(file)

            input_dict = {}
            input_dict['file_v'] = file
            input_dict['horizontals'] = horizontals
            input_dict['output_directory'] = output_directory
            input_dict['functions'] = functions
            input_dict['do_parallel'] = do_parallel
            input_dict['jobs'] = jobs
            input_dict['channel'] = channel
            
            input_dicts.append(input_dict)

        print('----------\n')

    return input_dicts


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def apply_serial(main_function, input_directory, functions, do_parallel, jobs, channel):
    '''
    Takes a list of functions, reads all input files, applies each function to all files
    IN SERIES...

    Parameters
    ----------
    main_function : list
        Length one list containing main function.
    input_directory: str
        input_directory name where all the data is.
    functions : list
        list of functions to apply to inputs
    do_parallel : bool
    jobs :  int
        number of paralle workers
    channel : str
        required output channel of seismograms

    Returns
    ----------
    Nothing : None
        Functions take care of returns/writing to file  
           
    The functions must return all input variables
    Functions can write to logfile & outfile.
    '''

    # get input dicts
    if main_function[0].__name__ == 'process_one_file':
        input_dicts = get_input_dicts(input_directory, functions, do_parallel, jobs, channel)

    elif main_function[0].__name__ == 'process_one_station':
        input_dicts = get_input_dicts_rotate(input_directory, functions, do_parallel, jobs, channel)

    if input_dicts:
        # Serial processing
        for input_dict in input_dicts:
    
            for main_f in main_function:
                res = main_f(input_dict)
    
        return
    else:
        pass

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def apply_parallel(main_function, input_directory, functions, do_parallel, jobs, channel):
    '''
    Takes a list of functions, reads all input files, applies each function to all of the files
    IN PARALLEL...

    Parameters
    ----------
    main_function : list
        Length one list containing main function.
    input_directory: str
        input_directory name where all the data is.
    functions : list
        list of functions to apply to inputs
    do_parallel : bool
    jobs :  int
        number of paralle workers
    channel : str
        required output channel of seismograms

    Returns
    ----------
    Nothing : None
        Functions take care of returns/writing to file  
    
    The functions must return all input variables
    Functions can write to logfile & outfile.
    '''

    # get input dicts
    if main_function[0].__name__ == 'process_one_file':
        input_dicts = get_input_dicts(input_directory, functions, do_parallel, jobs, channel)

    elif main_function[0].__name__ == 'process_one_station':
        input_dicts = get_input_dicts_rotate(input_directory, functions, do_parallel, jobs, channel)

    if input_dicts:
        # Parallel processing
        for main_f in main_function:

            with concurrent.futures.ProcessPoolExecutor(max_workers = jobs) as executor:
                results = executor.map(main_f, input_dicts)

            return
    else:
        pass

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def execute(main_function, input_directory, functions, do_parallel, jobs, channel):
    '''
    Takes main function inputs and passes to serial or parallel executer, see below.
    '''
    if do_parallel:
        # Execute in parallel
        apply_parallel(main_function, input_directory, functions, do_parallel, jobs, channel)
    else:
        # Execute in series
        apply_serial(main_function, input_directory, functions, do_parallel, jobs, channel)
    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def process_one_file(input_dict):
    '''
    Processes one file reading the file, applying the list functions to them

    Parameters
    ----------
    input_dicts : list (of dictionaries as below)
        input_dict : dict : specific to each file
            input dictionary containing:
                file : file name to process.
                sacfile_dict : dict of parameters describing file
                output_directory : formatted output dir
                cmt_df : pd df of CMT CATALOG
                stations_df : STATION file from Specfem execution as pd df
                functions : list of functions to execute on sacfile
                do_parallel : boolean
                jobs : number of jobs
                channel : desired output data channel

    The functions in the list must take the following inputs: 

    input_dict : dict
    file : file location string
    seis : obspy stream object onwhich operations are performed.
    fail : int - fail flag = 1 if function fails, else 0.

    Returns
    ----------
    Nothing: None
    '''

    file = input_dict['file']
    functions = input_dict['functions']

    # read sac file
    try:
        seis = read(file,'sac', debug_headers=True)
        fail = 0
        print(f'Read seis {file}')
    except:
        #  FAILED_LOADING...
        fail = 1
        return
    
    if len(functions) ==  0: 
        # NO_FUNCTIONS_TO_APPLY.....
        return

    else:
        # apply functions, only executed when fail = 0
        for function in functions:
            input_dict, file, seis, fail = function(input_dict, file, seis, fail)

    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def process_one_station(input_dict):
    '''
    Processes one station reading the verticals and horizontals and applying the list functions to them

    Parameters
    ----------
    input_dicts : list (of dictionaries as below)
        input_dict : dict : specific to each file
            input dictionary containing:
                file_v : file name of vertical component to process.
                horizontals : list of horizontal file names to process
                output_directory : formatted output dir
                functions : list of functions to execute on sacfile
                do_parallel : boolean
                jobs : number of jobs
                channel : desired output data channel

    The functions in the list must take the following inputs: 

    input_dict : dict
    file_vertical, file_east, file_north : file location strings
    vert_component, east_component, north_component : obspy stream object onwhich operations are performed.
    fail : int - fail flag = 1 if function fails, else 0.

    Returns
    ----------
    Nothing: None
    '''

    file_vertical = input_dict['file_v']
    functions = input_dict['functions']
    horizontals = input_dict['horizontals']
    
    if len(horizontals) >= 2: 
        # read sac files
        try:
            # Read components
            vert_component = read(file_vertical,'sac')

            if str(file_vertical[:-1]) + 'E' in horizontals:
                file_east = str(file_vertical[:-1]) + 'E'
                east_component = read(file_east, format='sac')
            else:
                file_east = str(file_vertical[:-1]) + '1'
                east_component = read(file_east, format='sac')

            if str(file_vertical[:-1]) + 'N' in horizontals:
                file_north = str(file_vertical[:-1]) + 'N'
                north_component = read(file_north, format='sac')
            else:
                file_north = str(file_vertical[:-1]) + '2'
                north_component = read(file_north, format='sac')

            fail = 0
            print(f'Read seis {file_vertical}, {file_east} and {file_north}...')
        except:
            #  FAILED_LOADING...
            fail = 1
            return

        if len(functions) ==  0: 
            # NO_FUNCTIONS_TO_APPLY.....
            return

        else:
            # apply functions, only executed when fail = 0
            for function in functions:
                input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail = function(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail)
    else:
        return

    return



# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def check_length(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
    '''
    Check trace length

    Inputs and returns:
        input_dict : dict
        file_vertical, file_east, file_north : seismogram file names
        vert_component, east_component, north_component : obspy streams
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        
        # Check trace length and sample rate
        try:
            if vert_component[0].stats.npts % 10 == 1 & east_component[0].stats.npts % 10 == 1 & north_component[0].stats.npts % 10 == 1:
                # print('Caught XXXXX1 sample points')
                vert_component[0].data=vert_component[0].data[:-1]
                east_component[0].data=east_component[0].data[:-1]
                north_component[0].data=north_component[0].data[:-1]
            
            # print(f'Passed length QC....')
            fail = 0
        except:
            fail = 1
            print(f'Failed length QC....')
        
    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail



# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def check_sample_rate(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
    '''
    Check sample rate

    Inputs and returns:
        input_dict : dict
        file_vertical, file_east, file_north : seismogram file names
        vert_component, east_component, north_component : obspy streams
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:

        try:
            if input_dict['channel'] == 'LH':
                if vert_component[0].stats.sampling_rate < 1.0 or east_component[0].stats.sampling_rate < 1.0 or north_component[0].stats.sampling_rate < 1.0:
                    vert_component.resample(1)
                    east_component.resample(1)
                    north_component.resample(1)
                    # print('Some issue with sample rate <1.0 for channel: ', input_dict['channel'])
            if input_dict['channel'] == 'BH':
                if vert_component[0].stats.sampling_rate < 0.02 or east_component[0].stats.sampling_rate < 0.02 or north_component[0].stats.sampling_rate < 0.02:
                    vert_component.resample(0.02)
                    east_component.resample(0.02)
                    north_component.resample(0.02)
                    # print('Some issue with sample rate <0.02 for channel: ', input_dict['channel'])
            if input_dict['channel'] == 'HH':
                if vert_component[0].stats.sampling_rate < 0.005 or east_component[0].stats.sampling_rate < 0.005 or north_component[0].stats.sampling_rate < 0.005:
                    vert_component.resample(0.005)
                    east_component.resample(0.005)
                    north_component.resample(0.005)
                    # print('Some issue with sample rate <0.005 for channel: ', input_dict['channel'])
            
            # print(f'Passed sample rate QC....')
            fail = 0     
        except:
            fail = 1
            print(f'Failed sample rate QC....')
        
    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail




# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def trim_components(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
    '''
    Check Component length

    Inputs and returns:
        input_dict : dict
        file_vertical, file_east, file_north : seismogram file names
        vert_component, east_component, north_component : obspy streams
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        # Trim files:

        # find streams with min/max start/end times
        start_cut = UTCDateTime(vert_component[0].stats.starttime)
        end_cut = UTCDateTime(vert_component[0].stats.endtime)

        if UTCDateTime(east_component[0].stats.starttime) > start_cut : 
            start_cut = UTCDateTime(east_component[0].stats.starttime)

        if UTCDateTime(north_component[0].stats.starttime) > start_cut : 
            start_cut = UTCDateTime(north_component[0].stats.starttime)

        if UTCDateTime(east_component[0].stats.endtime) > end_cut : 
            end_cut = UTCDateTime(east_component[0].stats.endtime)

        if UTCDateTime(north_component[0].stats.endtime) > end_cut : 
            end_cut = UTCDateTime(north_component[0].stats.endtime)

        # Cut components to the same length
        # If this cant happen now its due to gappy waveforms of varying lengths.
        try:
            vert_component.trim(starttime=start_cut, endtime=end_cut)
            east_component.trim(starttime=start_cut, endtime=end_cut)
            north_component.trim(starttime=start_cut, endtime=end_cut)
            # print('Trimmed to :' + str(vert_component[0].stats['npts']))

            # print(f'Passed Trimming....')
            fail = 0
        except:
            print('Trace error in seisN/seisE/seisZ')
            print('Station: ' + str(vertical[:-1]) + ', Event: ' + str(s) + '\n')
            print('Issues with starttime/endtime/npts' + '\n')
            print(f'Failed Trimming....')
            fail = 1

    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail



# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def rotate_components(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
    '''
   Rotate components

    Inputs and returns:
        input_dict : dict
        file_vertical, file_east, file_north : seismogram file names
        vert_component, east_component, north_component : obspy streams
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        
        # Check trace length and sample rate
        try:

            loc_Z=vert_component[0].stats.location
            seisZ = vert_component.select(channel='*HZ',location=loc_Z)
            seisN_channel  = north_component[0].stats['channel']
            seisE_channel  = east_component[0].stats['channel']
            seisN = north_component.select(channel=seisN_channel,location=loc_Z)
            seisE = east_component.select(channel=seisE_channel,location=loc_Z)


            ori_seisN_temp=seisN[0].stats['sac']['cmpaz'] # .stats['orientation']
            ori_seisE_temp=seisE[0].stats['sac']['cmpaz'] # .stats['orientation']
            dip_seisE_temp=seisE[0].stats['sac']['cmpinc'] # .stats['dip']
            dip_seisN_temp=seisN[0].stats['sac']['cmpinc'] # .stats['dip']
            ori_seisZ_temp=seisZ[0].stats['sac']['cmpaz'] # .stats['orientation']
            dip_seisZ_temp=seisZ[0].stats['sac']['cmpinc'] # .stats['dip']

            # make sure the orientations are different and dip is not vertical
            if ori_seisN_temp == ori_seisE_temp:            
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisE_temp == -90:
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisN_temp == -90:
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisZ_temp != -90.0:
                print('Trace error in seisZ')
                print('Failed on Z component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


            EVLA = seisZ[0].stats['sac']['evla']
            EVLO = seisZ[0].stats['sac']['evlo']
            EVDP = seisZ[0].stats['sac']['evdp']
            STLA = seisZ[0].stats['sac']['stla']
            STLO = seisZ[0].stats['sac']['stlo']
            STEL = seisZ[0].stats['sac']['stel']
            STATION = seisZ[0].stats['station']
            NETWORK = seisZ[0].stats['network']

            DISTM, AZ, BAZ = obspy.geodetics.base.gps2dist_azimuth(EVLA, EVLO, STLA, STLO)
            DISTDG = DISTM / (6371.e3 * np.pi / 180.)

            # Check North and East exist
            if len(seisN) == 1 and len(seisE) == 1:
                
                # Make sure components are of the same length
                if len(seisN[0].data) > len(seisE[0].data):
                    seisN[0].data = seisN[0].data[:-1]
                if len(seisN[0].data) < len(seisE[0].data):
                    seisE[0].data = seisE[0].data[:-1]
                if len(seisZ[0].data) > len(seisE[0].data) and len(seisZ[0].data) > len(seisN[0].data):
                    seisZ[0].data = seisZ[0].data[:-1]  

                # Final check of length - if not remove.
                if len(seisZ[0].data) != len(seisE[0].data) or len(seisZ[0].data) != len(seisN[0].data):
                    print('Trace length error in seisZ/seisN/seisE')
                    print('Length seisZ[0].data: '+str(len(seisZ[0].data)))
                    print('Length seisN[0].data: '+str(len(seisN[0].data)))
                    print('Length seisE[0].data: '+str(len(seisE[0].data)))

                    fail = 1
                    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail

                # Make sure North component is North and East component is East
                if 45 < ori_seisN_temp < 135 or 225 < ori_seisN_temp < 315:
                    seisN1 = seisE
                    seisE1 = seisN
                    seisN = seisN1
                    seisE = seisE1

                # print('Orientation N: ' + str(ori_seisN_temp))
                # print('Orientation E: ' + str(ori_seisE_temp))

                # print('Rotating to true N/E')
                [seisZtmp, seisNtmp, seisEtmp] = obspy.signal.rotate.rotate2zne(
                    seisZ[0].data, ori_seisZ_temp, dip_seisZ_temp,
                    seisN[0].data, ori_seisN_temp, dip_seisN_temp,
                    seisE[0].data, ori_seisE_temp, dip_seisE_temp)
                seisN = seisN[0].copy() 
                seisN.stats['channel'] = str(input_dict['channel']) + 'N'
                seisN.data = seisNtmp
                seisE = seisE[0].copy()
                seisE.stats['channel'] = str(input_dict['channel']) + 'E'
                seisE.data = seisEtmp    

                # should be setting the new orientation and dip to 0,90,-90 etc

                # rotate components to from North and East to Radial and Transverse
            
                # print('Rotating to RT using BAZ: ',BAZ)
                [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(
                    seisN.data, seisE.data, BAZ)
                seisR = seisN.copy()
                seisT = seisN.copy()

                seisR.stats['channel'] = str(input_dict['channel']) + 'R'
                seisT.stats['channel'] = str(input_dict['channel']) + 'T'
            
                seisR.data = seisRtmp
                seisT.data = seisTtmp


                ########## Sort out SAC headers #################

                # Radial
                # seisR.stats['channel'] = str(input_dict['channel']) + 'R'
                seisR.stats['sac']['cmpaz']=AZ
                seisR.stats['sac']['az']=AZ
                seisR.stats['sac']['baz']=BAZ
                seisR.stats['sac']['gcarc']=DISTDG
                seisR.stats['sac']['kcmpnm']=str(input_dict['channel'])+'R'

                # Tangential
                # seisT.stats['channel'] = str(input_dict['channel']) + 'T'
                seisT.stats['sac']['cmpaz']=AZ+90.0
                seisT.stats['sac']['az']=AZ
                seisT.stats['sac']['baz']=BAZ
                seisT.stats['sac']['gcarc']=DISTDG
                seisT.stats['sac']['kcmpnm']=str(input_dict['channel'])+'T'

                # East
                # seisE.stats['channel']=str(channel)+'E'
                seisE.stats['sac']['cmpaz']=90.0
                seisE.stats['sac']['az']=AZ
                seisE.stats['sac']['baz']=BAZ
                seisE.stats['sac']['gcarc']=DISTDG
                seisE.stats['sac']['kcmpnm']=str(input_dict['channel'])+'E'

                # North
                # seisN.stats['channel']=str(channel)+'N'
                seisN.stats['sac']['cmpaz']=0.0
                seisN.stats['sac']['az']=AZ
                seisN.stats['sac']['baz']=BAZ
                seisN.stats['sac']['gcarc']=DISTDG
                seisN.stats['sac']['kcmpnm']=str(input_dict['channel'])+'N'

                # Vertical
                seisZ[0].stats['sac']['az']=AZ
                seisZ[0].stats['sac']['baz']=BAZ
                seisZ[0].stats['sac']['gcarc']=DISTDG


                filename_E = str(file_vertical[:-1]) + 'E'
                filename_N = str(file_vertical[:-1]) + 'N'
                filename_V = str(file_vertical[:-1]) + 'Z'
                filename_R = str(file_vertical[:-1]) + 'R'
                filename_T = str(file_vertical[:-1]) + 'T'

                # Write over sac files as we now have everything sorted :)
                seisE.write(filename_E, 'sac')
                seisN.write(filename_N, 'sac')
                seisZ.write(filename_V, 'sac')
                seisR.write(filename_R, 'sac')
                seisT.write(filename_T, 'sac')

            
            # print(f'Passed Rotation....')
            fail = 0
        except:
            print(f'Failed Rotation....')
            fail = 1

    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail



def cleanup_files(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
    '''
   Clean up unused files

    Inputs and returns:
        input_dict : dict
        file_vertical, file_east, file_north : seismogram file names
        vert_component, east_component, north_component : obspy streams
        fail : 1/0
    '''

    # Remove all that are unnecessary - so ignore the fail flag
    try:
                    
        # Make a directory for unused files...
        Not_useable = input_dict['output_directory'] + '/../Not_useable/'
        if not os.path.exists(Not_useable):
            os.makedirs(Not_useable)

        stat = file_vertical[:-3]

        # Keep ZRT components
        to_delete = sorted(glob.glob(stat+'*'))

        if file_vertical[:-1]+'Z' in to_delete:
            to_delete.remove(file_vertical[:-1]+'Z')
        if file_vertical[:-1]+'R' in to_delete:
            to_delete.remove(file_vertical[:-1]+'R')
        if file_vertical[:-1]+'T' in to_delete:
            to_delete.remove(file_vertical[:-1]+'T')

        if len(to_delete) > 1:
            for file in to_delete:
                try:
                    # print('Moving :' +str(file)+' -> '+str(Not_useable))
                    shutil.move(file, Not_useable)
                except:
                    pass

            # print(f'Passed Cleanup....')
            fail = 0
    except:
        print(f'Failed Cleanup of {stat}* ....')
        fail = 1

    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def update_headers_seis(input_dict, file, seis, fail):
    '''
    Updates all headers of the seismogram

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''

    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        ###
        # print(f'Update headers seis: {file}')

        seis[0].stats['sac']['o'] = UTCDateTime(input_dict['cmt_df']['datetime'].to_list()[0]) - UTCDateTime(seis[0].stats['starttime'])

        seis[0].stats['sac']['evla'] =  float(input_dict['cmt_df']['latitude'].to_numpy()[0])
        seis[0].stats['sac']['evlo'] =  float(input_dict['cmt_df']['longitude'].to_numpy()[0])
        seis[0].stats['sac']['evdp'] =  float(input_dict['cmt_df']['depth'].to_numpy()[0])
        seis[0].stats['sac']['mag'] =  float(input_dict['cmt_df']['magnitude'].to_numpy()[0])
        seis[0].stats['sac']['kevnm'] =  str(input_dict['cmt_df']['event_id'].to_list()[0])
        seis[0].stats['sac']['lovrok'] =  1
        seis[0].stats['sac']['lcalda'] =  1

        seis[0].stats['sac']['stla'] = input_dict['sacfile_dict']['stlat']
        seis[0].stats['sac']['stlo'] = input_dict['sacfile_dict']['stlon']
        seis[0].stats['sac']['stel'] = input_dict['sacfile_dict']['stel']
        seis[0].stats['sac']['stdp'] = input_dict['sacfile_dict']['stdp']
        seis[0].stats['sac']['kstnm'] = input_dict['sacfile_dict']['kstnm']
        seis[0].stats['sac']['knetwk'] = input_dict['sacfile_dict']['ntwrk']
        seis[0].stats['sac']['kcmpnm'] = input_dict['sacfile_dict']['kcmpnm']
        seis[0].stats['sac']['khole'] = input_dict['sacfile_dict']['kstloc']

        seis[0].stats['network'] = input_dict['sacfile_dict']['ntwrk']
        seis[0].stats['station'] = input_dict['sacfile_dict']['kstnm']
        seis[0].stats['location'] = input_dict['sacfile_dict']['kstloc']
        seis[0].stats['channel'] = str(input_dict['sacfile_dict']['kcmpnm'])


    return input_dict, file, seis, fail


def rem_inst_resp(input_dict, file, seis, fail):
    '''
    Removes instrument response from seismogram
    using appropriate pre-filter and station xml file
    outputs to displacement

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''

    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:

        if input_dict['sacfile_dict']['chan_in'] == "LH":
            # For LH data:
            pre_filt = [0.01, 0.02, 0.1, 0.4]

        if input_dict['sacfile_dict']['chan_in'] == "BH":
            # For BH data:
            pre_filt = [0.01, 0.02, 10, 20]

        if input_dict['sacfile_dict']['chan_in'] == "HH":
            # For HH data:
            pre_filt = [0.01, 0.02, 40, 60]

        inv_st = read_inventory(input_dict['sacfile_dict']['resp_file'], format="stationxml")

        #### OUTPUT TO DISPLACEMENT ###
        try:
            seis.remove_response(inventory=inv_st, pre_filt=pre_filt, output="DISP", water_level=None)  
        except:
            print(f'Failed to remove response for: {file}')
            fail = 1

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def remove_trend_taper(input_dict, file, seis, fail):
    '''
    Removes trend and tapers seismogram

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''

    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        ###
        # print(f'Remove trend seis: {file}')
        
        # merge gappy waveforms and overlap
        try:
            seis.merge(fill_value='interpolate')
        except:
            seis.resample(40)
            seis.merge(fill_value='interpolate')

        seis.detrend()
        seis.taper(max_percentage=0.015, type='cosine')

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def filter_seis(input_dict, file, seis, fail):
    '''
    Filter seismogram

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''

    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        ###
        # print(f'Filter seis: {file}')

        if input_dict['channel'] == 'LH':
            fmin = 0.01
            fmax = 0.4
        else:
            fmin = 0.01
            fmax = 20
        seis.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def interpolate_seis(input_dict, file, seis, fail):
    '''
    Interpolate seismogram

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''

    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        # print(f'Interpolate seis: {file}')
        ###
        if input_dict['channel'] == 'LH':
            delta = 1
        else:
            delta = 0.05
        seis.resample(1 / delta)

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def save_seis(input_dict, file, seis, fail):
    '''
    Saves processed sac file

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        out_filename = input_dict['output_directory'] + '/' + input_dict['sacfile_dict']['output_sacfile']
        
        try:
            # print(f'Function Save seismogram: {file}')

            print(f'Saving {out_filename}')
            seis.write(out_filename, format='sac')
            fail = 0
        except:
            fail = 1
            print(f'FAILED to save: {out_filename}')

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def check_no_file_exists(input_dict, file, seis, fail):
    '''
    Checks if adequate file exists, else we skip to save time...

    Inputs and returns:
        input_dict : dict
        file : seismogram file name
        seis : obspy stream
        fail : 1/0
    '''
    if fail:
        ###
        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
        fail = 1
    else:
        out_filename = input_dict['output_directory'] + '/' + input_dict['sacfile_dict']['output_sacfile']
        
        if not os.path.isfile(out_filename):
            # print(f'Lets process the seismogram: {file}')
            fail = 0
        else:
            print(f'Suitable seismogram exists at: {out_filename}')
            # print(f'Return Fail = 1 to skip processing.')
            fail = 1

    return input_dict, file, seis, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():

    start_time = time.time()

    if len(sys.argv) != 4:
        print("Please specify Start date and years")
        print("USAGE:   python parallel_proc_obspyDMT_seis.py <data_loc> <jobs> <channel>")
        print("EXAMPLE: python parallel_proc_obspyDMT_seis.py ./20081230_194956.a 96 LH")
        sys.exit()
    else:
        data_loc = str(sys.argv[1])
        if not os.path.isdir(data_loc):
            print("Directory given is not present.... exiting....")
            sys.exit()

        jobs = int(sys.argv[2])
        channel = str(sys.argv[3])
        print(f"Processing {data_loc}, using {jobs} processes, assuming data channel: {channel}")


    if jobs > 0  and jobs <= 96 :
        do_parallel = True
    else:
        do_parallel = False

    # Format the data
    main_function = [process_one_file]

    function_list = [check_no_file_exists, rem_inst_resp, update_headers_seis, remove_trend_taper, filter_seis, interpolate_seis, save_seis]
    # Executeabove functions to produce formatted data
    execute(main_function, data_loc, function_list, do_parallel, jobs, channel)

    # rotate_components
    main_function = [process_one_station]

    function_list = [check_length, check_sample_rate, trim_components, rotate_components, cleanup_files]
    # Execute functions to Rotate the data to ZRT
    execute(main_function, data_loc, function_list, do_parallel, jobs, channel)

    # Sometimes there remains some N & E files, since they have no associated vertical component.
    del_list = glob.glob(data_loc + '/py_formatted/??????????????/??????????????_??_*E') + glob.glob(data_loc + '/py_formatted/??????????????/??????????????_??_*N')
    # print(del_list)
    for file in del_list:
        os.remove(file)

    
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()