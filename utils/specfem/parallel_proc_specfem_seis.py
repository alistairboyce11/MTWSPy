# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:46:02 2024
@author: alistairboyce
if this fails: >> conda activate env3.12
"""

import time, sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace, UTCDateTime, read, Stream
from scipy.signal import convolve
import pandas as pd
import concurrent.futures

class ParallelProcSpecfemSeis:

    """
    Class to handle serial/parallel post processing of SPECFEM3D simulations
    to format the data ready for MTWSPy travel time picking
    
    First we obtain dictionaries of CMT solutions and station details.
    Then apply a common series of functions to all data files produced by 
    SPECFEM and save output to a formatted directory.

    """

    def __init__(self):
        pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def generate_triangular_source_time_function(self, half_duration, dt):
        """
        return triangular source time function using 
        half duration and dt

        :param half_duration: half_duration of source time function
        :type half_duration: float
        :param dt: sample interval (s)
        :type dt: float
        :return stf: stf normalised to unit area
        :type stf: np.array
        """
        t = np.arange(0, half_duration * 2, dt) 
        tri_wave = np.maximum(0, 1 - np.abs(t - half_duration) / half_duration)
        return tri_wave / np.sum(tri_wave)  # Normalize to have unit area


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def read_CMTSOLUTION(self, file):
        """
        Read CMTSOLUTION file into dict

        :param file: file name of CMT file to be read
        :type file: str
        :return cmt_dict: dictionary of CMT parameters
        :type cmt_dict: dict

        Returns
        ----------
        cmt_dict : dict 
            Original timing params: orig_time,orig_year,orig_month,orig_day,orig_hour,orig_min,orig_sec
            Shifted timing params: s_time,s_year,s_month,s_day,s_jday,s_hour,s_min,s_sec,s_msec
            ev_id : formatted event id string
            evname : GCMT event name
            tshift : time shift
            hdur : stf half duration
            evlat : event latitude
            evlon : event longitude
            evdep : event depth
            evmag : event magnitude
        """
        cmt_dict = {}

        f=open(file,'r')
        lines=f.readlines()
        f.close()

        orig_year = int(lines[0].split()[1])
        orig_month = int(lines[0].split()[2])
        orig_day = int(lines[0].split()[3])
        orig_hour = int(lines[0].split()[4])
        orig_min = int(lines[0].split()[5])
        orig_sec = float(lines[0].split()[6])
        orig_evmag = float(lines[0].split()[11])

        evname = str(lines[1].split()[2])
        tshift = float(lines[2].split()[2])
        hdur = float(lines[3].split()[2])
        evlat = float(lines[4].split()[1])
        evlon = float(lines[5].split()[1])
        evdep = float(lines[6].split()[1])

        time = UTCDateTime( str(orig_year) + ',' + str(orig_month) + ',' + str(orig_day) + ',' + str(orig_hour) + ',' + str(orig_min) + ',' + str(orig_sec) )

        time_shifted = UTCDateTime(time) + tshift
        
        cmt_dict['orig_time'] = time
        cmt_dict['orig_year'] = orig_year
        cmt_dict['orig_month'] = orig_month
        cmt_dict['orig_day'] = orig_day
        cmt_dict['orig_hour'] = orig_hour
        cmt_dict['orig_min'] = orig_min
        cmt_dict['orig_sec'] = orig_sec

        cmt_dict['s_time'] = time_shifted
        cmt_dict['s_year'] = time_shifted.year
        cmt_dict['s_month'] = time_shifted.month
        cmt_dict['s_day'] = time_shifted.day
        cmt_dict['s_jday'] = time_shifted.julday
        cmt_dict['s_hour'] = time_shifted.hour
        cmt_dict['s_min'] = time_shifted.minute
        cmt_dict['s_sec'] = time_shifted.second
        cmt_dict['s_msec'] = time_shifted.microsecond

        # Guard against rounding error in time shift gcmt catalog (1 s.f.) versus CMT solution (2.s.f.)
        # It affects the ev_id on ~5% of events.
        if UTCDateTime(UTCDateTime(time) + tshift).second != UTCDateTime(UTCDateTime(time) + np.round(tshift,1)).second:
            time_shifted = UTCDateTime(UTCDateTime(time) + np.round(tshift,1))

        cmt_dict['ev_id'] = '{0:4d}{1:02d}{2:02d}{3:02d}{4:02d}{5:02d}'.format(time_shifted.year, time_shifted.month, time_shifted.day, time_shifted.hour, time_shifted.minute, time_shifted.second)
        cmt_dict['evname'] = evname
        cmt_dict['tshift'] = tshift
        cmt_dict['hdur'] = hdur
        cmt_dict['evlat'] = evlat
        cmt_dict['evlon'] = evlon
        cmt_dict['evdep'] = evdep
        cmt_dict['evmag'] = orig_evmag

        return cmt_dict


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def read_STATION_file(self, file):
        """
        Read STATION file into pandas dataframe

        :param file: file name of station file to be read
        :type file: str
        :return stations_df: dictionary of station parameters
        :type stations_df: pd.df
        """

        with open(file, 'r') as f:
            # Skip header lines
            try:
                # Define column headers
                column_headers = ['sta_name', 'net_name', 'sta_lat', 'sta_lon', 'sta_elev', 'sta_depth']

                # Read the data into a pandas DataFrame
                stations_df = pd.read_csv(f, delim_whitespace=True, header=None, names= column_headers)

                # Network Name NA is read as nan... correct.
                stations_df.replace(np.nan, "NA", inplace=True)
            except:
                #Return empy dataframe if not possible.
                stations_df = pd.DataFrame()
        return stations_df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_sacfile_dict(self, input_sacfile, cmt_dict, stations_df, channel):
        """
        Read input_sacfilename get parameter for processing.

        :param input_sacfile: filename of input sac file
        :type input_sacfile: str
        :param cmt_dict: CMTSOLUTION dataframe
        :type cmt_dict: dict
        :param stations_df : Stations dataframe
        :type stations_df: pd.df 
        :param channel : desired data channel after processing
        :type channel: str
                   
        :return sacfile_dict: dict of sac file details
        :type sacfile_dict: dict

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
        """

        sacfile_dict = {}
        sacfile_dict['input_sacfile'] = input_sacfile

        network_name = input_sacfile.split('/')[-1].split('.')[0]
        station_name = input_sacfile.split('/')[-1].split('.')[1]

        sacfile_dict['ntwrk'] = network_name 
        sacfile_dict['kstnm'] = station_name 
        sacfile_dict['chan_in'] = input_sacfile.split('/')[-1].split('.')[2][0:2]
        sacfile_dict['comp_in'] = input_sacfile.split('/')[-1].split('.')[2][2]
        sacfile_dict['kstloc'] = str('00')
        sacfile_dict['kcmpnm'] = str(channel) + sacfile_dict['comp_in']

        # $oyear1$omonth1$oday1$ohr1$omin1$osec1"_"$ntwrk"_"$kstnm"."$kstloc"."$kcmpnm
        sacfile_dict['output_sacfile'] = cmt_dict['ev_id'] + '_' + sacfile_dict['ntwrk'] + '_' + sacfile_dict['kstnm'] + '.' + sacfile_dict['kstloc'] + '.' + sacfile_dict['kcmpnm']

        net_df = stations_df.query("net_name == @network_name")
        sta_df = net_df.query("sta_name == @station_name")

        if len(sta_df) != 1: 
            print(f'Failed searching for : {network_name}, {station_name} in STATIONS file')
            print('Exiting....')
            sys.exit()
        else:
            sacfile_dict['stlat'] = sta_df['sta_lat'].to_numpy()[0]
            sacfile_dict['stlon'] = sta_df['sta_lon'].to_numpy()[0]
            sacfile_dict['stel'] = sta_df['sta_elev'].to_numpy()[0]
            sacfile_dict['stdp'] = sta_df['sta_depth'].to_numpy()[0]

        return sacfile_dict


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_input_dicts(self, input_directory, functions, do_parallel, jobs, channel):
        """
        Parallel processing needs a list of inputs for each independent iteration
        So put all input into a dict and create a list of dicts

        :param input_directory: input_directory name where the data is
        :type input_directory: str
        :param functions : list of functions to apply to inputs
        :type functions: list
        :param do_parallel: Execute in parallel or not
        :type do_parallel : bool
        :param jobs: number of parallel workers
        :type jobs: int
        :param channel: sesimogram channel used in processing
        :type channel: str

        :return input_dicts: list of dicts with info on files to be processed
        :type input_dicts: list
            
        Returns
        ----------
        input_dicts : list (of dictionaries as below)
            input_dict : dict : specific to each file
                input dictionary containing:
                    file : file name to process.
                    sacfile_dict : dict of parameters describing file
                    output_directory : formatted output dir
                    cmt_dict : dictionary of CMTSOLUTION parameters from Specfem execution
                    stations_df : STATION file from Specfem execution as pd df
                    functions : list of functions to execute on sacfile
                    do_parallel : boolean
                    jobs : number of jobs
                    channel : desired output data channel
        """

        # Get CMT solutions
        cmt_filepath = input_directory + '/CMTSOLUTION'
        if  os.path.isfile(cmt_filepath):
            cmt_dict = self.read_CMTSOLUTION(cmt_filepath)
            # print(f'Got CMT dict from {cmt_filepath}')

        else:
            print(f'**WARNING: cannot find CMTSOLUTION at: {str(cmt_filepath):s}\n')
            print('Exiting...')
            sys.exit()

        # Get station file dataframe
        station_filepath = input_directory + '/STATIONS'
        if  os.path.isfile(station_filepath):
            stations_df = self.read_STATION_file(station_filepath)
            # print(f'Got STATIONS df from {station_filepath}')

        else:
            print(f'**WARNING: cannot find STATIONS at: {str(station_filepath):s}\n')
            print('Exiting...')
            sys.exit()

        # Make the output formatted directory
        output_directory =  input_directory + '/py_formatted/' + cmt_dict['ev_id']
        if not os.path.exists(output_directory):

            # print(f'Making formatted dir: {output_directory}')
            os.makedirs(output_directory)

        
        # Finds files, makes input dictionary of files, functions, inputs
        files = sorted(glob.glob(input_directory+'/*.sem.sac'))
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

                sacfile_dict = self.get_sacfile_dict(file, cmt_dict, stations_df, channel)

                input_dict = {}
                input_dict['file'] = file
                input_dict['sacfile_dict'] = sacfile_dict
                input_dict['output_directory'] = output_directory
                input_dict['cmt_dict'] = cmt_dict
                input_dict['stations_df'] =  stations_df
                input_dict['functions'] = functions
                input_dict['do_parallel'] = do_parallel
                input_dict['jobs'] = jobs
                input_dict['channel'] = channel

                input_dicts.append(input_dict)

            print('----------\n')

        return input_dicts



    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def apply_serial(self, main_function, input_directory, functions, do_parallel, jobs, channel):
        """
        Takes a list of functions, reads all input files, applies each function to all files
        IN SERIES...

        :param main_function : list of the 1 main function applied to inputs
        :type main_function: list
        :param input_directory: input_directory name where the data is
        :type input_directory: str
        :param functions : list of functions to apply to inputs
        :type functions: list
        :param do_parallel: Execute in parallel or not
        :type do_parallel : bool
        :param jobs: number of parallel workers
        :type jobs: int
        :param channel: sesimogram channel used in processing
        :type channel: str
            
        The functions in the list must take the following inputs: 

        input_dict : dict
        file : file location string
        seis : obspy stream object read into memory onwhich operations are performed.
        fail : int - fail flag = 1 if function fails, else 0.
        
        The functions must return all input variables
        Functions can write to logfile & outfile.
        """

        # get input dicts
        input_dicts = self.get_input_dicts(input_directory, functions, do_parallel, jobs, channel)
        
        if input_dicts:
            # Serial processing
            for input_dict in input_dicts:
        
                for main_f in main_function:
                    res = main_f(input_dict)
        
            return
        else:
            pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def apply_parallel(self, main_function, input_directory, functions, do_parallel, jobs, channel):
        """
        Takes a list of functions, reads all input files, applies each function to all of the files
        IN PARALLEL...

        :param main_function : list of the 1 main function applied to inputs
        :type main_function: list
        :param input_directory: input_directory name where the data is
        :type input_directory: str
        :param functions : list of functions to apply to inputs
        :type functions: list
        :param do_parallel: Execute in parallel or not
        :type do_parallel : bool
        :param jobs: number of parallel workers
        :type jobs: int
        :param channel: sesimogram channel used in processing
        :type channel: str
            
        The functions in the list must take the following inputs: 

        input_dict : dict
        file : file location string
        seis : obspy stream object read into memory onwhich operations are performed.
        fail : int - fail flag = 1 if function fails, else 0.
        
        The functions must return all input variables
        Functions can write to logfile & outfile.
        """

        # get input dicts
        input_dicts = self.get_input_dicts(input_directory, functions, do_parallel, jobs, channel)

        if input_dicts:
            # Parallel processing
            for main_f in main_function:

                with concurrent.futures.ProcessPoolExecutor(max_workers = jobs) as executor:
                    results = executor.map(main_f, input_dicts)

                return
        else:
            pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def execute(self, main_function, input_directory, functions, do_parallel, 
                jobs, channel):
        """
        Takes main function inputs and passes to serial or parallel 
        executer, see below.

        :param main_function : list of the 1 main function applied to inputs
        :type main_function: list
        :param input_directory: input_directory name where the data is
        :type input_directory: str
        :param functions : list of functions to apply to inputs
        :type functions: list
        :param do_parallel: Execute in parallel or not
        :type do_parallel : bool
        :param jobs: number of parallel workers
        :type jobs: int
        :param channel: sesimogram channel used in processing
        :type channel: str
        """
        if do_parallel:
            # Execute in parallel
            self.apply_parallel(main_function, input_directory, functions, do_parallel, jobs, channel)
        else:
            # Execute in series
            self.apply_serial(main_function, input_directory, functions, do_parallel, jobs, channel)
        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def process_one_file(self, input_dict):
        """
        Processes one file reading the file, applying the list functions 
        to them

        :return input_dict: dictionary of info used in processing
        :type input_dict: dict

        Parameters
        ----------
        input_dict : dict : specific to each file
            input dictionary containing:
                file : file name to process.
                sacfile_dict : dict of parameters describing file
                output_directory : formatted output dir
                cmt_dict : dictionary of CMTSOLUTION parameters from Specfem execution
                stations_df : STATION file from Specfem execution as pd df
                functions : list of functions to execute on sacfile
                do_parallel : boolean
                jobs : number of jobs
                channel : desired output data channel
        """

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


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def update_headers_seis(self, input_dict, file, seis, fail):
        """
        Updates all headers of the seismogram

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            DUMMY = 1
        else:
            ###
            # print(f'Update headers seis: {file}')

            # seis[0].stats['starttime'] -= input_dict['cmt_dict']['tshift'] 
            # seis[0].stats['sac']['o'] = 0

            seis[0].stats['sac']['evla'] =  input_dict['cmt_dict']['evlat']
            seis[0].stats['sac']['evlo'] =  input_dict['cmt_dict']['evlon']
            seis[0].stats['sac']['evdp'] =  input_dict['cmt_dict']['evdep']
            seis[0].stats['sac']['mag'] =  input_dict['cmt_dict']['evmag']
            seis[0].stats['sac']['kevnm'] =  str(input_dict['cmt_dict']['evname'])
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


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def convolve_seis(self, input_dict, file, seis, fail):
        """
        Convolve seismogram with appropriate sourcetime function

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            DUMMY = 1
        else:
            ###
            # print(f'Convolve seis: {file}')

            dt = seis[0].stats['delta']
            hdur = input_dict['cmt_dict']['hdur']
            
            # Get triangular source time function of unit area
            stf = self.generate_triangular_source_time_function(hdur, dt)

            # Convolve data with stf
            convolved_data = convolve(seis[0].data, stf, mode='full', method='auto')
            
            # Center the convolution result
            center_idx = len(stf) // 2
            convolved_data = convolved_data[center_idx : center_idx + len(seis[0].data)]

            # Create a new ObsPy Stream and add the convolved trace, and return as seis
            convolved_trace = Trace(data = convolved_data, header = seis[0].stats)
            convolved_seis = Stream(traces = [convolved_trace])

        return input_dict, file, convolved_seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def remove_trend_taper(self, input_dict, file, seis, fail):
        """
        Removes trend and tapers seismogram

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            DUMMY = 1
        else:
            ###
            # print(f'Remove trend seis: {file}')
            seis.detrend()
            seis.taper(max_percentage=0.015, type='cosine')

        return input_dict, file, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def filter_seis(self, input_dict, file, seis, fail):
        """
        Filter seismogram

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            DUMMY = 1
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


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def interpolate_seis(self, input_dict, file, seis, fail):
        """
        Interpolate seismogram

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            DUMMY = 1
        else:
            # print(f'Interpolate seis: {file}')
            ###
            if input_dict['channel'] == 'LH':
                delta = 1
            else:
                delta = 0.05
            seis.resample(1 / delta)

        return input_dict, file, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def save_seis(self, input_dict, file, seis, fail):
        """
        Saves processed sac file

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param file: seismogram file name
        :type file: str
        :param seis: seismic data in python format in memory
        :type seis: obspy stream
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

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


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def main():

    parallel_proc_specfem_seis = ParallelProcSpecfemSeis()

    start_time = time.time()

    if len(sys.argv) != 4:
        print("Please specify Start date and years")
        print("USAGE:   python parallel_proc_specfem_seis.py <data_loc> <jobs> <channel>")
        print("EXAMPLE: python parallel_proc_specfem_seis.py ./OUTPUT_FILES_3D_seis_hd0 96 LH")
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

    main_function = [parallel_proc_specfem_seis.process_one_file]

    function_list = [parallel_proc_specfem_seis.update_headers_seis,
                    parallel_proc_specfem_seis.convolve_seis,
                    parallel_proc_specfem_seis.remove_trend_taper,
                    parallel_proc_specfem_seis.filter_seis,
                    parallel_proc_specfem_seis.remove_trend_taper,
                    parallel_proc_specfem_seis.interpolate_seis,
                    parallel_proc_specfem_seis.save_seis]

    parallel_proc_specfem_seis.execute(main_function, data_loc, function_list, do_parallel, jobs, channel)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()