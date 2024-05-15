import glob, sys, os
import yaml
from yaml.loader import SafeLoader
import pandas as pd
from datetime import datetime, timezone
import numpy as np
import concurrent.futures

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
#                                                                                     #
#                                     toolkit.py                                      #
#                                                                                     # 
# ---                                                                             --- # 

# A bunch of functions utilised throughout the code.

class Toolkit:
    def __init__(self):
        pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def check_files(self, loc, year, component):
        '''
        Check files exist for given obs and synth data location for
        specified year and component
        Otherwise quit the code...
        '''
        check_name = f'{loc}/{str(year)}*/{str(year)}*{component}'

        filelist = glob.glob(check_name)
        
        numfiles = len(filelist)

        if numfiles > 0:
            print(f'Found: {str(numfiles):s} files\n')
            print(f'At: {check_name:s}\n')
            print('Continuing...\n')
        else:
            print(f'ERROR....\n')
            print(f'Nothing found at: {check_name:s}\n')
            print(f'Retry...\n')
            sys.exit()

        return numfiles

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_params(self, params_filename):
        '''
        Get input parameters from yaml file
        '''
        with open(params_filename) as f:
            params = yaml.load(f, Loader = SafeLoader)
        # add the starttime of the code
        now = datetime.now(timezone.utc)
        params['code_start_time'] = now.strftime("%Y%m%d%H%M%S")

        cha_obs, cha_syn = f'{params['twin_obs_out_loc'][-2:]}{params['component']}', f'{params['twin_syn_out_loc'][-2:]}{params['component']}'
        params['mtf_outfilename'] = f'{str(params['mtf_prefix'])}_{str(params['year'])}.{cha_obs}-{cha_syn}'

        return params

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_event_id_table(self, params):
        '''
        Get event id table as pandas data frame
        '''
        cmt_outfile_name = f'{params['cmt_outfile']}_{str( params['year'])}.csv'

        tab = pd.read_csv(cmt_outfile_name)

        return tab

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_log_statement(self, event_id, filename):
        '''
        Return statement for logfile 
        '''
        filename_p = filename.split('/')[-1]
        statement = f'{event_id}, {filename_p}'

        return statement

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def subvec(self, t, x, interval):
        '''
        Return indices and values within a given interval [t1,t2] for timeseries x(t) 
        '''
        t_new = (t >= interval[0]) & (t <= interval[1])
        x_new = x[t_new]

        return t_new, x_new

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def flatten_concatenation(self, matrix):
        '''
        Return a flattened list
        '''
        flat_list = []
        for row in matrix:
            flat_list  +=  row
        return flat_list

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def print_log(self, params, logfile, statement):
        '''
        Print statement to terminal (depending on input params) and open logfile 
        '''
        if params['verbose']:
            print(statement)
        
        try:
            logfile.write(f'{statement}\n')
        except:
            # Probably something that cannot be written just printed.
            pass

        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_io_log_file(self, params_in):
        '''
        Return an open input/output log file in log_loc/'code_start_time'/steps_results.log
        '''

        lf_loc = f'{params_in['home']}/{params_in['log_loc']}/{str(params_in['code_start_time'])}'
        if not os.path.exists(lf_loc):
            os.makedirs(lf_loc, exist_ok=True)

        lf_name = f'{lf_loc}/io_results.log'

        if not os.path.isfile(lf_name):
            print_header = 1
        else:
            print_header = 0

        logfile = open(lf_name, 'a')
        
        if print_header: 
            logfile.write(f'----------////      TRACKING INPUT / OUTPUT           ////----------\n\n')
        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_io_statement(self, step_name, num_files_in, num_obj_in, num_files_out, num_obj_out):
        '''
        Return an statement for the input/output log file
        '''
        
        statement = f'----------////               {step_name}                ////----------\n'

        try:
            statement += f'Files   in/out: {num_files_in} / {num_files_out}, percentage gain/loss (+/-): {-1 * np.round(((num_files_in - num_files_out) / num_files_in) * 100, 2)}%\n'
        except:
            statement += f'Files   in/out: {num_files_in} / {num_files_out}\n'

        try: 
            statement += f'Objects in/out: {num_obj_in} / {num_obj_out}, percentage gain/loss (+/-): {-1 * np.round(((num_obj_in - num_obj_out) / num_obj_in) * 100, 2)}%\n'
        except: 
            statement += f'Objects in/out: {num_obj_in} / {num_obj_out}\n'
        statement += f'\n'

        return statement


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_input_dicts(self, input_directory, evt_id_tab, functions, params_in, phases):
        '''
        Parallel processing needs a list of inputs for each independent iteration
        So put all input into a dict and create a list of dicts

        Parameters
        ----------
        input_directory: str
            input_directory name where all the data is.
        evt_id_tab: pd dataframe
            pandas dataframe containing cmt event details
        functions : list
            list of functions to apply to inputs
        params_in : dict

        phases : dict
            phase dictionary for picks
            
        Returns
        ----------
        input_dicts : list (of dictionaries as below)
            input_dict : dict : specific to each event
                input dictionary containing:
                    event : Path to event
                    event_id : cmt event id
                    evtm : event time full format
                    evla : ev latitude
                    evlo : ev longitude
                    evdp : ev depth
                    mag : ev magnitude
                    id_cmt : cmt catalogue id
                    id_ctm : centroid id
                    id_fmt_ctm : formatted centroid id
                    phases : phases : dict
                    functions : list of functions to execute
                    params_in : input paramters
        '''

        # Finds files, makes input dictionary of files, functions, inputs
        events = sorted(glob.glob(f'{input_directory}/{str(params_in['year'])}*'))
        input_dicts = []

        if len(events) == 0:
            # No events found...
            print('----------')
            print(f'----------////   NO-EVENTS-IN: {input_directory}    ////----------')
        else:

            print('----------')
            print('Processing the following events list:   ')
            print('----------')

            # If the match file exists, read it and add to dictionaries.
            try:
                if os.path.isfile(params_in['mtf_outfilename']):            
                    with open(params_in['mtf_outfilename'], 'r') as file:
                        match_twin_files = [line.strip() for line in file]
            except:
                pass

            for event in events:
                # First check we have the event in the cmt-event-table
                id_fmt_ctm = int(event.split('/')[-1])
                ind = np.where( np.array(evt_id_tab['id_fmt_ctm']) == id_fmt_ctm )[0]

                # Should only be one unique event...
                if len(ind) == 1: 
                    print(f'{str(ind[0])}, {str(event)}')
                    input_dict = {}
                    input_dict['event'] = event
                    input_dict['event_id'] = ind[0]
                    input_dict['evtm'] = evt_id_tab['evtm'][ind[0]]
                    input_dict['evla'] = evt_id_tab['evla'][ind[0]]
                    input_dict['evlo'] = evt_id_tab['evlo'][ind[0]]
                    input_dict['evdp'] = evt_id_tab['evdp'][ind[0]]
                    input_dict['mag'] = evt_id_tab['mag'][ind[0]]
                    input_dict['id_cmt'] = evt_id_tab['id_cmt'][ind[0]]
                    input_dict['id_ctm'] = evt_id_tab['id_ctm'][ind[0]]
                    input_dict['id_fmt_ctm'] = id_fmt_ctm

                    input_dict['phases'] = phases
                    input_dict['functions'] = functions
                    input_dict['params_in'] = params_in

                    # Counting stats
                    input_dict['num_files_in'] =  0 # Files
                    input_dict['num_obj_in'] =    0 # Twin or picks
                    input_dict['num_files_out'] = 0 # Files 
                    input_dict['num_obj_out'] =   0 # Twin or picks

                    try: 
                        if len(match_twin_files) > 0:
                            input_dict['match_twin_files'] = match_twin_files
                    except:
                        pass
                    input_dicts.append(input_dict)
                elif len(ind) == 0 :
                    print(f'**WARNING: cannot find event {str(id_fmt_ctm):s} in catalog\n')
                else:
                    print(f'**ERROR Something wrong event {str(id_fmt_ctm):s} duplcated in catalog?\n')
            
            print('----------\n')

        return input_dicts


    def execute(self, main_function, input_directory, evt_id_tab, functions, params_in, phases):
        '''
        Takes main function inputs and passes to serial or parallel executer, see below.
        '''
        if params_in['parallel']:
            # Execute in parallel
            self.apply_parallel(main_function, input_directory, evt_id_tab, functions, params_in, phases)
        else:
            # Execute in series
            self.apply_serial(main_function, input_directory, evt_id_tab, functions, params_in, phases)
        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def apply_serial(self, main_function, input_directory, evt_id_tab, functions, params_in, phases):
        '''
        Takes a list of functions, reads all input files, applies each function to all files
        IN SERIES...

        Parameters
        ----------
        main_function : list
            Length one list containing main function.
        input_directory: str
            input_directory name where all the data is.
        evt_id_tab: pd dataframe
            pandas dataframe containing cmt event details
        functions : list
            list of functions to apply to inputs
        params_in : dict

        phases : dict
            phase dictionary for picks
            
        Returns
        ----------
        Nothing : None
            Functions take care of returns/writing to file  
            
        The functions in the list must take the following inputs: 

        input_dict : dict
        file : file location string
        seis : obspy stream object read into memory onwhich operations are performed.
        logfile : open logfile to report progress
        outfile : open outfile to report results
        fail : int - fail flag = 1 if function fails, else 0.
        
        The functions must return all input variables
        Functions can write to logfile & outfile.
        '''

        # Counting stats
        num_files_in =  0 # Files
        num_obj_in =    0 # Twin or picks
        num_files_out = 0 # Files 
        num_obj_out =   0 # Twin or picks

        # get input dicts
        input_dicts = self.get_input_dicts(input_directory, evt_id_tab, functions, params_in, phases)
        
        io_logfile = self.open_io_log_file(params_in)


        if input_dicts:
            # Serial processing
            for input_dict in input_dicts:
        
                for main_f in main_function:
                    results = main_f(input_dict)
                    num_files_in  += results['num_files_in'] 
                    num_obj_in    += results['num_obj_in']
                    num_files_out += results['num_files_out']
                    num_obj_out   += results['num_obj_out']

            statement  = self.get_io_statement(results['step_name'], num_files_in, num_obj_in, num_files_out, num_obj_out)
            ###
            self.print_log(params_in, io_logfile, statement)
            ###
            io_logfile.close()

            return
        else:
            pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def apply_parallel(self, main_function, input_directory, evt_id_tab, functions, params_in, phases):
        '''
        Takes a list of functions, reads all input files, applies each function to all of the files
        IN PARALLEL...

        Parameters
        ----------
        main_function : list
            Length one list containing main function.
        input_directory: str
            input_directory name where all the data is.
        evt_id_tab: pd dataframe
            pandas dataframe containing cmt event details
        functions : list
            list of functions to apply to inputs
        params_in : dict

        phases : dict
            phase dictionary for picks
            
        Returns
        ----------
        Nothing : None
            Functions take care of returns/writing to file    
            
        The functions in the list must take the following inputs: 

        input_dict : dict
        file : file location string
        seis : obspy stream object read into memory onwhich operations are performed.
        logfile : open logfile to report progress
        outfile : open outfile to report results
        fail : int - fail flag = 1 if function fails, else 0.
        
        The functions must return all input variables
        Functions can write to logfile & outfile.
        '''
        
        # Counting stats
        num_files_in =  0 # Files
        num_obj_in =    0 # Twin or picks
        num_files_out = 0 # Files 
        num_obj_out =   0 # Twin or picks

        # get input dicts
        input_dicts = self.get_input_dicts(input_directory, evt_id_tab, functions, params_in, phases)

        io_logfile = self.open_io_log_file(params_in)

        if input_dicts:
            # Parallel processing
            for main_f in main_function:

                with concurrent.futures.ProcessPoolExecutor(max_workers = params_in['cores']) as executor:
                    results = executor.map(main_f, input_dicts)
                    for res in results:
                        num_files_in  += res['num_files_in'] 
                        num_obj_in    += res['num_obj_in']
                        num_files_out += res['num_files_out']
                        num_obj_out   += res['num_obj_out']

            statement  = self.get_io_statement(res['step_name'], num_files_in, num_obj_in, num_files_out, num_obj_out)
            ###
            self.print_log(params_in, io_logfile, statement)
            ###
            io_logfile.close()


            return
        else:
            pass


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def read_twin_file(self, file_path):
        '''
        Read twin file from file_path into pandas dataframe
        Ignore header lines before column headers
        '''
        with open(file_path, 'r') as file:
            # Skip header lines
            try:
                header = 1
                while header:
                    fline = file.readline()
                    if 'columns format' in fline:
                        header = 0

                # Read data into DataFrame
                df = pd.read_csv(file, delimiter = ',')
            except:
                #Return empy dataframe if not possible.
                df = pd.DataFrame()

        return df

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def read_tdelay_file(self, file_path):
        '''
        Read tdelay file from file_path into pandas dataframe
        Add header lines of Earthquake details
        '''
        with open(file_path, 'r') as file:
            # Skip header lines
            try:
                header = 1
                while header:
                    fline = file.readline()
                    if 'event_name' in fline:
                        evid = str(fline.split()[-1])
                    if 'date_time' in fline:
                        date_time = str(fline.split()[-1])
                    if 'latitude' in fline:
                        evt_lat = str(fline.split()[-1])
                    if 'longitude' in fline:
                        evt_lon = str(fline.split()[-1])
                    if 'depth' in fline:
                        evt_dep = str(fline.split()[-1])
                    if 'Mw' in fline:
                        evt_mag = str(fline.split()[-1])
                    if 'Channels' in fline or 'channels' in fline:
                        channel = fline.split()[-1][-1]
                    if 'columns format' in fline:
                        header = 0

                convert_dict = {'nslc': str,
                            'stla': float,
                            'stlo': float,
                            'stel': float,
                            'phase': str,
                            'tdelay': float,
                            'tderr': float,
                            'ccmx': float,
                            'ttaup': float,
                            'tp_obs': float,
                            'tp_syn': float,
                            'Ap_obs': float,
                            'Ap_syn': float}

                # Read data into DataFrame
                s_df = pd.read_csv(file, delimiter = ',', dtype = convert_dict)

                eq_df = pd.DataFrame(columns = ['evid','date_time','evt_lat','evt_lon','evt_dep','evt_mag','channel'], index = range(0,len(s_df)))
                for index, row in eq_df.iterrows():
                    eq_df.loc[index, "evid"] =  str(evid)
                    eq_df.loc[index, "date_time"] =  str(date_time)
                    eq_df.loc[index, "evt_lat"] =  float(evt_lat)
                    eq_df.loc[index, "evt_lon"] =  float(evt_lon)
                    eq_df.loc[index, "evt_dep"] =  float(evt_dep)
                    eq_df.loc[index, "evt_mag"] =  float(evt_mag)
                    eq_df.loc[index, "channel"] =  channel

                # Merge eq_df, s_df and write out where no nans in line. of merged df.
                merged_df = pd.merge(eq_df, s_df, left_index = True, right_index = True)
                # Drop Nans non_nan_rows = df.dropna()
                df = merged_df.dropna()

            except:
                #Return empy dataframe if not possible.
                df = pd.DataFrame()

        return df
