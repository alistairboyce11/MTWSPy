import glob, sys, os
import yaml
from yaml.loader import SafeLoader
import pandas as pd
from datetime import datetime
import numpy as np
import concurrent.futures

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
#                                                                                     #
#                                     toolkit.py                                      #
#                                                                                     # 
# ---                                                                             --- # 

# A bunch of functions utilised throughout the code.

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def check_files(loc, year, component):
    '''
    Check files exist for given obs and synth data location for
    specified year and component
    Otherwise quit the code...
    '''
    check_name=str(loc)+'/'+str(year)+'*/'+str(year)+'*'+str(component)
    filelist=glob.glob(check_name)
    
    numfiles=len(filelist)

    if numfiles > 0:
        print(f'Found: {str(numfiles):s} files\n')
        print(f'At: {check_name:s}\n')
        print('Continuing...\n')
    else:
        print(f'ERROR....\n')
        print(f'Nothing found at: {check_name:s}\n')
        print(f'Retry...\n')
        sys.exit()

    return

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_params(params_filename):
    '''
    Get input parameters from yaml file
    '''
    with open(params_filename) as f:
        params = yaml.load(f, Loader=SafeLoader)
    # add the starttime of the code
    now = datetime.now()
    params['code_start_time']=now.strftime("%Y%m%d%H%M%S")

    cha_obs, cha_syn = str(params['twin_obs_out_loc'][-2:])+str(params['component']), str(params['twin_syn_out_loc'][-2:])+str(params['component'])
    params['mtf_outfilename'] = params['home']+'/'+params['code_loc']+'/'+str(params['mtf_prefix']) + '.' + cha_obs + '-' + cha_syn

    return params

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_event_id_table(cmt_outfile):
    '''
    Get event id table as pandas data frame
    '''
    tab=pd.read_csv(cmt_outfile)

    return tab

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_log_statement(event_id, filename):
    '''
    Return statement for logfile 
    '''
    filename_p=filename.split('/')[-1]
    statement=str(event_id)+', '+str(filename_p)

    return statement

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def subvec(t, x, interval):
    '''
    Return indices and values within a given interval [t1,t2] for timeseries x(t) 
    '''
    t_new = (t >= interval[0]) & (t <= interval[1])
    x_new = x[t_new]

    return t_new, x_new

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def print_log(params,logfile,statement):
    '''
    Print statement to terminal (depending on input params) and open logfile 
    '''
    if params['verbose']:
        print(statement)
    
    try:
        logfile.write(statement+'\n')
    except:
        # Probably something that cannot be written just printed.
        pass

    return

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_input_dicts(input_directory, evt_id_tab, functions, params_in, phases):
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
    events=sorted(glob.glob(input_directory+'/'+str(params_in['year'])+'*'))
    input_dicts=[]

    if len(events) == 0:
        # No events found...
        print('----------')
        print('----------////   NO-EVENTS-IN: '+str(input_directory)+'    ////----------')
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
            id_fmt_ctm=int(event.split('/')[-1])
            ind=np.where( np.array(evt_id_tab['id_fmt_ctm']) == id_fmt_ctm )[0]

            # Should only be one unique event...
            if len(ind) == 1: 
                print(str(ind[0])+', '+str(event))
                input_dict={}
                input_dict['event']=event
                input_dict['event_id']=ind[0]
                input_dict['evtm']=evt_id_tab['evtm'][ind[0]]
                input_dict['evla']=evt_id_tab['evla'][ind[0]]
                input_dict['evlo']=evt_id_tab['evlo'][ind[0]]
                input_dict['evdp']=evt_id_tab['evdp'][ind[0]]
                input_dict['mag']=evt_id_tab['mag'][ind[0]]
                input_dict['id_cmt']=evt_id_tab['id_cmt'][ind[0]]
                input_dict['id_ctm']=evt_id_tab['id_ctm'][ind[0]]
                input_dict['id_fmt_ctm']=id_fmt_ctm
                input_dict['phases']=phases
                input_dict['functions']=functions
                input_dict['params_in']=params_in
                try: 
                    if len(match_twin_files) > 0:
                        input_dict['match_twin_files']=match_twin_files
                except:
                    pass
                input_dicts.append(input_dict)
            elif len(ind) == 0 :
                print(f'**WARNING: cannot find event {str(id_fmt_ctm):s} in catalog\n')
            else:
                print(f'**ERROR Something wrong event {str(id_fmt_ctm):s} duplcated in catalog?\n')
        
        print('----------\n')

    return input_dicts


def execute(main_function, input_directory, evt_id_tab, functions, params_in, phases):
    '''
    Takes main function inputs and passes to serial or parallel executer, see below.
    '''
    if params_in['parallel']:
        # Execute in parallel
        apply_parallel(main_function,input_directory,evt_id_tab,functions,params_in,phases)
    else:
        # Execute in series
        apply_serial(main_function,input_directory,evt_id_tab,functions,params_in,phases)
    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def apply_serial(main_function, input_directory, evt_id_tab, functions, params_in, phases):
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

    # get input dicts
    input_dicts=get_input_dicts(input_directory, evt_id_tab, functions, params_in, phases)
    
    if input_dicts:
        # Serial processing
        for input_dict in input_dicts:
    
            for main_f in main_function:
                res=main_f(input_dict)
    
        return
    else:
        pass

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def apply_parallel(main_function, input_directory, evt_id_tab, functions, params_in, phases):
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

    # get input dicts
    input_dicts=get_input_dicts(input_directory, evt_id_tab, functions, params_in, phases)

    if input_dicts:
        # Parallel processing
        for main_f in main_function:

            with concurrent.futures.ProcessPoolExecutor(max_workers=params_in['cores']) as executor:
                results = executor.map(main_f, input_dicts)

            return
    else:
        pass


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def read_twin_file(file_path):
    '''
    Read twin file from file_path into pandas dataframe
    Ingore header lines before column headers
    '''
    with open(file_path, 'r') as file:
        # Skip header lines
        try:
            header = 1
            while header:
                fline=file.readline()
                if 'columns format' in fline:
                    header = 0

            # Read data into DataFrame
            df = pd.read_csv(file, delimiter=',')
        except:
            #Return empy dataframe if not possible.
            df = pd.DataFrame()

    return df



