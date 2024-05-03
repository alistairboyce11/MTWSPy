import sys,glob,time,os
import pandas as pd
import numpy as np
import scipy
from obspy import read, UTCDateTime
from obspy.taup import TauPyModel
import inspect
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 8
from matplotlib.ticker import (MultipleLocator)
import toolkit

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def process_one_event(input_dict):
    '''
    Processes one event (earthquake), reading the file, applying the list functions to them

    Parameters
    ----------
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

    Returns : input_dict
    ----------
    Nothing: None
    '''
    event = input_dict['event']
    event_id = input_dict['event_id']
    functions = input_dict['functions']
    params_in = input_dict['params_in']

    # Initiate logfile
    logfile = open_log_file(input_dict)
    logfile = write_params_logfile(input_dict, logfile)

    # Initiate outfile
    outfile = open_outfile_file(input_dict)
    outfile = write_params_outfile(input_dict, outfile)

    # List of files in the event
    files = sorted(glob.glob(f'{event}/{str(params_in['year'])}*{str(params_in['component'])}'))

    num_files = len(files)

    # Counting stats
    input_dict['step_name'] = str(os.path.basename(__file__).split('.')[0])
    input_dict['num_files_in'] = num_files
    input_dict['num_obj_in'] = 0 # No timewindows yet

    if num_files ==  0:
        # No files found...
        ###
        toolkit.print_log(params_in, logfile, f'----------////   NO-FILES-IN: {str(event)}    ////----------')
        ###

    else:
        ###
        toolkit.print_log(params_in, logfile, f'----------////    WORKING ON: {str(event)}    ////----------')
        ###

        # Loop through file list and apply functions to each
        for k, file in enumerate(files):
            log_statement = f'{toolkit.get_log_statement(event_id, file)}, [{str(k + 1)}/{str(num_files)}]'

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  PROCESSING.....')
            ###

            # read sac file
            try:
                seis = read(file,'sac')
                fail = 0
            except:
                fail = 1
                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  FAILED_LOADING...')
                ###

                continue
            
            if len(functions) ==  0: 
                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  NO_FUNCTIONS_TO_APPLY.....')
                ###

            else:
                # apply functions, only executed when fail = 0
                for function in functions:
                    input_dict, file, seis, logfile, outfile, fail = function(input_dict, file, seis, logfile, outfile, fail)
            
            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  , FINISHED.....')
            ###

    logfile.close()
    outfile.close()
    return input_dict

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def write_params_logfile(input_dict, logfile):
    '''
    write key params to logfile. To be changed for each script
    '''
    justify = 30

    params_in = input_dict['params_in']
    phases = input_dict['phases']

    logfile.write(' ')
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('id_fmt_ctm',' : ',str(input_dict['id_fmt_ctm']), x = justify) )

    params_list = ['obs_loc','component', 'twin_plot_pic', 'twin_save_pic', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'min_snr_P', 'min_snr_A', 'npow']
    for k, param in enumerate(params_list):
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x = justify) )

    # Complex parameters where multiplication factors are used...  
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min signal window length [s]',' : ',str(params_in['min_sig_win_f']*params_in['T2']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min noise window length [s]',' : ',str(params_in['min_nois_win_f']*params_in['T2']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output loc',' : ', f'{params_in['home']}/{params_in['twin_obs_out_loc']}{params_in['component']}', x = justify) )
    logfile.write('')
    if '/syn/' in input_dict['event']:
        min_amp = float(params_in['min_amp_syn'])
        min_snr = float(params_in['min_snr_syn'])
        wsz_lim = params_in['wsz_lim_syn']
        bnd2pk_t = params_in['bnd2pk_t_syn']
        bnd2pk_A = params_in['bnd2pk_A_syn']

    else:
        min_amp = float(params_in['min_amp_obs'])
        min_snr = params_in['min_snr']
        wsz_lim = params_in['wsz_lim_obs']
        bnd2pk_t = params_in['bnd2pk_t_obs']
        bnd2pk_A = params_in['bnd2pk_A_obs']

    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min peak amp',' : ',str(min_amp), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min snr',' : ',str(min_snr), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('wsz_lim [s]',' : ',str(wsz_lim), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('bnd2pk_t [s]',' : ',str(bnd2pk_t), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('bnd2pk_A',' : ',str(bnd2pk_A), x = justify) )
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

    return logfile

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def write_params_outfile(input_dict, outfile):
    '''
    write key params to outfile. To be changed for each script
    '''
    justify = 30
    params_in = input_dict['params_in']

    outfile.write(' ')
    outfile.write('----------////               EVENT PARAMETERS                ////----------\n')
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('event_name',' : ',str(input_dict['id_cmt']), x = justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('date_time',' : ',str(input_dict['evtm']), x = justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('latitude',' : ',str(input_dict['evla']), x = justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('longitude',' : ',str(input_dict['evlo']), x = justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('depth',' : ',str(input_dict['evdp']), x = justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Mw',' : ',str(input_dict['mag']), x = justify) )
    outfile.write('----------////               EVENT PARAMETERS                ////----------\n')
    outfile.write('----------\n')
    if params_in['filttwin']:
        if '/syn/' in input_dict['event']:
            min_amp = float(params_in['min_amp_syn'])
            min_snr = float(params_in['min_snr_syn'])
            wsz_lim = params_in['wsz_lim_syn']
            bnd2pk_t = params_in['bnd2pk_t_syn']
            bnd2pk_A = params_in['bnd2pk_A_syn']

        else:
            min_amp = float(params_in['min_amp_obs'])
            min_snr = params_in['min_snr']
            wsz_lim = params_in['wsz_lim_obs']
            bnd2pk_t = params_in['bnd2pk_t_obs']
            bnd2pk_A = params_in['bnd2pk_A_obs']

        outfile.write('----------////              FILTTWIN PARAMETERS               ////----------\n')
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min peak amp',' : ',str(min_amp), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min snr [peak/noise]',' : ',str(min_snr), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('window size limits [s]',' : ',str(wsz_lim), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('edge to peak tshift limits [s]',' : ',str(bnd2pk_t), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('edge to peak Aratio limits',' : ',str(bnd2pk_A), x = justify) )
        outfile.write('----------\n')

    outfile.write('{0:<20s} {1:s} {2:s}\n'.format('columns format',' : ',params_in['twin_outfmt']))
    for h,header in enumerate(params_in['twin_obs_outcols']):
        if h < len(params_in['twin_obs_outcols']) - 1:
            outfile.write('{0:s}'.format(f'{header},'))
        else:
            outfile.write('{0:s}'.format(f'{header}\n'))
    
    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_tt_times(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Adds all phase traveltimes requested to seis[0].stats.traveltimes dictionary
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''

    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        phases = input_dict['phases']

        #### Determine Taup Travel times ######
        model = TauPyModel(model = params_in['taup_model_name'])
        phs = phases[str(params_in['phases_key'])][str(params_in['component'])]
        evdp = float(seis[0].stats['sac']['evdp']) # Event depth
        gcarc = float(seis[0].stats['sac']['gcarc']) # Great circle arc distance
        stdp = float(seis[0].stats['sac']['stdp'])/1000 # Station depth in Km

        ttimes = []; ttt = {}
        if not hasattr(seis[0].stats, 'traveltimes'):
            seis[0].stats.traveltimes = dict()

        for ph in range(len(phs)):
            # extract only first time value, else NaN
            try:
                arrivals = model.get_travel_times(source_depth_in_km = evdp,distance_in_degree = gcarc,receiver_depth_in_km = stdp,phase_list = [phs[ph]])
                ttime = arrivals[0].time
            except:
                ttime  = np.NaN
            ttimes.append(ttime)
            seis[0].stats.traveltimes[phs[ph]] = ttime
            ttt[phs[ph]] = ttime

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -- {str(phs[ph]):s}: {ttime:.3f}s')
            ###

        tlim = [np.nanmin(ttimes),np.nanmax(ttimes)] # taup predicted sigwin

        # Save to input_dict: tlim,ttt
        write_to_dict(input_dict,['tlim','ttt'],[tlim,ttt])

    return input_dict, filename, seis, logfile, outfile, fail

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def primary_qc(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Primary quality control for seismogram in memory
    based on SNR_amplitude, SNR_power
    Using signal and noise windows
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''
        
    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        dt = seis[0].stats['delta'] # 1/sample-rate
        npts = seis[0].stats['npts']
           
        # Recover from input_dict values for the specified keys: tlim
        tlim = read_from_dict(input_dict,['tlim'])[0]

        ### Determine time axis for obspy trace, time axis relative to ev time
        t = seis[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))

        # Run check on minimum signal window length
        min_sig_win = params_in['min_sig_win_f']*params_in['T2']
        if t[0] > tlim[0] or t[-1] - tlim[0] < min_sig_win:
            fail = 1

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **too short [t = {t[0]:.1f}-{t[-1]:.1f}s; sigwin = {tlim[0]:.1f}-{tlim[1]:.1f}s -FAILED- ]')
            ###

            return input_dict, filename, seis, logfile, outfile, fail

        # slice time axis to get signal window indices and array
        xb = [tlim[0],np.max([tlim[0] + 1800,min_sig_win])]
        id_sig = np.arange(np.argmax(t >=  xb[0]), np.argmin(t <=  xb[1]) + 1)
         
        if params_in['sig_win_type'] > 1:
            tlim[1] = t[-1] # extend to the end

        # determine noise window (20min length) for snr evaluation
        min_nois_win = params_in['min_nois_win_f']*params_in['T2']
        xb = [max(t[round(params_in['T2']/dt)], tlim[0] - params_in['sig_win_ext'] - 1200), tlim[0] - params_in['sig_win_ext']]
        id_noise = np.arange(np.argmax(t >=  xb[0]), np.argmin(t <=  xb[1]) + 1)

        if len(id_noise)*dt <=  min_nois_win: # use last quarter for snr
            # evaluation if pre-signal noise is insufficient
            id_noise = np.arange(int(npts * 0.75), npts)

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **use last quarter ({t[id_noise][0]:.1f}-{t[id_noise][-1]:.1f}s) as noiswin [t = {t[0]:.1f}-{t[-1]:.1f}s, first arrival = {t[0]:.1f}s]')
            ###

        # Pre-process seismograms
        seis.detrend()
        seis.taper(max_percentage = 0.05, type = 'cosine')
        seis.filter('bandpass',freqmin = 1/params_in['T2']*2*dt,freqmax = 1/params_in['T1']*2*dt,corners = 4,zerophase = True)


        if len(id_sig) ==  0 or len(id_noise) ==  0:
            fail = 1
            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **poor sig/noise win selection: skip -FAILED- ]')
            ###
            return input_dict, filename, seis, logfile, outfile, fail

        # Get Signal-to-noise properties
        sig_pow = np.var(seis[0].data[id_sig]) 
        nois_pow = np.var(seis[0].data[id_noise])

        if nois_pow ==  0 or nois_pow ==  np.nan:
            fail = 1
            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **poor noise variance: skip -FAILED- ]')
            ###
            return input_dict, filename, seis, logfile, outfile, fail
        
        snr_pow = sig_pow / nois_pow
        sig_max = np.max(np.abs(seis[0].data[id_sig]))
        snr_amp = sig_max / np.sqrt(nois_pow) / 2

        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  nois_std = {np.sqrt(nois_pow):.2e}, sig_max = {sig_max:.2e}')
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  snr_P = {snr_pow:.2f}, snr_A = {snr_amp:.2f}')
        ###

        if sig_pow <=  params_in['min_snr_P'] * nois_pow or snr_amp <=  params_in['min_snr_A']:
            fail = 1
            
            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **poor quality: skip -FAILED- ]')
            ###

        # Save to input_dict: id_sig,id_noise
        write_to_dict(input_dict,['id_sig','id_noise'],[id_sig,id_noise])

    return input_dict, filename, seis, logfile, outfile, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def filter_seis(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Simple function to filter synth seis to match obs primary_qc
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''
        
    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        # Pre-process seismograms
        dt = seis[0].stats['delta']
        seis.detrend()
        seis.taper(max_percentage = 0.05, type = 'cosine')
        seis.filter('bandpass',freqmin = 1/params_in['T2']*2*dt,freqmax = 1/params_in['T1']*2*dt,corners = 4,zerophase = True)

        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  Filtered synthetic to match observed')
        ###

    return input_dict, filename, seis, logfile, outfile, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def detect_window_peaks(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Morphological Time Window Selection (MTWS) method of Li et al., (2023)
    Uses a number of SNR constraints to propose many candidate time windows
    based on a power envelope around maxima at id_pks and left and right
    edges at i1 and i2.
    Wrties all candidate windows to outfile
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''

    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        dt = seis[0].stats['delta']
        stel = float(seis[0].stats['sac']['stel']) # Station elevation in meters
        stla = float(seis[0].stats['sac']['stla']) # Station latitude
        stlo = float(seis[0].stats['sac']['stlo']) # Station longitude

        nslc = f'{str(seis[0].stats['network'])}_{str(seis[0].stats['station'])}.{str(seis[0].stats['location'])}.{str(seis[0].stats['channel'])}'

        t = seis[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))

        ####### Detect window peaks on envelope ######
        e1_h = scipy.signal.hilbert(seis[0].data) # Hilbert transform
        e1 = np.abs(e1_h) # envelope of signal
        e = e1**params_in['npow'] #  n-power envelope

        mxf_win = params_in['mxf_win_f']*params_in['Tc']

        # Use Pandas dataframe.rolling to compute a rolling maximum over the n-power envelope array, e : Mirrors MATLAB MOVMAX
        b = np.array(pd.DataFrame(e).rolling(window = int(np.floor(mxf_win/dt)), min_periods = 1, center = True).max().values.flatten())
        b[[0, -1]] = e[[0, -1]] + np.finfo(float).eps  # fix ends, avoid divide by zero errors
        id = np.where(e ==  b)[0]  # max extrema
        
        ### Execute for observed data
        if '/obs/' in input_dict['event']:
            
            # Recover from input_dict values for the specified keys: tlim,id_sig,id_noise
            tlim, id_sig, id_noise = read_from_dict(input_dict,['tlim', 'id_sig', 'id_noise'])

            e_noi_m = np.mean(e[id_noise]) # mean of noise window 
            e_noi_sd = np.std(e[id_noise]) # std of noise window 
            e_sig_m = np.mean(e[id_sig]) # mean of signal window
            e_sig_sd = np.std(e[id_sig]) # std of signal window

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  e_noi_m = {e_noi_m:.2e}, e_noi_sd = {e_noi_sd:.2e}')
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  e_sig_m = {e_sig_m:.2e}, e_sig_sd = {e_sig_sd:.2e}')
            ###

            id_low = e[id] <=  e_noi_m * params_in['min_snr']

            bool_1 = params_in['min_snr'] is not None
            bool_2 = params_in['min_snr'] > 0
        
        ### Execute for synthetic data
        if '/syn/' in input_dict['event']:

            # Recover from input_dict values for the specified keys: tlim,id_sig
            tlim, id_sig = read_from_dict(input_dict,['tlim', 'id_sig'])

            e_sig_m = np.median(e[id_sig]) # median of signal window
            e_sig_sd = np.std(e[id_sig]) # std of signal window

            logfile.write(f'{log_statement:s}  ,  e_sig_m = {e_sig_m:.2e}, e_sig_sd = {e_sig_sd:.2e}\n')

            id_low = e[id] < e_sig_m * params_in['k_em']

            bool_1 = params_in['k_em'] is not None
            bool_2 = params_in['k_em']  > 0

            # to write np.median(e[id_sig]) to final column of twin file instead of e_noi_m in observed case
            e_noi_m = e_sig_m

        if bool_1 and bool_2:  # a hard threshold
            if np.any(id_low):
                
                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_low)} / {len(id)} low peaks')
                ###

            id = id[~id_low]

        if params_in['sig_win_type'] ==  0:  # use full trace
            id_pks = id
        else:  # reject pre- and post-signal peaks
            id_in = (t[id] > tlim[0] - params_in['sig_win_ext']) & (t[id] < tlim[1] + params_in['sig_win_ext'])
            id_pks = id[id_in]
            if np.sum(id_in) < len(id):

                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {len(id) - np.sum(id_in)} / {len(id)} outside tlim [{tlim[0]:.2f},{tlim[1]:.2f}] s')
                ###

        # Remove a peak at the first/last index...
        id_pks = id_pks[id_pks !=  0]
        id_pks = id_pks[id_pks !=  len(t)-1]

        if not id_pks.any():
            fail = 1

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **no prominent peaks [skip]')
            ###

            return input_dict, filename, seis, logfile, outfile, fail

        # determine the window around each peak
        walkaway =  params_in['walkaway_f']*params_in['Tc']
        i1,i2 = tw_det(e,t,id_pks,walkaway)

        twin_span = (i2 - i1) * dt
        
        min_win_span = np.max([float(params_in['min_win_span_f'][0]*params_in['T1']),float(params_in['min_win_span_f'][1])])

        if not params_in['filttwin']:
            tw_info = ''
            for ii in range(len(twin_span)):
                # if twin_span[ii] < min_win_span:
                tw_info  +=  params_in['twin_outfmt'].format(nslc, stla, stlo, stel, t[i1[ii]], t[id_pks[ii]], t[i2[ii]], e1[i1[ii]], e1[id_pks[ii]], e1[i2[ii]], e_noi_m) + '\n' # need e_sig_m as final parameter when using synth data.

            outfile.write(tw_info)

            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  --- time windows --- \n' + tw_info + '---')
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  Save all {len(i1)} windows to {outfile:s}')
            ###

            # Counting stats
            input_dict['num_obj_out']   += len(twin_span)
            input_dict['num_files_out'] += 1


        # Save to input_dict: t,e_noi_m,e_sig_m,min_win_span, nslc, mxf_win, e1, e, i1, id_pks, i2
        write_to_dict(input_dict,['t','dt','e_noi_m','e_sig_m','min_win_span', 'nslc', 'stel', 'stla', 'stlo', 'mxf_win', 'e1', 'e', 'i1', 'id_pks', 'i2', 'twin_span'],[t,dt,e_noi_m,e_sig_m,min_win_span,nslc,stel,stla,stlo,mxf_win, e1, e, i1, id_pks, i2,twin_span])

    return input_dict, filename, seis, logfile, outfile, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def filter_window_peaks(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Uses a number of SNR constraints to filter candidate time windows
    based on a XXX around maxima at id_pks and left and right
    edges at i1 and i2.
    Wrties remaining candidate windows to outfile
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''

    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        if params_in['filttwin']:

            #### Recover from input_dict values for the specified keys:
            t,e_noi_m,nslc,stel,stla,stlo,e1, i1, id_pks, i2,twin_span = read_from_dict(input_dict,['t','e_noi_m', 'nslc', 'stel', 'stla', 'stlo', 'e1', 'i1', 'id_pks', 'i2', 'twin_span'])

            if '/syn/' in input_dict['event']:
                min_amp = float(params_in['min_amp_syn'])
                min_snr = float(params_in['min_snr_syn'])
                wsz_lim = params_in['wsz_lim_syn']
                bnd2pk_t = params_in['bnd2pk_t_syn']
                bnd2pk_A = params_in['bnd2pk_A_syn']
            else:
                min_amp = float(params_in['min_amp_obs'])
                min_snr = params_in['min_snr']
                wsz_lim = params_in['wsz_lim_obs']
                bnd2pk_t = params_in['bnd2pk_t_obs']
                bnd2pk_A = params_in['bnd2pk_A_obs']

            assert all(value <=  1 for value in bnd2pk_A), "Assertion failed: One or more values are not less than 1."

            #### reject twin too narrow or too wide
            if wsz_lim is not None:
                if np.isscalar(wsz_lim):
                    id_sz = twin_span < wsz_lim
                else:
                    id_sz = (twin_span < wsz_lim[0]) | (twin_span > wsz_lim[1])

                if np.any(id_sz):
                    ###
                    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_sz)} / {len(id_pks)} windows with span outside [{wsz_lim[0]}, {wsz_lim[1]}]s')
                    ###

                i1, id_pks, i2, twin_span = filter_arrays(i1, id_pks, i2, twin_span, id_sz)

            #### Check if bnd2pk_t is not empty
            if bnd2pk_t is not None:
                p2e_span = np.column_stack([t[id_pks] - t[i1], t[i2] - t[id_pks]])
                id_sz = (np.min(p2e_span, axis = 1) < bnd2pk_t[0]) | (np.max(p2e_span, axis = 1) > bnd2pk_t[1])

                if np.any(id_sz):
                    ###
                    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_sz)} / {len(id_pks)} windows with edge-to-peak outside [{bnd2pk_t[0]}, {bnd2pk_t[1]}]s')
                    ###

                i1, id_pks, i2, twin_span = filter_arrays(i1, id_pks, i2, twin_span, id_sz)

            #### reject non-prominent twin based on minimum amplitude
            if min_amp is not None and min_amp > 0:
                id_lowamp = e1[id_pks] < min_amp
                
                if np.any(id_lowamp):
                    ###
                    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_lowamp)} / {len(id_pks)} windows of amp < {min_amp}')
                    ###

                i1, id_pks, i2, twin_span = filter_arrays(i1, id_pks, i2, twin_span, id_lowamp)

            #### reject non-prominent twin based on minimum SNR
            if min_snr is not None and min_snr > 0:
                id_lowsnr = e1[id_pks] < min_snr * e_noi_m # This params is e_sig_m in synth case

                if np.any(id_lowsnr):
                    ###
                    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_lowsnr)} / {len(id_pks)} windows of snr < {min_snr:.2f}')
                    ###

                i1, id_pks, i2, twin_span = filter_arrays(i1, id_pks, i2, twin_span, id_lowsnr)

            #### reject twin with high-valued boundary
            if bnd2pk_A is not None:                
                id_hb2p = (np.maximum(e1[i1], e1[i2]) > e1[id_pks] * bnd2pk_A[0]) | \
                        ((e1[i1] + e1[i2]) > e1[id_pks] * bnd2pk_A[1])

                if np.any(id_hb2p):
                    ###
                    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject {np.sum(id_hb2p)} / {len(id_pks)} windows with high edges')
                    ###

                i1, id_pks, i2, twin_span = filter_arrays(i1, id_pks, i2, twin_span, id_hb2p)

            if not np.any(twin_span):
                fail = 1
                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **no prominent twin remains....\n')
                ###
                
                return input_dict, filename, seis, logfile, outfile, fail

            else:
                tw_info = ''
                for ii in range(len(twin_span)):
                    # if twin_span[ii] < min_win_span:
                        # ###
                        # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **reject narrow window: {t[i1[ii]]} -- {t[i2[ii]]} = {twin_span[ii]:.1f}')
                        # ###
                    tw_info  +=  params_in['twin_outfmt'].format(nslc, stla, stlo, stel, t[i1[ii]], t[id_pks[ii]], t[i2[ii]], e1[i1[ii]], e1[id_pks[ii]], e1[i2[ii]], e_noi_m) + '\n' # need e_sig_m as final parameter when using synth data.

                outfile.write(tw_info)

                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  --- time windows --- \n' + tw_info + '---\n')
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  Save filtered {len(i1)} windows to {outfile}')
                ###

                # Counting stats
                input_dict['num_obj_out']   += len(twin_span)
                input_dict['num_files_out'] += 1

                # Save to input_dict: 
                write_to_dict(input_dict,['i1', 'id_pks', 'i2', 'twin_span'],[i1, id_pks, i2,twin_span])

    return input_dict, filename, seis, logfile, outfile, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def plot_waveform_envelope_peaks_windows(input_dict, filename, seis, logfile, outfile, fail):
    '''
    Plots trace and envelope (upper subfigure) and power envelope 
    with time windows and predicted arrival times
    Saves or shows plots
    
    Reads & Returns
    ----------
    input_dict : dict
    file : file location string
    seis : obspy stream object read into memory onwhich operations are performed.
    logfile : open logfile to report progress
    outfile : open outfile to report results
    fail : int - fail flag = 1 if function fails, else 0.
    '''
    
    # setup logfile statement using event_id, filename and function name (with inspect)
    log_statement = f'{toolkit.get_log_statement(input_dict['event_id'], filename)}, {str(inspect.stack()[0][3])}'
    params_in = input_dict['params_in']

    if fail:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
        ###

        if params_in['twin_plot_pic']:
            # Plot pic, else, skip...
            evdp = float(seis[0].stats['sac']['evdp'])
            gcarc = float(seis[0].stats['sac']['gcarc']) 
            baz = float(seis[0].stats['sac']['baz']) # Backazimuth

            # Recover from input_dict: e_noi_m,e_sig_m,min_win_span,nslc, e1, e, i1, id_pks, i2,ttt,t
            e_noi_m,e_sig_m,nslc, e1, e, i1, id_pks, i2,ttt,t,tlim = read_from_dict(input_dict,['e_noi_m','e_sig_m', 'nslc', 'e1', 'e', 'i1', 'id_pks', 'i2','ttt','t','tlim'])
            dt = t[1] - t[0]
            twin_span = (i2 - i1) * dt

            # Only windows longer than `min_win_span`
            id_spanfilt = twin_span >=  input_dict['min_win_span']
            i1, id_pks, i2 = i1[id_spanfilt], id_pks[id_spanfilt], i2[id_spanfilt]
            t_pks, twin_span = t[id_pks], (i2 - i1) * dt

            # Initialise Figure and add axes with constraints
            fig = plt.figure(figsize  = (params_in['twin_plot_width']*params_in['cm'],params_in['twin_plot_height']*params_in['cm']))

            ax0 = fig.add_axes([0.1, params_in['twin_upper_pos'], params_in['twin_sf_width'], params_in['twin_upper_height']], projection = None, polar = False,facecolor = 'white',frame_on = True)
            ax1 = fig.add_axes([0.1, params_in['twin_lower_pos'], params_in['twin_sf_width'], params_in['twin_lower_height']], projection = None, polar = False,facecolor = 'white',frame_on = True)

            axes_list_all = [ax0, ax1]
            for j, ax in enumerate(axes_list_all):
                ax.set_xticks([])
                # ax.set_yticks([])
                ax.set_xlim([max(t[0], tlim[0] - 100), min(t[-1], tlim[-1] + 300)])

                ax.spines["right"].set_linewidth(1.5)
                ax.spines["left"].set_linewidth(1.5)
                ax.spines["top"].set_linewidth(1.5)
                ax.spines["bottom"].set_linewidth(1.5)
                ax.tick_params(labelsize = 12)

            ax0.set_title(f'Event {input_dict['id_cmt']}; Station {nslc}\nDis {gcarc:.2f} deg; Baz {baz:.1f} deg; Depth {evdp:.1f} km',fontsize = 14)
            ax0.spines["bottom"].set_linewidth(0.2) 
            ax1.spines["top"].set_linewidth(0.2)
            max_d = np.max(seis[0].data)
            ax0.set_ylim([-1*(1.2*max_d), max_d*1.2])

            ax1.xaxis.set_minor_locator(MultipleLocator(params_in['twin_x_minor']))
            ax1.xaxis.set_major_locator(MultipleLocator(params_in['twin_x_major']))
            ax1.set_xlabel('Time [s]',fontsize = 14)
            max_e = np.max(e)
            ax1.set_ylim([-1*(0.05*max_e), max_e*1.02])
            
            # Add Waveform to panel (ax0) upper subfigure
            ax0.plot(t, seis[0].data, color = [0.2, 0.2, 0.2], lw = 0.5, label = f'{params_in['component']} waveform')
            ax0.plot(t, e1, color = [0.7, 0.3, 0.2], lw = 0.5, label = f'{params_in['component']} envelope')
            ax0.axhline(y = e_sig_m, xmin = t[0], xmax = t[-1], ls = '--', lw = 0.25, color = [0.2, 0.6, 0.8], label = 'mean(sigwin)')
            ax0.legend(loc = 'upper right', fontsize = 10)

            # Add Envelope to panel (ax1) lower subfigure
            ax1.plot(t, e, ls = '-', lw = 0.5, color = [0.2, 0.2, 0.2], label = f'ord-{params_in['npow']} envelope')
            ax1.scatter(t[id_pks], e[id_pks], s = 16, color = 'r', label = 'candidate peaks')
            b = np.array(pd.DataFrame(e).rolling(window = int(np.floor(input_dict['mxf_win']/dt)), min_periods = 1, center = True).max().values.flatten())
            ax1.plot(t, b, ls = '-', lw = 0.5, color = [0.1, 0.6, 0.3], label = f'{input_dict['mxf_win']:.1f} s movmax')

            # background SNR level
            if '/obs/' in input_dict['event']:
                ax1.axhline(y = e_noi_m * params_in['min_snr'], xmin = t[0], xmax = t[-1], ls = '--', lw = 0.25, color = 'c',label = f'mean(env(noise)) * [1, {params_in['min_snr']}]')
                pic_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['twin_obs_pic_loc']}{params_in['component']}'

            if '/syn/' in input_dict['event']:
                ax1.axhline(y = e_sig_m, xmin = t[0], xmax = t[-1], ls = '--', lw = 0.25, color = 'c',label = f'median(env(sigwin))')
                pic_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['twin_syn_pic_loc']}{params_in['component']}'

            # Add predicted times to both subfigures
            for key, value in ttt.items():
                if  not np.isnan(value):
                    ax0.axvline(x = value, ymin = -1.2*max_d, ymax = 1.2*max_d, ls = '--', lw = 0.25, color = [0.5, 0.5, 0.5], clip_on = True)
                    ax1.axvline(x = value, ymin = 0, ymax = 1.5*max_e, ls = '--', lw = 0.25, color = [0.5, 0.5, 0.5], clip_on = True)
                    ax1.text(value, max_e*0.7, str(key), color = [0.3, 0, 0.3], fontsize = 14, va = 'center', ha = 'left', rotation = 45,rotation_mode = 'anchor')

            # Add the time windows    
            for k in range(len(t_pks)):
                ax1.axvline(x = t[i1[k]], ymin = 0, ymax = 1.02 - 0.15,color = generate_colors(k), lw = 0.5, clip_on = True)
                ax1.axvline(x = t[i2[k]], ymin = 0, ymax = 1.02 - 0.15,color = generate_colors(k), lw = 0.5, clip_on = True)
            
            # Use a dummy variable for the legend
            ax1.axvline(x = np.NaN, ymin = np.NaN, ymax = np.NaN, color = 'r', lw = 0.5, label = 'time windows')
            ax1.legend(loc = 'upper right', fontsize = 10)
            
            # Label subfigures
            ax0.annotate('a',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',fontsize = 12,textcoords = 'offset points', color = 'k', backgroundcolor = 'none',ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))
            ax1.annotate('b',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',fontsize = 12,textcoords = 'offset points', color = 'k', backgroundcolor = 'none',ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))

            # Save or show 
            if params_in['twin_save_pic']:
                pic_filename =  f'{input_dict['id_cmt']}_{nslc}.twin.pdf'
                ###
                toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  Save {pic_loc}/{pic_filename}')
                ###
                
                if not os.path.exists(pic_loc):
                    os.makedirs(pic_loc, exist_ok=True)
                plt.savefig(f'{pic_loc}/{pic_filename}', format = 'pdf')
                plt.close()
            else:
                plt.show()
        else:
            
            ###
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -NO PLOT REQUESTED- ')
            ###
    return input_dict, filename, seis, logfile, outfile, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def filter_arrays(arr1, arr2, arr3, arr4, condition):
    '''
    This function takes four arrays and a condition as input, then returns the filtered arrays based on the provided condition using the ~ (logical NOT) operator.
    '''
    return arr1[~condition], arr2[~condition], arr3[~condition], arr4[~condition]


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def generate_colors(input_variable):
    '''
    Generate colours for time windowing plotting repeating every 20 windows.
    '''
    num_colors = 20
    var = input_variable % num_colors
    cmap = plt.get_cmap('gnuplot2')  # You can choose a different colormap if desired
    norm = plt.Normalize(0, 20)
    colors = cmap(norm(var))
    return colors


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def tw_det(e, t, ipeak, walkaway, balance = False):
    """
    Determine the left and right locations of windows around peaks.

    Parameters:
        e: Signal envelope
        t: Time axis or time interval
        ipeak: Peak locations [ipk1,...] or {wsz, base, k_med}
        walkaway: The dip line across [t(ipk), e(ipk)] and [t(ipk)  + /- walkaway, 0]
        balance: Scale peak to walkaway (default, False)

    Returns:
        Locations of the left and right edges of windows
    """

    if len(e.shape) > 1:
        e = e.flatten()

    t = t.flatten()

    dt = t[1] - t[0]
    assert walkaway > 10 * dt, 'Wrong walkaway'

    ipeak = ipeak.flatten()
    ileft, iright = np.full_like(ipeak, np.nan, dtype = float), np.full_like(ipeak, np.nan, dtype = float)

    for ii in range(len(ipeak)): #len(ipeak)
        ipk = ipeak[ii]
        y = e.copy()

        if balance:
            y = y * walkaway / e[ipk]  # scale y axis

        # Left side
        sid, x = toolkit.subvec(t, t, [t[ipk] + dt - walkaway, t[ipk] - dt])
        dis = point_line_distance([x, y[sid]], [t[ipk], y[ipk]], [t[ipk] - walkaway, 0])
        mid = np.argmax(dis)  # find max dis to the dip line
        ileft[ii] = np.where(sid)[0][mid]  # left bound at t(i1)
        # Right side
        sid, x = toolkit.subvec(t, t, [t[ipk] + dt, t[ipk] + walkaway - dt])
        dis = point_line_distance([x, y[sid]], [t[ipk], y[ipk]], [t[ipk] + walkaway, 0])
        mid = np.argmax(dis)  # max dis to the dip line
        iright[ii] = np.where(sid)[0][mid]  # right bound at t(i2)


        if ileft[ii] < ipk and ipk < iright[ii]:
            # print('ileft[ii] < ipk and ipk < iright[ii]')
            dummy = 1
        else:
            # Make modifications to ileft[ii] && iright[ii] so they are not included in the output
            # Should not happen really....
            ileft[ii]  = ipk
            iright[ii] = ipk

    return ileft.astype(int), iright.astype(int)


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def point_line_distance(*args):
    """
    Calculate distances from points to a line defined by two points.

    Parameters:
        args: Either (x, y, x1, y1, x2, y2) or ([x, y], [x1, y1], [x2, y2])

    Returns:
        distances: Distances to the line
        projections: Projection points (optional)
    """

    if len(args) ==  3:  # (pt = [x, y], pt1 = [x1, y1], pt2 = [x2, y2])
        points, pt1, pt2 = args
        x, y = points[0], points[1]
        x1, y1 = pt1[0], pt1[1]
        x2, y2 = pt2[0], pt2[1]
    else:  # (x, y, x1, y1, x2, y2)
        x, y, x1, y1, x2, y2 = args

    x, y = np.array(x), np.array(y)

    # Switch points on the line to ensure dis > 0 for points below the line
    if x1 > x2:
        x1, y1, x2, y2 = x2, y2, x1, y1

    # Distance = cross(pt2-pt1, pt-pt1) / norm(pt2-pt1)
    #        = (a*x + b*y + c) / sqrt(a*a + b*b)
    numerator = (x2 - x1) * (y1 - y) - (x1 - x) * (y2 - y1)
    denominator = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    distances = numerator / denominator

    # Projections of points onto the line (optional)
    if len(args) > 3:
        if x1 ==  x2:  # Vertical line
            xp = x1 * np.ones_like(y)
            yp = y
        elif y1 ==  y2:  # Horizontal line
            xp = x
            yp = y1 * np.ones_like(x)
        else:
            k = (y2 - y1) / (x2 - x1)
            b = y1 - k * x1
            kv = -1 / k
            bv = -kv * x + y
            xp = (bv - b) / (k - kv)
            yp = kv * xp + bv

        projections = np.column_stack((xp, yp))
        return distances, projections

    return distances


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_log_file(input_dict):
    '''
    Return an open log file in log_loc/'code_start_time'/'filename'/event_name.log
    '''
    params_in = input_dict['params_in']

    lf_loc = f'{params_in['home']}/{params_in['log_loc']}/{str(params_in['code_start_time'])}/{os.path.basename(__file__).split('.')[0]}'

    if not os.path.exists(lf_loc):
        os.makedirs(lf_loc, exist_ok=True)

    lf_name = f'{lf_loc}/{str(input_dict['id_fmt_ctm'])}.log'

    logfile = open(lf_name,'w')
    
    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_outfile_file(input_dict):
    '''
    Return an open twin file in twin_loc/filename_out_loc + component/event_name."OS/MX" + component.twin
    '''
    params_in = input_dict['params_in']

    of_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['twin_obs_out_loc']}{params_in['component']}'

    if not os.path.exists(of_loc):
        os.makedirs(of_loc, exist_ok=True)

    of_name = f'{of_loc}/{str(input_dict['id_ctm'])}.{params_in['twin_obs_out_loc'][-2:]}{params_in['component']}.twin'

    outfile = open(of_name,'w')

    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def write_to_dict(dict, key_list, vals_list):
    '''
    Return the input dict with the key-value parts from the input lists added
    '''
    if len(key_list) !=  len(vals_list):
        # Problem
        sys.exit('Issue with writing to dict')
    for i in range(len(key_list)):
        dict[key_list[i]] = vals_list[i]

    return dict


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def read_from_dict(dict, keys_to_extract):
    '''
    Return the variables in the list keys_to_extract from the input dict
    '''

    return [dict.get(key, None) for key in keys_to_extract]
    

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():
    start_time = time.time()

    # Params go here.
    params_in = toolkit.get_params('params_in.yaml')

    # Define input data directory and function list.
    input_directory = f'{str(params_in['obs_loc'])}/e{str(params_in['year'])}{str(params_in['fmt_data_loc'])}'

    functions = [get_tt_times, primary_qc, detect_window_peaks, filter_window_peaks, plot_waveform_envelope_peaks_windows]

    # Get event id table as pandas data frame
    evt_id_tab = toolkit.get_event_id_table(params_in)

    # Get phases names dictionary 
    import v01_phasenames
    phases = v01_phasenames.phases()

    main_function = [process_one_event]
    
    toolkit.execute(main_function, input_directory, evt_id_tab, functions, params_in, phases)
    
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ ==  '__main__':
    main()




