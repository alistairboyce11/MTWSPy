import glob,time,os 
from obspy import read
import toolkit, find_twin_obs

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
    files = sorted(glob.glob(event + '/' + str(params_in['year']) + '*' + str(params_in['component'])))
    num_files = len(files)

    # Counting stats
    input_dict['step_name'] = str(os.path.basename(__file__).split('.')[0])
    input_dict['num_files_in'] = num_files
    input_dict['num_obj_in'] = 0 # No timewindows yet

    if num_files ==  0:
        # No files found...
        ###
        toolkit.print_log(params_in, logfile, f'----------////   NO-FILES-IN: ' + str(event) + '    ////----------')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'----------////    WORKING ON: ' + str(event) + '    ////----------')
        ###

        # Loop through file list and apply functions to each
        for k, file in enumerate(files):
            log_statement = toolkit.get_log_statement(event_id, file) + ', [' + str(k + 1) + '/' + str(num_files) + ']'
            
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

    params_list = ['data_loc','component', 'twin_plot_pic', 'twin_save_pic', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'k_em', 'npow']
    for k, param in enumerate(params_list):
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x = justify) )

    # Complex parameters where multiplication factors are used...  
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output loc',' : ',str(params_in['home'] + '/' + params_in['twin_syn_out_loc'] + params_in['component']), x = justify) )
    logfile.write('')
    if '/specfem/' in input_dict['event']:
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
    write key params to outfile # To be changed for each script
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
        if '/specfem/' in input_dict['event']:
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
    for h,header in enumerate(params_in['twin_syn_outcols']):
        if h < len(params_in['twin_syn_outcols']) - 1:
            outfile.write('{0:s}'.format(header + ','))
        else:
            outfile.write('{0:s}'.format(header + '\n'))


    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_log_file(input_dict):
    '''
    Return an open log file in log_loc/'code_start_time'/'filename'/event_name.log
    '''
    params_in = input_dict['params_in']

    lf_loc = params_in['home'] + '/' + params_in['log_loc'] + '/' + str(params_in['code_start_time']) + '/' + os.path.basename(__file__).split('.')[0]
    if not os.path.exists(lf_loc):
        os.makedirs(lf_loc)

    lf_name = lf_loc + '/' + str(input_dict['id_fmt_ctm']) + '.log'

    logfile = open(lf_name,'w')
    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_outfile_file(input_dict):
    '''
    Return an open twin file in twin_loc/filename_out_loc + component/event_name."OS/MX" + component.twin
    '''
    params_in = input_dict['params_in']

    of_loc = params_in['home'] + '/' + params_in['twin_loc'] + '/' + params_in['twin_syn_out_loc'] + params_in['component']
    if not os.path.exists(of_loc):
        os.makedirs(of_loc)

    of_name = of_loc + '/' + str(input_dict['id_ctm']) + '.' + params_in['twin_syn_out_loc'][-2:] + params_in['component'] + '.twin'

    outfile = open(of_name,'w')
    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():

    # Params go here.
    params_in = toolkit.get_params('params_in.yaml')
  
    # Define input data directory and function list.
    input_directory = str(params_in['synth_loc']) + '/e' + str(params_in['year']) + str(params_in['fmt_data_loc'])
    functions = [find_twin_obs.get_tt_times, find_twin_obs.filter_seis, find_twin_obs.detect_window_peaks, find_twin_obs.filter_window_peaks, find_twin_obs.plot_waveform_envelope_peaks_windows]

    # Get event id table as pandas data frame
    evt_id_tab = toolkit.get_event_id_table(params_in['cmt_outfile'])

    # Get phases names dictionary 
    import v01_phasenames
    phases = v01_phasenames.phases()

    start_time = time.time()
    main_function = [process_one_event]

    toolkit.execute(main_function, input_directory, evt_id_tab, functions, params_in, phases)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ ==  '__main__':
    main()


