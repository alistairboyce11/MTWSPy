import time,os
import pandas as pd
import numpy as np
# from obspy.taup import TauPyModel
# import obspy.geodetics.base
import inspect
import toolkit, phase_association_obs

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

    Returns
    ----------
    Nothing: None
    '''
    # event=input_dict['event']
    # event_id=input_dict['event_id']
    id_ctm=input_dict['id_ctm']
    functions=input_dict['functions']
    params_in=input_dict['params_in']

    # Initiate logfile
    logfile=open_log_file(input_dict)
    logfile=write_params_logfile(input_dict, logfile)

    # Proceed with list of matched events in input_dict:
    if len(input_dict['match_twin_files']) > 0:
        if id_ctm in input_dict['match_twin_files']:
            # This event is matched... proceed.

            tw_loc=params_in['home']+'/'+params_in['twin_loc']+'/'+params_in['twin_syn_out_loc']+params_in['component']
            twin_in=tw_loc+'/'+str(id_ctm)+'.'+params_in['twin_syn_out_loc'][-2:]+params_in['component']+'.twin'
            
            if os.path.isfile(twin_in):
                # read twin file
                twin_df = toolkit.read_twin_file(twin_in)
                # Check df is not empty and then proceed:
                if not twin_df.empty:
                    
                    # Initiate outfile
                    outfile=open_outfile_file(input_dict)
                    outfile=write_params_outfile(input_dict, outfile)

                    # Synth rename A_median (if any) to A_noise for syn twin
                    try:
                        twin_df.rename(columns={'A_median': 'A_noise'}, inplace=True)
                    except:
                        pass

                    ###
                    toolkit.print_log(params_in, logfile, f'----------////    WORKING ON: '+str(id_ctm)+'    ////----------')
                    ###

                    fail = 0

                    if len(functions) == 0: 
                        ###
                        toolkit.print_log(params_in, logfile, f'----------////    NO_FUNCTIONS_TO_APPLY    ////----------')
                        ###
                        pass

                    else:
                        # apply functions, only executed when fail = 0
                        for function in functions:
                            input_dict,twin_in,twin_df,logfile,outfile,fail = function(input_dict,twin_in,twin_df,logfile,outfile,fail)
                        
                    ###    
                    toolkit.print_log(params_in, logfile, f'----------////    FINISHED    ////----------')
                    ###

                    outfile.close()

                else:
                    # Returned empty dataframe
                    ###
                    toolkit.print_log(params_in, logfile, f'----------////   NO-TWINS-IN: '+str(twin_in)+'    ////----------')
                    ###
                    pass

            else:
                # No corresponding twin file found at twin_in
                ###
                toolkit.print_log(params_in, logfile, f'----------////   NO-TWIN-FILE-AT: '+str(twin_in)+'    ////----------')
                ###
                pass

    else:
        # No file found at params_in['mtf_outfilename']
        ###
        toolkit.print_log(params_in, logfile, f'----------////   NO-MATCHED_TWIN-FILE-AT: '+str(params_in['mtf_outfilename'])+'    ////----------')
        ###
        pass

    logfile.close()
    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def write_params_logfile(input_dict, logfile):
    '''
    write key params to logfile. To be changed for each script
    '''
    justify=30

    params_in=input_dict['params_in']
    phases=input_dict['phases']

    logfile.write(' ')
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('id_fmt_ctm',' : ',str(input_dict['id_fmt_ctm']), x=justify) )

    params_list=['data_loc','component', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'min_snr_P', 'min_snr_A', 'npow']
    for k, param in enumerate(params_list):
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x=justify) )

    # Complex parameters where multiplication factors are used...  
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min signal window length [s]',' : ',str(params_in['min_sig_win_f']*params_in['T2']), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min noise window length [s]',' : ',str(params_in['min_nois_win_f']*params_in['T2']), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output loc',' : ',str(params_in['home']+'/'+params_in['phase_a_syn_out_loc']+params_in['component']), x=justify) )
    logfile.write('')
    if '/specfem/' in input_dict['event']:
        min_amp=float(params_in['min_amp_syn'])
        min_snr = float(params_in['min_snr_syn'])
        wsz_lim = params_in['wsz_lim_syn']
        bnd2pk_t = params_in['bnd2pk_t_syn']
        bnd2pk_A = params_in['bnd2pk_A_syn']

    else:
        min_amp=float(params_in['min_amp_obs'])
        min_snr = params_in['min_snr']
        wsz_lim = params_in['wsz_lim_obs']
        bnd2pk_t = params_in['bnd2pk_t_obs']
        bnd2pk_A = params_in['bnd2pk_A_obs']

    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min peak amp',' : ',str(min_amp), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min snr',' : ',str(min_snr), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('wsz_lim [s]',' : ',str(wsz_lim), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('bnd2pk_t [s]',' : ',str(bnd2pk_t), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('bnd2pk_A',' : ',str(bnd2pk_A), x=justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x=justify) )

    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def write_params_outfile(input_dict, outfile):
    '''
    write key params to outfile. To be changed for each script
    '''
    justify = 30
    params_in=input_dict['params_in']
    phases=input_dict['phases']

    outfile.write(' ')
    outfile.write('----------////               EVENT PARAMETERS                ////----------\n')
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('event_name',' : ',str(input_dict['id_cmt']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('date_time',' : ',str(input_dict['evtm']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('latitude',' : ',str(input_dict['evla']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('longitude',' : ',str(input_dict['evlo']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('depth',' : ',str(input_dict['evdp']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Mw',' : ',str(input_dict['mag']), x=justify) )
    outfile.write('----------////               EVENT PARAMETERS                ////----------\n')
    outfile.write('----------\n')
    if params_in['filttwin']:
        if '/specfem/' in input_dict['event']:
            min_amp=float(params_in['min_amp_syn'])
            min_snr = float(params_in['min_snr_syn'])
            wsz_lim = params_in['wsz_lim_syn']
            bnd2pk_t = params_in['bnd2pk_t_syn']
            bnd2pk_A = params_in['bnd2pk_A_syn']

        else:
            min_amp=float(params_in['min_amp_obs'])
            min_snr = params_in['min_snr']
            wsz_lim = params_in['wsz_lim_obs']
            bnd2pk_t = params_in['bnd2pk_t_obs']
            bnd2pk_A = params_in['bnd2pk_A_obs']

        outfile.write('----------////              FILTTWIN PARAMETERS               ////----------\n')
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min peak amp',' : ',str(min_amp), x=justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min snr [peak/noise]',' : ',str(min_snr), x=justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('window size limits [s]',' : ',str(wsz_lim), x=justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('edge to peak tshift limits [s]',' : ',str(bnd2pk_t), x=justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('edge to peak Aratio limits',' : ',str(bnd2pk_A), x=justify) )
        outfile.write('----------\n')

    outfile.write('----------////              PHASE ASSOCIATION PARAMETERS               ////----------\n')
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output channel',' : ',str(params_in['phase_a_syn_out_loc']+params_in['component']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x=justify) )
    outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x=justify) )
    outfile.write('----------\n')

    outfile.write('{0:<20s} {1:s} {2:s}\n'.format('columns format',' : ',params_in['phase_a_outfmt']))
    for h,header in enumerate(params_in['phase_a_outcols']):
        if h < len(params_in['phase_a_outcols']) - 1:
            outfile.write('{0:s}'.format(header+','))
        else:
            outfile.write('{0:s}'.format(header+'\n'))
    
    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_log_file(input_dict):
    '''
    Return an open log file in log_loc/'filename'/'code_start_time'/event_name.log
    '''
    params_in=input_dict['params_in']

    lf_loc=params_in['home']+'/'+params_in['log_loc']+'/'+str(params_in['code_start_time'])+'/'+os.path.basename(__file__).split('.')[0]
    if not os.path.exists(lf_loc):
        os.makedirs(lf_loc)

    lf_name=lf_loc+'/'+str(input_dict['id_fmt_ctm'])+'.log'

    logfile=open(lf_name,'w')
    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_outfile_file(input_dict):
    '''
    Return an open twin file in twin_loc/filename_out_loc+component/event_name."OS/MX"+component.twin
    '''
    params_in=input_dict['params_in']

    of_loc=params_in['home']+'/'+params_in['twin_loc']+'/'+params_in['phase_a_syn_out_loc']+params_in['component']
    if not os.path.exists(of_loc):
        os.makedirs(of_loc)

    of_name=of_loc+'/'+str(input_dict['id_ctm'])+'.'+params_in['phase_a_syn_out_loc'][-2:]+params_in['component']+'.twin'

    outfile=open(of_name,'w')
    return outfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():
    # Params go here.
    params_in=toolkit.get_params('params_in.yaml')

    # Define input data directory and function list.
    input_directory=str(params_in['synth_loc'])+'/e'+str(params_in['year'])+str(params_in['fmt_data_loc'])

    functions=[phase_association_obs.associate_twin_phase]

    # Get event id table as pandas data frame
    evt_id_tab = toolkit.get_event_id_table(params_in['cmt_outfile'])

    # Get phases names dictionary 
    import v01_phasenames
    phases = v01_phasenames.phases()

    start_time = time.time()
    main_function=[process_one_event]

    toolkit.execute(main_function, input_directory, evt_id_tab, functions, params_in, phases)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()




