import time,os
import pandas as pd
import numpy as np
from obspy.taup import TauPyModel
import obspy.geodetics.base
import inspect
from toolkit import Toolkit


class PhaseAssociation:
    tk = Toolkit()

    def __init__(self):
        pass

    def common_function(self):
        # Common function used by both phase_association_obs.py and phase_association_syn.py
        pass


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def associate_twin_phase(self, input_dict, twin_in, twin_df, logfile, outfile, fail):
        '''
        Takes dataframe of twin file, iterates through unique stations
        Collects predicted travel times
        Attempts to locate phases within max_tshift of each time window
        rejects twins with multiple phases in a window, that are not depth phases.
        Saves associated time windows to file.
        
        Reads & Returns
        ----------
        input_dict : dict
        twin_in : twin file location string
        twin_df : twin dataframe object
        logfile : open logfile to report progress
        outfile : open outfile to report results
        fail : int - fail flag = 1 if function fails, else 0.
        '''
        params_in = input_dict['params_in']
        phases = input_dict['phases']

        # setup logfile statement using event_id, filename and function name (with inspect)
        log_statement = f'{self.tk.get_log_statement(input_dict['event_id'], twin_in)}, {str(inspect.stack()[0][3])}'
        if fail:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING- \n')
            ###
        else:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING- \n')
            ###

            # Get event parameters
            evla = input_dict['evla']
            evlo = input_dict['evlo']
            evdp = input_dict['evdp']

            # find unique stations in df.
            unique_stations = list(twin_df['nslc'].unique())
            nst = len(unique_stations)

            # Make empty df for new data:
            df = pd.DataFrame(columns = ['phase','n_depth_phase','t_taup'], index = range(0,len(twin_df)))

            for s, station in enumerate(unique_stations):
                # Find the indices where nslc ==  station
                stat_df = twin_df.query("nslc ==  @station")
                len_stat_df = len(stat_df)
                inds = stat_df.index.tolist()

                ###
                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  {len_stat_df} windows for station {station} [{str(s + 1)}/{str(nst)}]')
                ###

                # Get station lat, lon, elev
                stla = stat_df.latitude.to_numpy()[0]
                stlo = stat_df.longitude.to_numpy()[0]
                # stel = stat_df.elevation.to_numpy()[0]

                model = TauPyModel(model = params_in['taup_model_name'])
                phs = phases[str(params_in['phases_key'])][str(params_in['component'])]

                distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(evla, evlo, stla, stlo)
                gcarc = distm / (6371.e3 * np.pi / 180.)

                # Get phase arrival times for station event pair
                ttimes = []; ttt = {}
                for ph in range(len(phs)):
                    # extract only first time value, else NaN
                    # Loop faster than an list comprehension.
                    try:
                        # removed ,receiver_depth_in_km = stdp
                        arrivals = model.get_travel_times(source_depth_in_km = evdp,distance_in_degree = gcarc,phase_list = [phs[ph]])
                        ttime = arrivals[0].time
                    except:
                        ttime = np.NaN
                    ttimes.append(ttime)
                    ttt[phs[ph]] = ttime

                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -- {str(phs[ph]):s}: {ttime:.3f}s')
                    ###

                ttimes = np.array(ttimes)

                phs = np.array(phs)

                # Iterate through rows of station dataframe
                for index, row in stat_df.iterrows():
                    ndph = 0  # number of depth phases other than the main phase
                    tshift = row['t_peak'] - ttimes

                    # Find phase index where less than max_tshift from the peak
                    iph = np.array(np.where(np.abs(tshift) <=  params_in['max_tshift'])[0])
                    
                    if iph.size ==  0:
                        # No nearby phases
                        j = np.nanargmin(np.abs(tshift))
                        ###
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  reject twin @ {row['t_peak']:.1f}s without any phase [nearest: {phs[j] } @ {ttimes[j]:.1f}s]')
                        ###
                        continue
                    elif iph.size ==  1:
                        # Unique phase found
                        phw = phs[iph[0]]
                    else:
                        uphs, nuph = np.unique([phs[i] for i in iph], return_counts = True)

                        if len(iph) > 3 or any(nuph > 1):  # multiphase window
                            ###
                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  reject multiphase twin @ {row['t_peak']:.1f}s: {', '.join([f'{phase} @ {ttimes[iph[p]]:.1f}s' for p,phase in enumerate([phs[i] for i in iph])])}\n')
                            ###
                            continue

                        # a phase X and possibly corresponding depth phases p|sX
                        uuphs = self.strip_first_letter(uphs) # phase names without leading 'p|s'

                        if len(uuphs) ==  1:
                            # Only have main phase and depth phase in window
                            phw = uuphs[0]
                            ndph = len(iph) - 1
                        else:
                            ###
                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  reject multiphase twin @ {row['t_peak']:.1f}s: {', '.join([f'{phase} @ {ttimes[iph[p]]:.1f}s' for p,phase in enumerate([phs[i] for i in iph])])}')
                            ###
                            continue

                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  got twin @ {row['t_peak']:.1f} for {phw} with {ndph} depth phases {', '.join([f'{phase} @ {ttimes[iph[p]]:.1f}s' for p,phase in enumerate([phs[i] for i in iph])])}')
                    ###

                    # # Store results in df
                    # df['phase'][index] = str(phw)
                    # df['n_depth_phase'][index] = ndph
                    # df['t_taup'][index] = ttimes[iph[0]]

                    # Change of style: df.loc[row_indexer, "col"] = values
                    df.loc[index, "phase"] = str(phw)
                    df.loc[index, "n_depth_phase"] = ndph
                    df.loc[index, "t_taup"] = ttimes[iph[0]]


            # Merge twin_df with df and write out where no nans in line. of merged df.
            merged_df = pd.merge(twin_df, df, left_index = True, right_index = True)

            # Drop Nans non_nan_rows = df.dropna()
            non_nan_merged_df = merged_df.dropna()
            
            # Find maximum "pronouced peak" window amplitude, normalised by window length - shorter windows
            max_pronouced_peak_amplitude = ((non_nan_merged_df['A_peak'] / non_nan_merged_df[['A_left', 'A_right']].mean(axis=1))   / non_nan_merged_df['t_right'].sub( non_nan_merged_df['t_left'])      ).max()

            # Remove rows of df that refer to the same station & phase & taup time
            # keep row that minimises absolute difference between the taup_time and the t_peak ()
            drop_list = []

            unique_stations = list(non_nan_merged_df['nslc'].unique())

            for s, station in enumerate(unique_stations):
                # Find the indices where nslc ==  station
                stat_df = non_nan_merged_df.query("nslc ==  @station")

                unique_phases = list(stat_df['phase'].unique())

                for u, unique_phase in enumerate(unique_phases):
                    phase_df = stat_df.query("phase == @unique_phase")
                    
                    # Find taup_time of phase
                    taup_time = phase_df['t_taup'].iloc[0]

                    # If len(phase_df == 0) dont do anything.
                    if len(phase_df) > 1:

                        phase_df_copy = phase_df.copy()
                        # For each row, divide A_peak by average amplitude of edges and normalise by most pronouced window in dataframe
                        phase_df_copy['A_peak_normalized'] = ((phase_df_copy['A_peak'] / phase_df_copy[['A_left', 'A_right']].mean(axis=1))    / phase_df_copy['t_right'].sub( phase_df_copy['t_left'])      ) / max_pronouced_peak_amplitude
                        # For each row normalise the absolute time difference between the t_peak and taup_tiime by the max_thsift and take 1 minus this value.
                        phase_df_copy['t_peak_normalized'] = 1 - (phase_df['t_peak'].sub(taup_time).abs() / params_in['max_tshift'])
                        # Sum the two normalised columns and keep maximum index
                        phase_df_copy['peak_normalized'] = phase_df_copy[['A_peak_normalized', 't_peak_normalized']].sum(axis=1)
                        idx_select = phase_df_copy['peak_normalized'].idxmax()

                        # # find index of row that minimises the absolute difference between the taup_time and the t_peak
                        # idx_select_time = phase_df['t_peak'].sub(taup_time).abs().idxmin()

                        # if idx_select != idx_select_time:
                        #     print(phase_df_copy)
                        #     sys.exit()

                        # Select all other entries
                        idx = [x for x in phase_df.index.tolist() if idx_select !=  x]
                        drop_list.append(idx)


            # Remove rows that fail QC above in dataframe
            drop_list = self.tk.flatten_concatenation(drop_list)
            if drop_list:
                # print('dropping indices: ', drop_list)
                ###
                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Removing {len(drop_list)} multi-window phases to avoid pick duplication...')
                ###
                non_nan_merged_df = non_nan_merged_df.drop(drop_list)

            if not non_nan_merged_df.empty:
                for index, row in non_nan_merged_df.iterrows():
                    # print(str(row['phase']),row['n_depth_phase'],row['t_taup'])
                    tw_info = params_in['phase_a_outfmt'].format(row['nslc'],row['latitude'],row['longitude'],row['elevation'],row['t_left'],row['t_peak'],row['t_right'],row['A_left'],row['A_peak'],row['A_right'],row['A_noise'],str(row['phase']),row['n_depth_phase'],row['t_taup']) + '\n'
                    outfile.write(tw_info)

                # Counting stats
                input_dict['num_obj_out']   += len(non_nan_merged_df)
                input_dict['num_files_out'] += 1

        return input_dict, twin_in, twin_df, logfile, outfile, fail

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def strip_first_letter(self, uphs):
        '''
        This function takes an array of phase names and strips any starting with 's' or 'p'
        Returns unique array without depth phases
        '''
        result = [phase[1:] if phase[0] in {'s', 'p'} else phase for phase in uphs]
        return np.unique(np.array(result))


class PhaseAssociationObs(PhaseAssociation):
    tk = Toolkit()

    def __init__(self):
        super().__init__()

    def obs_specific_function(self):
        # Function specific to phase_association_obs.py
        pass


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def process_one_event(self, input_dict):
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
        # event = input_dict['event']
        # event_id = input_dict['event_id']
        id_ctm = input_dict['id_ctm']
        functions = input_dict['functions']
        params_in = input_dict['params_in']

        # Initiate logfile
        logfile = self.open_log_file(input_dict)
        logfile = self.write_params_logfile(input_dict, logfile)
        input_dict['step_name'] = str(os.path.basename(__file__).split('.')[0])

        # Proceed with list of matched events in input_dict:
        if len(input_dict['match_twin_files']) > 0:
            if id_ctm in input_dict['match_twin_files']:
                # This event is matched... proceed.
                tw_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['twin_obs_out_loc']}{params_in['component']}'
                twin_in = f'{tw_loc}/{str(id_ctm)}.{params_in['twin_obs_out_loc'][-2:]}{params_in['component']}.twin'

                if os.path.isfile(twin_in):
                    # read twin file
                    twin_df = self.tk.read_twin_file(twin_in)
                    # Check df is not empty and then proceed:
                    if not twin_df.empty:
                        
                        # Initiate outfile
                        outfile = self.open_outfile_file(input_dict)
                        outfile = self.write_params_outfile(input_dict, outfile)
                        
                        # Counting stats
                        input_dict['num_files_in'] += 1
                        input_dict['num_obj_in'] = len(twin_df)

                        ###
                        self.tk.print_log(params_in, logfile, f'----------////    WORKING ON: {str(id_ctm)}    ////----------')
                        ###

                        fail = 0

                        if len(functions) ==  0: 
                            ###
                            self.tk.print_log(params_in, logfile, f'----------////    NO_FUNCTIONS_TO_APPLY    ////----------')
                            ###
                            pass

                        else:
                            # apply functions, only executed when fail = 0
                            for function in functions:
                                input_dict, twin_in, twin_df, logfile, outfile, fail = function(input_dict, twin_in, twin_df, logfile, outfile, fail)

                        ###    
                        self.tk.print_log(params_in, logfile, f'----------////    FINISHED    ////----------')
                        ###

                        outfile.close()

                    else:
                        # Returned empty dataframe
                        ###
                        self.tk.print_log(params_in, logfile, f'----------////   NO-TWINS-IN: {str(twin_in)}    ////----------')
                        ###
                        pass

                else:
                    # No corresponding twin file found at twin_in
                    ###
                    self.tk.print_log(params_in, logfile, f'----------////   NO-TWIN-FILE-AT: {str(twin_in)}    ////----------')
                    ###
                    pass
                
        else:
            # No file found at params_in['mtf_outfilename']
            ###
            self.tk.print_log(params_in, logfile, f'----------////   NO-MATCHED_TWIN-FILE-AT: {str(params_in['mtf_outfilename'])}    ////----------')
            ###
            pass

        logfile.close()
        
        return input_dict


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def write_params_logfile(self, input_dict, logfile):
        '''
        write key params to logfile. To be changed for each script
        '''
        justify = 30

        params_in = input_dict['params_in']
        phases = input_dict['phases']

        logfile.write(' ')
        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('id_fmt_ctm',' : ',str(input_dict['id_fmt_ctm']), x = justify) )

        params_list = ['obs_loc','component', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'min_snr_P', 'min_snr_A', 'npow']
        for k, param in enumerate(params_list):
            logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x = justify) )

        # Complex parameters where multiplication factors are used...  
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min signal window length [s]',' : ',str(params_in['min_sig_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min noise window length [s]',' : ',str(params_in['min_nois_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output loc', ' : ', f'{params_in['home']}/{params_in['phase_a_obs_out_loc']}{params_in['component']}', x = justify) )
        logfile.write('')
        if 'syn' in input_dict['event']:
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
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x = justify) )

        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_params_outfile(self, input_dict, outfile):
        '''
        write key params to outfile. To be changed for each script
        '''
        justify = 30
        params_in = input_dict['params_in']
        phases = input_dict['phases']

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
            if 'syn' in input_dict['event']:
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

        outfile.write('----------////              PHASE ASSOCIATION PARAMETERS               ////----------\n')
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output channel',' : ', f'{params_in['phase_a_obs_out_loc']}{params_in['component']}', x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ', str(params_in['max_tshift']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ', str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
        outfile.write('----------\n')

        outfile.write('{0:<20s} {1:s} {2:s}\n'.format('columns format',' : ',params_in['phase_a_outfmt']))
        for h,header in enumerate(params_in['phase_a_outcols']):
            if h < len(params_in['phase_a_outcols']) - 1:
                outfile.write('{0:s}'.format(f'{header},'))
            else:
                outfile.write('{0:s}'.format(f'{header}\n'))
        
        return outfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_log_file(self, input_dict):
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
    def open_outfile_file(self, input_dict):
        '''
        Return an open twin file in twin_loc/filename_out_loc + component/event_name."OS/MX" + component.twin
        '''
        params_in = input_dict['params_in']

        of_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['phase_a_obs_out_loc']}{params_in['component']}'
        
        if not os.path.exists(of_loc):
            os.makedirs(of_loc, exist_ok=True)

        of_name = f'{of_loc}/{str(input_dict['id_ctm'])}.{params_in['phase_a_obs_out_loc'][-2:]}{params_in['component']}.twin'

        outfile = open(of_name,'w')
        return outfile





class PhaseAssociationSyn(PhaseAssociation):
    tk = Toolkit()

    def __init__(self):
        super().__init__()

    def syn_specific_function(self):
        # Function specific to phase_association_syn.py
        pass



    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def process_one_event(self, input_dict):
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
        # event = input_dict['event']
        # event_id = input_dict['event_id']
        id_ctm = input_dict['id_ctm']
        functions = input_dict['functions']
        params_in = input_dict['params_in']

        # Initiate logfile
        logfile = self.open_log_file(input_dict)
        logfile = self.write_params_logfile(input_dict, logfile)
        input_dict['step_name'] = str(os.path.basename(__file__).split('.')[0])

        # Proceed with list of matched events in input_dict:
        if len(input_dict['match_twin_files']) > 0:
            if id_ctm in input_dict['match_twin_files']:
                # This event is matched... proceed.

                tw_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['twin_syn_out_loc']}{params_in['component']}'
                twin_in = f'{tw_loc}/{str(id_ctm)}.{params_in['twin_syn_out_loc'][-2:]}{params_in['component']}.twin'
                
                if os.path.isfile(twin_in):
                    # read twin file
                    twin_df = self.tk.read_twin_file(twin_in)
                    # Check df is not empty and then proceed:
                    if not twin_df.empty:
                        
                        # Initiate outfile
                        outfile = self.open_outfile_file(input_dict)
                        outfile = self.write_params_outfile(input_dict, outfile)
        
                        # Counting stats
                        input_dict['num_files_in'] += 1
                        input_dict['num_obj_in'] = len(twin_df)


                        # Synth rename A_median (if any) to A_noise for syn twin
                        try:
                            twin_df.rename(columns = {'A_median': 'A_noise'}, inplace = True)
                        except:
                            pass

                        ###
                        self.tk.print_log(params_in, logfile, f'----------////    WORKING ON: '+str(id_ctm)+'    ////----------')
                        ###

                        fail = 0

                        if len(functions) == 0: 
                            ###
                            self.tk.print_log(params_in, logfile, f'----------////    NO_FUNCTIONS_TO_APPLY    ////----------')
                            ###
                            pass

                        else:
                            # apply functions, only executed when fail = 0
                            for function in functions:
                                input_dict, twin_in, twin_df, logfile, outfile, fail = function(input_dict, twin_in, twin_df, logfile, outfile, fail)
                            
                        ###    
                        self.tk.print_log(params_in, logfile, f'----------////    FINISHED    ////----------')
                        ###

                        outfile.close()

                    else:
                        # Returned empty dataframe
                        ###
                        self.tk.print_log(params_in, logfile, f'----------////   NO-TWINS-IN: '+str(twin_in)+'    ////----------')
                        ###
                        pass

                else:
                    # No corresponding twin file found at twin_in
                    ###
                    self.tk.print_log(params_in, logfile, f'----------////   NO-TWIN-FILE-AT: '+str(twin_in)+'    ////----------')
                    ###
                    pass

        else:
            # No file found at params_in['mtf_outfilename']
            ###
            self.tk.print_log(params_in, logfile, f'----------////   NO-MATCHED_TWIN-FILE-AT: '+str(params_in['mtf_outfilename'])+'    ////----------')
            ###
            pass

        logfile.close()
        return input_dict


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def write_params_logfile(self, input_dict, logfile):
        '''
        write key params to logfile. To be changed for each script
        '''
        justify = 30

        params_in = input_dict['params_in']
        phases = input_dict['phases']

        logfile.write(' ')
        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('id_fmt_ctm',' : ',str(input_dict['id_fmt_ctm']), x = justify) )

        params_list = ['syn_loc','component', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'min_snr_P', 'min_snr_A', 'npow']
        for k, param in enumerate(params_list):
            logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x = justify) )

        # Complex parameters where multiplication factors are used...  
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min signal window length [s]',' : ',str(params_in['min_sig_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min noise window length [s]',' : ',str(params_in['min_nois_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output loc',' : ',str(params_in['home']+'/'+params_in['phase_a_syn_out_loc']+params_in['component']), x = justify) )
        logfile.write('')
        if 'syn' in input_dict['event']:
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
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x = justify) )

        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_params_outfile(self, input_dict, outfile):
        '''
        write key params to outfile. To be changed for each script
        '''
        justify = 30
        params_in = input_dict['params_in']
        phases = input_dict['phases']

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
            if 'syn' in input_dict['event']:
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

        outfile.write('----------////              PHASE ASSOCIATION PARAMETERS               ////----------\n')
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('output channel',' : ',str(params_in['phase_a_syn_out_loc']+params_in['component']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
        outfile.write('----------\n')

        outfile.write('{0:<20s} {1:s} {2:s}\n'.format('columns format',' : ',params_in['phase_a_outfmt']))
        for h,header in enumerate(params_in['phase_a_outcols']):
            if h < len(params_in['phase_a_outcols']) - 1:
                outfile.write('{0:s}'.format(f'{header},'))
            else:
                outfile.write('{0:s}'.format(f'{header}\n'))
        
        return outfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_log_file(self, input_dict):
        '''
        Return an open log file in log_loc/'code_start_time'/'filename'/event_name.log
        '''
        params_in = input_dict['params_in']

        lf_loc = f'{params_in['home']}/{params_in['log_loc']}/{str(params_in['code_start_time'])}/{os.path.basename(__file__).split('.')[0]}'

        if not os.path.exists(lf_loc):
            os.makedirs(lf_loc, exist_ok=True)

        lf_name = f'{lf_loc}/{str(input_dict['id_fmt_ctm'])}.log'

        logfile = open(lf_name, 'w')
        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_outfile_file(self, input_dict):
        '''
        Return an open twin file in twin_loc/filename_out_loc+component/event_name."OS/MX"+component.twin
        '''
        params_in = input_dict['params_in']

        of_loc = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['phase_a_syn_out_loc']}{params_in['component']}'
        if not os.path.exists(of_loc):
            os.makedirs(of_loc, exist_ok=True)

        of_name = f'{of_loc}/{str(input_dict['id_ctm'])}.{params_in['phase_a_syn_out_loc'][-2:]}{params_in['component']}.twin'

        outfile = open(of_name, 'w')
        return outfile


def main():


    # Create instance of toolkit as tk
    tk = Toolkit()

    start_time = time.time()
    # Params go here.
    params_in = tk.get_params('params_in.yaml')
    # Get event id table as pandas data frame
    evt_id_tab = tk.get_event_id_table(params_in)

    from v01_phasenames import Phases
    phases = Phases().get_phase_dictionary()

    #################### Do the observed. ############################

    # Define observed data input directory
    obs_input_directory = f'{str(params_in['obs_loc'])}/e{str(params_in['year'])}{str(params_in['fmt_data_loc'])}'
    
    # Instantiate the classes and call their methods as needed
    phase_association_obs = PhaseAssociationObs()

    functions = [phase_association_obs.associate_twin_phase]

    main_function = [phase_association_obs.process_one_event]

    tk.execute(main_function, obs_input_directory, evt_id_tab, functions, params_in, phases)

    #################### Do the synthetics. ############################

    # Define input data directory and function list.
    syn_input_directory = f'{str(params_in['syn_loc'])}/e{str(params_in['year'])}{str(params_in['fmt_data_loc'])}'

    phase_association_syn = PhaseAssociationSyn()

    functions = [phase_association_syn.associate_twin_phase]

    main_function = [phase_association_syn.process_one_event]

    tk.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)

    # phase_association_syn.common_function()
    # phase_association_syn.syn_specific_function()


if __name__ == '__main__':
    main()
