import time,os
import pandas as pd
import numpy as np
from obspy import read, UTCDateTime
import inspect
from toolkit import Toolkit
from scipy.interpolate import interp1d
from scipy import signal
from scipy.signal.windows import tukey

import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 10
from matplotlib.ticker import (MultipleLocator)

from obspy.taup import TauPyModel
from ellipticipy import ellipticity_correction
import obspy.geodetics

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
class CorrelateTwin:
    """
    Class to handle calculation of travel time delays either by
    cross correlation or Zaroli et al., (2010) method - which is much slower
    and perhaps not much more accurate
    """
    tk = Toolkit()

    def __init__(self):
        pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def process_one_event(self, input_dict):
        """
        Processes one event (earthquake), reading the file, 
        applying the list of functions to them

        :param input_dict: dictionary as below
        :type  input_dict: dict
        :return input_dict: dictionary
        :type  input_dict: dict

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
        """
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

                # Find location of both data and twin files for obs and syn.

                tw_loc_obs = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['phase_a_obs_out_loc']}{params_in['component']}'
                twin_in_obs = f'{tw_loc_obs}/{str(id_ctm)}.{params_in['phase_a_obs_out_loc'][-2:]}{params_in['component']}.twin'

                tw_loc_syn = f'{params_in['home']}/{params_in['twin_loc']}/{params_in['phase_a_syn_out_loc']}{params_in['component']}'
                twin_in_syn = f'{tw_loc_syn}/{str(id_ctm)}.{params_in['phase_a_syn_out_loc'][-2:]}{params_in['component']}.twin'

                if os.path.isfile(twin_in_obs) and os.path.isfile(twin_in_syn):
                    # read twin file
                    twin_df_obs = self.tk.read_twin_file(twin_in_obs)
                    twin_df_syn = self.tk.read_twin_file(twin_in_syn)

                    # Check df is not empty and then proceed:
                    if not twin_df_obs.empty and not twin_df_syn.empty:
                    
                        # Initiate outfile
                        outfile = self.open_outfile_file(input_dict)
                        outfile = self.write_params_outfile(input_dict, outfile)

                        # Counting stats
                        input_dict['num_files_in'] += 1
                        input_dict['num_obj_in'] = np.min([len(twin_df_obs), len(twin_df_syn)])

                        # Initiate empty df of correlated time windows
                        corr_df = pd.DataFrame(columns = params_in['correlate_outcols'])

                        ###
                        self.tk.print_log(params_in, logfile, f'----------////    WORKING ON: {str(id_ctm)}    ////----------')
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
                                input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail = function(input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail)
                            
                        ###    
                        self.tk.print_log(params_in, logfile, f'----------////    FINISHED    ////----------')
                        ###

                        outfile.close()

                    else:
                        # Returned empty dataframe
                        ###
                        self.tk.print_log(params_in, logfile, f'----------////            NO-TWINS-IN:            ////----------')
                        self.tk.print_log(params_in, logfile, f'----------////      {str(twin_df_obs)}          ////----------')
                        self.tk.print_log(params_in, logfile, f'----------////                  OR:               ////----------')
                        self.tk.print_log(params_in, logfile, f'----------////      {str(twin_df_syn)}          ////----------')
                        ###
                        pass

                else:
                    # No corresponding twin file found at twin_in
                    ###
                    self.tk.print_log(params_in, logfile, f'----------////        NO-TWIN-FILE-AT:            ////----------')
                    self.tk.print_log(params_in, logfile, f'----------////      {str(twin_in_obs)}          ////----------')
                    self.tk.print_log(params_in, logfile, f'----------////                  OR:               ////----------')
                    self.tk.print_log(params_in, logfile, f'----------////      {str(twin_in_syn)}          ////----------')
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
        """
        write key params to logfile

        :param input_dict: loaded input parameters
        :type input_dict: dict
        :param logfile: open logfile to write params
        :type logfile: open txt file
        :return logfile: open logfile
        :type logfile: open txt file
        """
        justify = 30

        params_in = input_dict['params_in']
        phases = input_dict['phases']

        logfile.write(' ')
        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('code start',' : ',str(params_in['code_start_time']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('id_fmt_ctm',' : ',str(input_dict['id_fmt_ctm']), x = justify) )

        params_list = ['obs_loc', 'syn_loc', 'component', 'T1', 'Tc', 'T2','sig_win_ext', 'sig_win_type', 'min_snr_P', 'min_snr_A', 'npow']
        for k, param in enumerate(params_list):
            logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params_in[param]), x = justify) )

        # Complex parameters where multiplication factors are used...  
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('movmax window [s]',' : ',str(params_in['mxf_win_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('walk away [s]',' : ',str(params_in['walkaway_f']*params_in['Tc']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min signal window length [s]',' : ',str(params_in['min_sig_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min noise window length [s]',' : ',str(params_in['min_nois_win_f']*params_in['T2']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min window size [s]',' : ', str(max(params_in['min_win_span_f'][0]*params_in['T1'],params_in['min_win_span_f'][1])), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taup phases',' : ',str(phases[str(params_in['phases_key'])][str(params_in['component'])]), x = justify) )
        logfile.write('')

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time shift [s]',' : ',str(params_in['max_tshift']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('interpolate delta',' : ',str(params_in['interp_delta']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('use velocity',' : ',str(params_in['use_velocity']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('normalise traces',' : ',str(params_in['normalise_traces']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taper windows',' : ',str(params_in['taper_windows']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('taper alpha',' : ',str(params_in['taper_alpha']), x = justify) )

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max tpeak delay',' : ',str(params_in['max_tpeak_delay']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max twin len diff',' : ',str(params_in['max_twin_len_diff']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('limit corr window',' : ',str(params_in['limit_corr_window']), x = justify) )

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max XC lag',' : ',str(params_in['XC_max_lag']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min XC peak amp',' : ',str(params_in['XC_min_peak_search']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min XC main peak amp',' : ',str(params_in['XC_min_main_peak']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('2nd XC peak max amp.',' : ',str(params_in['XC_max_secondary_peak_percentage']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min XC peak separation',' : ',str(params_in['XC_min_peak_distance']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('f spec max diff',' : ',str(params_in['tw_max_f_diff']), x = justify) )

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('use Zaroli 2010 EQ8',' : ',str(params_in['Zaroli']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('XC Zaroli TT max diff [s]',' : ',str(params_in['XC_Zaroli_diff']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min F3 peak amp',' : ',str(params_in['F3_min_peak_search']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min F3 main peak amp',' : ',str(params_in['F3_min_main_peak']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('2nd F3 peak max amp.',' : ',str(params_in['F3_max_secondary_peak_percentage']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min F3 peak separation',' : ',str(params_in['F3_min_peak_distance']), x = justify) )

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Use Ellipticity corr.',' : ',str(params_in['ellip_corr']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Use Station elev corr.',' : ',str(params_in['stel_corr']), x = justify) )

        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('plot correlation pic',' : ',str(params_in['corr_plot_pic']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('save correlation pic',' : ',str(params_in['corr_save_pic']), x = justify) )
        
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time delay [s]',' : ',str(params_in['max_tdelay']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time delay error [s]',' : ',str(params_in['max_tderr']), x = justify) )
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min Cross Correl. Coeff',' : ',str(params_in['min_xcc']), x = justify) )

        logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_params_outfile(self, input_dict, outfile):
        """
        write key params to outfile

        :param input_dict: loaded input parameters
        :type input_dict: dict
        :param outfile: open outfile to write code output
        :type outfile: open txt file
        :return outfile: open outfile
        :type outfile: open txt file  
        """
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
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('code start',' : ',str(params_in['code_start_time']), x = justify) )
        outfile.write('----------\n')
        outfile.write('----------////              CORRELATE TWIN PARAMETERS               ////----------\n')
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Channels',' : ', f'{params_in['phase_a_obs_out_loc'][-2:]}{params_in['component']}-{params_in['phase_a_syn_out_loc'][-2:]}{params_in['component']}', x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('period band [s]',' : ', f'[{str(params_in['T1'])}, {str(params_in['T2'])}]', x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('interpolate delta',' : ',str(params_in['interp_delta']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Use Zaroli F3',' : ',str(params_in['Zaroli']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Use Ellipticity corr.',' : ',str(params_in['ellip_corr']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('Use Station elev corr.',' : ',str(params_in['stel_corr']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time delay [s]',' : ',str(params_in['max_tdelay']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('max time delay error [s]',' : ',str(params_in['max_tderr']), x = justify) )
        outfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('min Cross Correl. Coeff',' : ',str(params_in['min_xcc']), x = justify) )
        outfile.write('----------\n')

        outfile.write('{0:<20s} {1:s} {2:s}\n'.format('columns format',' : ',params_in['correlate_outfmt']))
        for h,header in enumerate(params_in['correlate_outcols']):
            if h < len(params_in['correlate_outcols']) - 1:
                outfile.write('{0:s}'.format(f'{header},'))
            else:
                outfile.write('{0:s}'.format(f'{header}\n'))
        
        return outfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def correlate_windows(self, input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail):
        """
        Correlate time windows for arrival time residual determination
        Based on Cross correlation (XC) or Zaroli et al (2010) GJI (F3)

        Approximate arrival time errors using Chevrot (2002)
        Optional ellipticity correction based on Russell et al., (2022)
        
        ----------
        Reads and returns same inputs/outputs
        ** param == return **

        :param input_dict: event details for processing 
        :type input_dict: dict
        :param twin_in_obs: obs twin file location string
        :type twin_in_obs: str
        :param twin_df_obs: obs twin dataframe object
        :type twin_df_obs: pd.df
        :param twin_in_syn: syn twin file location string
        :type twin_in_syn: str
        :param twin_df_syn: syn twin dataframe object
        :type twin_df_syn: pd.df
        :param corr_df: correlated twin dataframe object
        :type corr_df: pd.df
        :param logfile: open logfile to report progress
        :type: logfile: open textfile for writing
        :param outfile: open outfile to report results
        :type: outfile: open textfile for writing
        :param fail: fail flag = 1 if function fails, else 0
        :type fail: int
        """
        params_in = input_dict['params_in']
        np.seterr(divide = 'ignore', invalid = 'ignore')

        # setup logfile statement using event_id, filename and function name (with inspect)
        log_statement = f'{self.tk.get_log_statement(input_dict['event_id'], twin_in_obs)}, {str(inspect.stack()[0][3])}'

        if fail:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
        else:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
            ###

            # data loc / synth loc
            input_directory_obs = f'{params_in['obs_loc']}/e{str(params_in['year'])}{params_in['fmt_data_loc']}'
            input_directory_syn = f'{params_in['syn_loc']}/e{str(params_in['year'])}{params_in['fmt_data_loc']}'

            # find unique stations in obs & syn df.
            unique_stations_obs = list(twin_df_obs['nslc'].unique())
            nst_obs = len(unique_stations_obs)

            unique_stations_syn = list(twin_df_syn['nslc'].unique())
            nst_syn = len(unique_stations_syn)

            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  obs: {twin_df_obs.shape[0]} windows for {nst_obs} stations')
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  syn: {twin_df_syn.shape[0]} windows for {nst_syn} stations')
            ###
            
            unique_stations = sorted(list(set(unique_stations_obs).intersection(unique_stations_syn)))
            nst = len(unique_stations)

            if nst > 0:
                # Loop through stations occuring in both dataframes
                for s, station in enumerate(unique_stations):
                    # Find the indices where nslc == station
                    stat_df_obs = twin_df_obs.query("nslc == @station")
                    len_stat_df_obs = len(stat_df_obs)
                    stat_df_syn = twin_df_syn.query("nslc == @station")
                    len_stat_df_syn = len(stat_df_syn)

                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  {len_stat_df_obs} obs windows for station {station} [{str(s + 1)}/{str(nst)}]')
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  {len_stat_df_syn} syn windows for station {station} [{str(s + 1)}/{str(nst)}]')
                    ###

                    # read sac files
                    file_obs = f'{input_directory_obs}/{str(input_dict['id_fmt_ctm'])}/{str(input_dict['id_fmt_ctm'])}_{str(station)}'
                    file_syn = f'{input_directory_syn}/{str(input_dict['id_fmt_ctm'])}/{str(input_dict['id_fmt_ctm'])}_{str(station)}'

                    try:
                        seis_obs = read(file_obs,'sac')
                        seis_syn = read(file_syn,'sac')
                    except:
                        ###
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------////     FAILED_LOADING SAC FILES    ////----------\n')
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------////       FAILED_LOADING.....       ////----------\n')
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------////          {twin_df_obs}          ////----------\n')
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------////               OR:               ////----------\n')
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------////          {twin_df_syn}          ////----------\n')
                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  CONTINUE...\n')
                        ###
                        continue

                    nslc = f'{str(seis_obs[0].stats['network'])}_{str(seis_obs[0].stats['station'])}.{str(seis_obs[0].stats['location'])}.{str(seis_obs[0].stats['channel'])}'

                    stel = float(seis_obs[0].stats['sac']['stel']) # in meters
                    stla = float(seis_obs[0].stats['sac']['stla'])
                    stlo = float(seis_obs[0].stats['sac']['stlo'])

                    if params_in['use_velocity']:  
                        seis_obs.differentiate()
                        seis_syn.differentiate()

                    # t_obs = seis_obs[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))
                    # dt_obs = seis_obs[0].stats['delta']
                    # seis_obs.filter('bandpass',freqmin = 1/params_in['T2']*2*dt_obs,freqmax = 1/params_in['T1']*2*dt_obs,corners = 4,zerophase = True)
                    seis_obs.filter('bandpass',freqmin = 1/params_in['T2'],freqmax = 1/params_in['T1'],corners = 4,zerophase = True)

                    seis_obs.detrend()
                    seis_obs.taper(max_percentage = 0.1)

                    # t_syn = seis_syn[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))
                    # dt_syn = seis_syn[0].stats['delta']
                    # seis_syn.filter('bandpass',freqmin = 1/params_in['T2']*2*dt_syn,freqmax = 1/params_in['T1']*2*dt_syn,corners = 4,zerophase = True)
                    seis_syn.filter('bandpass',freqmin = 1/params_in['T2'],freqmax = 1/params_in['T1'],corners = 4,zerophase = True)

                    seis_syn.detrend()
                    seis_syn.taper(max_percentage = 0.1)

                    if params_in['normalise_traces']:
                        # Normalising traces.. Shouldnt affect normalised cross-correlation co-efficients

                        seis_obs.normalize()            
                        seis_syn.normalize()
                        
                    # Interpolate to common fs
                    
                    seis_obs_start = seis_obs[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))[0]
                    seis_obs_end = seis_obs[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))[-1]
                    interp_function = interp1d(seis_obs[0].times(reftime = UTCDateTime(str(input_dict['evtm']))), seis_obs[0].data, kind = 'linear', fill_value = 0)
                    try:
                        t_obs = np.arange(0, seis_obs_end, params_in['interp_delta'])
                        s_obs = interp_function(t_obs)
                    except:
                        t_obs = np.arange(seis_obs_start, seis_obs_end, params_in['interp_delta'])
                        s_obs = interp_function(t_obs)


                    seis_syn_start = seis_syn[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))[0]
                    seis_syn_end = seis_syn[0].times(reftime = UTCDateTime(str(input_dict['evtm'])))[-1]
                    interp_function = interp1d(seis_syn[0].times(reftime = UTCDateTime(str(input_dict['evtm']))), seis_syn[0].data, kind = 'linear', fill_value = 0)
                    try:
                        t_syn = np.arange(0, seis_syn_end, params_in['interp_delta'])
                        s_syn = interp_function(t_syn)
                    except:
                        t_syn = np.arange(seis_syn_start, seis_syn_end, params_in['interp_delta'])
                        s_syn = interp_function(t_syn)

                    for index_obs, row_obs in stat_df_obs.iterrows():
                        r_phase = row_obs['phase']
                        r_t_taup = row_obs['t_taup']
                        itw_syn_1 = stat_df_syn.query("phase == @r_phase")
                        itw_syn_2 = stat_df_syn.query("t_taup == @r_t_taup")

                        if not itw_syn_1.equals(itw_syn_2):
                            # Some issue with the dfs - not looking at same phase and arrival time
                            continue
                        else:
                        
                            if len(itw_syn_1) == 0:
                                ###
                                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  No correspondence for {row_obs['phase']} @ {station}')
                                ###
                                continue
                            else:

                                # index for current obs window and take autocorrelation
                                i_obs, x_obs = self.tk.subvec(t_obs, s_obs, [row_obs['t_left'], row_obs['t_right']])


                                for index_syn, row_syn in itw_syn_1.iterrows():
                                    if np.abs(row_obs['t_peak'] - row_syn['t_peak']) > params_in['max_tpeak_delay']:
                                        # greater than max_tpeak_delay between envelope peaks - skip
                                            # DOES NOT ensure that correlated windows have delay times less than max_tpeak_delay
                                            
                                        ###
                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - distant correspondence {row_obs['phase']} @ {station}')
                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  abs( obs_t_peak - syn_t_peak ) > max_tpeak_delay; np.abs({row_obs['t_peak']} - {row_syn['t_peak']}) > {params_in['max_tpeak_delay']}')
                                        ###

                                        continue

                                    elif np.abs((row_syn['t_right'] - row_syn['t_left']) - (row_obs['t_right'] - row_obs['t_left'])) > params_in['max_twin_len_diff'] :
                                        # Time window lengths are too different...

                                        ###
                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - time window length difference greater than {params_in['max_twin_len_diff']}s {row_obs['phase']} @ {station}')
                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Obs twin len {row_obs['t_right'] - row_obs['t_left']}s >> OR << Syn twin len: {row_syn['t_right'] - row_syn['t_left']}s')
                                        ###

                                        continue

                                    else:
                                        # index for current syn window and take autocorrelation
                                        i_syn, x_syn = self.tk.subvec(t_syn, s_syn, [row_syn['t_left'], row_syn['t_right']])
                                        lags_auto_syn = signal.correlation_lags(x_syn.size, x_syn.size, mode = "full")  * params_in['interp_delta']

                                        if params_in['normalise_traces']:
                                            #  Normalising time windows... Shouldnt affect normalised cross-correlation co-efficients

                                            x_syn = x_syn / np.max(np.abs(x_syn))
                                            x_obs = x_obs / np.max(np.abs(x_obs))

                                        # Taper the windows prior to cross correlation
                                        if params_in['taper_windows']:
                                            # Tapering time windows using params_in['taper_alpha']... Shouldnt affect normalised cross-correlation co-efficients

                                            x_syn = x_syn * tukey(len(x_syn), params_in['taper_alpha'])
                                            x_obs = x_obs * tukey(len(x_obs), params_in['taper_alpha'])

                                        # Get autocorrelations
                                        auto_obs = signal.correlate(x_obs, x_obs, mode = "full", method = "fft")
                                        auto_syn = signal.correlate(x_syn, x_syn, mode = "full", method = "fft")
                                        # auto_syn_max_pos = np.argmax(auto_syn)

                                        # corrdelay estimates
                                        correlation = signal.correlate(x_obs, x_syn, mode = "full", method = "fft")
                                        lags = signal.correlation_lags(x_obs.size, x_syn.size, mode = "full") * params_in['interp_delta']
                                        
                                        # if params_in['limit_corr_window']:

                                        #     i_corr_lim, corr_lim = self.tk.subvec(lags, correlation, [-1 * params_in['XC_max_lag'], params_in['XC_max_lag']])
                                        #     lags_lim = lags[i_corr_lim]
                                        #     correlation = corr_lim
                                        #     lags = lags_lim
                                        #     # If max of auto correlation is always at zero no need to cut autocorrelation windows

                                        i_ccmx = np.argmax(correlation)
                                        lag = lags[i_ccmx]

                                        if i_ccmx == 0 or i_ccmx == len(correlation):
                                            # We have picked somehting at the end of the window - so we shouldnt trust it...
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - XC maximum position error (lags extremities) {row_obs['phase']} @ {station}')
                                            ###
                                            continue

                                        if lag < -1*params_in['XC_max_lag'] or lag > params_in['XC_max_lag']:
                                            # Got some problem where the max corr-coeff lag is very large
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - XC maximum position error (XC_max_lag) {row_obs['phase']} @ {station}')
                                            ###
                                            continue


                                        # Normalise correlation function
                                        norm_xc = np.round(np.sqrt((np.maximum(correlation, 0) / np.max(auto_syn)) * (np.maximum(correlation, 0) /  np.max(auto_obs))),4)
                                        i_ccmx = np.argmax(norm_xc)

                                        # Check peaks in Normed XC function
                                        XC_peaks_ind, XC_peak_dicts = signal.find_peaks(norm_xc, height = params_in['XC_min_peak_search'], distance = params_in['XC_min_peak_distance']) # Find peaks in XC
                                        XC_sort_args = np.argsort(XC_peak_dicts['peak_heights'])[::-1] # Decreasing order.

                                        if len(XC_peaks_ind) == 0:
                                            # Reject - main peak too weak
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - All XC peaks below {params_in['XC_min_peak_search']} threshold {row_obs['phase']} @ {station}')
                                            ###
                                            continue
                                        else:
                                            if XC_peak_dicts['peak_heights'][XC_sort_args[0]] <  params_in['XC_min_main_peak']: 
                                                # Reject - main peak too weak
                                                ###
                                                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - XC Main peak amplitude {XC_peak_dicts['peak_heights'][XC_sort_args[0]]:.3f} below {params_in['XC_min_main_peak']} threshold {row_obs['phase']} @ {station}')
                                                ###
                                                continue
                                            else:
                                                # Central peak is strong enough, now check secondaries
                                                if len(XC_peaks_ind) > 1:
                                                    if XC_peak_dicts['peak_heights'][XC_sort_args[1]] > XC_peak_dicts['peak_heights'][XC_sort_args[0]] * params_in['XC_max_secondary_peak_percentage']:
                                                        # Reject - main peak indistinct from secondary
                                                        ###
                                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - XC main peak amplitude {XC_peak_dicts['peak_heights'][XC_sort_args[0]]:.3f} indistinct from secondary {XC_peak_dicts['peak_heights'][XC_sort_args[1]]:.3f} {row_obs['phase']} @ {station}')
                                                        ###
                                                        continue
                                                    else:
                                                        # Strong main peak, weaker secondary peak
                                                        ###
                                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Accept XC {XC_peak_dicts['peak_heights'][XC_sort_args[0]]:.3f} peak at {lags[XC_peaks_ind[XC_sort_args[0]]] + t_obs[i_obs][0] - t_syn[i_syn][0]:.3f}s, with weaker secondaries {row_obs['phase']} @ {station}')
                                                        ###
                                                else:
                                                    # Accept 1 peak
                                                    ###
                                                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Accept XC {XC_peak_dicts['peak_heights'][XC_sort_args[0]]:.3f} peak at {lags[XC_peaks_ind[XC_sort_args[0]]] + t_obs[i_obs][0] - t_syn[i_syn][0]:.3f}s {row_obs['phase']} @ {station}')
                                                    ###


                                        i_ccmx_norm = XC_peaks_ind[XC_sort_args[0]] # or np.argmax(norm_xc) Index of max XC coefficient
                                        ccmx = norm_xc[i_ccmx_norm] # Maximum cross correlation coefficent

                                        if i_ccmx !=  i_ccmx_norm:
                                            # Somehow normalization went wrong
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - Norm XC maximum position error {row_obs['phase']} @ {station}')
                                            ###
                                            continue

                                        # Calculate frequency spectra of time windows.
                                        # obs
                                        ft_obs = np.fft.rfft(x_obs)
                                        freqs_obs = np.fft.rfftfreq(len(x_obs), t_obs[i_obs][1]-t_obs[i_obs][0]) # Get frequency axis from the time axis
                                        mags_obs = (np.abs(ft_obs)**2) * (2 / (len(x_obs)**2)) # We don't care about the phase information here
                                        max_mags_obs = freqs_obs[np.argmax(mags_obs)]
                                        # Syn
                                        ft_syn = np.fft.rfft(x_syn)
                                        freqs_syn = np.fft.rfftfreq(len(x_syn), t_syn[i_syn][1]-t_syn[i_syn][0]) # Get frequency axis from the time axis
                                        mags_syn = (np.abs(ft_syn)**2) * (2 / (len(x_syn)**2)) # We don't care about the phase information here
                                        max_mags_syn = freqs_syn[np.argmax(mags_syn)]

                                        if max_mags_obs > max_mags_syn * (1 + params_in['tw_max_f_diff']) or max_mags_obs < max_mags_syn * (1 - params_in['tw_max_f_diff']):
                                            # Somehow wildly different peak frequency contents in the windows.
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - Frequency spectra max diff too great {row_obs['phase']} @ {station}')
                                            ###
                                            continue

                                        # auto_obs_max_pos = np.argmax(auto_obs) # COMMENT
                                        # lags_auto_obs = signal.correlation_lags(x_obs.size, x_obs.size, mode = "full")  * params_in['interp_delta'] # COMMENT
                                        
                                        # print(f'len(i_obs): {len(t_obs[i_obs])}, len(t_syn[i_syn]): {len(t_syn[i_syn])}: np.abs(len(t_syn[i_syn]) - len(t_obs[i_obs])): {np.abs(len(t_syn[i_syn]) - len(t_obs[i_obs]))}')  # COMMENT


                                        # print(f'len(auto_syn): {len(auto_syn)}, np.argmax(auto_syn): {auto_syn_max_pos}, max auto_syn: {auto_syn[auto_syn_max_pos]:.3f} at {lags_auto_syn[auto_syn_max_pos]:.3f}s')  # COMMENT
                                        # print(f'len(auto_obs): {len(auto_obs)}, np.argmax(auto_obs): {auto_obs_max_pos}, max auto_obs: {auto_obs[auto_obs_max_pos]:.3f} at {lags_auto_obs[auto_obs_max_pos]:.3f}s') # COMMENT


                                        start_ind_shift = int((t_obs[i_obs][0] - t_syn[i_syn][0])/params_in['interp_delta'])
                                        # print(f'start_ind_shift = int((t_obs[0] - t_syn[0])/delta_t): {start_ind_shift} = int (({t_obs[i_obs][0]} - {t_syn[i_syn][0]}) / {params_in['interp_delta']}))') # COMMENT

                                        # Time delay based on XC
                                        XC_tdl = np.round(lag + t_obs[i_obs][0] - t_syn[i_syn][0], 4)
                                        # XC_ind_shift = auto_syn_max_pos - i_ccmx - start_ind_shift
                                        # print(f'XC_ind_shift = auto_syn_max_pos - i_ccmx - start_ind_shift: {XC_ind_shift:d} = {auto_syn_max_pos:d} - {i_ccmx:d} - {start_ind_shift:d}') # COMMENT

                                        t_min = np.min([t_syn[0], t_obs[0]])
                                        t_max = np.max([t_syn[-1], t_obs[-1]])
                                        ax_time = np.arange(t_min, t_max, params_in['interp_delta'])

                                        
                                        if params_in['limit_corr_window']:
                                            # This helps to speed up Zaroli.
                                            t_max = lags[-1] - start_ind_shift * params_in['interp_delta']
                                            t_min = lags[0] - start_ind_shift * params_in['interp_delta']
                                            ax_time_F3 = np.arange(t_min, t_max, params_in['interp_delta'])
                                        else:
                                            ax_time_F3 = ax_time


                                        # Zaroli time delay calculation
                                        if params_in['Zaroli']:
                                            ###
                                            # self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Using Zaroli Delay time calculation...')
                                            ###
                                            # stz = time.time()
                                            F3_output = self.Func_3(x_obs, x_syn, t_obs[i_obs], t_syn[i_syn], ax_time, ax_time_F3, params_in['interp_delta']) # Roughly 1s for the calculation
                                            # print("Got Zaroli F3 in --- %s seconds ---" % (time.time() - stz))

                                            # Check peaks in F3 function
                                            F3_peaks_ind, F3_peak_dicts = signal.find_peaks(F3_output, height = params_in['F3_min_peak_search'], distance = params_in['F3_min_peak_distance']) # Find peaks in F3
                                            F3_sort_args = np.argsort(F3_peak_dicts['peak_heights'])[::-1] # Decreasing order.

                                            if len(F3_peaks_ind) == 0:
                                                # Reject - main peak too weak
                                                ###
                                                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - All F3 peaks below {params_in['F3_min_peak_search']} threshold {row_obs['phase']} @ {station}')
                                                ###
                                                continue
                                            else:
                                                if F3_peak_dicts['peak_heights'][F3_sort_args[0]] <  params_in['F3_min_main_peak']: 
                                                    # Reject - main peak too weak
                                                    ###
                                                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - F3 Main peak amplitude {F3_peak_dicts['peak_heights'][F3_sort_args[0]]:.3f} below {params_in['F3_min_main_peak']} threshold {row_obs['phase']} @ {station}')
                                                    ###
                                                    continue
                                                else:
                                                    # Central peak is strong enough, now check secondaries
                                                    if len(F3_peaks_ind) > 1:
                                                        if F3_peak_dicts['peak_heights'][F3_sort_args[1]] > F3_peak_dicts['peak_heights'][F3_sort_args[0]] * params_in['F3_max_secondary_peak_percentage']:
                                                            # Reject - main peak indistinct from secondary
                                                            ###
                                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - F3 main peak amplitude {F3_peak_dicts['peak_heights'][F3_sort_args[0]]:.3f} indistinct from secondary {F3_peak_dicts['peak_heights'][F3_sort_args[1]]:.3f} {row_obs['phase']} @ {station}')
                                                            ###
                                                            continue
                                                        else:
                                                            # Strong main peak, weaker secondary peak
                                                            ###
                                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Accept F3 {F3_peak_dicts['peak_heights'][F3_sort_args[0]]:.3f} peak at {ax_time_F3[F3_peaks_ind[F3_sort_args[0]]]:.3f}s, with weaker secondaries {row_obs['phase']} @ {station}')
                                                            ###
                                                    else:
                                                        # Accept 1 peak
                                                        ###
                                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Accept F3 {F3_peak_dicts['peak_heights'][F3_sort_args[0]]:.3f} peak at {ax_time_F3[F3_peaks_ind[F3_sort_args[0]]]:.3f}s {row_obs['phase']} @ {station}')
                                                        ###

                                            # F3_peak_ind = np.argmax(F3_output)        # Index of highest peak in F3
                                            F3_peak_ind = F3_peaks_ind[F3_sort_args[0]] # Index of highest peak in F3

                                            # Get time delay (tdl) and ind_shift (for plotting) and equivalent of Normed Cross correlation coeffecient from F3
                                            tdl = np.round(ax_time_F3[F3_peak_ind], 4)
                                            # ax_time_zero = np.argmin(np.abs(ax_time_F3))
                                            # ind_shift = ax_time_zero - F3_peak_ind
                                            # print(f'ind_shift_Z = ax_time_zero - F3_peak_ind: {ind_shift:d} = {ax_time_zero:d} - {F3_peak_ind:d}') # COMMENT

                                            # print(f'ZR ALIGNED TRACES IND_SHIFT: {ind_shift}') # COMMENT
                                            # print(f'XC ALIGNED TRACES IND_SHIFT: {XC_ind_shift}') # COMMENT

                                            ccmx = F3_output[F3_peak_ind]
                                        else:
                                            # Just use cross correlation output.
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Using XC based Delay time calculation...')
                                            ###
                                                                                    
                                            tdl = XC_tdl
                                            Z_tdl = []
                                            F3_output = []
                                            F3_peaks_ind = []
                                            F3_peak_dicts = {}
                                            F3_sort_args = []
                                            F3_peak_ind = []


                                        if np.abs(XC_tdl - tdl) > params_in['XC_Zaroli_diff']:
                                            # Large difference between XC and Zaroli delay time picks...
                                            # tdl == XC_tdl when using XC calculation
                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject - Delay time pick discrepancy {row_obs['phase']} @ {station}')
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  tdl = {tdl:.2f}, XC tdl = {XC_tdl:.2f}')
                                            ###
                                            continue


                                        #  Use Chevrot (2002) to find approx error for measrured time residual using atuo correlation of synthetic
                                        try:
                                            if np.max(auto_syn) > np.max(correlation):
                                                # # Minimum positive lag where auto_syn < max(correlation)
                                                lags_p = lags_auto_syn[lags_auto_syn>0]
                                                tdl_err = np.min(lags_p[np.where(auto_syn[lags_auto_syn>0] < np.max(correlation))[0]]) #  only accurate to params_in['interp_delta']
                                            else:
                                                ###
                                                self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Error np.max(auto_syn) < np.max(correlation) {np.max(auto_syn):.3f} < {np.max(correlation):.3f}')
                                                ###
                                                tdl_err = -999.0
                                        except:
                                            tdl_err = -999.0

                                        if params_in['ellip_corr'] or params_in['stel_corr']:
                                            # compute distances azimuth and backazimuth
                                            # and arrival information
                                        
                                            distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(float(input_dict['evla']), 
                                                                                                   float(input_dict['evlo']), 
                                                                                                   float(stla), float(stlo))
                                            distdg = distm / (6371.e3 * np.pi / 180.)
                                            model = TauPyModel(model = params_in['taup_model_name'])

                                            arrival = model.get_ray_paths(
                                                source_depth_in_km=input_dict['evdp'], 
                                                distance_in_degree=distdg, 
                                                phase_list=[row_obs['phase']])                   

                                        # Ellipticity correction Russell et al., (2022)
                                        if params_in['ellip_corr']:
                                            # Often not required if specfem option turned on.
                                            correction = ellipticity_correction(arrival, 
                                                                                az, 
                                                                                input_dict['evla'])
                                            corr_tdl = tdl - correction[0]

                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Use ellipticity correction: {corr_tdl:.3f} = {tdl:.3f} - {correction[0]:.3f}')
                                            ###
                                            tdl = corr_tdl

                                        # Compute an elevation correction
                                        if params_in['stel_corr']:
                                            # Often not required if specfem option turned on.
                                            # Define the parameters from crust2.0
                                            r = 6371  # Radius of the Earth in km
                                            d1 = 1
                                            d2 = 1
                                            V1 = 3.5
                                            V2 = 4.1
                                            V3 = 6.1
                                            V4 = 5.8

                                            if stel > 0 and stel <= 1000:
                                                tt = stel / 1000 / V1
                                            elif stel > 1000 and stel <= 2000:
                                                tt = 1 / V1 + (stel - 1000) / 1000 / V2
                                            elif stel > 2000:
                                                tt = 1 / V1 + 1 / V2 + (stel - 2000) / 1000 / V3
                                            elif stel < -10:
                                                tt = stel / 1000 / V4
                                            else:
                                                tt = 0
                                            
                                            theta = np.arcsin(arrival[0].ray_param * V1 / r)

                                            correction = tt / np.cos(theta)

                                            corr_tdl = tdl - correction

                                            ###
                                            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Use station elevation correction at {stel:.3f}m: {corr_tdl:.3f} = {tdl:.3f} - {correction:.3f}')
                                            ###
                                            tdl = corr_tdl


                                        # Execute plotting for correlated twin
                                        if params_in['corr_plot_pic']:
                                            self.corr_twin_plot(input_dict, params_in, logfile, log_statement, nslc, 
                                                row_obs, t_obs, x_obs, i_obs, freqs_obs, mags_obs, max_mags_obs, 
                                                row_syn, t_syn, x_syn, i_syn, freqs_syn, mags_syn, max_mags_syn,
                                                auto_syn, lags_auto_syn, correlation, XC_tdl, tdl,
                                                ax_time, ax_time_F3, start_ind_shift, norm_xc, lags, i_ccmx_norm,
                                                F3_output, F3_peaks_ind, F3_peak_dicts, F3_sort_args, F3_peak_ind,
                                                tdl_err)

                                        # Log final result
                                        self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  nslc,stla,stlo,stel,row_obs[phase],tdl,tdl_err,ccmx,row_obs[t_taup],row_obs[t_peak],row_syn[t_peak]')
                                        self.tk.print_log(params_in, logfile, '{0:s}  ,  {1:s},{2:.4f},{3:.4f},{4:.4f},{5:s},{6:.2f},{7:.2f},{8:.2f},{9:.2f},{10:.2f},{11:.2f}'.format(log_statement,nslc,stla,stlo,stel,row_obs['phase'],tdl,tdl_err,ccmx,row_obs['t_taup'],row_obs['t_peak'],row_syn['t_peak']))
                                        
                                        # Add result to dataframe
                                        # nslc,stla,stlo,stel,phase,tdelay,tderr,ccmx,ttaup,tp_obs,tp_syn,Ap_obs,Ap_syn
                                        lt = [nslc,stla,stlo,stel,row_obs['phase'],tdl,tdl_err,ccmx,row_obs['t_taup'],row_obs['t_peak'],row_syn['t_peak'],row_obs['A_peak'],row_syn['A_peak']]

                                        corr_df.loc[len(corr_df)] = lt

        return input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def corr_twin_plot(self, input_dict, params_in, logfile, log_statement, nslc, 
                row_obs, t_obs, x_obs, i_obs, freqs_obs, mags_obs, max_mags_obs, 
                row_syn, t_syn, x_syn, i_syn, freqs_syn, mags_syn, max_mags_syn,
                auto_syn, lags_auto_syn, correlation, XC_tdl, Z_tdl,
                ax_time, ax_time_F3, start_ind_shift, norm_xc, lags, i_ccmx_norm,
                F3_output, F3_peaks_ind, F3_peak_dicts, F3_sort_args, F3_peak_ind,
                tdl_err):
        """
        Produce 5 panel figure of correlated time windows
        (a) Raw obs & syn time windows
        (b) Aligned time windows
        (c) Frequency spectrum comparison
        (d) Normalised XC Function and Zaroli et al (2010) F3
        (e) Autocorrelation derived error (Chevrot, 2002)
        (f) displays of saves plot as PDF

        :param input_dict: event details for processing 
        :type input_dict: dict
        :param params_in: loaded parameters from yaml file
        :type params_in: dict
        :param logfile: logfile for writing
        :type logfile: open txt file
        :param log_statement: text string to write to logfile
        :type log_statement: str
        :param nslc: station name
        :type nslc: str
        :param row_obs: row from obs data twin dataframe
        :type row_obs: pd.df
        :param t_obs: time axis for obs twin
        :type t_obs: np.array
        :param x_obs: variable axis for obs twin
        :type x_obs: np.array
        :param i_obs: indicies of obs twin
        :type i_obs: np.array
        :param freqs_obs: obs twin spectrum x variable
        :type freqs_obs: np.array
        :param mags_obs: obs twin spectrum y variable
        :type mags_obs: np.array
        :param max_mags_obs: obs twin maximum frequency 
        :type max_mags_obs: float
        :param row_syn: row from syn data twin dataframe
        :type row_syn: pd.df
        :param t_syn: time axis for syn twin
        :type t_syn: np.array
        :param x_syn: variable axis for syn twin
        :type x_syn: np.array
        :param i_syn: indicies of syn twin
        :type i_syn: np.array
        :param freqs_syn: syn twin spectrum x variable
        :type freqs_syn: np.array
        :param mags_syn: syn twin spectrum y variable
        :type mags_syn: np.array
        :param max_mags_syn: syn twin maximum frequency
        :type max_mags_syn: float
        :param auto_syn: syn twin auto correlation function
        :type auto_syn: np.array
        :param lags_auto_syn: syn twin correlation lags
        :type lags_auto_syn: np.array
        :param correlation: obs-syn correlation function
        :type correlation: np.array
        :param XC_tdl: time delay from cross correlation
        :type XC_tdl: float
        :param Z_tdl: time delay from Zaroli et al., (2010) method
        :type Z_tdl: float
        :param ax_time: time axis for twin plotting
        :type ax_time: np.array
        :param ax_time_F3: time axis for calculation of Zaroli et al (2010) F3
        :type ax_time_F3: np.array
        :param start_ind_shift: initial trace time shift for cross correlation
        :type start_ind_shift: int
        :param norm_xc: normalized obs-syn cross correlation function
        :type norm_xc: np.array
        :param lags: correlation lags array
        :type lags: np.array
        :param i_ccmx_norm: Index of max XC coefficient
        :type i_ccmx_norm: int
        :param F3_output: Zaroli et al., (2010) F3 for obs-syn twins
        :type F3_output: np.array
        :param F3_peaks_ind: indicies fo peaks in F3
        :type F3_peaks_ind: np.array
        :param F3_peak_dicts: dict of peak properties in F3
        :type F3_peak_dicts: dict
        :param F3_sort_args: sorted peaks by magntiude
        :type F3_sort_args: np.array
        :param F3_peak_ind: index of peak in F3
        :type F3_peak_ind: int
        :param tdl_err: travel time delay error
        :type tdl_err: float
        """

        # Set up the plot
        fig, (ax_tw, ax_tw_aligned, ax_freq, ax_XC_Zr, ax_auto) = plt.subplots(5, 1, figsize = (8, 12))
        fig_title = f'Event {input_dict['id_cmt']}, Station {nslc}\n' \
                    f'{row_obs['phase']:s} @ {row_obs['t_taup']:.2f}s, tdl: {Z_tdl:.2f}s\n'
        fig.suptitle(fig_title, fontsize = 16)

        axes_list_all = [ax_tw, ax_tw_aligned, ax_freq, ax_XC_Zr, ax_auto]
        for j, ax in enumerate(axes_list_all):
            ax.spines["right"].set_linewidth(1.5)
            ax.spines["left"].set_linewidth(1.5)
            ax.spines["top"].set_linewidth(1.5)
            ax.spines["bottom"].set_linewidth(1.5)
            ax.tick_params(labelsize = 12)

        ############## Raw time windows. ##############

        obs, = ax_tw.plot(t_obs[i_obs],x_obs,'b')
        tp_obs, = ax_tw.plot([row_obs['t_peak'],row_obs['t_peak']],[-1*np.max(x_obs),np.max(x_obs)],'b--')
        syn, = ax_tw.plot(t_syn[i_syn],x_syn,'r')
        tp_syn, = ax_tw.plot([row_syn['t_peak'],row_syn['t_peak']],[-1*np.max(x_syn),np.max(x_syn)],'r--')
        ax_tw.legend([obs, tp_obs, syn, tp_syn], ['Obs', 'Tp Obs', 'Syn', 'Tp Syn'], loc = 'upper right', fontsize = 8)

        ax_tw.set_title('Synth & Obs windows', fontsize = 14)
        ax_tw.set_ylabel('Amplitude', fontsize = 14)
        ax_tw.set_xlabel('Time [s]', fontsize = 14)
        ax_tw.xaxis.set_minor_locator(MultipleLocator(5))
        ax_tw.xaxis.set_major_locator(MultipleLocator(20))
        ax_tw.set_xlim(np.min([t_obs[i_obs][0],t_syn[i_syn][0]])-5,np.max([t_obs[i_obs][-1],t_syn[i_syn][-1]]) + 5)
        ax_tw.annotate('a', (0, 1), xytext = (5,-5), xycoords = 'axes fraction', fontsize = 12, textcoords = 'offset points', color = 'k', backgroundcolor = 'none', ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))

        ############### Aligned traces ##############

        synth, = ax_tw_aligned.plot(t_syn[i_syn],x_syn,'r')

        # print(f'XC_tdl: {XC_tdl}')
        XC_obs_shifted = self.shift_and_truncate(x_obs, t_obs[i_obs], -XC_tdl/params_in['interp_delta'], ax_time, 1/params_in['interp_delta'])
        XC_obs, = ax_tw_aligned.plot(ax_time, XC_obs_shifted , 'b')
        
        if params_in['Zaroli']:
            # print(f'Z_tdl: {Z_tdl}')
            ZR_obs_shifted = self.shift_and_truncate(x_obs, t_obs[i_obs], -Z_tdl/params_in['interp_delta'], ax_time, 1/params_in['interp_delta'])
            ZR_obs, = ax_tw_aligned.plot(ax_time, ZR_obs_shifted , 'g')

            ax_tw_aligned.legend([ZR_obs, XC_obs, synth], ['ZR_obs','XC_obs','Syn'],loc = 'upper right',fontsize = 8)
        else:
            ax_tw_aligned.legend([obs, synth], ['Obs','Syn'],loc = 'upper right',fontsize = 8)

        ax_tw_aligned.set_title('Aligned Synth & Obs windows', fontsize = 14)
        ax_tw_aligned.set_ylabel('Amplitude', fontsize = 14)
        ax_tw_aligned.set_xlabel('Time [s]', fontsize = 14)
        ax_tw_aligned.set_xlim(np.min([t_obs[i_obs][0],t_syn[i_syn][0]])-5,np.max([t_obs[i_obs][-1],t_syn[i_syn][-1]]) + 5)
        ax_tw_aligned.xaxis.set_minor_locator(MultipleLocator(5))
        ax_tw_aligned.xaxis.set_major_locator(MultipleLocator(20))
        ax_tw_aligned.annotate('b', (0, 1), xytext = (5,-5), xycoords = 'axes fraction', fontsize = 12, textcoords = 'offset points', color = 'k', backgroundcolor = 'none', ha = 'left', va = 'top', bbox = dict(facecolor = 'white', edgecolor = 'black', pad = 2.0))

        ############# Frequency spectrums ##############

        # Obs
        obs, = ax_freq.plot(freqs_obs, mags_obs, 'b')
        ax_freq.plot([max_mags_obs, max_mags_obs], [0,np.max(mags_obs)], 'k--')
        # syn 
        syn, = ax_freq.plot(freqs_syn, mags_syn, 'r')
        ax_freq.plot([max_mags_syn, max_mags_syn], [0,np.max(mags_syn)], 'k--')

        ax_freq.set_xlim([0,0.3])
        ax_freq.set_title('Frequency Spectrum', fontsize = 14)
        ax_freq.set_ylabel('Amplitude', fontsize = 14)
        ax_freq.set_xlabel('Frequency [Hz]', fontsize = 14)
        ax_freq.xaxis.set_minor_locator(MultipleLocator(0.01))
        ax_freq.xaxis.set_major_locator(MultipleLocator(0.05))
        ax_freq.legend([obs, synth], ['Obs','Syn'],loc = 'upper right',fontsize = 8)
        ax_freq.annotate('c', (0, 1), xytext = (5,-5), xycoords = 'axes fraction', fontsize = 12, textcoords = 'offset points', color = 'k', backgroundcolor = 'none', ha = 'left', va = 'top', bbox = dict(facecolor = 'white', edgecolor = 'black', pad = 2.0))

        ################## Norm XC and Zaroli F3 ##################

        xc_norm_shifted = self.shift_and_truncate(norm_xc, lags, start_ind_shift, ax_time_F3, 1/params_in['interp_delta'])

        n_xc, = ax_XC_Zr.plot(ax_time_F3, xc_norm_shifted, 'k')
        n_xc_peak_time = lags[i_ccmx_norm]+(start_ind_shift * params_in['interp_delta'])
        n_xc_peak, = ax_XC_Zr.plot([n_xc_peak_time, n_xc_peak_time], [0,norm_xc[i_ccmx_norm]], 'k--')

        if params_in['Zaroli']:
            F3_l, = ax_XC_Zr.plot(ax_time_F3, F3_output, 'r')
            F3_peak, = ax_XC_Zr.plot([ax_time_F3[F3_peak_ind],ax_time_F3[F3_peak_ind]], [0,F3_output[F3_peak_ind]], 'r--', markersize = 10)
            ax_XC_Zr.set_title('Norm. Cross Corr. & Zaroli et al 2010 F3', fontsize = 14)
            for arg in F3_sort_args[1:]:
                F3_peaks, = ax_XC_Zr.plot(ax_time_F3[F3_peaks_ind[arg]], F3_peak_dicts['peak_heights'][arg], 'r+', markersize = 8)

            ax_XC_Zr.legend([n_xc, n_xc_peak, F3_l, F3_peak, F3_peaks], ['XC', 'XC peak', 'F3', 'F3 peak', 'F3 2nd peaks'],loc = 'upper right',fontsize = 8)

        else:
            ax_XC_Zr.legend([n_xc, n_xc_peak], ['XC', 'XC peak'],loc = 'upper right',fontsize = 8)
            ax_XC_Zr.set_title('Norm. Cross Corr. Function', fontsize = 14)

        # ax_XC_Zr.set_xlim([lags[0] + (start_ind_shift * params_in['interp_delta']), lags[-1] + (start_ind_shift * params_in['interp_delta'])])
        ax_XC_Zr.set_xlim([ax_time_F3[0], ax_time_F3[-1]])

        ax_XC_Zr.set_ylabel('Amplitude', fontsize = 14)
        ax_XC_Zr.set_xlabel('Lag [s]', fontsize = 14)
        ax_XC_Zr.xaxis.set_minor_locator(MultipleLocator(5))
        ax_XC_Zr.xaxis.set_major_locator(MultipleLocator(20))
        ax_XC_Zr.annotate('d', (0, 1), xytext = (5,-5), xycoords = 'axes fraction', fontsize = 12, textcoords = 'offset points', color = 'k', backgroundcolor = 'none', ha = 'left', va = 'top', bbox = dict(facecolor = 'white', edgecolor = 'black', pad = 2.0))

        ################## Auto Correlation based error ##################

        syn_auto, = ax_auto.plot(lags_auto_syn[lags_auto_syn>0], auto_syn[lags_auto_syn>0],'k')
        corr_max, = ax_auto.plot([lags_auto_syn[lags_auto_syn>0][0] , lags_auto_syn[lags_auto_syn>0][-1] ], [np.max(correlation) , np.max(correlation)],'r')
        # err_tdl, = ax_auto.plot(tdl_err, np.max(correlation), 'rx', markersize = 10)
        err_tdl, = ax_auto.plot([tdl_err, tdl_err], [np.min(auto_syn), np.max(auto_syn)], 'r--')

        ax_auto.set_title('Chevrot 2002 Error Estimate', fontsize = 14)
        ax_auto.set_ylabel('Amplitude', fontsize = 14)
        ax_auto.set_xlabel('Lag [s]', fontsize = 14)
        ax_auto.set_xlim([0, params_in['XC_max_lag']])
        ax_auto.xaxis.set_minor_locator(MultipleLocator(1))
        ax_auto.xaxis.set_major_locator(MultipleLocator(5))
        ax_auto.legend([syn_auto, corr_max, err_tdl], ['synth autocorrelation','correlation max', 'error estimate'],loc = 'upper right',fontsize = 8)
        ax_auto.annotate('e', (0, 1), xytext = (5,-5), xycoords = 'axes fraction', fontsize = 12, textcoords = 'offset points', color = 'k', backgroundcolor = 'none', ha = 'left', va = 'top', bbox = dict(facecolor = 'white', edgecolor = 'black', pad = 2.0))

        fig.tight_layout()

        ################## Save or show  ##################

        if params_in['corr_save_pic']:

            pic_loc = f'{params_in['home']}/{params_in['log_loc']}/{str(params_in['code_start_time'])}/{os.path.basename(__file__).split('.')[0]}'

            pic_filename = f'Event_{input_dict['id_cmt']}_{nslc}_{row_obs['phase']:s}_{row_obs['t_taup']:.2f}s.pdf'

            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Save {pic_loc}/{pic_filename}')
            ###

            if not os.path.exists(pic_loc):
                os.makedirs(pic_loc, exist_ok=True)

            plt.savefig(f'{pic_loc}/{pic_filename}', format = 'pdf')

            plt.close()
        else:
            plt.show()

        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def select_windows(self, input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail):
        """
        Select time windows in corr_df based on some set criteria: max_tdelay, min_xcc
        
        ----------
        Reads and returns same inputs/outputs
        ** param == return **

        :param input_dict: event details for processing 
        :type input_dict: dict
        :param twin_in_obs: obs twin file location string
        :type twin_in_obs: str
        :param twin_df_obs: obs twin dataframe object
        :type twin_df_obs: pd.df
        :param twin_in_syn: syn twin file location string
        :type twin_in_syn: str
        :param twin_df_syn: syn twin dataframe object
        :type twin_df_syn: pd.df
        :param corr_df: correlated twin dataframe object
        :type corr_df: pd.df
        :param logfile: open logfile to report progress
        :type: logfile: open textfile for writing
        :param outfile: open outfile to report results
        :type: outfile: open textfile for writing
        :param fail: fail flag = 1 if function fails, else 0
        :type fail: int
        """
        params_in = input_dict['params_in']

        # setup logfile statement using event_id, filename and function name (with inspect)
        log_statement = f'{self.tk.get_log_statement(input_dict['event_id'], twin_in_obs)}, {str(inspect.stack()[0][3])}'

        if fail:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
        else:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
            ###

            drop_list = []    
            num_rows = len(corr_df)   
            for index, row in corr_df.iterrows():

                # Quality control time delay:
                if np.abs(row['tdelay']) > params_in['max_tdelay']:
                    drop_list.append(index)
                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject row {index + 1:d}/{num_rows:d} - Time delay too large {row['phase']} @ {row['nslc']}: {row['tdelay']:.2f}s > {params_in['max_tdelay']}s limit')
                    ###
                    continue

                # Quality control time delay error:
                if row['tderr'] > params_in['max_tderr'] or row['tderr'] == -999.0:
                    drop_list.append(index)

                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject row {index + 1:d}/{num_rows:d} - Time delay error too large {row['phase']} @ {row['nslc']}: {row['tderr']:.2f}s > {params_in['max_tderr']}s limit')
                    ###
                    continue
                
                # Quality control XC coefficient
                if row['ccmx'] < params_in['min_xcc']:
                    drop_list.append(index)
                    
                    ###
                    self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  Reject row {index + 1:d}/{num_rows:d} - XC coefficient too small {row['phase']} @ {row['nslc']}: {row['ccmx']:.2f} < {params_in['min_xcc']} limit')
                    ###
                    continue

            # Remove rows that fail QC above in dataframe
            corr_df = corr_df.drop(drop_list)

            # Check for repeat picks... stations, locations, phases etc.

        return input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def save_tdelay_files(self, input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail):
        """
        Use this to save tdelay files whether with or without the sort for high cross correlation co-efficients.
        
        ----------
        Reads and returns same inputs/outputs
        ** param == return **

        :param input_dict: event details for processing 
        :type input_dict: dict
        :param twin_in_obs: obs twin file location string
        :type twin_in_obs: str
        :param twin_df_obs: obs twin dataframe object
        :type twin_df_obs: pd.df
        :param twin_in_syn: syn twin file location string
        :type twin_in_syn: str
        :param twin_df_syn: syn twin dataframe object
        :type twin_df_syn: pd.df
        :param corr_df: correlated twin dataframe object
        :type corr_df: pd.df
        :param logfile: open logfile to report progress
        :type: logfile: open textfile for writing
        :param outfile: open outfile to report results
        :type: outfile: open textfile for writing
        :param fail: fail flag = 1 if function fails, else 0
        :type fail: int
        """
        params_in = input_dict['params_in']

        # setup logfile statement using event_id, filename and function name (with inspect)
        log_statement =  f'{self.tk.get_log_statement(input_dict['event_id'], twin_in_obs)}, {str(inspect.stack()[0][3])}'

        if fail:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
        else:
            ###
            self.tk.print_log(params_in, logfile, f'{log_statement:s}  ,  -EXECUTING-')
            ###

            # Drop any Nans (should be zero)
            non_nan_merged_df = corr_df.dropna()
            if not non_nan_merged_df.empty:
                for index, row in non_nan_merged_df.iterrows():
                    # Print each row of corr_df to outfile
                    tw_info = params_in['correlate_outfmt'].format(row['nslc'],row['stla'],row['stlo'],row['stel'],row['phase'],row['tdelay'],row['tderr'],row['ccmx'],row['ttaup'],row['tp_obs'],row['tp_syn'],row['Ap_obs'],row['Ap_syn']) + '\n'
                    outfile.write(tw_info)

                # Counting stats
                input_dict['num_obj_out']   += len(non_nan_merged_df)
                input_dict['num_files_out'] += 1


        return input_dict, twin_in_obs, twin_df_obs, twin_in_syn, twin_df_syn, corr_df, logfile, outfile, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_log_file(self, input_dict):
        """
        Return an open log file in 
        log_loc/'code_start_time'/'filename'/event_name.log

        :param input_dict: loaded input parameters
        :type input_dict: dict
        :return logfile: open logfile to write code output
        :type logfile: open txt file
        """
        params_in = input_dict['params_in']

        lf_loc = f'{params_in['home']}/{params_in['log_loc']}/{str(params_in['code_start_time'])}/{os.path.basename(__file__).split('.')[0]}'

        if not os.path.exists(lf_loc):
            os.makedirs(lf_loc, exist_ok=True)

        lf_name = f'{lf_loc}/{str(input_dict['id_fmt_ctm'])}.log'

        logfile = open(lf_name,'w')
        return logfile


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def open_outfile_file(self, input_dict):
        """
        Return an open twin file in 
        twin_loc/filename_out_loc+component/event_name."OS/MX"+component.twin
        
        :param input_dict: loaded input parameters
        :type input_dict: dict
        :return outfile: open outfile to write code output
        :type outfile: open txt file
        """
        params_in = input_dict['params_in']

        of_loc = f'{params_in['home']}/{params_in['tdelay_loc']}/{params_in['phase_a_obs_out_loc'][-2:]}{params_in['component']}-{params_in['phase_a_syn_out_loc'][-2:]}{params_in['component']}'

        if not os.path.exists(of_loc):
            os.makedirs(of_loc, exist_ok=True)

        of_name = f'{of_loc}/{str(input_dict['id_ctm'])}.T.tdl'

        outfile = open(of_name,'w')
        return outfile


    ##############################################################################
    def A1(self, obs, syn):
        """
        A1(tau) EQ 7 Zaroli et al 2010 GJI
        Input: observed, synth windows
        Returns A1 output timeseries (length of correlation)

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :return output: function output of A1
        :type output: np.array
        """
        auto_syn = signal.correlate(syn, syn, mode = "full", method = "fft")

        result = signal.correlate(obs, syn, mode = "full", method = "fft")

        output = result / np.max(auto_syn)
        return output

    ##############################################################################
    def A2(self, obs, syn): 
        """
        A2(tau) EQ 7 Zaroli et al 2010 GJI
        Input: observed, synth windows
        Returns A2 output timeseries (length of correlation)

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :return output: function output of A2
        :type output: np.array
        """
        np.seterr(divide = 'ignore', invalid = 'ignore')

        auto_obs = signal.correlate(obs, obs, mode = "full", method = "fft")

        result = signal.correlate(obs, syn, mode = "full", method = "fft")

        output = np.max(auto_obs) / (result + np.finfo(float).eps) # Avoid div by zero errors

        return output

    ##############################################################################
    def Func_2(self, obs, syn):
        """
        F2(tau) EQ 5 Zaroli et al 2010 GJI
        Input: observed, synth windows
        Returns F2 output timeseries (length of correlation)

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :return output: function output of F2
        :type output: np.array
        """
        np.seterr(divide = 'ignore', invalid = 'ignore')

        R1 = self.A1(obs, syn)
        R2 = self.A2(obs, syn)

        F2 = np.minimum(R1,R2) / (np.maximum(R1,R2)  + np.finfo(float).eps) # Avoid div by zero errors

        # EQ 4 Zaroli et al 2010 GJI constraint 
        F2[F2 < 0] = 0
        F2[F2 > 1] = 0

        return F2


    ##############################################################################
    def Func_3(self, obs, syn, t_obs, t_syn, ax_time, ax_time_F3, delta):
        """
        F3(tau) EQ 8 Zaroli et al 2010 GJI
        Input: observed, synth windows, obs and syn time axes, 
        output time axes and delta
        Returns F3 output timeseries (length of ax_time)

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :param t_obs: obs twin independent variable (time axis)
        :type t_obs: np.array
        :param t_syn: syn twin independent variable (time axis)
        :type t_syn: np.array
        :param ax_time: time axis for twin plotting
        :type ax_time: np.array
        :param ax_time_F3: time axis for calculation of Zaroli et al (2010) F3
        :type ax_time_F3: np.array
        :param delta: 1/sample rate - time between points
        :type delta: float
        :return output: function output of F3
        :type output: np.array

        """
        y_syn_pad = np.interp(ax_time, t_syn, syn, left = 0, right = 0)
        y_obs_pad = np.interp(ax_time, t_obs, obs, left = 0, right = 0)
        
        F1_res = self.Func_1(y_obs_pad, y_syn_pad, ax_time, ax_time_F3, delta)
        F1_times = ax_time_F3
        F1_res_pad = np.interp(ax_time_F3, F1_times, F1_res, left = 0, right = 0) 
        
        F2_res = self.Func_2(y_obs_pad, y_syn_pad)
        F2_times =  np.arange(-1 * len(y_obs_pad), len(y_syn_pad)-1, 1) * delta
        F2_res_pad = np.interp(ax_time_F3, F2_times, F2_res, left = 0, right = 0) 

        output = (F1_res_pad + F2_res_pad) / 2

        return output


    # ##############################################################################
    def shift_and_truncate(self, time_series, orig_time, shift_indices, interp_time, sample_rate = 128):
        """
        Shift the input time series by a specified number of indices
        and truncate the output to the range between -7 and 7 seconds.

        :param time_series: Input time series to shift and truncate
        :type time_series: np.array
        :param orig_time: Input time series independent variable vector (time)
        :type orig_time: np.array
        :param shift_indices: number of indicies to shift time series
        :type shift_indices: int
        :param interp_time: time array onto which to interpolate data
        :type interp_time: np.array
        :params sample_rate: time series sampling rate, samples per sec
        :type params: float

        :return shifted_time_series: Shifted and truncated time series
        :type shifted_time_series: np.array
        """
        # Calculate the time shift in seconds
        time_shift = shift_indices / sample_rate

        # Calculate the time array for the shifted time series
        shifted_time = orig_time + time_shift

        # Use linear interpolation to get the values at the desired time points
        shifted_time_series = np.interp(interp_time, shifted_time, time_series, left = 0, right = 0)

        return shifted_time_series


    # ##############################################################################
    def res_trap_d(self, obs, syn, delta):
        """
        Trapezoidal integration for numerator of F1(tau)
        EQ 4 Zaroli et al 2010 GJI
        Inputs: obs and shifted synth time window
        Output: result of integral

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :param delta: 1/sample rate - time between points
        :type delta: float

        :return res_trap_d: function output of integral
        :type res_trap_d: np.array
        """
        squared_diff = (obs[:-1] - syn[:-1])**2 + (obs[1:] - syn[1:])**2
        res_trap_d = np.sum(squared_diff) * delta / 2
        return res_trap_d


    # ##############################################################################
    def res_trap(self, obs, delta):
        """
        Trapezoidal integration for denominator of F1(tau)
        EQ 4 Zaroli et al 2010 GJI
        Inputs: obs time window
        Output: result of integral

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param delta: 1/sample rate - time between points
        :type delta: float

        :return res_trap: function output of integral
        :type res_trap: np.array
        """
        res_trap = np.sum(delta * (obs[:-1]**2 + obs[1:]**2) / 2)
        return res_trap


    # ##############################################################################
    def Func_1(self, obs, syn, ax_time, ax_time_F3, delta):
        """
        F1(tau) EQ 4 Zaroli et al 2010 GJI
        Input: observed, synth windows, output time axes and delta
        Returns F1 output timeseries (length of ax_time)

        :param obs: obs twin dependent variable
        :type obs: np.array
        :param syn: syn twin dependent variable
        :type syn: np.array
        :param ax_time: time axis for twin plotting
        :type ax_time: np.array
        :param ax_time_F3: time axis for calculation of Zaroli et al (2010) F3
        :type ax_time_F3: np.array
        :param delta: 1/sample rate - time between points
        :type delta: float
        :return F1: function output of F1
        :type F1: np.array
        """

        F1_n = []
        for i in range(len(ax_time_F3)):
            shifted_t_series = self.shift_and_truncate(syn, ax_time, ax_time_F3[i]/delta, ax_time, sample_rate = 1/delta)
            res = self.res_trap_d(obs, shifted_t_series , delta)
            F1_n.append(res)

        F1 = 1 - np.array(F1_n)/ self.res_trap(obs, delta)

        # EQ 4 Zaroli et al 2010 GJI constraint 
        F1[F1 < 0] = 0
        F1[F1 > 1] = 0

        return F1


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 

def main():

    # Create instance of toolkit as tk
    tk = Toolkit()

    # Params go here.
    params_in = tk.get_params('params_in.yaml')

    # Define input data directory and function list.
    input_directory = f'{params_in['syn_loc']}/e{str(params_in['year'])}{params_in['fmt_data_loc']}'

    # Get event id table as pandas data frame
    evt_id_tab = tk.get_event_id_table(params_in)

    # Get phases names dictionary 
    from v01_phasenames import Phases
    phases = Phases().get_phase_dictionary()

    start_time = time.time()

    correlate_twin = CorrelateTwin()

    main_function = [correlate_twin.process_one_event]

    functions = [correlate_twin.correlate_windows, correlate_twin.select_windows, correlate_twin.save_tdelay_files]

    tk.execute(main_function, input_directory, evt_id_tab, functions, params_in, phases)
    
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()




