# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:10:02 2024
@author: alistairboyce
if this fails: >> conda activate env3.12
"""

import time
# Import individual codes
import mk_events_csv
import toolkit
import v01_phasenames, find_twin_obs, find_twin_syn, match_twin_files
import phase_association_obs, phase_association_syn, correlate_twin


def main():

    start_time = time.time()
    # Params go here.
    params_in = toolkit.get_params('params_in.yaml')

    # Define observed data input directory
    obs_input_directory = str(params_in['data_loc']) + '/e'+str(params_in['year']) + str(params_in['fmt_data_loc'])
    # Define synthetic data input directory
    syn_input_directory = str(params_in['synth_loc']) + '/e'+str(params_in['year']) + str(params_in['fmt_data_loc'])


    ######### Check for presence of Data files #########
    toolkit.check_files(obs_input_directory, params_in['year'], params_in['component'])
    
    ######### Check for presence of Synth files #########
    toolkit.check_files(syn_input_directory, params_in['year'], params_in['component'])

    #########          mk_events_csv.py        #########
    # print(mk_events_csv)
    # Can choose between below
    # evt_id_tab = mk_events_csv.execute(params_in)
    evt_id_tab = toolkit.get_event_id_table(params_in['cmt_outfile'])

    #########     v01_phasenames.py           ##########
    # Get phases
    phases = v01_phasenames.phases()



    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     find_twin_obs.py           ##########
    # Data: get travel times, primary quality control (snr_pow, snr_amp), detect time windows, filter time windows, plot time windows
    main_function = [find_twin_obs.process_one_event]

    # Functions to perform on observed data 
    functions = [find_twin_obs.get_tt_times, find_twin_obs.primary_qc, find_twin_obs.detect_window_peaks, find_twin_obs.filter_window_peaks, find_twin_obs.plot_waveform_envelope_peaks_windows]
    
    toolkit.execute(main_function, obs_input_directory, evt_id_tab, functions, params_in, phases)



    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     find_twin_syn.py           ##########
    # Synth: get travel times, filter seismograms, detect time windows, filter time windows, plot time windows (NO NEED FOR primary_qc)
    # Manually update two main parameters:
    params_in['sig_win_ext'] = params_in['Tc'] 
    # params_in['sig_win_type'] = 0 
    
    main_function = [find_twin_syn.process_one_event]
    # Functions to perform on synthetic data 
    functions = [find_twin_obs.get_tt_times, find_twin_obs.filter_seis, find_twin_obs.detect_window_peaks, find_twin_obs.filter_window_peaks, find_twin_obs.plot_waveform_envelope_peaks_windows]
        
    toolkit.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)
    

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     match_twin_files.py           ##########
    # find common twin files for obs and syn  s05_twin_files.m

    match_twin_files.find_common_events(params_in)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     phase_association_obs.py           ##########
    # Data: Phase association s|w40_phase_association.m 
    # Associate time window with a given phase

    main_function = [phase_association_obs.process_one_event]
    # Functions to perform on observed data 
    functions = [phase_association_obs.associate_twin_phase] # 

    toolkit.execute(main_function, obs_input_directory, evt_id_tab, functions, params_in, phases)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     phase_association_syn.py           ##########
    # Synth: Phase association s|w40_phase_association.m
    # Associate time window with a given phase

    main_function = [phase_association_syn.process_one_event]
    # Functions to perform on synthetic data 
    functions = [phase_association_obs.associate_twin_phase] # 

    toolkit.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     correlate_twin.py           ##########
    # Measure time delay for available syn-obs time window pairs 
    # Select high quality results....
    # (s|w10_twin_corr.m)

    main_function = [correlate_twin.process_one_event]
    functions = [correlate_twin.correlate_windows, correlate_twin.select_windows, correlate_twin.save_tdelay_files]

    toolkit.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)

  
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()