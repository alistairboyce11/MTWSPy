# tests/test_MTWSPy_main.py

import time
# Import individual codes
from MTWSPy.toolkit import Toolkit
from MTWSPy.match_catalog import MatchCatalog
from MTWSPy.v01_phasenames import Phases
from MTWSPy.find_twin import FindTwinObs, FindTwinSyn
from MTWSPy.match_twin_files import MatchTwinFiles
from MTWSPy.phase_association import PhaseAssociationObs, PhaseAssociationSyn
from MTWSPy.correlate_twin import CorrelateTwin
from MTWSPy.post_processing.process_tdl_files import ProcessTdlFiles
from MTWSPy.post_processing.create_inv_files import CreateInvFiles

import glob, os 
import numpy as np
import platform

def test_MTWSPy_main():
    # Create instance of toolkit as tk
    tk = Toolkit()

    # Params go here.

    params_in = tk.get_params('params_in.yaml')

    # Modify a few paths so that we dont put test data where real data will go

    params_in['cmt_outfile'] = 'test_output/cmt-events' # Prefix for match_catalog events file
    params_in['mtf_outfilename'] = f'test_output/{params_in['mtf_outfilename']}' # Matched twin files prefix
    params_in['log_loc'] = 'test_output/log' # Log location
    params_in['twin_loc'] = 'test_output/twin' # Time Window (twin) file location
    params_in['tdelay_loc'] = 'test_output/tdelay' # Time delay (tdelay) file location
    params_in['proc_tdl_loc'] = 'test_output/proc_tdelay' # Processed dataset location
    params_in['inv_out_loc'] = 'test_output/inv_out_files' # Inversion ready files location
    params_in['verbose'] = False # Hide terminal output during test

    # Define observed data input directory
    obs_input_directory = f'{str(params_in['obs_loc'])}/e{str(params_in['year'])}{str(params_in['fmt_data_loc'])}'
    # Define synthetic data input directory
    syn_input_directory = f'{str(params_in['syn_loc'])}/e{str(params_in['year'])}{str(params_in['fmt_data_loc'])}'


    # Create an instance of the match catalog as mc
    mc = MatchCatalog(params_in)
    evt_id_tab = mc.execute()
    # evt_id_tab = tk.get_event_id_table(params_in)

    #########     v01_phasenames.py           ##########
    # Get phases
    phases = Phases().get_phase_dictionary()

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     find_twin_obs.py           ##########
    # Data: get travel times, primary quality control (snr_pow, snr_amp), detect time windows, filter time windows, plot time windows
    find_twin_obs = FindTwinObs()

    main_function = [find_twin_obs.process_one_event]

    # Functions to perform on observed data 
    functions = [find_twin_obs.get_tt_times, find_twin_obs.primary_qc, find_twin_obs.detect_window_peaks, find_twin_obs.filter_window_peaks, find_twin_obs.plot_waveform_envelope_peaks_windows]
    
    tk.execute(main_function, obs_input_directory, evt_id_tab, functions, params_in, phases)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     find_twin_syn.py           ##########
    # Synth: get travel times, filter seismograms, detect time windows, filter time windows, plot time windows (NO NEED FOR primary_qc)
    # Manually update two main parameters:
    params_in['sig_win_ext'] = params_in['Tc'] 
    # params_in['sig_win_type'] = 0 
    
    find_twin_syn = FindTwinSyn()

    main_function = [find_twin_syn.process_one_event]
    # Functions to perform on synthetic data 
    functions = [find_twin_syn.get_tt_times, find_twin_syn.filter_seis, find_twin_syn.detect_window_peaks, find_twin_syn.filter_window_peaks, find_twin_syn.plot_waveform_envelope_peaks_windows]
        
    tk.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)
    

    # # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    # #########     match_twin_files.py           ##########
    # # find common twin files for obs and syn  s05_twin_files.m

    match_twin_files = MatchTwinFiles(params_in)
    # Make the list of common evids and write to file.
    match_twin_files.find_common_events(params_in)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     phase_association_obs.py           ##########
    # Data: Phase association s|w40_phase_association.m 
    # Associate time window with a given phase

    phase_association_obs = PhaseAssociationObs()

    main_function = [phase_association_obs.process_one_event]
    # Functions to perform on observed data 
    functions = [phase_association_obs.associate_twin_phase] # 

    tk.execute(main_function, obs_input_directory, evt_id_tab, functions, params_in, phases)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     phase_association_syn.py           ##########
    # Synth: Phase association s|w40_phase_association.m
    # Associate time window with a given phase

    phase_association_syn = PhaseAssociationSyn()

    main_function = [phase_association_syn.process_one_event]
    # Functions to perform on synthetic data 
    functions = [phase_association_syn.associate_twin_phase] # 

    tk.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########     correlate_twin.py           ##########
    # Measure time delay for available syn-obs time window pairs 
    # Select high quality results....
    # (s|w10_twin_corr.m)

    correlate_twin = CorrelateTwin()

    main_function = [correlate_twin.process_one_event]

    functions = [correlate_twin.correlate_windows, correlate_twin.select_windows, correlate_twin.save_tdelay_files]

    tk.execute(main_function, syn_input_directory, evt_id_tab, functions, params_in, phases)





    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    #########           process_tdl_files.py     #########


    # Instantiate the classes and call their methods as needed
    process_tdl_files = ProcessTdlFiles()

    filter_functions = [process_tdl_files.filt_date_time_max,
                        process_tdl_files.filt_date_time_min,
                        process_tdl_files.filt_networks,
                        process_tdl_files.filt_components,
                        process_tdl_files.filt_evt_lat_max,
                        process_tdl_files.filt_evt_lat_min,
                        process_tdl_files.filt_evt_lon_max,
                        process_tdl_files.filt_evt_lon_min,
                        process_tdl_files.filt_evt_dep_max,
                        process_tdl_files.filt_evt_dep_min,
                        process_tdl_files.filt_evt_mag_max,
                        process_tdl_files.filt_evt_mag_min,
                        process_tdl_files.filt_sta_lat_max,
                        process_tdl_files.filt_sta_lat_min,
                        process_tdl_files.filt_sta_lon_max,
                        process_tdl_files.filt_sta_lon_min,
                        process_tdl_files.filt_phases,
                        process_tdl_files.filt_tdl_max,
                        process_tdl_files.filt_ccmx_min,
                        process_tdl_files.filt_tderr_max,
                        process_tdl_files.filt_dist_max,
                        process_tdl_files.filt_dist_min]

    filtered_df = process_tdl_files.apply_stuff(params_in, filter_functions)
    
    print(filtered_df)

    # Write outfile
    # process_tdl_files.write_tdl_summary_file(params, filtered_df)

    output_directory = f"{params_in['home']}/{params_in['proc_tdl_loc']}/"

    filename = f"{params_in['tt_out_f_name']}"
    process_tdl_files.write_dataframe(output_directory, filename, filtered_df)

    # # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    # #########           create_inv_files.py     #########


    # Instantiate the classes and call their methods as needed
    create_inv_files = CreateInvFiles()

    filter_functions = [create_inv_files.filt_date_time_max,
                        create_inv_files.filt_date_time_min,
                        create_inv_files.filt_networks,
                        create_inv_files.filt_components,
                        create_inv_files.filt_phases,
                        create_inv_files.filt_tdl_max,
                        create_inv_files.filt_ccmx_min]    


    filename = f"{params_in['tt_out_f_name']}"

    ######################### AL data - XC - T comp #########################
    
    # The create_inv_files only works reliable on Linux as Mac is often not case-sensitive
    if not "Darwin" in platform.uname():
        print('Checking create inv')


        output_directory = f"{params_in['home']}/{params_in['proc_tdl_loc']}/"

        python_filtered_XC_df = create_inv_files.load_dataframe(params_in, 
                                                                filter_functions, 
                                                                output_directory, 
                                                                filename)

        python_filtered_XC_df = create_inv_files.filter_dataframe(params_in, 
                                                                filter_functions, 
                                                                python_filtered_XC_df)
        
        print(python_filtered_XC_df)

        create_inv_files.write_dt_files(params_in, python_filtered_XC_df)

        create_inv_files.write_path_corr_files(params_in, python_filtered_XC_df)

    else:
        print('Not checking create_inv')
        


    ##################### TESTING CODES ##########################

    # First check that len (EVID_2008.OST-MXT) is 35

    params_in['mtf_outfilename']
    with open(params_in['mtf_outfilename'], 'r') as f:
        lines = f.readlines()
    
    assert len(lines) == 35

    # Check we have 35 twin files for obs and syn and 35 tdelay output files.
    obs_twin_files = glob.glob(f'{params_in["home"]}/{params_in['twin_loc']}/phs/obs/OST/*twin')
    syn_twin_files = glob.glob(f'{params_in["home"]}/{params_in['twin_loc']}/phs/syn/MXT/*twin')
    tdl_files = glob.glob(f'{params_in['home']}/{params_in['tdelay_loc']}/OST-MXT/*tdl')

    assert len(obs_twin_files) == 35
    assert len(syn_twin_files) == 35
    assert len(tdl_files) == 35

    # Now check proc_tdelay file:

    tdl_df_test = create_inv_files.load_dataframe(params_in, 
                                                    filter_functions, 
                                                    f"{params_in['home']}/{params_in['proc_tdl_loc']}/", 
                                                    f"{params_in["tt_out_f_name"]}")

    tdl_df_orig = create_inv_files.load_dataframe(params_in, 
                                                    filter_functions, 
                                                    './MTWSPy/data/output_df/', 
                                                    '200801_IU_T_df')

    assert len(tdl_df_test) == len(tdl_df_orig)

    # Compare all values in dataframes
    for index, row in tdl_df_test.iterrows():
        true_vals = np.where((row == tdl_df_orig.iloc[index]) == True)[0]
        assert len(true_vals) == 25


    if not "Darwin" in platform.uname():
        print('Checking create inv')

        # Check data_error.dat
        data_error_file = glob.glob(f'{params_in["home"]}/{params_in['inv_out_loc']}/list/data_error.dat')
        assert os.path.isfile(data_error_file[0]) == True
        with open(data_error_file[0], 'r') as f:
            lines = f.readlines()
        assert len(lines) == 496

        # Check listS_evt_sta.dat
        listS_file = glob.glob(f'{params_in["home"]}/{params_in['inv_out_loc']}/list/listS_evt_sta.dat')
        assert os.path.isfile(listS_file[0]) == True
        with open(listS_file[0], 'r') as f:
            lines = f.readlines()
        assert len(lines) == 142
        
        # Check correction files
        corr_files = glob.glob(f'{params_in["home"]}/{params_in['inv_out_loc']}/correction/*dat')
        assert len(corr_files) == 14

        Scorr_files = glob.glob(f'{params_in["home"]}/{params_in['inv_out_loc']}/correction/Sdt_correction.dat')
        with open(Scorr_files[0], 'r') as f:
            lines = f.readlines()
        assert len(lines) == 142

        # Number of path files
        path_files = glob.glob(f'{params_in["home"]}/{params_in['inv_out_loc']}/paths/*/*dat')
        assert len(path_files) == 496
    
    else:
        print('Not checking create_inv')
        
def main():
    test_MTWSPy_main()

    return


if __name__ == '__main__':
    main()
