import os, sys
import pandas as pd
import numpy as np
from toolkit import Toolkit
from post_processing.compare_tdl_files import CompareTdlFiles
import concurrent.futures
from subprocess import call
import subprocess

from obspy.taup import TauPyModel
import obspy.geodetics

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
class CreateInvFiles(CompareTdlFiles):
    """
    Class to handle the creation of inversion ready files from
    time delay (.tdl) outputs by the MTWSPy algorithm
    in SEISGLOB2 formats
    May operate in parallel using parameters in param_in.yaml file

    Applies certain filters to the output dataset
    Creates formatted *evt_sta.dat, *data_error.dat files
    Creates formatted *Path/*Path??????.dat files
    Creates formatted correction/*dt_correction.dat files
    """
    
    tk = Toolkit()

    def __init__(self):
        super().__init__()


    def write_dt_files(self, params, input_df):
        """
        Write dt files for inversion
        Creates formatted *evt_sta.dat, *data_error.dat files

        :param params: loaded parameter file
        :type params: dict
        :param input_df: dataframe of input tdls
        :type input_df: pd.df
        """

        of_loc = f'{params['home']}/inv_out_files/list/'
        # print(of_loc)

        if not os.path.exists(of_loc):
            os.makedirs(of_loc, exist_ok=True)

        of_name_error = f'{of_loc}data_error.dat'
        outfile_error = open(of_name_error,'w')

        outfile_error_info = ''

        # Determine list of phases to print.
        if params['phases'][0] == 'All':
            from v01_phasenames import Phases
            phases_dict = Phases().get_phase_dictionary()
            phases = phases_dict[str(params['phases_key'])][str(params['component'])]
        else:
            phases = params['phases']
        
        # print(phases)

        for phase in phases:

            of_name_phase = f'{of_loc}list{str(phase)}_evt_sta.dat'
            outfile_phase = open(of_name_phase,'w')
            outfile_phase_info = ''

            phase_df = input_df.query("phase == @phase")
            
            for index, row in phase_df.iterrows():
                # list{phase}_evt_sta.dat: elat(k),elon(k),edep(k),slat(k),slon(k),dt(k),sig(k)
                # data_error.dat: dt(k),sig(k)
                outfile_phase_info  +=  params['outfile_phase_outfmt'].format(row['evt_lat'], 
                                                                              row['evt_lon'], 
                                                                              row['evt_dep'], 
                                                                              row['stla'], 
                                                                              row['stlo'], 
                                                                              row['tdelay'], 
                                                                              row['tderr']) + '\n' 
                
                outfile_error_info += params['data_error_outfmt'].format(row['tdelay'], 
                                                                         row['tderr']) + '\n'

            outfile_phase.write(outfile_phase_info)
            outfile_phase.close()

        outfile_error.write(outfile_error_info)
        outfile_error.close()


    #########################################################
    # Process the files in parallel.
    #########################################################

    def write_path_corr_files(self, params, input_df):
        """
        Takes a input tdl dataframe and writes path and correction
        outputs for each phase.

        :param params: loaded parameter file
        :type params: dict
        :param input_df: dataframe of input tdls
        :type input_df: pd.df
        """

        of_loc = f'{params['home']}/inv_out_files/paths/'
        # print(of_loc)

        if not os.path.exists(of_loc):
            os.makedirs(of_loc, exist_ok=True)

        of_loc_corr = f'{params['home']}/inv_out_files/correction/'
        # print(of_loc_corr)

        if not os.path.exists(of_loc_corr):
            os.makedirs(of_loc_corr, exist_ok=True)

        #parallel processing needs a list of single inputs, so we put all input into a dict and create a list of dicts

        # Determine list of phases to print.
        if params['phases'][0] == 'All':
            from v01_phasenames import Phases
            phases_dict = Phases().get_phase_dictionary()
            phases = phases_dict[str(params['phases_key'])][str(params['component'])]
        else:
            phases = params['phases']
        
        # print(phases)

        for phase in phases:
            input_dicts = []

            path_name = f'{phase}Path'

            path_loc = f'{of_loc}{path_name}'
            if not os.path.exists(path_loc):
                os.makedirs(path_loc, exist_ok=True)

            corr_name = f'{phase}dt_correction.dat'


            phase_df = input_df.query("phase == @phase")

            phase_df = phase_df.reset_index()

            for index, row in phase_df.iterrows():
                if 1 :
                    input_dict = {}
                    input_dict['params'] = params
                    input_dict['phase'] = phase
                    input_dict['phase_num'] = "{0:0=8d}".format(index + 1) 
                    input_dict['path_loc'] = f'{path_loc}'
                    input_dict['path_file'] = f'{path_name}_{input_dict['phase_num']}'

                    input_dict['corr_loc'] = f'{of_loc_corr}'
                    input_dict['corr_file'] = f'{corr_name}'

                    input_dict['model'] = params['taup_model_name']
                    input_dict['evt_lat'] = row['evt_lat']
                    input_dict['evt_lon'] = row['evt_lon']
                    input_dict['evt_dep'] = row['evt_dep']
                    input_dict['stla'] = row['stla']
                    input_dict['stlo'] = row['stlo']
                    input_dict['stel'] = row['stel']
                    # print(input_dict)
                    input_dicts.append(input_dict)


            if input_dicts:

                if params['parallel'] and params['cores'] > 1:
                    # Parallel processing
                    with concurrent.futures.ProcessPoolExecutor(max_workers = params['cores']) as executor:
                        r1 = executor.map(self.process_one_file_path, input_dicts)
                        r2 = executor.map(self.process_one_file_corr, input_dicts)

                        # Write path to file
                        outfile_path = open(f'{of_loc_corr}{corr_name}','w')
                        outfile_path_info = ''

                        for res in r2:
                            # print(phase, params['corr_outfmt'].format(res[0], res[1]))
                            outfile_path_info  +=  params['corr_outfmt'].format(res[0], res[1]) + '\n' 

                        outfile_path.write(outfile_path_info)
                        outfile_path.close()
                else:
                    # Serial processing
                    # Write path to file
                    outfile_path = open(f'{of_loc_corr}{corr_name}','w')
                    outfile_path_info = ''

                    for input_dict in input_dicts:
                        r1 = self.process_one_file_path(input_dict)
                        r2 = self.process_one_file_corr(input_dict)
                        # print(phase, params['corr_outfmt'].format(r2[0], r2[1]))
                        outfile_path_info  +=  params['corr_outfmt'].format(r2[0], r2[1]) + '\n' 

                    outfile_path.write(outfile_path_info)
                    outfile_path.close()

        return 


    #########################################################
    # Function to process one file for a path.
    #########################################################

    def process_one_file_path(self, input_dict):
        """
        Processes one line of the input dictionary to give the taup path
        Creates formatted *Path/*Path??????.dat files

        :param input_dict: input dict containing the params necessary for 
        taup path
        :type input_dict: dict
        """

        params = input_dict['params']
        model = input_dict['model']
        phase = input_dict['phase']
        STLA = input_dict['stla']
        STLO = input_dict['stlo']
        EVLA = input_dict['evt_lat']
        EVLO = input_dict['evt_lon']
        EVDP = input_dict['evt_dep']
        temp_outfile = f'{input_dict['path_file']}'
        outfile = f'{input_dict['path_loc']}/{input_dict['path_file']}.dat'

        test = ['taup_path -mod ' + str(model) +' -h ' + str(
                EVDP) + ' -ph ' + str(phase) + ' -sta ' + str(
                STLA) + ' ' + str(
                    STLO) + ' -evt ' + str(
                        EVLA) + ' ' + str(
                            EVLO) + ' -o ' + str(temp_outfile)]

        # print(test)
        # Run test in terminal
        out = subprocess.check_output(
            test,
            shell=True,
            universal_newlines=True)

        # Read the file called f'{temp_outfile}.gmt'
        try: 
            path_df = pd.read_csv(f'./{temp_outfile}.gmt', 
                                skiprows=1, 
                                sep='\\s+', 
                                names=['dist', 'radius', 'lat', 'lon'])
        except:
            # Need something more complex if multipathing occurs
            path_df = pd.DataFrame(columns = ['dist', 'radius', 'lat', 'lon'])

            # read file
            f=open(f'./{temp_outfile}.gmt','r')
            lines=f.readlines()
            f.close()     

            # Skip header:
            for line in lines [1:]:
                data=line.split()
                if len(data) == 4:
                    # Read:
                    dist = float(data[0])
                    radius = float(data[1])
                    lat = float(data[2])
                    lon = float(data[3])


                    lt = [dist,radius,lat,lon]
                    # print(lt)
                    path_df.loc[len(path_df)] = lt

                else:
                    # Found the multipathed phase
                    # So continue.
                    continue

        # Write path to file
        outfile_path = open(outfile,'w')
        outfile_path_info = ''
        outfile_path_info  +=  '{0: 12.7f}'.format(len(path_df)) + '\n'
        for index, row in path_df.iterrows():
            outfile_path_info  +=  params['path_outfmt'].format(row['lat'], 
                                                                row['lon'], 
                                                                6371.0 - row['radius'], 
                                                                row['dist']) + '\n' 

        outfile_path.write(outfile_path_info)
        outfile_path.close()
        
        # Remove temp taup files
        try:
            # print(f'deleting: {temp_outfile}')
            os.remove(f'./{temp_outfile}.gmt')
        except: 
            Dummy=1
            
        return 


    #########################################################
    # Function to process one file for a correction
    #########################################################

    def process_one_file_corr(self, input_dict):
        """
        Processes one entry of the input dict for the taup corrections
        Creates formatted correction/*dt_correction.dat files

        :param input_dict: input dict containing the params necessary for taup
        :type input_dict: dict

        :return arr_mod1[0].time: Arrival time of phase for taup model 1, PREM
        :type arr_mod1[0].time: float
        :return arr_mod2[0].time: Arrival time of phase for taup model 2
        :type arr_mod2[0].time: float
        """

        params = input_dict['params']
        model = input_dict['model']
        phase = input_dict['phase']
        STLA = input_dict['stla']
        STLO = input_dict['stlo']
        EVLA = input_dict['evt_lat']
        EVLO = input_dict['evt_lon']
        EVDP = input_dict['evt_dep']

        model1 = TauPyModel(model="prem")
        model2 = TauPyModel(model=model)
        
        arr_mod1 = model1.get_travel_times_geo(source_depth_in_km = EVDP,
                                    source_latitude_in_deg = EVLA,
                                    source_longitude_in_deg = EVLO, 
                                    receiver_latitude_in_deg = STLA,
                                    receiver_longitude_in_deg = STLO,
                                    phase_list = [phase])

        if len(arr_mod1) == 0:
            # Possibly phase that is at limit of the distance range
            # check upgoing s instead of S.
            # Others may be required

            distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(EVLA, 
                                                                    EVLO, 
                                                                    STLA, 
                                                                    STLO)
            distdg = distm / (6371.e3 * np.pi / 180.)            
            
            if phase == 'S' and distdg <= 22.0:
                #likely an upgoing s in Prem
                # print('Probably upgoing s phase...')
                arr_mod1 = model1.get_travel_times_geo(source_depth_in_km = EVDP,
                                            source_latitude_in_deg = EVLA,
                                            source_longitude_in_deg = EVLO, 
                                            receiver_latitude_in_deg = STLA,
                                            receiver_longitude_in_deg = STLO,
                                            phase_list = ['s'])
                # print(f's: {np.round(arr_mod1[0].time,4)}s')
        
        arr_mod2 = model2.get_travel_times_geo(source_depth_in_km = EVDP,
                                    source_latitude_in_deg = EVLA,
                                    source_longitude_in_deg = EVLO, 
                                    receiver_latitude_in_deg = STLA,
                                    receiver_longitude_in_deg = STLO,
                                    phase_list = [phase])
        
        T2 = arr_mod2[0].time

        # If we cant find the correct time make T1 equal to T2
        try:
            T1 = arr_mod1[0].time
        except:
            T1 = T2

        #if arr_mod1[0].time and arr_mod2[0].time:
        if T1 and T2:
            # print(phase, params['corr_outfmt'].format(T1, T2))
            return T1, T2

        else:
            print('MASSIVE ISSUE - arrivals not found')
            print(input_dict)
            print(f'mod source_depth_in_km = {EVDP},\
                                source_latitude_in_deg = {EVLA},\
                                source_longitude_in_deg = {EVLO}, \
                                receiver_latitude_in_deg = {STLA},\
                                receiver_longitude_in_deg = {STLO},\
                                phase_list = [{phase}]')
            return 0.0, 0.0


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def main():

    # Create instance of toolkit as tk
    tk = Toolkit()
    params = tk.get_params('params_in.yaml')

    # Instantiate the classes and call their methods as needed
    create_inv_files = CreateInvFiles()

    filter_functions = [create_inv_files.filt_date_time_max,
                        create_inv_files.filt_date_time_min,
                        create_inv_files.filt_networks,
                        create_inv_files.filt_components,
                        create_inv_files.filt_phases,
                        create_inv_files.filt_tdl_max,
                        create_inv_files.filt_ccmx_min]    


    filename = f"{params['tt_out_f_name']}"

    ######################### AL data - XC - T comp #########################
    # Location where MTWSPy code is installed
    # params['home'] = '/Users/alistair/Google_Drive/GITHUB_AB/MTWSPy'

    output_directory = f"{params['home']}/{params['proc_tdl_loc']}/" # /FINAL/"

    python_filtered_XC_df = create_inv_files.load_dataframe(params, 
                                                            filter_functions, 
                                                            output_directory, 
                                                            filename)

    python_filtered_XC_df = create_inv_files.filter_dataframe(params, 
                                                              filter_functions, 
                                                              python_filtered_XC_df)
    
    print(python_filtered_XC_df)

    create_inv_files.write_dt_files(params, python_filtered_XC_df)

    create_inv_files.write_path_corr_files(params, python_filtered_XC_df)



if __name__ == '__main__':
    main()
