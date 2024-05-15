#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 01:39 2023

@author: alistair
"""
import numpy as np
import glob
from obspy import UTCDateTime
import concurrent.futures
import os
from geographiclib.geodesic import Geodesic
import pandas as pd
from toolkit import Toolkit
import obspy.geodetics


class ProcessTdlFiles:
    tk = Toolkit()

    def __init__(self):
        pass

    #########################################################
    # Process the files in parallel.
    #########################################################

    def apply_stuff(self, params, sort_functions):
        '''
        Takes a list of functions, reads all tdl files, applies each sort function to all of the tdl files 
        returns result as concatenated dataframe

        Parameters
        ----------
        input_directory: str
            input_directory name where all the data is.
        sort functions : list
            list of functions for sorting tdl measurements
            
        Returns
        -------
        outputs_all : dataframe
            contains the dataframe output of all the sort functions.
        '''
        input_directory = f'{params['home']}/{params['tdelay_loc']}/{params['phase_a_obs_out_loc'][-2:]}{params['component']}-{params['phase_a_syn_out_loc'][-2:]}{params['component']}'

        cores = params['cores']

        print(f'Processing tdl results for: {str(input_directory)} using {str(cores)} cores...')

        files=glob.glob(f'{input_directory}/*.tdl')

        files.sort()
        
        output_df = pd.DataFrame()

        #parallel processing needs a list of single inputs, so we put all input into a dict and create a list of dicts
        input_dicts = []
        for file in files:
            input_dict = {}
            input_dict['file'] = file
            input_dict['sort_functions'] = sort_functions
            input_dict['params'] = params

            input_dicts.append(input_dict)
                
        # Parallel processing
        with concurrent.futures.ProcessPoolExecutor(max_workers = cores) as executor:
            
            results = executor.map(self.process_one_file, input_dicts)
            
            # collapsing of results
            for res in results:
                if not res.empty:
                    cat_df = pd.concat([output_df, res], ignore_index=True)
                    output_df = cat_df

        return output_df
        
    #########################################################
    # Function to process one file.
    #########################################################

    def process_one_file(self, input_dict):
        '''
        Processes one file, reading the models in the file, applying functions
        Also adds, dist, az, baz, mid-lat, mid-lon

        Parameters
        ----------
        input_dict : dict
            input dictionary containing the file name (file), the list of functions (functions)

        Returns
        -------
        outputs : dataframe
            tdl measurements from a given file that pass all sort criteria

        '''
        
        file = input_dict['file']
        functions = input_dict['sort_functions']
        params = input_dict['params']

        tdl_df = self.tk.read_tdelay_file(file)

        if not tdl_df.empty:

            if len(tdl_df) > 0:

                mp_df = pd.DataFrame(columns = ['dist', 'az', 'baz', 'mid_lat','mid_lon'], index = range(0,len(tdl_df)))
                for index, row in tdl_df.iterrows():
                    # compute distances azimuth and backazimuth
                    distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(
                        float(row['evt_lat']), float(row['evt_lon']), float(row['stla']), float(row['stlo']))
                    distdg = distm / (6371.e3 * np.pi / 180.)

                    # Calculate midpoints for phases if necessary:
                    # Define the path from 1 to 2
                    l = Geodesic.WGS84.InverseLine(float(row['evt_lat']), float(row['evt_lon']), float(row['stla']), float(row['stlo']))

                    # Compute the midpoint
                    m = l.Position(0.5 * l.s13)

                    mp_df.loc[index, "dist"] = np.round(distdg,3)
                    mp_df.loc[index, "az"] = np.round(az,3)
                    mp_df.loc[index, "baz"] = np.round(baz,3)

                    mp_df.loc[index, "mid_lat"] = np.round(m['lat2'],3)
                    mp_df.loc[index, "mid_lon"] = np.round(m['lon2'],3)
                
                # Merge the dataframes
                tdl_df = pd.merge(tdl_df, mp_df, left_index = True, right_index = True)

            # apply functions
            for k, function in enumerate(functions):
                output_df = function(params, tdl_df)
                tdl_df = output_df

        else:
            print(f'Failed to read: {file}, or empty...')

        return tdl_df


    #################################################################
    def filt_date_time_max(self, params, input_df):
        '''
        limit pick date/time
        '''
        drop_list = []
        for index, row in input_df.iterrows():

            if UTCDateTime(str(row['date_time'])) > UTCDateTime(str(params['date_time_max'])):
                drop_list.append(index)

        # Remove rows that fail QC above in dataframe
        if drop_list:
            try: 
                drop_list = self.tk.flatten_concatenation(drop_list)
            except:
                pass
            output_df = input_df.drop(drop_list)
        else:
            output_df = input_df 
        
        return output_df

    #################################################################
    def filt_date_time_min(self, params, input_df):
        '''
        limit pick date/time
        '''
        drop_list = []
        for index, row in input_df.iterrows():

            if UTCDateTime(str(row['date_time'])) < UTCDateTime(str(params['date_time_min'])):
                drop_list.append(index)

        # Remove rows that fail QC above in dataframe
        if drop_list:
            try: 
                drop_list = self.tk.flatten_concatenation(drop_list)
            except:
                pass
            output_df = input_df.drop(drop_list)
        else:
            output_df = input_df 
        
        return output_df

    #################################################################
    def filt_evt_lat_max(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_lat_max']
        output_df = input_df.query("evt_lat <= @limit")

        return output_df

    #################################################################
    def filt_evt_lat_min(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_lat_min']
        output_df = input_df.query("evt_lat >= @limit")
    
        return output_df

    #################################################################
    def filt_evt_lon_max(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_lon_max']
        output_df = input_df.query("evt_lon <= @limit")
    
        return output_df

    #################################################################
    def filt_evt_lon_min(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_lon_min']
        output_df = input_df.query("evt_lon >= @limit")
    
        return output_df

    #################################################################
    def filt_evt_dep_max(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_dep_max']
        output_df = input_df.query("evt_dep <= @limit")
    
        return output_df

    #################################################################
    def filt_evt_dep_min(self, params, input_df):
        '''
        limit pick location
        '''
        limit = params['evt_dep_min']
        output_df = input_df.query("evt_dep >= @limit")
    
        return output_df

    #################################################################
    def filt_evt_mag_max(self, params, input_df):
        '''
        limit pick magnitude
        '''
        limit = params['evt_mag_max']
        output_df = input_df.query("evt_mag <= @limit")
    
        return output_df

    #################################################################
    def filt_evt_mag_min(self, params, input_df):
        '''
        limit pick magnitude
        '''
        limit = params['evt_mag_min']
        output_df = input_df.query("evt_mag >= @limit")
    
        return output_df

    #################################################################
    def filt_networks(self, params, input_df):
        '''
        limit network name
        '''
        if not params['networks'][0] == 'All':
            output_df = pd.DataFrame()

            for network in params['networks']:

                string = f'network_'

                res = input_df[input_df['nslc'].str.contains(string)]
                cat_df = pd.concat([output_df, res], ignore_index=True)
                output_df = cat_df

        else:
            # Do nothing
            output_df = input_df
        return output_df
        
    #################################################################
    def filt_components(self, params, input_df):
        '''
        limit Component
        '''
        if not params['components'][0] == 'All':
            output_df = pd.DataFrame()

            for comp in params['components']:

                res = input_df.query("channel == @comp")
                cat_df = pd.concat([output_df, res], ignore_index=True)
                output_df = cat_df

        else:
            # Do nothing
            output_df = input_df
        return output_df


    #################################################################
    def filt_phases(self, params, input_df):
        '''
        limit phases
        '''
        if not params['phases'][0] == 'All':
            output_df = pd.DataFrame()

            for phase in params['phases']:

                res = input_df.query("phase == @phase")
                cat_df = pd.concat([output_df, res], ignore_index=True)
                output_df = cat_df

        else:
            # Do nothing
            output_df = input_df
        return output_df

    #################################################################
    def filt_sta_lat_max(self, params, input_df):
        '''
        limit station location
        '''
        limit = params['sta_lat_max']
        output_df = input_df.query("stla <= @limit")
    
        return output_df

    #################################################################
    def filt_sta_lat_min(self, params, input_df):
        '''
        limit station location
        '''
        limit = params['sta_lat_min']
        output_df = input_df.query("stla >= @limit")
    
        return output_df

    #################################################################
    def filt_sta_lon_max(self, params, input_df):
        '''
        limit station location
        '''
        limit = params['sta_lon_max']
        output_df = input_df.query("stlo <= @limit")
    
        return output_df

    #################################################################
    def filt_sta_lon_min(self, params, input_df):
        '''
        limit station location
        '''
        limit = params['sta_lon_min']
        output_df = input_df.query("stlo >= @limit")
    
        return output_df

    #################################################################
    def filt_tdl_max(self, params, input_df):
        '''
        limit time delay
        '''
        limit = params['tdl_max']
        output_df = input_df.query("tdelay <= @limit")
        limit = -1 * limit
        output_df = output_df.query("tdelay >= @limit")
        return output_df

    #################################################################
    def filt_ccmx_min(self, params, input_df):
        '''
        limit cross correlation minimum value
        '''
        limit = params['ccmx_min']
        output_df = input_df.query("ccmx >= @limit")
    
        return output_df




    #################################################################
    def filt_tderr_max(self, params, input_df):
        '''
        limit time delay error maximum value
        '''
        limit = params['max_tderr']
        output_df = input_df.query("tderr <= @limit")
    
        return output_df

    # #################################################################
    def filt_dist_max(self, params, input_df):
        '''
        limit epicentral distance
        '''
        limit = params['dist_max']
        output_df = input_df.query("dist <= @limit")
    
        return output_df

    # #################################################################
    def filt_dist_min(self, params, input_df):
        '''
        limit epicentral distance
        '''
        limit = params['dist_min']
        output_df = input_df.query("dist >= @limit")
    
        return output_df


    ############################################################################
    def get_station_means(self, params, input_df):
        '''
        Get station means by phase for given delay time dataframe
        '''

        if len(params['phases']) == 1 and params['phases'][0] == 'All':
            # Make unique_phases list.
            unique_phases = sorted(list(input_df['phase'].unique()))
        else:
            unique_phases = params['phases']
        
        n_phases = len(unique_phases)

        print(f'Producing station means for {n_phases} phases....')
        print(f'For the following: {unique_phases}')

        column_names = ['nslc','stla','stlo','stel','channel','phase','av_tdelay','av_tderr','av_ccmx','n_picks']
        
        out_df = []
        
        for p, phase in enumerate(unique_phases):



            # Make phase data frame
            # ph_df = input_df.query("phase == @phase")  

            # ph_df = input_df[input_df['phase'] == str(phase)]

            ph_df = input_df[input_df['phase'].str.fullmatch(phase, case=True)]

            if len(ph_df) > 0:
                # Find unique stations in ev_df
                unique_stations = list(ph_df['nslc'].unique())
                n_stations = len(unique_stations)

                # Loop through all unique stations
                if n_stations > 0:
                    for s, station in enumerate(unique_stations):
                        # Make station data frame
                        st_df = ph_df.query("nslc == @station")
                        n_picks = len(st_df)

                        # Get columnn values.
                        nslc = str(st_df['nslc'].iloc[0])
                        stla = float(st_df['stla'].iloc[0])
                        stlo = float(st_df['stlo'].iloc[0])
                        stel = float(st_df['stel'].iloc[0])
                        channel = str(st_df['channel'].iloc[0])
                        av_tdelay = np.round(float(st_df['tdelay'].mean()),3)
                        av_tderr = np.round(float(st_df['tderr'].mean()),3)
                        av_ccmx = np.round(float(st_df['ccmx'].mean()),3)

                        # Maybe only add if we have actually taken a mean....
                        if n_picks > 1: 
                            out_df.append([nslc,stla,stlo,stel,channel,str(phase),av_tdelay,av_tderr,av_ccmx,n_picks])
                        

        output_df = pd.DataFrame(sorted(out_df), columns = column_names)
        print(output_df)

        return output_df

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_station_means(self, params, input_df):
        '''
        Writes travel time delay file for station means.
        
        Parameters
        ----------
        params : dict parameter file
        input_df : dataframe of tdl measurements to save

        Returns
        -------
        outputs : saved output file: filename_out
        '''

        output_directory = f'{params['home']}/{params['proc_tdl_loc']}/'

        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)

        filename_out = f'{output_directory}{params['tt_out_f_name']}_station_means.out'

        input_df.to_csv(filename_out, index=False) 

        return


    ############################################################################
    def compute_tt_diffs(self, input_df, diff_phases):
        '''
        Compute travel time differences between phases in diff_phases for dataframe input_df
        '''

        if len(diff_phases) != 2:
            print('Can only difference two phases')
            return pd.DataFrame()


        # Get unique event ids
        unique_evid = list(input_df['evid'].unique())
        n_evid = len(unique_evid)
        
        
        out_df = []
        # Loop through all unique evids
        if n_evid > 0:
            for ev, ev_id in enumerate(unique_evid):
                # Make ev data frame
                ev_df = input_df.query("evid == @ev_id")
                # Find unique stations in ev_df
                unique_stations = list(ev_df['nslc'].unique())
                n_stations = len(unique_stations)

                # Loop through all unique stations
                if n_stations > 0:
                    for s, station in enumerate(unique_stations):
                        # Make station data frame
                        st_df = ev_df.query("nslc == @station")

                        # Now find the phases from diff_phases.
                        ph1_df = st_df.query("phase == @diff_phases[0]")
                        ph2_df = st_df.query("phase == @diff_phases[1]")

                        if len(ph1_df) != 1 or len(ph2_df) != 1:
                            # Phase not found or Multipathing detected
                            # print('Phase not found or Multipathing...')
                            # sys.exit()
                            continue
                        
                        evid = str(ph1_df['evid'].iloc[0])
                        date_time = str(ph1_df['date_time'].iloc[0])
                        evt_lat = float(ph1_df['evt_lat'].iloc[0])
                        evt_lon = float(ph1_df['evt_lon'].iloc[0])
                        evt_dep = float(ph1_df['evt_dep'].iloc[0])
                        evt_mag = float(ph1_df['evt_mag'].iloc[0])
                        channel = str(ph1_df['channel'].iloc[0])
                        nslc = str(ph1_df['nslc'].iloc[0])
                        stla = float(ph1_df['stla'].iloc[0])
                        stlo = float(ph1_df['stlo'].iloc[0])
                        stel = float(ph1_df['stel'].iloc[0])
                        phase = f'{str(diff_phases[0])}-{str(diff_phases[1])}'
                        tderr = np.round(np.mean([ph1_df['tderr'].iloc[0], ph2_df['tderr'].iloc[0]]),3)
                        ccmx = np.round(np.mean([ph1_df['ccmx'].iloc[0], ph2_df['ccmx'].iloc[0]]),3)
                        tdelay_p1 = float(ph1_df['tdelay'].iloc[0])
                        ttaup_p1 = float(ph1_df['ttaup'].iloc[0])
                        tdelay_p2 = float(ph2_df['tdelay'].iloc[0])
                        ttaup_p2 = float(ph2_df['ttaup'].iloc[0])
                        tdelay = np.round(tdelay_p1 - tdelay_p2,3)
                        dist = float(ph1_df['dist'].iloc[0])
                        az = float(ph1_df['az'].iloc[0])
                        baz = float(ph1_df['baz'].iloc[0])
                        mid_lat = float(ph1_df['mid_lat'].iloc[0])
                        mid_lon = float(ph1_df['mid_lon'].iloc[0])

                        # Add values to list for dataframe
                        out_df.append([evid,date_time,evt_lat,evt_lon,evt_dep,evt_mag,channel,nslc,stla,stlo,stel,phase,tdelay,tderr,ccmx,tdelay_p1,ttaup_p1,tdelay_p2,ttaup_p2,dist,az,baz,mid_lat,mid_lon])

                else:
                    print('No unique stations found')
        
            # Make the dataframe
            column_names = ['evid','date_time','evt_lat','evt_lon','evt_dep','evt_mag','channel','nslc','stla','stlo','stel','phase','tdelay','tderr','ccmx','tdelay_p1','ttaup_p1','tdelay_p2','ttaup_p2','dist','az','baz','mid_lat','mid_lon']
            output_df = pd.DataFrame(sorted(out_df), columns = column_names)

        else:
            print('No unique evid found')
            return pd.DataFrame()

        return output_df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_tdl_summary_file(self, params, input_df):
        '''
        Writes travel time delay summary file.
        
        Parameters
        ----------
        params : dict parameter file
        input_df : dataframe of tdl measurements to save

        Returns
        -------
        outputs : saved output file: filename_out
        '''

        output_directory = f'{params['home']}/{params['proc_tdl_loc']}/'

        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)

        print(f'Sending output to: {str(output_directory)}')

        filename_out = f'{output_directory}{params['tt_out_f_name']}.out'

        input_df.to_csv(filename_out, index=False) 

        print(f'Written travel time delay summary file to: {str(filename_out)}')

        return

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_tt_diff_summary_file(self, params, input_df):
        '''
        Writes travel time difference delay summary file.
        
        Parameters
        ----------
        params : dict parameter file
        input_df : dataframe of tdl measurements to save

        Returns
        -------
        outputs : saved output file: filename_out
        '''

        output_directory = f'{params['home']}/{params['proc_tdl_loc']}/'

        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)
        filename_out = f'{output_directory}{params['tt_out_f_name']}_{str(params['diff_phases'][0])}-{str(params['diff_phases'][1])}.out'

        input_df.to_csv(filename_out, index=False) 

        print(f'Written travel time diffs delay summary file to: {str(filename_out)}')
        return

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def read_dataframe(self, file_path):
        '''
        Read csv file from filepath into dataframe
        '''
        try:
            with open(file_path, 'r') as file:
                output_df = pd.read_csv(file, delimiter = ',')
        except:
            #Return empy dataframe if not possible.
            output_df = pd.DataFrame()
        return output_df
        
############################################################################
#        main()
############################################################################

def main():
    from toolkit import Toolkit

    # Open the file and load the file
    tk = Toolkit()
    params = tk.get_params("params_in.yaml")

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

    filtered_df = process_tdl_files.apply_stuff(params, filter_functions)
    
    print(filtered_df)

    print('applied sort functions')
    # Write outfile
    process_tdl_files.write_tdl_summary_file(params, filtered_df)


    # Get station means per phase:

    stat_means_df = process_tdl_files.get_station_means(params, filtered_df)

    process_tdl_files.write_station_means(params, stat_means_df)


    # Find travel time differences
    print('Computing travel time differences for: '+str(params['diff_phases'][0])+' - '+str(params['diff_phases'][1]))
    tt_diff_output = process_tdl_files.compute_tt_diffs(filtered_df, params['diff_phases'])

    print(tt_diff_output)

    # # Write out file...
    process_tdl_files.write_tt_diff_summary_file(params, tt_diff_output)


if __name__ == '__main__':
    main()
    
