# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:25:02 2024
@author: alistairboyce
if this fails: >> conda activate env3.12
"""

import time, sys, os, glob, shutil
import numpy as np
from obspy import UTCDateTime, read, read_inventory, Stream
import obspy.signal
import obspy.signal.rotate
import obspy.geodetics.base
import obspy.geodetics
import pandas as pd
import concurrent.futures
import subprocess
import random
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

class FetchData:

    '''
    Class to handle serial/parallel download of seimic data and
    to format the data ready for MTWSPy travel time picking
    
    First we obtain dictionaries of CMT solutions and station details.
    Then apply a common series of functions to all data files downloaded by 
    the iris-FetchData scripts and save output to a formatted directory.
    Added complexity with rotation of components compared to 
    Synthetics that are already rotated in ZRT.

    '''

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_event_metadata(self, data_loc, year):
        '''
        Use Iris Event server to find all events for given year using 
        FetchEvent perl script. Try to read this into pandas df
        Add ev_id column -> formatted event name for saving purposes.

        Parameters:
        -------
        data_loc : str
            Location to save data
        year : int/str
            Year of data investigated.    

        Returns:
        -------
        cmt_table : df : event details in pandas dataframe
        '''

        # Create output directory
        if not os.path.exists(f'{data_loc}/e{year}'):
            os.makedirs(f'{data_loc}/e{year}', exist_ok=True)

        # Syntax:
        # ./FetchEvent -s {year}-01-01,01:01:00.000 -e {year}-12-31,23:59:59.999 --depth 0:700 --mag 5.5:7.5 --cat GCMT -X {data_loc}/e{year}/{year}_events.xml -o {data_loc}/e{year}/{year}_events.out

        # Execute Perl Scipt using subprocess
        service = f'export SERVICEBASE=http://service.iris.edu;'
        event = f'export EVENTWS=http://service.iris.edu/fdsnws/event/1;'
        event_metadata = [f'{service} {event} /home/aboyce/bin/FetchEvent -s {year}-01-01,00:00:00.000 -e {year}-01-01,23:59:59.999 --depth 0:700 --mag 5.5:7.5 --cat GCMT -o {data_loc}/e{year}/{year}_events.out']

        out = subprocess.check_output(event_metadata, shell=True)

        try:
            # Read cmt table, sort and re-index.
            cmt_table = pd.read_csv(f'{data_loc}/e{year}/{year}_events.out', delimiter="|", header=None, names = ['ev_num','ev_time','lat', 'lon', 'depth', 'cat', 'temp_cat', 'cmt_id', 'mag', 'location'])
            cmt_table = cmt_table.sort_values(by='ev_time').reset_index(drop=True)

            # Add ev_ids
            ev_df = pd.DataFrame(columns = ['ev_id'], index = range(0,len(cmt_table)))

            for index, row in cmt_table.iterrows(): 
                eq_time = UTCDateTime(row['ev_time'])
                # ev_df['ev_id'][index] = f'{eq_time.year:4d}{eq_time.month:02d}{eq_time.day:02d}{eq_time.hour:02d}{eq_time.minute:02d}{eq_time.second:02d}'
                ev_df.loc[index, 'ev_id'] = f'{eq_time.year:4d}{eq_time.month:02d}{eq_time.day:02d}{eq_time.hour:02d}{eq_time.minute:02d}{eq_time.second:02d}'

            cmt_table = pd.merge(cmt_table, ev_df, left_index = True, right_index = True)

        except:
            cmt_table = pd.DataFrame()
        
        return cmt_table


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_datacenters(self):
        '''
        Produce dataframe of datacenter names and their locations online.

        Returns:
        -------
        df : pd dataframe
        '''


        data = {
            "datacenter": ["AUSPASS", "GEOFON", "GEONET", "INGV", "IPGP", "IRIS", "KOERI", "LMU", "NIEP", "NOA", "ORFEUS", "RESIF", "UIB-NORSAR", "USP"],
            "loc": ["http://auspass.edu.au", "http://geofon.gfz-potsdam.de", "http://service.geonet.org.nz", "http://webservices.ingv.it", "http://ws.ipgp.fr", "http://service.iris.edu", "http://eida.koeri.boun.edu.tr", "http://erde.geophysik.uni-muenchen.de", "http://eida-sc3.infp.ro", "http://eida.gein.noa.gr", "http://www.orfeus-eu.org", "http://ws.resif.fr", "http://eida.geo.uib.no", "http://sismo.iag.usp.br"]
        }

        df = pd.DataFrame(data)

        return df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_station_metadata(self, data_loc, year, channel, dc_name, dc_loc_string):
        '''
        Use Given datacenter server to find all available stations for given
         year using FetchMetadata perl script. Try to read this into pandas df
        Keep most suitable channel for request based on heirachy/availability
        This reduces redundnacy of requests
        Download station XML file for use in instrument reponse removal

        Parameters:
        -------
        data_loc : str
            Location to save data
        year : int/str
            Year of data investigated.    
        channel : str
            Search channel
        dc_name : str
            datacenter name
        dc_loc_string : str
            web address of datacenter

        Returns:
        -------
        out_stat_df : df 
            station details in pandas dataframe
        inv_loc : str
            Location of station XML file for inst. resp.
        '''

        # Syntax:
        # FetchMetadata -C HH?,BH?,LH? -s 2011-01-01,00:00:00 -e 2011-01-01,01:00:00 -o IRIS_2011_stations.out
        # Export environment variables and Execute perl script:

        service = f'export SERVICEBASE={dc_loc_string};'
        meta = f'export METADATAWS={dc_loc_string}/fdsnws/station/1;'

        if channel == 'LH':
            channels = 'HH?,BH?,LH?'
        elif channel == 'BH':
            channels = 'HH?,BH?'
        elif channel == 'HH':
            channels = 'HH?'


        station_metadata = [f'{service} {meta} /home/aboyce/bin/FetchMetadata -C {channels} -s {year}-01-01,00:00:00 -e {year}-12-31,23:59:59 -o {data_loc}/e{year}/{dc_name}_{year}_stations.out -X {data_loc}/e{year}/{dc_name}_{year}.xml -resp']
        out = subprocess.check_output(station_metadata, shell=True)

        try:
            # Read station metadata into df
            sta_table = pd.read_csv(f'{data_loc}/e{year}/{dc_name}_{year}_stations.out', delimiter="|", dtype={'loc': 'str'})
            sta_table.rename(columns={"#net": "net"}, inplace=True)

            # Correct '00' errors
            sta_table['loc'] = sta_table['loc'].replace(0.0, '00')

            unique_net_sta = sta_table[['net', 'sta']].drop_duplicates()

            # Deal with the heirachy of requesting channels
            out_stat_df = pd.DataFrame(columns=sta_table.columns)

            for iter, row in unique_net_sta.iterrows():
                network = row['net']
                station = row['sta']
                
                tab = sta_table.query("(net == @network) and (sta == @station)")
                if channel == 'LH':
                    keep = tab[tab['chan'].str.contains('LH')]
                    if len(keep) == 0:
                        keep = tab[tab['chan'].str.contains('BH')]
                        if len(keep) == 0:
                            keep = tab[tab['chan'].str.contains('HH')]
                elif channel == 'BH':
                    keep = tab[tab['chan'].str.contains('BH')]
                    if len(keep) == 0:
                        keep = tab[tab['chan'].str.contains('HH')]
                elif channel == 'HH':
                    keep = tab[tab['chan'].str.contains('HH')]

                # now append keep to overall df
                if len(keep) > 0:
                    # Can have multiple start times, locs, channels BH1,2,E,N,Z
                    # Sort and keep first three.
                    if len(keep) > 3:
                        keep = keep.sort_values(['start', 'loc', 'chan'], ascending=[True, True, True])[0:3]

                    if len(out_stat_df) == 0:
                        out_stat_df = keep.copy()
                    else:
                        out_stat_df = pd.concat([out_stat_df, keep], ignore_index=True)

            # Remove nan locations and exchange for empty string.
            out_stat_df['loc'] = out_stat_df['loc'].replace(np.nan, '*', regex=True)

            # Read downloaded inventory for inst. resp. removal:
            inv_loc = f'{data_loc}/e{year}/{dc_name}_{year}.xml'
        
        except:
            out_stat_df = pd.DataFrame()
            inv_loc = f'{data_loc}/e{year}/temp.out'

        if not os.path.exists(inv_loc):
            print(f'Somehow we did not download a station xml file for {dc_name}, {year} at: {inv_loc}')

        return out_stat_df, inv_loc


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def split_dataframe_by_ev_id(self, df, packet_request_size):
        '''
        Takes large download request df and splits into chunks
        based on packet size request

        Parameters:
        -------
        df : df
            data frame of all requested files for given data center
        packet_request_size : int
            Size of each request sent to datacenter (traces)

        Returns:
        -------
        chunks : list
            List of lists with chunks for downloading.
        ev_ids : list
            ev_ids used (for directory names etc)
        ev_id_counts : list
            numerical count of ev_ids (for labelling of chunks)
        
        '''

        # Initialize an empty list to store chunks
        chunks = []
        ev_ids = []
        ev_id_counts = []

        # Initialize variables for chunk creation
        current_chunk = []
        current_ev_id = df['ev_id'][0]
        current_chunk_size = 0
        ev_id_count = 1

        # Iterate over rows in the sorted DataFrame
        for index, row in df.iterrows():

            # Check if current row has the same 'ev_id' as the current chunk
            if row['ev_id'] == current_ev_id and current_chunk_size < packet_request_size:
                # current_chunk.append(row)
                current_chunk.append(f'{row['net']} {row['sta']} {row['loc']} {row['chan']} {row['start']} {row['end']}')
                current_chunk_size += 1
            else:
                # Check if the current chunk has reached the packet_request_size
                if current_chunk_size >= packet_request_size or row['ev_id'] != current_ev_id:
                    # chunks.append(pd.DataFrame(current_chunk))
                    chunks.append(current_chunk)
                    ev_ids.append(current_ev_id)
                    ev_id_counts.append(ev_id_count)

                    if current_chunk_size >= packet_request_size:
                        ev_id_count += 1

                    if row['ev_id'] != current_ev_id:
                        ev_id_count = 1

                    current_chunk = []
                    current_chunk_size = 0


                # Start a new chunk with the current row
                # current_chunk.append(row)
                current_chunk.append(f'{row['net']} {row['sta']} {row['loc']} {row['chan']} {row['start']} {row['end']}')
                current_chunk_size += 1
                current_ev_id = row['ev_id']
        
        # Append the last chunk
        if current_chunk:
            # chunks.append(pd.DataFrame(current_chunk))
            chunks.append(current_chunk)
            ev_ids.append(current_ev_id)
            ev_id_counts.append(ev_id_count)

        return chunks, ev_ids, ev_id_counts 


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def write_chunk_file(self, data_loc, year, output_directory, ev_id_count, chunk):
        '''
        Write out the chunk request file to the correct location

        Parameters:
        -------
        data_loc : str : location of download
        year : int/str : year of download
        output_directory : str : ev_id formatted name
        ev_id_count : int : unique number of request for given earthquake
        chunk : list : request file contents

        Returns:
        ------
        chunk_f_name : str
        
        '''

        # mkdir
        if not os.path.exists(f'{data_loc}/e{year}/py_formatted/{output_directory}'):
            os.makedirs(f'{data_loc}/e{year}/py_formatted/{output_directory}', exist_ok=True)

        chunk_f_name = f'{data_loc}/e{year}/py_formatted/{output_directory}/{output_directory}_chunk_request_{ev_id_count:03d}.txt'

        # Write chunk to file
        outfile = open(chunk_f_name,'w')
        for j in range(len(chunk)):
            outfile.write(chunk[j]+'\n')
        outfile.close()

        return chunk_f_name


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def process_one_chunk(self, input_dict):
        '''
        Attempts processing for one chunk file
        Attempts download of mseed file, and reads into obspy if available
        Processes one file reading the file, applying the list functions
        Reads station inventory into memory for instrument response removal
        Applies list of functions including instrument response removal
        Moves Mseed and request file to tidy location


        Parameters
        ----------
        input_dicts : list (of dictionaries as below)
            input_dict : dict : specific to each file
                input dictionary containing:
                    file : file name to process.
                    sacfile_dict : dict of parameters describing file
                    output_directory : formatted output dir
                    cmt_dict : dict CMTSOLUTION params from Specfem execution
                    stations_df : STATION file from Specfem execution as pd df
                    functions : list of functions to execute on sacfile
                    do_parallel : boolean
                    jobs : number of jobs
                    channel : desired output data channel

        Returns
        ----------
        Nothing: None
        '''

        chunk_f_name = input_dict['chunk_f_name']
        dc_loc_string = input_dict['dc_loc_string']
        output_directory = input_dict['output_directory']
        data_loc = input_dict['data_loc']
        year = input_dict['year']
        inv_loc = input_dict['inv_loc']
        functions = input_dict['functions']

        mseed_f_name = chunk_f_name.split('.')[0] + '.mseed'

        service = f'export SERVICEBASE={dc_loc_string};'
        meta = f'export METADATAWS={dc_loc_string}/fdsnws/station/1;'
        time = f'export TIMESERIESWS={dc_loc_string}/fdsnws/dataselect/1;'

        # Syntax:
        # FetchData -l myselection.txt -o mydata.mseed
        # Execute perl script using subprocess:
        station_data = [f'{service} {time} {meta} /home/aboyce/bin/FetchData -l {chunk_f_name} -o {mseed_f_name}']
        out = subprocess.check_output(station_data, shell=True)

        # Add inventory to the dict
        input_dict['inventory'] = read_inventory(inv_loc)

        if os.path.exists(mseed_f_name):
            # Read downloaded mseed file
            try:
                sta_inv = read(mseed_f_name)
            except:
                print(f'Cannot read .mseed file at {mseed_f_name}')
                return

            if sta_inv:
                print('station inventory read')
                for seis in sta_inv:
                    fail = 0     
                    if len(functions) ==  0: 
                        # NO_FUNCTIONS_TO_APPLY.....
                        return

                    else:
                        # apply functions, only executed when fail = 0
                        for function in functions:
                            input_dict, seis, fail = function(input_dict, seis, fail)


            # Tidy up the mseed files.
            tidy_loc = f'{data_loc}/e{year}/py_formatted/{output_directory}/request_files/'
            if not os.path.exists(tidy_loc):
                os.makedirs(tidy_loc, exist_ok=True)      
            try:
                # Move mseed file to tidy loc
                shutil.move(mseed_f_name, tidy_loc)
            except:
                pass

            try:
                # Move chunk file to tidy loc
                shutil.move(chunk_f_name, tidy_loc)
            except:
                pass

        else:
            print(f'Cannot find .mseed file at {mseed_f_name}')

        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_data(self, data_loc, year, channel, main_function, function_list, packet_request_size, do_parallel, jobs, processing_channel):
        '''
        Main function for data downloadn and post processing.
        Collects, Events, Datacenters, Stations, chunks
        Collates in input_dicts
        Sends to main_function for processing.

        Parameters:
        -------
        data_loc : str
            Location to save data
        year : int/str
            Year of data investigated.    
        channel : str
            Search channel
        main_function : list
            Name of main function operating on each member on input_dicts
        function_list : list
            List of functions to apply to data following download
        packet_request_size : int
            Size of each request sent to datacenter (traces)    
        do_parallel : bool
            Do parallel processing or not
        jobs : int
            Number of parallel workers
        processing_channel : str
            Output channel of processing


        Returns:
        -------
        None
        '''

        # Get events
        cmt_table = self.get_event_metadata(data_loc, year)
        
        if cmt_table.empty:
            print(f'No events found for {year}... Return')
            return

        # Get datacenter details for download
        datacenters = self.get_datacenters()

        input_dicts = []

        # Loop through each datacenter.
        for i, dc_row in datacenters.iterrows():
            
            # Data center parameters:
            dc_name = dc_row['datacenter']
            dc_loc_string = dc_row['loc']

            # Get station metadata
            out_stat_df, inv_loc = self.get_station_metadata(data_loc, year, channel, dc_name, dc_loc_string)

            if not out_stat_df.empty:
                # Loop through CMT_table and then out_stat_df and create a combined request_df.
                column_headers = ['ev_time','ev_id','net','sta','loc','chan','start','end']
                request_list = []
                for j, cmt_row in cmt_table.iterrows():
                    ev_time = UTCDateTime(cmt_row['ev_time'])
                    ev_id = cmt_row['ev_id']
                    s_time = ev_time - 60
                    e_time = ev_time + 7200
                    start = f'{s_time.year:4d}-{s_time.month:02d}-{s_time.day:02d}T{s_time.hour:02d}:{s_time.minute:02d}:{s_time.second:02d}'
                    end = f'{e_time.year:4d}-{e_time.month:02d}-{e_time.day:02d}T{e_time.hour:02d}:{e_time.minute:02d}:{e_time.second:02d}'

                    for k, sta_row in out_stat_df.iterrows():
                        request_list.append([ev_time, ev_id, sta_row['net'], sta_row['sta'], sta_row['loc'], sta_row['chan'], start, end])


                request_df = pd.DataFrame(request_list, columns = column_headers)

                # Split request_df into chunks 
                chunks, ev_ids, ev_id_counts = self.split_dataframe_by_ev_id(request_df, packet_request_size)

                # Collate chunks and other info into input_dicts
                if len(chunks) == len(ev_ids) and len(ev_ids) == len(ev_id_counts):

                    for i in range(len(chunks)):
                        input_dict = {}
                        chunk = chunks[i]
                        output_directory = ev_ids[i]
                        ev_id_count = ev_id_counts[i]

                        chunk_f_name = self.write_chunk_file(data_loc, year, output_directory, ev_id_count, chunk)

                        input_dict['data_loc'] = data_loc
                        input_dict['year'] = year
                        input_dict['channel'] = channel
                        input_dict['processing_channel'] = processing_channel
                        input_dict['packet_request_size'] = packet_request_size
                        input_dict['functions'] = function_list
                        input_dict['dc_name'] = dc_name
                        input_dict['dc_loc_string'] = dc_loc_string
                        input_dict['inv_loc'] = inv_loc
                        input_dict['output_directory'] = output_directory
                        input_dict['chunk_f_name'] = chunk_f_name
                        input_dict['cmt_line'] = cmt_table.query('ev_id == @output_directory')
                        input_dict['out_stat_df'] = out_stat_df
                        input_dicts.append(input_dict)

        # Send input dicts to main_function for processing.
        if input_dicts:

            # Shuffle the input_dicts so we dont overload data centers
            random.shuffle(input_dicts)

            ############## PARALLEL ################
            if do_parallel:

                # Parallel processing
                for main_f in main_function:
                    with concurrent.futures.ProcessPoolExecutor(max_workers = jobs) as executor:
                        results = executor.map(main_f, input_dicts)

                    return
        
            ############## SERIES ################
            else:
                # Serial processing
                for input_dict in input_dicts:
        
                    for main_f in main_function:
                        res = main_f(input_dict)
        
            return
        
        else:
            pass


        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def update_headers_seis(self, input_dict, seis, fail):
        '''
        Updates all headers of the seismogram in obspy and sac formats

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            ###
            # print(f'Update headers seis: {file}')

            ev_time = UTCDateTime(input_dict['cmt_line']['ev_time'].values[0])

            seis.stats['sac'] = {}
            seis.stats['sac']['o'] = UTCDateTime(ev_time) - UTCDateTime(seis.stats['starttime'])

            seis.stats['sac']['evla'] =  input_dict['cmt_line']['lat'].values[0]
            seis.stats['sac']['evlo'] =  input_dict['cmt_line']['lon'].values[0]
            seis.stats['sac']['evdp'] =  input_dict['cmt_line']['depth'].values[0]
            seis.stats['sac']['mag'] =  input_dict['cmt_line']['mag'].values[0].split(',')[-1]
            seis.stats['sac']['kevnm'] =  input_dict['cmt_line']['cmt_id'].values[0].split(',')[-1]
            seis.stats['sac']['lovrok'] =  1
            seis.stats['sac']['lcalda'] =  1

            station = seis.stats['station']
            network = seis.stats['network']
            channel = seis.stats['channel']
            location = seis.stats['location']

            seis.stats['sac']['kstnm'] = station
            seis.stats['sac']['knetwk'] = network
            seis.stats['sac']['kcmpnm'] = channel
            seis.stats['sac']['khole'] = location

            st_df = input_dict['out_stat_df'].query("(net == @network) and (sta == @station) and (chan  == @channel)")

            seis.stats['sac']['stla'] = st_df['lat'].values[0]
            seis.stats['sac']['stlo'] = st_df['lon'].values[0]
            seis.stats['sac']['stel'] = st_df['elev'].values[0]
            seis.stats['sac']['stdp'] = st_df['depth'].values[0]

            seis.stats['sac']['cmpaz'] = st_df['azimuth'].values[0]
            seis.stats['sac']['cmpinc'] = st_df['dip'].values[0]


        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def rem_inst_resp(self, input_dict, seis, fail):
        '''
        Removes instrument response from seismogram
        using appropriate pre-filter and station xml file
        outputs to displacement

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:

            if "LH" in seis.stats['channel']:
                # For LH data:
                pre_filt = [0.01, 0.02, 0.1, 0.4]

            if "BH" in seis.stats['channel']:
                # For BH data:
                pre_filt = [0.01, 0.02, 10, 20]

            if "HH" in seis.stats['channel']:
                # For HH data:
                pre_filt = [0.01, 0.02, 40, 60]

            #### OUTPUT TO DISPLACEMENT ###
            try:
                seis.remove_response(inventory=input_dict['inventory'], pre_filt=pre_filt, output="DISP", water_level=None)
                # print(f'rem_inst_resp for: {seis}')
                fail = 0
            except:
                print(f'Failed to remove response for: {seis}')
                fail = 1

        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def remove_trend_taper(self, input_dict, seis, fail):
        '''
        Removes trend and tapers seismogram

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            ###
            # print(f'Remove trend seis: {file}')
            
            seis.detrend()
            seis.taper(max_percentage=0.015, type='cosine')

            # print(f'remove_trend_taper for: {seis}')
        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def filter_seis(self, input_dict, seis, fail):
        '''
        Filter seismogram

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            ###
            # print(f'Filter seis: {file}')

            if input_dict['processing_channel'] == 'LH':
                fmin = 0.01
                fmax = 0.4
            else:
                fmin = 0.01
                fmax = 20
            seis.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=2,zerophase=True)
            # print(f'filter_seis for: {seis}')
        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def interpolate_seis(self, input_dict, seis, fail):
        '''
        Interpolate seismogram

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''

        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            # print(f'Interpolate seis: {file}')
            ###
            if input_dict['processing_channel'] == 'LH':
                delta = 1
            else:
                delta = 0.05
            seis.resample(1 / delta)
            # print(f'interpolate_seis for: {seis}')
        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def save_seis(self, input_dict, seis, fail):
        '''
        Saves processed sac file

        Inputs and returns:
            input_dict : dict
            file : seismogram file name
            seis : obspy stream
            fail : 1/0
        '''
        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:

            output_directory = input_dict['output_directory']
            data_loc = input_dict['data_loc']
            year = input_dict['year']
            processing_channel = input_dict['processing_channel']

            net = seis.stats['network']
            sta = seis.stats['station']
            # loc = seis.stats['location']
            # if not loc:
            #     loc = '00'

            # Change all Loc t0 '00' to match specfem:
            loc = '00'

            chan = processing_channel + seis.stats['channel'][-1]

            save_dir = f'{data_loc}/e{year}/py_formatted/{output_directory}'

            if not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)

            out_filename = f'{save_dir}/{output_directory}_{net}_{sta}.{loc}.{chan}'

            try:
                # print(f'Function Save seismogram: {file}')


                seis.stats['channel'] = chan
                seis.stats['sac']['kcmpnm'] = chan

                seis.stats['location'] = loc
                seis.stats['sac']['khole'] = loc

                print(f'Saving {out_filename}')
                seis.write(out_filename, format='sac')
                fail = 0
                # print(f'save_seis for: {seis}')
            except:
                fail = 1
                print(f'FAILED to save: {out_filename}')

        return input_dict, seis, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    # \\\\\\\\\\\\\\\\\\\\\      ROTATION         \\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_input_dicts_rotate(self, year, channel, data_loc, function_list, do_parallel, jobs, processing_channel):
        '''
        Parallel processing needs a list of inputs for each iteration
        So put all input into a dict and create a list of dicts

        Parameters
        ----------
        input_directory: str
            input_directory name where all the data is.
        functions : list
            list of functions to apply to inputs
        do_parallel : bool
        jobs :  int
            number of paralle workers
        channel : str
            required output channel of seismograms
            
        Returns
        ----------
        input_dicts : list (of dictionaries as below)
            input_dict : dict : specific to each file
                input dictionary containing:
                    file_v : file name of vertical component to process.
                    horizontals : list of horizontal file names to process
                    output_directory : formatted output dir
                    functions : list of functions to execute on sacfile
                    do_parallel : boolean
                    jobs : number of jobs
                    channel : desired output data channel
        '''


        ev_list = glob.glob(f'{data_loc}/e{year}/py_formatted/{year}??????????')
        input_dicts = []
        
        for output_directory in ev_list:

            
            # Finds files, makes input dictionary of files, functions, inputs
                
            vert_files = sorted(glob.glob(f'{output_directory}/{year}??????????_*_*.{processing_channel}Z'))
            num_files = len(vert_files)


            if num_files == 0:
                # No files found...
                print('----------')
                print('----------////   NO-VERT-FILES-IN: '+str(output_directory)+'    ////----------')
            else:

                print(f'----------')
                print(f'Processing list of {num_files} files...')
                print(f'----------')

                for file in vert_files:

                    horizontals = sorted(glob.glob(file[:-1] + '*'))
                    horizontals.remove(file)

                    input_dict = {}
                    input_dict['file_v'] = file
                    input_dict['horizontals'] = horizontals
                    input_dict['output_directory'] = output_directory
                    input_dict['functions'] = function_list
                    input_dict['do_parallel'] = do_parallel
                    input_dict['jobs'] = jobs
                    input_dict['channel'] = channel
                    input_dict['processing_channel'] = processing_channel
                    
                    input_dicts.append(input_dict)

                # print('----------\n')

        return input_dicts


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def apply_serial(self, data_loc, year, channel, main_function, function_list, do_parallel, jobs, processing_channel):
        '''
        Takes a list of functions, reads all input files, applies each 
        function to all files.
        IN SERIES...

        Parameters
        ----------
        main_function : list
            Length one list containing main function.
        input_directory: str
            input_directory name where all the data is.
        functions : list
            list of functions to apply to inputs
        do_parallel : bool
        jobs :  int
            number of paralle workers
        channel : str
            required output channel of seismograms

        Returns
        ----------
        Nothing : None
            Functions take care of returns/writing to file  
            
        The functions must return all input variables
        Functions can write to logfile & outfile.
        '''

        input_dicts = self.get_input_dicts_rotate(year, channel, data_loc, function_list, do_parallel, jobs, processing_channel)

        if input_dicts:

            # Serial processing
            for input_dict in input_dicts:
                
                for main_f in main_function:
                    res = main_f(input_dict)
        
            return
        else:
            pass


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def apply_parallel(self, data_loc, year, channel, main_function, function_list, do_parallel, jobs, processing_channel):
        '''
        Takes a list of functions, reads all input files, applies 
        each function to all of the files.
        IN PARALLEL...

        Parameters
        ----------
        main_function : list
            Length one list containing main function.
        input_directory: str
            input_directory name where all the data is.
        functions : list
            list of functions to apply to inputs
        do_parallel : bool
        jobs :  int
            number of paralle workers
        channel : str
            required output channel of seismograms

        Returns
        ----------
        Nothing : None
            Functions take care of returns/writing to file  
        
        The functions must return all input variables
        Functions can write to logfile & outfile.
        '''
        input_dicts = self.get_input_dicts_rotate(year, channel, data_loc, function_list, do_parallel, jobs, processing_channel)

        if input_dicts:
            # Parallel processing
            for main_f in main_function:

                with concurrent.futures.ProcessPoolExecutor(max_workers = jobs) as executor:
                    results = executor.map(main_f, input_dicts)

                return
        else:
            pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def execute(self, data_loc, year, channel, main_function, function_list, do_parallel, jobs, processing_channel):
        '''
        Takes main function inputs and passes to serial or parallel 
        executer, see below.
        '''
        if do_parallel:
            # Execute in parallel
            self.apply_parallel(data_loc, year, channel, main_function, function_list, do_parallel, jobs, processing_channel)
        else:
            # Execute in series
            self.apply_serial(data_loc, year, channel, main_function, function_list, do_parallel, jobs, processing_channel)
        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def process_one_station(self, input_dict):
        '''
        Processes one station reading the verticals and horizontals and 
        applying the list functions to them

        Parameters
        ----------
        input_dicts : list (of dictionaries as below)
            input_dict : dict : specific to each file
                input dictionary containing:
                    file_v : file name of vertical component to process.
                    horizontals : list of horizontal file names to process
                    output_directory : formatted output dir
                    functions : list of functions to execute on sacfile
                    do_parallel : boolean
                    jobs : number of jobs
                    channel : desired output data channel

        The functions in the list must take the following inputs: 

        input_dict : dict
        file_vertical, file_east, file_north : file location strings
        vert_component, east_component, north_component : obspy stream object
            onwhich operations are performed.
        fail : int - fail flag = 1 if function fails, else 0.

        Returns
        ----------
        Nothing: None
        '''

        file_vertical = input_dict['file_v']
        functions = input_dict['functions']
        horizontals = input_dict['horizontals']

        if len(horizontals) >= 2: 
            # read sac files
            try:
                # Read components
                vert_component = read(file_vertical,'sac')

                if str(file_vertical[:-1]) + 'E' in horizontals:
                    file_east = str(file_vertical[:-1]) + 'E'
                    east_component = read(file_east, format='sac')
                else:
                    file_east = str(file_vertical[:-1]) + '1'
                    east_component = read(file_east, format='sac')

                if str(file_vertical[:-1]) + 'N' in horizontals:
                    file_north = str(file_vertical[:-1]) + 'N'
                    north_component = read(file_north, format='sac')
                else:
                    file_north = str(file_vertical[:-1]) + '2'
                    north_component = read(file_north, format='sac')

                fail = 0
                print(f'Read seis {file_vertical}, {file_east} and {file_north}...')
            except:
                #  FAILED_LOADING...
                fail = 1
                return

            if len(functions) ==  0: 
                # NO_FUNCTIONS_TO_APPLY.....
                return

            else:
                # apply functions, only executed when fail = 0
                for function in functions:
                    input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail = function(input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail)

        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def check_length(self, input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
        '''
        Check trace length
        Often traces are one sample longer than the others
        Test this trimming.

        Inputs and returns:
            input_dict : dict
            file_vertical, file_east, file_north : seismogram file names
            vert_component, east_component, north_component : obspy streams
            fail : 1/0
        '''
        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            
            # Check trace length and sample rate
            try:
                if vert_component[0].stats.npts % 10 == 1 & east_component[0].stats.npts % 10 == 1 & north_component[0].stats.npts % 10 == 1:
                    # print('Caught XXXXX1 sample points')
                    vert_component[0].data=vert_component[0].data[:-1]
                    east_component[0].data=east_component[0].data[:-1]
                    north_component[0].data=north_component[0].data[:-1]
                
                # print(f'Passed length QC....')
                fail = 0
            except:
                fail = 1
                print(f'Failed length QC....')
            
        return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def check_sample_rate(self, input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
        '''
        Check sample rate - required for component rotation

        Inputs and returns:
            input_dict : dict
            file_vertical, file_east, file_north : seismogram file names
            vert_component, east_component, north_component : obspy streams
            fail : 1/0
        '''
        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:

            try:
                if input_dict['processing_channel'] == 'LH':
                    if vert_component[0].stats.sampling_rate < 1.0 or east_component[0].stats.sampling_rate < 1.0 or north_component[0].stats.sampling_rate < 1.0:
                        vert_component.resample(1)
                        east_component.resample(1)
                        north_component.resample(1)
                        # print('Some issue with sample rate <1.0 for channel: ', input_dict['processing_channel'])
                if input_dict['processing_channel'] == 'BH':
                    if vert_component[0].stats.sampling_rate < 0.02 or east_component[0].stats.sampling_rate < 0.02 or north_component[0].stats.sampling_rate < 0.02:
                        vert_component.resample(0.02)
                        east_component.resample(0.02)
                        north_component.resample(0.02)
                        # print('Some issue with sample rate <0.02 for channel: ', input_dict['processing_channel'])
                if input_dict['processing_channel'] == 'HH':
                    if vert_component[0].stats.sampling_rate < 0.005 or east_component[0].stats.sampling_rate < 0.005 or north_component[0].stats.sampling_rate < 0.005:
                        vert_component.resample(0.005)
                        east_component.resample(0.005)
                        north_component.resample(0.005)
                        # print('Some issue with sample rate <0.005 for channel: ', input_dict['processing_channel'])
                
                # print(f'Passed sample rate QC....')
                fail = 0     
            except:
                fail = 1
                print(f'Failed sample rate QC....')
            
        return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def trim_components(self, input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
        '''
        Check Component length prior to rotation - need same length

        Inputs and returns:
            input_dict : dict
            file_vertical, file_east, file_north : seismogram file names
            vert_component, east_component, north_component : obspy streams
            fail : 1/0
        '''
        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            # Trim files:



            # find streams with min/max start/end times
            start_cut = UTCDateTime(vert_component[0].stats.starttime)
            end_cut = UTCDateTime(vert_component[0].stats.endtime)

            if UTCDateTime(east_component[0].stats.starttime) > start_cut : 
                start_cut = UTCDateTime(east_component[0].stats.starttime)

            if UTCDateTime(north_component[0].stats.starttime) > start_cut : 
                start_cut = UTCDateTime(north_component[0].stats.starttime)

            if UTCDateTime(east_component[0].stats.endtime) > end_cut : 
                end_cut = UTCDateTime(east_component[0].stats.endtime)

            if UTCDateTime(north_component[0].stats.endtime) > end_cut : 
                end_cut = UTCDateTime(north_component[0].stats.endtime)

            # Cut components to the same length
            # If this cant happen now its due to gappy waveforms of varying lengths.
            try:
                vert_component.trim(starttime=start_cut, endtime=end_cut, pad=True, fill_value=0)
                east_component.trim(starttime=start_cut, endtime=end_cut, pad=True, fill_value=0)
                north_component.trim(starttime=start_cut, endtime=end_cut, pad=True, fill_value=0)
                # print('Trimmed to :' + str(vert_component[0].stats['npts']))

                # vert_component.merge(fill_value='interpolate')
                # east_component.merge(fill_value='interpolate')
                # north_component.merge(fill_value='interpolate')

                # print(f'Passed Trimming....')
                fail = 0
            except:
                print('Trace error in seisN/seisE/seisZ')
                print('Station: ' + str(file_vertical[:-1]) + '\n')
                print('Issues with starttime/endtime/npts' + '\n')
                print(f'Failed Trimming....')
                fail = 1

        return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def rotate_components(self, input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
        '''
    Rotate components to ZRT & save

        Inputs and returns:
            input_dict : dict
            file_vertical, file_east, file_north : seismogram file names
            vert_component, east_component, north_component : obspy streams
            fail : 1/0
        '''
        if fail:
            ###
            # toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  -SKIPPING-')
            ###
            fail = 1
        else:
            
            # Check trace length and sample rate
            # try:

            loc_Z=vert_component[0].stats.location
            seisZ = vert_component.select(channel='*HZ',location=loc_Z)
            seisN_channel  = north_component[0].stats['channel']
            seisE_channel  = east_component[0].stats['channel']
            seisN = north_component.select(channel=seisN_channel,location=loc_Z)
            seisE = east_component.select(channel=seisE_channel,location=loc_Z)


            ori_seisN_temp=seisN[0].stats['sac']['cmpaz'] # .stats['orientation']
            ori_seisE_temp=seisE[0].stats['sac']['cmpaz'] # .stats['orientation']
            dip_seisE_temp=seisE[0].stats['sac']['cmpinc'] # .stats['dip']
            dip_seisN_temp=seisN[0].stats['sac']['cmpinc'] # .stats['dip']
            ori_seisZ_temp=seisZ[0].stats['sac']['cmpaz'] # .stats['orientation']
            dip_seisZ_temp=seisZ[0].stats['sac']['cmpinc'] # .stats['dip']

            # make sure the orientations are different and dip is not vertical
            if ori_seisN_temp == ori_seisE_temp:            
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisE_temp == -90:
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisN_temp == -90:
                print('Trace error in seisN/seisE')
                print('Failed on N/E component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail
            elif dip_seisZ_temp != -90.0:
                print('Trace error in seisZ')
                print('Failed on Z component')
                fail = 1
                return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


            EVLA = seisZ[0].stats['sac']['evla']
            EVLO = seisZ[0].stats['sac']['evlo']
            EVDP = seisZ[0].stats['sac']['evdp']
            STLA = seisZ[0].stats['sac']['stla']
            STLO = seisZ[0].stats['sac']['stlo']
            STEL = seisZ[0].stats['sac']['stel']
            STATION = seisZ[0].stats['station']
            NETWORK = seisZ[0].stats['network']

            DISTM, AZ, BAZ = obspy.geodetics.base.gps2dist_azimuth(EVLA, EVLO, STLA, STLO)
            DISTDG = DISTM / (6371.e3 * np.pi / 180.)

            # Check North and East exist
            if len(seisN) == 1 and len(seisE) == 1:
                
                # Make sure components are of the same length
                if len(seisN[0].data) > len(seisE[0].data):
                    seisN[0].data = seisN[0].data[:-1]
                if len(seisN[0].data) < len(seisE[0].data):
                    seisE[0].data = seisE[0].data[:-1]
                if len(seisZ[0].data) > len(seisE[0].data) and len(seisZ[0].data) > len(seisN[0].data):
                    seisZ[0].data = seisZ[0].data[:-1]  

                # Final check of length - if not remove.
                if len(seisZ[0].data) != len(seisE[0].data) or len(seisZ[0].data) != len(seisN[0].data):
                    print('Trace length error in seisZ/seisN/seisE')
                    print('Length seisZ[0].data: '+str(len(seisZ[0].data)))
                    print('Length seisN[0].data: '+str(len(seisN[0].data)))
                    print('Length seisE[0].data: '+str(len(seisE[0].data)))

                    fail = 1
                    return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail

                # Make sure North component is North and East component is East
                if 45 < ori_seisN_temp < 135 or 225 < ori_seisN_temp < 315:
                    seisN1 = seisE
                    seisE1 = seisN
                    seisN = seisN1
                    seisE = seisE1

                # print('Orientation N: ' + str(ori_seisN_temp))
                # print('Orientation E: ' + str(ori_seisE_temp))

                # print('Rotating to true N/E')
                [seisZtmp, seisNtmp, seisEtmp] = obspy.signal.rotate.rotate2zne(
                    seisZ[0].data, ori_seisZ_temp, dip_seisZ_temp,
                    seisN[0].data, ori_seisN_temp, dip_seisN_temp,
                    seisE[0].data, ori_seisE_temp, dip_seisE_temp)
                seisN = seisN[0].copy() 
                seisN.stats['channel'] = str(input_dict['processing_channel']) + 'N'
                seisN.data = seisNtmp
                seisE = seisE[0].copy()
                seisE.stats['channel'] = str(input_dict['processing_channel']) + 'E'
                seisE.data = seisEtmp    

                # should be setting the new orientation and dip to 0,90,-90 etc

                # rotate components to from North and East to Radial and Transverse
            
                # print('Rotating to RT using BAZ: ',BAZ)
                [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(
                    seisN.data, seisE.data, BAZ)
                seisR = seisN.copy()
                seisT = seisN.copy()

                seisR.stats['channel'] = str(input_dict['processing_channel']) + 'R'
                seisT.stats['channel'] = str(input_dict['processing_channel']) + 'T'
            
                seisR.data = seisRtmp
                seisT.data = seisTtmp


                ########## Sort out SAC headers #################

                # Radial
                # seisR.stats['channel'] = str(input_dict['processing_channel']) + 'R'
                seisR.stats['sac']['cmpaz']=AZ
                seisR.stats['sac']['az']=AZ
                seisR.stats['sac']['baz']=BAZ
                seisR.stats['sac']['gcarc']=DISTDG
                seisR.stats['sac']['kcmpnm']=str(input_dict['processing_channel'])+'R'

                # Tangential
                # seisT.stats['channel'] = str(input_dict['processing_channel']) + 'T'
                seisT.stats['sac']['cmpaz']=AZ+90.0
                seisT.stats['sac']['az']=AZ
                seisT.stats['sac']['baz']=BAZ
                seisT.stats['sac']['gcarc']=DISTDG
                seisT.stats['sac']['kcmpnm']=str(input_dict['processing_channel'])+'T'

                # East
                # seisE.stats['channel']=str(channel)+'E'
                seisE.stats['sac']['cmpaz']=90.0
                seisE.stats['sac']['az']=AZ
                seisE.stats['sac']['baz']=BAZ
                seisE.stats['sac']['gcarc']=DISTDG
                seisE.stats['sac']['kcmpnm']=str(input_dict['processing_channel'])+'E'

                # North
                # seisN.stats['channel']=str(channel)+'N'
                seisN.stats['sac']['cmpaz']=0.0
                seisN.stats['sac']['az']=AZ
                seisN.stats['sac']['baz']=BAZ
                seisN.stats['sac']['gcarc']=DISTDG
                seisN.stats['sac']['kcmpnm']=str(input_dict['processing_channel'])+'N'

                # Vertical
                seisZ[0].stats['sac']['az']=AZ
                seisZ[0].stats['sac']['baz']=BAZ
                seisZ[0].stats['sac']['gcarc']=DISTDG


                filename_E = str(file_vertical[:-1]) + 'E'
                filename_N = str(file_vertical[:-1]) + 'N'
                filename_V = str(file_vertical[:-1]) + 'Z'
                filename_R = str(file_vertical[:-1]) + 'R'
                filename_T = str(file_vertical[:-1]) + 'T'

                # Write over sac files as we now have everything sorted :)
                seisE.write(filename_E, 'sac')
                seisN.write(filename_N, 'sac')
                seisZ.write(filename_V, 'sac')
                seisR.write(filename_R, 'sac')
                seisT.write(filename_T, 'sac')

            
                # print(f'Passed Rotation....')
                fail = 0
            # except:
            #     print(f'Failed Rotation....')
            #     fail = 1

        return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def cleanup_files(self, input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail):
        '''
        Clean up unused files

        Inputs and returns:
            input_dict : dict
            file_vertical, file_east, file_north : seismogram file names
            vert_component, east_component, north_component : obspy streams
            fail : 1/0
        '''

        # Remove all that are unnecessary - so ignore the fail flag
        try:
                        
            # Make a directory for unused files...
            Not_useable = input_dict['output_directory'] + '/Not_useable/'
            if not os.path.exists(Not_useable):
                os.makedirs(Not_useable, exist_ok=True)

            stat = file_vertical[:-3]

            # Keep ZRT components
            to_delete = sorted(glob.glob(f'{stat}*'))

            if file_vertical[:-1]+'Z' in to_delete:
                to_delete.remove(file_vertical[:-1]+'Z')
            if file_vertical[:-1]+'R' in to_delete:
                to_delete.remove(file_vertical[:-1]+'R')
            if file_vertical[:-1]+'T' in to_delete:
                to_delete.remove(file_vertical[:-1]+'T')

            if len(to_delete) > 1:
                for file in to_delete:
                    try:
                        # print('Moving :' +str(file)+' -> '+str(Not_useable))
                        shutil.move(file, Not_useable)
                    except:
                        pass

                # print(f'Passed Cleanup....')
                fail = 0
        except:
            print(f'Failed Cleanup of {stat}* ....')
            fail = 1

        return input_dict, file_vertical, file_east, file_north, vert_component, east_component, north_component, fail


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def main():

    fetch_data = FetchData()

    start_time = time.time()

    if len(sys.argv) != 3:
        print("Please specify year and search channel required")
        print("USAGE:   python fetch_data.py <year> <channel>")
        print("EXAMPLE: python fetch_data.py 2008 LH")
        sys.exit()
    else:
        year = int(sys.argv[1])
        if not os.path.isdir('e' + str(year)):
            print("Directory given is not present.... ")
            os.makedirs('e' + str(year))

        channel = str(sys.argv[2])
        print(f"Downloading {year}, for search data channel: {channel}")

    data_loc = '/home/aboyce/d_data_obs/iris_FetchData'
    processing_channel = 'LH' # output of processing, channel -> processing_channel: e.g., LH,BH,HH -> LH
    packet_request_size = 100 # Number of files requested in each packet
    jobs = 24 # Number of independent/parallel requests executed 
    if jobs > 1 and jobs <= 96:
        do_parallel = True

    # Find earthquakes, stations, make most appropriate request files. 
    # (Use FetchEvent, FetchMetaData)
    # Request using FetchData, process directly in obspy 
    # e.g., remove inst. resp save to py_formatted.

    main_function = [fetch_data.process_one_chunk]

    function_list = [fetch_data.rem_inst_resp, 
                    fetch_data.update_headers_seis, 
                    fetch_data.remove_trend_taper, 
                    fetch_data.filter_seis, 
                    fetch_data.interpolate_seis, 
                    fetch_data.save_seis]

    fetch_data.get_data(data_loc, year, channel, main_function, 
                        function_list, packet_request_size, do_parallel, 
                        jobs, processing_channel)


    # rotate available components in py_formatted to ZRT
    # Clean up remaining files
    main_function = [fetch_data.process_one_station]

    function_list = [fetch_data.check_length, 
                    fetch_data.check_sample_rate, 
                    fetch_data.trim_components, 
                    fetch_data.rotate_components, 
                    fetch_data.cleanup_files]

    # Execute functions to Rotate the data to ZRT
    fetch_data.execute(data_loc, year, channel, main_function, 
                        function_list, do_parallel, 96, processing_channel)

   
    # Sometimes there remains some N & E files, since they have no associated vertical component.
    del_list = glob.glob(f'{data_loc}/e{year}/py_formatted/{year}??????????/{year}??????????_??_*.{processing_channel}E') + 
                glob.glob(f'{data_loc}/e{year}/py_formatted/{year}??????????/{year}??????????_??_*.{processing_channel}N') +
                glob.glob(f'{data_loc}/e{year}/py_formatted/{year}??????????/{year}??????????_?_*.{processing_channel}E') +
                glob.glob(f'{data_loc}/e{year}/py_formatted/{year}??????????/{year}??????????_?_*.{processing_channel}N')

    # print(del_list)
    for file in del_list:
        os.remove(file)


    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()