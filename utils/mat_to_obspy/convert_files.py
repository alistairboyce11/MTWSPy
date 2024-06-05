
import sys,glob,time,os
import pandas as pd
import numpy as np
import scipy
from obspy import UTCDateTime, read, read_inventory, Stream, Trace
from datetime import datetime, timedelta
import obspy.signal
import obspy.signal.rotate
import obspy.geodetics.base
import obspy.geodetics
import concurrent.futures

from toolkit import Toolkit
import pprint

class ConvertMatFiles:
    """
    Class to handle conversion of old Matlab data files to python format
    for use with this code.
    """
    def __init__(self):
        pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def extract_values(self, item):
        """
        Helper function to extract values from the numpy object array.
        
        :param item: some item that has undetermined type
        :type item: unknown

        :return item: output with determined type
        :type item: variable
        """
        if isinstance(item, np.ndarray):
            if item.dtype.names:
                return {name: self.extract_values(item[name]) for name in item.dtype.names}
            else:
                if item.size == 1:
                    return item.item()  # Return a scalar
                else:
                    return item.tolist()  # Return a list
        return item

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def datenum_to_datetime(self, datenum):
        """
        Convert MATLAB datenum to Python datetime
        
        :param datenum: Matlab datenum stamp
        :type datenum: Matlab datenum format string

        :return output: translated time stamp for python
        :type output: datetime object
        """

        # MATLAB datenum starts from the year 0000-01-01
        # Python's datetime starts from 0001-01-01
        # So, we need to adjust by adding the difference in days between these two start points
        MATLAB_DATENUM_OFFSET = 367  # Days from 0000-01-01 to 0001-01-01

        days = datenum - MATLAB_DATENUM_OFFSET
        return datetime(1, 1, 1) + timedelta(days=days)
    

    def read_mat_file(self, input_dict, seis, data, fail):
        """
        Read matlab .mat file using scipy into python dict
 
        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param seis: seismic data in python format to be populated
        :type seis: obspy stream
        :param data: matlab data in scipy dict format to populate obspy stream
        :type seis: dict
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail != 1:
            file = input_dict['filename']
            try: 
                data = scipy.io.loadmat(file)
            except:
                fail = 1
                print(f'Failed to read {file}')
        else:
            print(f'Skipped: read_mat_file')

        return input_dict, seis, data, fail
    

    def load_headers(self, input_dict, seis, data, fail):
        """
        Interpret scipy dict object to populate headers in obspy stream

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param seis: seismic data in python format to be populated
        :type seis: obspy stream
        :param data: matlab data in scipy dict format to populate obspy stream
        :type seis: dict
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """
        if fail != 1:
            file = input_dict['filename']
            try: 

                # Create a new dictionary with readable format
                readable_dict = {}

                # Convert 'station' data
                station_data = data['station']
                station_list = []

                for station in station_data:
                    station_dict = {}
                    for field in station.dtype.names:
                        station_dict[field] = self.extract_values(station[field][0])
                    station_list.append(station_dict)

                readable_dict['station'] = station_list

                # Convert 'event' data
                event_data = data['event']
                event_list = []

                for event in event_data:
                    event_dict = {}
                    for field in event.dtype.names:
                        event_dict[field] = self.extract_values(event[field][0])
                    event_list.append(event_dict)

                readable_dict['event'] = event_list

                # # Add other fields
                converted_times = [self.datenum_to_datetime(t[0]) for t in data['time'].tolist()]
                readable_dict['time'] = [UTCDateTime(t) for t in converted_times]
                dt = round(readable_dict['time'][1]-readable_dict['time'][0],4)


                readable_dict['ZRT'] = data['ZRT']

                # # Add header, version, and globals
                # readable_dict['__header__'] = data['__header__']
                # readable_dict['__version__'] = data['__version__']
                # readable_dict['__globals__'] = data['__globals__']

                # Print the readable dictionary
                # pprint.pprint(readable_dict)

                sta = readable_dict['station'][0]['nslc'].split('.')[1]
                network = readable_dict['station'][0]['nslc'].split('.')[0]
                location = readable_dict['station'][0]['nslc'].split('.')[2]

                evid = readable_dict['event'][0]['centroid']['event_name'][0]
                time_stamp = readable_dict['event'][0]['centroid']['time'][0]
                year = int(time_stamp[0])
                month = int(time_stamp[1])
                day  = int(time_stamp[2])
                hour = int(time_stamp[3])
                minute = int(time_stamp[4])
                second = int(time_stamp[5])

                evt_string = f'{year:04d}{month:02d}{day:02d}{hour:02d}{minute:02d}{second:02d}'

                input_dict['outfile_name'] = f'{evt_string}_{network}_{sta}.{location}.'

                for c, component in enumerate(["Z", "R", "T"]):
                    channel = "LH" + component
                    tr = Trace(header = {'station': sta, 
                                       'network': network,
                                       'location': location,
                                       'channel': channel,
                                       'delta': dt,
                                       'npts': len(readable_dict['time']),
                                       'starttime': readable_dict['time'][0]},
                                data = readable_dict['ZRT'][:,c])
                    
                    tr.stats['sac']={}
                    tr.stats['sac']['stla'] = float(readable_dict['station'][0]['stla'])
                    tr.stats['sac']['stlo'] = float(readable_dict['station'][0]['stlo'])
                    tr.stats['sac']['stel'] = float(readable_dict['station'][0]['stel']) * 1000 # For Meters
                    tr.stats['sac']['stdp'] = 5.0
                    tr.stats['sac']['gcarc'] = float(readable_dict['station'][0]['dis'])
                    tr.stats['sac']['baz'] = float(readable_dict['station'][0]['baz'])

                    tr.stats['sac']['kstnm'] = sta
                    tr.stats['sac']['kevnm'] = evid
                    tr.stats['sac']['khole'] = location
                    tr.stats['sac']['kcmpnm'] = channel
                    tr.stats['sac']['knetwk'] = network

                    tr.stats['sac']['evla'] = float(readable_dict['event'][0]['centroid']['latitude'][0][0])
                    tr.stats['sac']['evlo'] = float(readable_dict['event'][0]['centroid']['longitude'][0][0])
                    tr.stats['sac']['evdp'] = float(readable_dict['event'][0]['centroid']['depth'][0][0])
                    tr.stats['sac']['o'] = UTCDateTime(year, month, day, hour, minute, time_stamp[5]) - UTCDateTime(tr.stats['starttime'])
                    print(tr.stats)

                    # The tangential seems to be flipped in polarity so multiply by -1.
                    if component == 'T':
                        tr.data = tr.data * -1
                    

                    seis.append(tr) 



            except:
                fail = 1
                print(f'Failed to load headers {file}')
        else:
            print(f'Skipped: load_headers')

        return input_dict, seis, data, fail

    
    def save_sac(self, input_dict, seis, data, fail):
        """
        Save newly populated obspy stream object to python compatible file

        Inputs and returns same/updated parameters
        :param input_dict : details for processing
        :type input_dict: dict
        :param seis: seismic data in python format to be populated
        :type seis: obspy stream
        :param data: matlab data in scipy dict format to populate obspy stream
        :type seis: dict
        :param fail: flag for whether processing step failed
        :type fail: 1/0 or bool
        """

        if fail != 1:
            file = input_dict['filename']
            try: 
                Dummy = 1

                for tr in seis:
                    channel = tr.stats['channel']
                    outfile = f'{input_dict['output_directory']}/{input_dict['outfile_name']}{channel}'

                    print(f'Saving {outfile}')
                    tr.write(outfile, format='sac')

            except:
                fail = 1
                print(f'Failed to save sac {file}')
        else:
            print(f'Skipped: save_sac')

        return input_dict, seis, data, fail
    
    def process_one_file(self, input_dict):
        """
        Processes one matlab file and applying the list functions to them

        :return input_dict: dictionary of info used in processing
        :type input_dict: dict
        """
        functions = input_dict['functions']

        seis = Stream()
        data = {}
        fail = 0
        # apply functions, only executed when fail = 0
        
        for function in functions:
            input_dict, seis, data, fail = function(input_dict, seis, data, fail)
        return
    
    
    def execute(self, main_function, input_directory, output_directory, 
                functions, do_parallel, jobs):
        """
        Takes main function inputs and passes to serial or parallel 
        executer, see below.

        :param main_function : list of the 1 main function applied to inputs
        :type main_function: list
        :param input_directory: input_directory name where the data is
        :type input_directory: str
        :param output_directory: output_directory name where the data is saved
        :type output_directory: str
        :param functions : list of functions to apply to inputs
        :type functions: list
        :param do_parallel: Execute in parallel or not
        :type do_parallel : bool
        :param jobs: number of parallel workers
        :type jobs: int
        """

        input_dicts = []

        files = sorted(glob.glob(f'{input_directory}/IU.FURI.00.ZRT.mat'))

        if len(files) == 0:
            # No events found...
            print('----------')
            print(f'----------////   NO-FILES-IN: {input_directory}    ////----------')
        else:

            print('----------')
            print('Processing the following files list:   ')
            print('----------')

        for filename in files:

            input_dict = {}
            input_dict['input_directory'] = input_directory
            input_dict['output_directory'] = output_directory
            input_dict['functions'] = functions
            input_dict['filename'] = filename

            input_dicts.append(input_dict)

        # Send input dicts to main_function for processing.
        if input_dicts:

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
    

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():

    start_time = time.time()

    # Create instance of toolkit as tk
    tk = Toolkit()
    params_in = tk.get_params('params_in.yaml')

    # Instantiate the classes and call their methods as needed
    convert_mat_files = ConvertMatFiles()

    input_directory = '/Volumes/INTENSO3/leli_home/d_data_and_docs/lyon/dmt/e2008/20081230_194956.a/obs'
    output_directory= '/Users/alistair/Lyon_Pdoc/Lei_Li_codes_data/d_data_obs/lyon/dmt/e2008/py_formatted/20081230194956'
    do_parallel=0
    jobs=4

    functions = [convert_mat_files.read_mat_file, convert_mat_files.load_headers, convert_mat_files.save_sac]

    main_function = [convert_mat_files.process_one_file]
    
    convert_mat_files.execute(main_function, input_directory, output_directory, functions, do_parallel, jobs)
    

    print("--- %s seconds ---" % (time.time() - start_time))



if __name__ == '__main__':
    main()