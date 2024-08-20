import os, glob
import pandas as pd
import numpy as np
from toolkit import Toolkit
from post_processing.process_tdl_files import ProcessTdlFiles
from obspy import UTCDateTime
from geographiclib.geodesic import Geodesic
import obspy.geodetics
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 8
from matplotlib.ticker import (MultipleLocator)


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
class CompareTdlFiles(ProcessTdlFiles):
    """
    Class to handle the comparison of time delay (.tdl)
    outputs by the MTWSPy algorithm when using differing parameters
    May operate in parallel using parameters in param_in.yaml file

    Applies certain filters to the output dataset
    Finds common time delays between datasets
    Plots time delay correlation, histogram of differences and 
    time delays against distance 

    Also saves output if process_tl_files not run previously.
    """
    
    tk = Toolkit()

    def __init__(self):
        super().__init__()

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def normalize_nslc(self, nslc):
        """
        Replace dots with underscores

        :param nslc: network-station string
        :type nslc: str
        :return nslc: modified network-station string
        :type nslc: str
        """
        nslc = nslc.replace("_", ".")
        nslc = nslc.replace("..", ".00.")
        nslc = nslc.replace(".10.", ".00.")
        nslc = nslc.replace(".20.", ".00.")

        nslc = nslc.replace(".LHT", ".LH")
        nslc = nslc.replace(".BHT", ".BH")
        nslc = nslc.replace(".HHT", ".HH")
        nslc = nslc.replace(".LHZ", ".LH")
        nslc = nslc.replace(".BHZ", ".BH")
        nslc = nslc.replace(".HHZ", ".HH")
        nslc = nslc.replace(".LHR", ".LH")
        nslc = nslc.replace(".BHR", ".BH")
        nslc = nslc.replace(".HHR", ".HH")

        return nslc


     # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def find_common_picks(self, df_1, df_2):
        """
        Try to find common arrivals between df_1 and df_2

        :param df_1: tdl dataframe from database 1
        :param df_1: pd.df
        :param df_2: tdl dataframe from database 1
        :param df_2: pd.df
        :return comp_df: tdl dataframe of pick comparison
        :type comp_df: pd.df
        """
        
        # Get unique event ids
        unique_evid_1 = list(df_1['evid'].unique())
        n_evid_1 = len(unique_evid_1)
        
        unique_evid_2 = list(df_2['evid'].unique())
        n_evid_2 = len(unique_evid_2)

        unique_evid = list(set(unique_evid_1) & set(unique_evid_2))
        n_evid = len(unique_evid)


        out_df = []
        # Loop through all unique evids
        if n_evid > 0:
            # Assuming df1 and df2 are your dataframes
            # Apply normalization function to "nslc" column in both dataframes
            df_1["nslc"] = df_1["nslc"].apply(self.normalize_nslc)
            df_2["nslc"] = df_2["nslc"].apply(self.normalize_nslc)

            # Perform merge operation on normalized "nslc" column
            comp_df = pd.merge(df_1, df_2, on=["evid", "phase", "nslc"],
                               suffixes=('_df1', '_df2'))

            # Drop the temporary normalized columns
            columns=["dist_df2", 
                        "az_df2", 
                        "baz_df2", 
                        "mid_lat_df2", 
                        "mid_lon_df2",
                        "date_time_df2", 
                        "evt_lat_df2", 
                        "evt_lon_df2",
                        "evt_dep_df2",
                        "evt_mag_df2",
                        "channel_df2", 
                        "stla_df2", 
                        "stlo_df2",
                        "stel_df2",
                        "ccmx_df1",
                        "ccmx_df2",
                        "tderr_df1",
                        "tderr_df2",
                        "ttaup_df2",
                        "tp_obs_df1",
                        "tp_obs_df2",
                        "tp_syn_df1",
                        "tp_syn_df2",
                        "Ap_obs_df1",
                        "Ap_syn_df1",
                        "Ap_obs_df2",
                        "Ap_syn_df2"]
            for column in columns:
                try: 
                    comp_df.drop(columns = column, inplace=True)
                except Exception as e:
                    print(f"Failed to drop column {column}: {e}")

            # Calculate differences between tdelay_df1 and tdelay_df2
            comp_df["tdelay_diff"] = comp_df["tdelay_df1"] - comp_df["tdelay_df2"]

            comp_df.rename(columns={"date_time_df1": "date_time"}, inplace=True)
            comp_df.rename(columns={"evt_lat_df1": "evt_lat"}, inplace=True)
            comp_df.rename(columns={"evt_lon_df1": "evt_lon"}, inplace=True)
            comp_df.rename(columns={"evt_dep_df1": "evt_dep"}, inplace=True)
            comp_df.rename(columns={"evt_mag_df1": "evt_mag"}, inplace=True)
            comp_df.rename(columns={"channel_df1": "channel"}, inplace=True)
            comp_df.rename(columns={"stla_df1": "stla"}, inplace=True)
            comp_df.rename(columns={"stlo_df1": "stlo"}, inplace=True)
            comp_df.rename(columns={"stel_df1": "stel"}, inplace=True)
            comp_df.rename(columns={"dist_df1": "dist"}, inplace=True)
            comp_df.rename(columns={"az_df1": "az"}, inplace=True)
            comp_df.rename(columns={"baz_df1": "baz"}, inplace=True)
            comp_df.rename(columns={"mid_lat_df1": "mid_lat"}, inplace=True)
            comp_df.rename(columns={"mid_lon_df1": "mid_lon"}, inplace=True)
            comp_df.rename(columns={"ttaup_df1": "ttaup"}, inplace=True)

        else:
            print('No unique evid found')
            print('So return empty dataframe...')
            
            comp_df = pd.DataFrame()
        
        return comp_df


     # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def remove_duplication(self, df):
        """
        The old Matlab code would often duplicate picks
        So lets average and remove them here to stop it influencing
        the comparison

        Causes loss of ccmx, tderr (non existent),
        tp_obs, tp_syn, Ap_obs, Ap_syn

        :param df: dataframe containing duplicate picks
                    (e.g., from two location codes)
        :type df: pd.df
        :return df: dataframe with ducplicates removed
        :type df: pd.df
        """

        group_by_columns = ['evid', 'date_time', 'evt_lat', 'evt_lon', 
                            'evt_dep', 'evt_mag', 'channel', 'nslc', 'stla', 
                            'stlo', 'stel', 'phase', 'ttaup', 'dist', 'baz', 'az', 
                            'mid_lat', 'mid_lon']

        tdl_df = df.groupby(group_by_columns, as_index=False)['tdelay'].mean()

        mp_df = pd.DataFrame(columns = ['ccmx', 'tderr', 'tp_obs', 'tp_syn', 
                                        'Ap_obs','Ap_syn'], 
                                        index = range(0,len(df)))

        mp_df.loc[:, "ccmx"] = 0
        mp_df.loc[:, "tderr"] = 0
        # mp_df.loc[:, "ttaup"] = 0
        mp_df.loc[:, "tp_obs"] = 0
        mp_df.loc[:, "tp_syn"] = 0
        mp_df.loc[:, "Ap_obs"] = 0
        mp_df.loc[:, "Ap_syn"] = 0

        df = pd.merge(tdl_df, mp_df, left_index = True, right_index = True)

        return df


     # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def plot_comparison(self, comp_df, output_directory, filename):
        """
        Try to plot the correlation between tdelay_df1 and tdelay_df2
        And a histogram of the differences
        And a travel time against distance plot

        :param comp_df: pick comparison dataframe to plot & save
        :type comp_df: pd.df
        :param output_directory: directory to save plot
        :type output_directory: str
        :param filename: filename string for plot
        :type filename: str
        """

        if comp_df.empty:
            # Comparison dataframe is empty so we do nothing...
            print('Empty dataframe given')
            print(f'So skipping plotting of {filename}...')

        else:

            # Initialise Figure and add axes with constraints
            fig = plt.figure(figsize  = (14,6))

            ax0 = fig.add_axes([0.05, 0.1, 0.25, 0.8], projection = None, 
                            polar = False,facecolor = 'white',frame_on = True)
            ax1 = fig.add_axes([0.39, 0.1, 0.25, 0.8], projection = None, 
                            polar = False,facecolor = 'white',frame_on = True)
            ax2 = fig.add_axes([0.73, 0.1, 0.25, 0.8], projection = None, 
                            polar = False,facecolor = 'white',frame_on = True)

            axes_list_all = [ax0, ax1, ax2]
            for j, ax in enumerate(axes_list_all):
                ax.set_xticks([])
                ax.set_yticks([])

                ax.spines["right"].set_linewidth(1.5)
                ax.spines["left"].set_linewidth(1.5)
                ax.spines["top"].set_linewidth(1.5)
                ax.spines["bottom"].set_linewidth(1.5)
                ax.tick_params(labelsize = 12)


            # Correlation plot 
            ax0.set_title(f'tdelay df1 v.s. df2',fontsize = 14)
            ax0.set_xlim([-1*(30), 30])
            ax0.set_ylim([-1*(30), 30])
            ax0.xaxis.set_minor_locator(MultipleLocator(5))
            ax0.xaxis.set_major_locator(MultipleLocator(25))
            ax0.yaxis.set_minor_locator(MultipleLocator(5))
            ax0.yaxis.set_major_locator(MultipleLocator(25))

            ax0.set_xlabel('tdelay df1 [s]',fontsize = 14)
            ax0.set_ylabel('tdelay df2 [s]',fontsize = 14)

            ax0.scatter(comp_df["tdelay_df1"], comp_df["tdelay_df2"], 
                        s = 2, color = 'k')

            ax0.plot([-30,30], [-30,30], ls = '-', lw = 0.5, color = 'g')

            # Histogram of differences 
            ax1.hist(comp_df["tdelay_diff"], bins=np.arange(-50, 50+1 , 2), 
                    color='blue', alpha=0.7)  # Adjust bins and color as needed

            # Add labels and title
            ax1.set_xlabel('tdelay Difference [s]',fontsize = 14)
            ax1.set_ylabel('Frequency',fontsize = 14)
            ax1.set_title('Histogram of tdelay Difference',fontsize = 14)
            ax1.set_xlim([-1*(50), 50])
            ax1.set_ylim([0, 5000])
            ax1.xaxis.set_minor_locator(MultipleLocator(5))
            ax1.xaxis.set_major_locator(MultipleLocator(15))
            ax1.yaxis.set_minor_locator(MultipleLocator(250))
            ax1.yaxis.set_major_locator(MultipleLocator(1000))

            # Dist tdelay plot:
            dist_t_df2 = ax2.scatter(comp_df["dist"], 
                                    comp_df["ttaup"] + comp_df["tdelay_df2"], 
                                    s = 2, color = 'b', label = f'df2')
            dist_t_df1 = ax2.scatter(comp_df["dist"], 
                                    comp_df["ttaup"] + comp_df["tdelay_df1"], 
                                    s = 2, color = 'r', label = f'df1') # 

            # Add labels and title
            ax2.set_xlabel('Epicentral Distance [°]',fontsize = 14)
            ax2.set_ylabel('tdelay [s]',fontsize = 14)
            ax2.set_title('tdelay against dist',fontsize = 14)
            ax2.set_xlim([0, 180])
            ax2.set_ylim([0, 3600])
            ax2.xaxis.set_minor_locator(MultipleLocator(15))
            ax2.xaxis.set_major_locator(MultipleLocator(45))
            ax2.yaxis.set_minor_locator(MultipleLocator(250))
            ax2.yaxis.set_major_locator(MultipleLocator(1000))

            # ax2.legend(loc = 'upper right', fontsize = 10)
            plt.sca(ax2)
            plt.legend([dist_t_df1, dist_t_df2], ['df1', 'df2'], 
                    loc = 'lower right', fontsize = 10)

            ax0.annotate('a',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k',
                        backgroundcolor = 'none',ha = 'left', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))
            num_picks = len(comp_df)
            ax0.annotate(f'Picks: {num_picks}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'right', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))
            
            ax1.annotate('b',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'left', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))
            

            mean_tdiff = np.round(np.mean(comp_df["tdelay_diff"]),2)
            ax1.annotate(f'Mean: {mean_tdiff}',(0.5, 1),xytext = (0,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'center', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))            

            std_tdiff = np.round(np.std(comp_df["tdelay_diff"]),2)
            ax1.annotate(f'Std: {std_tdiff}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'right', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))
            

            ax2.annotate('c',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'left', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))
            
            phase_name = filename.split('_')[-1]
            ax2.annotate(f'Phase: {phase_name}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                        fontsize = 12,textcoords = 'offset points', color = 'k', 
                        backgroundcolor = 'none',ha = 'right', va = 'top', 
                        bbox = dict(facecolor = 'white',edgecolor = 'black', 
                        pad = 2.0))

            # Save plot:
            if not os.path.exists(output_directory):
                os.makedirs(output_directory, exist_ok=True)

            filename_out = f'{output_directory}{filename}.png'
            print(f'Sending output to: {str(filename_out)}')

            plt.savefig(filename_out, format = 'png')
            plt.close()

        return


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def load_dataframe(self, params, filter_functions, 
                       output_directory, filename):
        """
        Load csv file from filepath if it exists
        Else use apply to make the dataframe:

        :param params: loaded parameter file 
        :type params: dict
        :param filter_functions: list of functions to filter tdl data
        :type filter_functions: list
        :param output_directory: location of already processed tdl data
        :type output_directory: str
        :param filename: name of tdl file to find/save
        :type filename: str
        :return input_df: loaded dataframe of tdl picks
        :type input_df: pd.df
        """

        if os.path.isfile(f'{output_directory}{filename}.out'):
            # Read:
            input_df = self.read_dataframe(f'{output_directory}{filename}.out')
        else:
            input_df = self.apply_stuff(params, filter_functions)
            self.write_dataframe(output_directory, filename, input_df)

        return input_df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def load_Zaroli_dataframe(self, output_directory, filename):
        """
        Load text file of Zaroli results into dataframe

        :param output_directory: location of already processed tdl data
        :type output_directory: str
        :param filename: name of tdl file to find/save
        :type filename: str
        :return output_df: loaded dataframe of tdl picks
        :type output_df: pd.df
        """
        
        path = f'{output_directory}{filename}'
        file_paths = sorted(glob.glob(path))
        output_df = pd.DataFrame()

        cols = ['year', 'day', 'hour', 'minute', 'second', 'evt_lat', 
                'evt_lon', 'evt_dep', 'stla', 'stlo', 'kntwrk', 
                'kstnm', 'phase', 'period', 'tdelay', 'tderr']

        for file_path in file_paths:
            print(file_path)
            try:
                # Read the file into a DataFrame, skipping lines starting with '#'
                temp = pd.read_csv(file_path, header=None, delim_whitespace=True, comment='#', names=cols)
                res = temp.dropna()

                cols_eq = ['evid', 'date_time', 'nslc', 'dist', 'az', 
                           'baz', 'mid_lat', 'mid_lon', 'evt_mag', 'stel',
                           'ttaup', 'tp_obs', 'tp_syn', 'Ap_obs', 'Ap_syn',
                           'channel', 'ccmx']
                eq_df = pd.DataFrame(columns=cols_eq, index=range(len(res)))

                for index, row in res.iterrows():
                    try:
                        year = int(row["year"])
                        julday = int(row["day"])
                        hour = int(row["hour"])
                        minute = int(row["minute"])
                        sec = float(row["second"])
                        second = int(np.floor(sec))
                        microsecond = int((sec - second) * 1e6)

                        date_time = UTCDateTime(year=year, 
                                                julday=julday, 
                                                hour=hour, 
                                                minute=minute, 
                                                second=second, 
                                                microsecond=microsecond)
                        
                        evid = f'C{date_time.year:04d}{date_time.month:02d}{date_time.day:02d}{date_time.hour:02d}{date_time.minute:02d}A'
                        eq_df.loc[index, "evid"] = evid
                        eq_df.loc[index, "date_time"] = date_time

                    except Exception as e:
                        print(f"Error processing row {index}: {e}")
                        eq_df.loc[index, ["evid", "date_time"]] = np.nan

                    eq_df.loc[index, "evt_mag"] = 6.5
                    eq_df.loc[index, "ccmx"] = 0.99
                    eq_df.loc[index, "channel"] = "T"
                    eq_df.loc[index, "nslc"] = f'{row["kntwrk"]}.{row["kstnm"]}.00.LHT'

                    #,Dummy Columns:

                    eq_df.loc[index, "stel"] = 0
                    eq_df.loc[index, "ttaup"] = 0
                    eq_df.loc[index, "tp_obs"] = 0
                    eq_df.loc[index, "tp_syn"] = 0
                    eq_df.loc[index, "Ap_obs"] = 0
                    eq_df.loc[index, "Ap_syn"] = 0

                    distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(float(row['evt_lat']), 
                                                                           float(row['evt_lon']), 
                                                                           float(row['stla']), 
                                                                           float(row['stlo']))
                    distdg = distm / (6371.e3 * np.pi / 180.)

                    l = Geodesic.WGS84.InverseLine(float(row['evt_lat']), 
                                                   float(row['evt_lon']), 
                                                   float(row['stla']), 
                                                   float(row['stlo']))
                    m = l.Position(0.5 * l.s13)

                    eq_df.loc[index, "dist"] = np.round(distdg, 3)
                    eq_df.loc[index, "az"] = np.round(az, 3)
                    eq_df.loc[index, "baz"] = np.round(baz, 3)
                    eq_df.loc[index, "mid_lat"] = np.round(m['lat2'], 3)
                    eq_df.loc[index, "mid_lon"] = np.round(m['lon2'], 3)

                merged_df = pd.merge(eq_df, res, left_index=True, 
                                     right_index=True)
                df = merged_df.dropna()
                df.drop(columns=["year", "day", "hour", "minute", "second", 
                                 "period", "kntwrk", "kstnm"], inplace=True)

                if not df.empty:
                    output_df = pd.concat([output_df, df], ignore_index=True)

                # Reindex the DataFrame to rearrange columns
                new_column_order = ['evid','date_time','evt_lat','evt_lon',
                                    'evt_dep','evt_mag','channel','nslc',
                                    'stla','stlo', 'stel', 'phase','tdelay',
                                    'tderr','ccmx','ttaup','tp_obs','tp_syn',
                                    'Ap_obs','Ap_syn','dist','az','baz',
                                    'mid_lat','mid_lon']
                output_df = output_df[new_column_order]

            except Exception as e:
                print(f"Failed to process file {file_path}: {e}")
        
        return output_df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def filter_dataframe(self, params, filter_functions, output_df):
        """
        Apply filter function to match other databases

        :param params: loaded parameter file 
        :type params: dict
        :param filter_functions: list of functions to filter tdl data
        :type filter_functions: list
        :param output_df: loaded dataframe of tdl picks
        :type output_df: pd.df       
        :return output_df: filtered dataframe of tdl picks
        :type output_df: pd.df       
        """
        # print(output_df)
        # apply functions
        for k, function in enumerate(filter_functions):
            tdl_df = function(params, output_df)
            output_df = tdl_df
        
        return output_df


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_Zaroli_dataframe(self, params, filter_functions, 
                              output_directory, filename, file_outname):
        """
        Get csv file from filepath if it exists
        Else load & save the dataframe:
        Filter the dataframe according to params

        :param params: loaded parameter file 
        :type params: dict
        :param filter_functions: list of functions to filter tdl data
        :type filter_functions: list
        :param output_directory: location of already processed tdl data
        :type output_directory: str
        :param filename: name of orig tdl file to find/save
        :type filename: str
        :param file_outname: name of saved proc tdl file to load
        :type file_outname: str
        :return input_df: loaded dataframe of tdl picks
        :type input_df: pd.df
        """

        if os.path.isfile(f'{output_directory}{file_outname}.out'):
            # Read:
            input_df = self.read_dataframe(f'{output_directory}{file_outname}.out')
        else:
            input_df = self.load_Zaroli_dataframe(output_directory, 
                                                        filename)
            self.write_dataframe(output_directory, file_outname, input_df)

                # Get station means per phase:

            stat_means_df = self.get_station_means(params, input_df)

            # process_tdl_files.write_station_means(params, stat_means_df)

            filename = f"{params['tt_out_f_name']}_station_means"
            self.write_dataframe(output_directory, filename, stat_means_df)
            

        print(input_df)
        input_filtered_df = self.filter_dataframe(params, filter_functions, input_df)

        return input_filtered_df


    def plot_eq_depths(self, df, output_directory):
        """
        Plot the variability in Earthquake depth
        as a function of the distance and time delay
        to understand the depth phases

        :param df: dataframe to plot & save
        :type df: pd.df
        :param output_directory: directory to save plot
        :type output_directory: str

        """

        # Ensure 'phase' column is of string type
        df['phase'] = df['phase'].astype(str)

        # Create the df_direct dataframe for rows where 'phase' starts with an uppercase 'S'
        df_direct = df[df['phase'].str.startswith('S')]

        # Create the df_depth dataframe for rows where 'phase' starts with a lowercase 's'
        df_depth = df[df['phase'].str.startswith('s')]


        # Initialise Figure and add axes with constraints
        fig = plt.figure(figsize  = (14,6))

        ax0 = fig.add_axes([0.08, 0.1, 0.22, 0.8], projection = None, 
                        polar = False,facecolor = 'white',frame_on = True)
        ax1 = fig.add_axes([0.4, 0.1, 0.22, 0.8], projection = None, 
                        polar = False,facecolor = 'white',frame_on = True)
        ax2 = fig.add_axes([0.73, 0.1, 0.22, 0.8], projection = None, 
                        polar = False,facecolor = 'white',frame_on = True)

        axes_list_all = [ax0, ax1, ax2]
        for j, ax in enumerate(axes_list_all):
            ax.set_xticks([])
            ax.set_yticks([])

            ax.spines["right"].set_linewidth(1.5)
            ax.spines["left"].set_linewidth(1.5)
            ax.spines["top"].set_linewidth(1.5)
            ax.spines["bottom"].set_linewidth(1.5)
            ax.tick_params(labelsize = 12)


        # Dist tdelay plot - all-phases:

        ax0.scatter(df["dist"], 
            df["ttaup"] + df["tdelay"], 
            s = 2, c = df["evt_dep"], cmap = 'viridis')

        ax0.set_xlabel('Epicentral Distance [°]',fontsize = 14)
        ax0.set_ylabel('tdelay [s]',fontsize = 14)
        ax0.set_title('tdelay v.s. dist - All',fontsize = 14)
        ax0.set_xlim([0, 180])
        ax0.set_ylim([0, 3600])
        ax0.xaxis.set_minor_locator(MultipleLocator(15))
        ax0.xaxis.set_major_locator(MultipleLocator(45))
        ax0.yaxis.set_minor_locator(MultipleLocator(250))
        ax0.yaxis.set_major_locator(MultipleLocator(1000))

        # Dist tdelay plot - direct-phases:
        ax1.scatter(df_direct["dist"], 
            df_direct["ttaup"] + df_direct["tdelay"], 
            s = 2, c = df_direct["evt_dep"], cmap = 'viridis')

        # Add labels and title
        ax1.set_xlabel('Epicentral Distance [°]',fontsize = 14)
        ax1.set_ylabel('tdelay [s]',fontsize = 14)
        ax1.set_title('tdelay v.s. dist - direct',fontsize = 14)
        ax1.set_xlim([0, 180])
        ax1.set_ylim([0, 3600])
        ax1.xaxis.set_minor_locator(MultipleLocator(15))
        ax1.xaxis.set_major_locator(MultipleLocator(45))
        ax1.yaxis.set_minor_locator(MultipleLocator(250))
        ax1.yaxis.set_major_locator(MultipleLocator(1000))

        # Dist tdelay plot - depth-phases:
        sc = ax2.scatter(df_depth["dist"], 
            df_depth["ttaup"] + df_depth["tdelay"], 
            s = 2, c = df_depth["evt_dep"], cmap = 'viridis')

        # Add a colorbar to the plot
        cbar = plt.colorbar(sc, ax=ax2)

        # Optionally, label the colorbar
        cbar.set_label('Event Depth')

        # Add labels and title
        ax2.set_xlabel('Epicentral Distance [°]',fontsize = 14)
        ax2.set_ylabel('tdelay [s]',fontsize = 14)
        ax2.set_title('tdelay v.s. dist - depth-phases',fontsize = 14)
        ax2.set_xlim([0, 180])
        ax2.set_ylim([0, 3600])
        ax2.xaxis.set_minor_locator(MultipleLocator(15))
        ax2.xaxis.set_major_locator(MultipleLocator(45))
        ax2.yaxis.set_minor_locator(MultipleLocator(250))
        ax2.yaxis.set_major_locator(MultipleLocator(1000))

        ax0.annotate('a',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k',
                    backgroundcolor = 'none',ha = 'left', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        num_picks = len(df)
        ax0.annotate(f'Picks: {num_picks}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k', 
                    backgroundcolor = 'none',ha = 'right', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        ax1.annotate('b',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k', 
                    backgroundcolor = 'none',ha = 'left', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        num_picks = len(df_direct)
        ax1.annotate(f'Picks: {num_picks}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k', 
                    backgroundcolor = 'none',ha = 'right', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        ax2.annotate('c',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k', 
                    backgroundcolor = 'none',ha = 'left', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        num_picks = len(df_depth)
        ax2.annotate(f'Picks: {num_picks}',(1, 1),xytext = (-5,-5),xycoords = 'axes fraction',
                    fontsize = 12,textcoords = 'offset points', color = 'k', 
                    backgroundcolor = 'none',ha = 'right', va = 'top', 
                    bbox = dict(facecolor = 'white',edgecolor = 'black', 
                    pad = 2.0))

        # Save plot:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)

        filename_out = f'{output_directory}All_direct_depth_phases.png'
        print(f'Sending output to: {str(filename_out)}')

        plt.savefig(filename_out, format = 'png')
        plt.close()

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
def main():
    # Create instance of toolkit as tk
    tk = Toolkit()
    params = tk.get_params('params_in.yaml')

    # Instantiate the classes and call their methods as needed
    compare_tdl_files = CompareTdlFiles()

    filter_functions = [compare_tdl_files.filt_date_time_max,
                        compare_tdl_files.filt_date_time_min,
                        compare_tdl_files.filt_networks,
                        compare_tdl_files.filt_components,
                        compare_tdl_files.filt_evt_lat_max,
                        compare_tdl_files.filt_evt_lat_min,
                        compare_tdl_files.filt_evt_lon_max,
                        compare_tdl_files.filt_evt_lon_min,
                        compare_tdl_files.filt_evt_dep_max,
                        compare_tdl_files.filt_evt_dep_min,
                        compare_tdl_files.filt_evt_mag_max,
                        compare_tdl_files.filt_evt_mag_min,
                        compare_tdl_files.filt_sta_lat_max,
                        compare_tdl_files.filt_sta_lat_min,
                        compare_tdl_files.filt_sta_lon_max,
                        compare_tdl_files.filt_sta_lon_min,
                        compare_tdl_files.filt_phases,
                        compare_tdl_files.filt_tdl_max,
                        compare_tdl_files.filt_ccmx_min,
                        compare_tdl_files.filt_tderr_max,
                        compare_tdl_files.filt_dist_max,
                        compare_tdl_files.filt_dist_min]    


    filename = f"{params['tt_out_f_name']}"

    ######################### AL JAN 2008 data - XC - T comp #########################

    output_directory_example = f"{params['home']}/MTWSPy/data/output_df/"

    al_python_filtered_XC_df = compare_tdl_files.load_dataframe(params, 
                                                                filter_functions, 
                                                                output_directory_example, 
                                                                "200801_IU_T_df")
    print(al_python_filtered_XC_df)



    #########################        USER PICKS             #########################

    output_directory = f"{params['home']}/{params['proc_tdl_loc']}/"

    user_python_filtered_XC_df = compare_tdl_files.load_dataframe(params, 
                                                                filter_functions, 
                                                                output_directory, 
                                                                filename)
    print(user_python_filtered_XC_df)


    #########################        PROCESSING            #########################

    # Check differences in means of S phase
    al_python_filtered_XC_df_S = al_python_filtered_XC_df.query("phase == 'S'")
    user_python_filtered_XC_df_S = user_python_filtered_XC_df.query("phase == 'S'")

    print(f'Mean S TCOMP AL: {np.mean(al_python_filtered_XC_df_S["tdelay"])}')
    print(f'Mean S TCOMP USER: {np.mean(user_python_filtered_XC_df_S["tdelay"])}')

    # Plot Earthquake depth analysis 
    compare_tdl_files.plot_eq_depths(al_python_filtered_XC_df, output_directory_example)
    compare_tdl_files.plot_eq_depths(user_python_filtered_XC_df, output_directory)

    # Find common picks between example data set (al_python_filtered_XC_df)
    # and that produced by the user (user_python_filtered_XC_df)
    # This will only work for Jan 2008 

    comp_al_user_py_XC_TCOMP_df_all = compare_tdl_files.find_common_picks(al_python_filtered_XC_df, user_python_filtered_XC_df)

    # Analyse differences between phases
    compare_tdl_files.plot_comparison(comp_al_user_py_XC_TCOMP_df_all, output_directory, 'comp_al_user_py_XC_TCOMP_df_all')

    # Analyse just S component
    comp_al_user_py_XC_TCOMP_df_all_S = comp_al_user_py_XC_TCOMP_df_all.query("phase == 'S'")
    compare_tdl_files.plot_comparison(comp_al_user_py_XC_TCOMP_df_all_S, output_directory, 'comp_al_user_py_XC_TCOMP_df_all_S')





if __name__ == '__main__':
    main()
