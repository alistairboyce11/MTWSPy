import pandas as pd
import numpy as np
from toolkit import Toolkit
from post_processing.process_tdl_files import ProcessTdlFiles
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 8
from matplotlib.ticker import (MultipleLocator)

class CompareTdlFiles(ProcessTdlFiles):
    tk = Toolkit()

    def __init__(self):
        super().__init__()

    def specific_function(self):
        # Function specific to CompareTdlFiles.py
        pass

    def normalize_nslc(self, nslc):
        # Replace dots with underscores
        nslc = nslc.replace("_", ".")
        nslc = nslc.replace("..", ".00.")
        nslc = nslc.replace(".10.", ".00.")
        nslc = nslc.replace(".20.", ".00.")

        return nslc


    def find_common_picks(self, df_1, df_2):
        '''
        Try to find common arrivals between df_1 and df_2
        '''
        
        # Get unique event ids
        unique_evid_1 = list(df_1['evid'].unique())
        n_evid_1 = len(unique_evid_1)
        
        unique_evid_2 = list(df_2['evid'].unique())
        n_evid_2 = len(unique_evid_2)

        unique_evid = list(set(unique_evid_1 + unique_evid_2))
        n_evid = len(unique_evid)


        
        out_df = []
        # Loop through all unique evids
        if n_evid > 0:
            # Assuming df1 and df2 are your dataframes
            # Apply normalization function to "nslc" column in both dataframes
            df_1["nslc"] = df_1["nslc"].apply(self.normalize_nslc)
            df_2["nslc"] = df_2["nslc"].apply(self.normalize_nslc)

            # Perform merge operation on normalized "nslc" column
            comp_df = pd.merge(df_1, df_2, on=["evid", "phase", "nslc"], suffixes=('_df1', '_df2'))

            # Drop the temporary normalized column
            comp_df.drop(columns=["dist_df2", 
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
                                "Ap_syn_df2"], inplace=True)

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
            comp_df = pd.DataFrame()
        
        return comp_df


    def plot_comparison(self, comp_df):
        '''
        Try to plot the correlation between tdelay_df1 and tdelay_df2
        And a histogram of the differences
        '''


        # Initialise Figure and add axes with constraints
        fig = plt.figure(figsize  = (14,6))

        ax0 = fig.add_axes([0.05, 0.1, 0.25, 0.8], projection = None, polar = False,facecolor = 'white',frame_on = True)
        ax1 = fig.add_axes([0.39, 0.1, 0.25, 0.8], projection = None, polar = False,facecolor = 'white',frame_on = True)
        ax2 = fig.add_axes([0.73, 0.1, 0.25, 0.8], projection = None, polar = False,facecolor = 'white',frame_on = True)

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
        ax0.set_xlim([-1*(60), 60])
        ax0.set_ylim([-1*(60), 60])
        ax0.xaxis.set_minor_locator(MultipleLocator(5))
        ax0.xaxis.set_major_locator(MultipleLocator(25))
        ax0.yaxis.set_minor_locator(MultipleLocator(5))
        ax0.yaxis.set_major_locator(MultipleLocator(25))

        ax0.set_xlabel('tdelay df1 [s]',fontsize = 14)
        ax0.set_ylabel('tdelay df2 [s]',fontsize = 14)

        ax0.scatter(comp_df["tdelay_df1"], comp_df["tdelay_df2"], s = 2, color = 'k') # , label = f'{params_in['component']} waveform'

        ax0.plot([-60,60], [-60,60], ls = '-', lw = 0.5, color = 'g')

        # Histogram of differences 
        ax1.hist(comp_df["tdelay_diff"], bins=np.arange(-50, 50+1 , 2), color='blue', alpha=0.7)  # Adjust bins and color as needed

        # Add labels and title
        ax1.set_xlabel('tdelay Difference',fontsize = 14)
        ax1.set_ylabel('Frequency',fontsize = 14)
        ax1.set_title('Histogram of tdelay Difference',fontsize = 14)
        ax1.set_xlim([-1*(50), 50])
        ax1.set_ylim([0, 5000])
        ax1.xaxis.set_minor_locator(MultipleLocator(5))
        ax1.xaxis.set_major_locator(MultipleLocator(15))
        ax1.yaxis.set_minor_locator(MultipleLocator(250))
        ax1.yaxis.set_major_locator(MultipleLocator(1000))

        # Dist tdelay plot:

        dist_t_df2 = ax2.scatter(comp_df["dist"], comp_df["ttaup"] + comp_df["tdelay_df2"], s = 2, color = 'b', label = f'df2') # , label = f'{params_in['component']} waveform'
        dist_t_df1 = ax2.scatter(comp_df["dist"], comp_df["ttaup"] + comp_df["tdelay_df1"], s = 2, color = 'r', label = f'df1') # 

        # Add labels and title
        ax2.set_xlabel('Epicentral Distance',fontsize = 14)
        ax2.set_ylabel('tdelay',fontsize = 14)
        ax2.set_title('tdelay against dist',fontsize = 14)
        ax2.set_xlim([0, 180])
        ax2.set_ylim([0, 3600])
        ax2.xaxis.set_minor_locator(MultipleLocator(15))
        ax2.xaxis.set_major_locator(MultipleLocator(45))
        ax2.yaxis.set_minor_locator(MultipleLocator(250))
        ax2.yaxis.set_major_locator(MultipleLocator(1000))



        # ax2.legend(loc = 'upper right', fontsize = 10)
        plt.sca(ax2)
        plt.legend([dist_t_df1, dist_t_df2], ['df1', 'df2'], loc = 'lower right', fontsize = 10)


        ax0.annotate('a',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',fontsize = 12,textcoords = 'offset points', color = 'k', backgroundcolor = 'none',ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))
        ax1.annotate('b',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',fontsize = 12,textcoords = 'offset points', color = 'k', backgroundcolor = 'none',ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))
        ax2.annotate('c',(0, 1),xytext = (5,-5),xycoords = 'axes fraction',fontsize = 12,textcoords = 'offset points', color = 'k', backgroundcolor = 'none',ha = 'left', va = 'top', bbox = dict(facecolor = 'white',edgecolor = 'black', pad = 2.0))

        # plt.savefig(f'{pic_loc}/{pic_filename}', format = 'pdf')
        plt.savefig(f'compare_codes.pdf', format = 'pdf')

        return


def main():
    # matlab_code/tdelay/OST-MXT/
    # Lei_Results/tdelay/OST-MXT/
    # python_code/tdelay/OST-MXT/

    # Create instance of toolkit as tk
    tk = Toolkit()
    params_in = tk.get_params('params_in.yaml')

    # Instantiate the classes and call their methods as needed
    compare_tdl_files = CompareTdlFiles()

    params_in['home'] = '/Users/alistair/Lyon_Pdoc/Lei_Li_codes_data/python_code' # Location where MTWSPy code is installed

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

    al_python_filtered_df = compare_tdl_files.apply_stuff(params_in, filter_functions)
    
    print(al_python_filtered_df)

    # Lei original results
    params_in['home'] = '/Users/alistair/Lyon_Pdoc/Lei_Li_codes_data/Lei_Results' # Location where MTWSPy code is installed

    lei_matlab_filtered_df = compare_tdl_files.apply_stuff(params_in, filter_functions)
    
    print(lei_matlab_filtered_df)

    # # # Lei original code, al results
    params_in['home'] = '/Users/alistair/Lyon_Pdoc/Lei_Li_codes_data/matlab_code' # Location where MTWSPy code is installed

    al_matlab_filtered_df = compare_tdl_files.apply_stuff(params_in, filter_functions)
    
    print(al_matlab_filtered_df)


    comp_df = compare_tdl_files.find_common_picks(al_python_filtered_df, lei_matlab_filtered_df)
    comp_mat_df = compare_tdl_files.find_common_picks(al_matlab_filtered_df, lei_matlab_filtered_df)

    print(comp_df)

    compare_tdl_files.plot_comparison(comp_df)


if __name__ == '__main__':
    main()