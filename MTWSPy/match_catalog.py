import pandas as pd
import os
from obspy import UTCDateTime
import numpy as np
import inspect

class MatchCatalog:
    """
    Class to handle matching observed and synthetic data earthquake catalogs
    """

    def __init__(self, params):
        """
        On initiation of class load params and write header of logfile
        """

        self.params = params
        self.logfile = self.open_log_file()
        self.write_params_logfile()


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def open_log_file(self):
        """
        Open a logfile for the match_catalog process
        Use params to define name and filepath
        Return open logfile
        """

        lf_loc = f"{self.params['home']}/{self.params['log_loc']}/{str(self.params['code_start_time'])}/{os.path.basename(__file__).split('.')[0]}"
        if not os.path.exists(lf_loc):
            os.makedirs(lf_loc, exist_ok=True)
        lf_name = f"{lf_loc}/match_catalog.log"
        return open(lf_name, "w")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def write_params_logfile(self):
        """
        Write header to match_catalog logfile
        """
        
        justify = 30
        self.logfile.write(" ")
        self.logfile.write("----------////               INPUT PARAMETERS                ////----------\n")
        params_list = ['cmt_infile', 'cmt_outfile', 'year', 'match_error_time', 'match_error_lat', 'match_error_lon', 'match_error_dep', 'match_error_mag']
        for k, param in enumerate(params_list):
            self.logfile.write("{0:>{x}s} {1:s} {2:s}\n".format(param, " : ", str(self.params[param]), x=justify))
        self.logfile.write("----------////               INPUT PARAMETERS                ////----------\n")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_cmt_catalog(self):
        """
        Read CMT catalog file into pandas dataframe
        It has the following columns:
        
        id,evid,evtm,evla,evlo,evdp,mag,tshift,hdur,m_rr,m_tt,m_pp,m_rt,m_rp,m_tp,
            strike_1,dip_1,rake_1,strike_2,dip_2,rake_2,azimuth_t,plunge_t,length_t,
                azimuth_p,plunge_p,length_p,azimuth_n,plunge_n,length_n

        Add the following formatted columns and return dataframe

        evtm: Centroid time in UTCDateTime format
        evhtm: Hypocentral time in UTCDateTime format
        id_ctm: centroid id
        id_htm: Hypocentral id
        id_fmt_ctm: format string for centroid time
        id_fmt_htm: format string for hypocentral time
        id_dmt_fmt_ctm: format string for dmt centroid time
        id_dmt_fmt_htm: format string for dmt hypocentral time

        :param self.params: dict inc loc of cmt_infile
        :type self.params: dict

        :return cat: dataframe of cmt catalog with added columns as above
        :type cat: pd.df
        """

        from toolkit import Toolkit
        tk = Toolkit()

        log_statement = f'* {os.path.basename(__file__).split('.')[0]}, {str(inspect.stack()[0][3])}'
        
        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Reading CMT catalog')
        ###
        
        cmt_fname = f'{self.params['home']}/{self.params['cmt_infile']}'
        cat = pd.read_csv(cmt_fname)

        ############## This recreates id from other field.
        evtm = []; evhtm = []; id_ctm = []; id_htm = [] 
        id_fmt_ctm = []; id_fmt_htm = []; id_dmt_fmt_ctm = []; id_dmt_fmt_htm = []

        # Iterate over rows in cat, add formatted strings to lists.
        for ind, row in cat.iterrows():

            evtm.append(UTCDateTime(row['evtm']))
            evhtm.append(evtm[ind]-row['tshift'])

            year = "{0:4d}".format(int(evtm[ind].year))
            month = "{0:0>2d}".format(int(evtm[ind].month))
            day = "{0:0>2d}".format(int(evtm[ind].day))
            days = year+month+day

            hour = "{0:0>2d}".format(int(evtm[ind].hour))
            minute = "{0:0>2d}".format(int(evtm[ind].minute))
            second = "{0:0>2d}".format(int(evtm[ind].second))
            times = hour+minute+second

            mag = "{0:0>3d}".format(int(np.round(row['mag']*100,0)))
            dep = "{0:0>3d}".format(int(np.round(row['evdp'],0)))

            id_ctm.append("{0:8s}T{1:6s}M{2:3s}Z{3:3s}".format(days, times, mag, dep))
            id_fmt_ctm.append("{0:8s}{1:6s}".format(days, times))
            id_dmt_fmt_ctm.append("{0:8s}_{1:6s}.a".format(days, times))

            year1 = "{0:4d}".format(int(evhtm[ind].year))
            month1 = "{0:0>2d}".format(int(evhtm[ind].month))
            day1 = "{0:0>2d}".format(int(evhtm[ind].day))
            days1 = year1+month1+day1

            hour1 = "{0:0>2d}".format(int(evhtm[ind].hour))
            minute1 = "{0:0>2d}".format(int(evhtm[ind].minute))
            second1 = "{0:0>2d}".format(int(evhtm[ind].second))
            times1 = hour1+minute1+second1

            id_htm.append("{0:8s}T{1:6s}M{2:3s}Z{3:3s}".format(days1, times1, mag, dep))
            id_fmt_htm.append("{0:8s}{1:6s}".format(days1, times1))
            id_dmt_fmt_htm.append("{0:8s}_{1:6s}.a".format(days1, times1))

        # Add lists to dataframe
        cat['evtm'] = evtm
        cat['evhtm'] = evhtm
        cat['id_ctm'] = id_ctm
        cat['id_htm'] = id_htm
        cat['id_fmt_ctm'] = id_fmt_ctm
        cat['id_fmt_htm'] = id_fmt_htm
        cat['id_dmt_fmt_ctm'] = id_dmt_fmt_ctm
        cat['id_dmt_fmt_htm'] = id_dmt_fmt_htm

        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  CMT catalog returned...')
        ###

        return cat


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_data_catalog_year(self):
        """
        Read obspy dmt catalog_table.txt file or
        IRIS FetchData ????_events.out file into pandas dataframe
        With the following columns:

        'N','LAT', 'LON', 'DEP', 'DATETIME', 'MAG', 'EV_ID'

        Return dataframe

        :param self.params: dict inc loc of data catalog
        :type self.params: dict

        :return data_cat: dataframe of data catalog with columns as above
        :type data_cat: pd.df
        """
        from toolkit import Toolkit
        tk = Toolkit()

        log_statement = f'* {os.path.basename(__file__).split('.')[0]}, {str(inspect.stack()[0][3])}'
        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Get sorted DMT catalog...')
        ###

        d_year = self.params['year']
        
        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Requesting DMT catalog year {d_year:d}...')
        ###

        f_loc_dmt = f'{self.params['obs_loc']}/e{str(d_year)}/EVENTS-INFO/catalog_table.txt'
        f_loc_iris= f'{self.params['obs_loc']}/e{str(d_year)}/{str(d_year)}_events.out'

        if os.path.isfile(f_loc_dmt): 
            # Data downloaded from IRIS dmt (old version) 
            
            ###
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Reading DMT file : {f_loc_dmt:s}')
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  For year : {d_year:d}')
            ###

            file = open(f_loc_dmt,'r')
            lines = file.readlines()
            file.close()

            out_cat = []
            N = []; LAT = []; LON = []; DEP = []
            DATETIME = []; MAG = []; EV_ID = []

            for line in lines:
                if 'Command' in line or 'pyDMT' in line or \
                    'DATETIME' in line or '--' in line or len(line) <= 1:
                    # Headers
                    continue
                else:
                    N.append(int(line.split()[0]))
                    LAT.append(float(line.split()[1]))
                    LON.append(float(line.split()[2]))
                    DEP.append(float(line.split()[3]))
                    DATETIME.append(str(UTCDateTime(str(line.split()[4]))))
                    MAG.append(float(line.split()[5]))
                    EV_ID.append(str(line.split()[7]))

                    out_cat.append([int(line.split()[0]), float(line.split()[1]), 
                                    float(line.split()[2]), float(line.split()[3]), 
                                    str(UTCDateTime(str(line.split()[4]))), 
                                    float(line.split()[5]),str(line.split()[7])])

            column_names = ['N','LAT', 'LON', 'DEP', 'DATETIME', 'MAG', 'EV_ID']
            data_cat = pd.DataFrame(sorted(out_cat), columns = column_names)
            
            # Drop duplicated rows
            data_cat = data_cat.drop_duplicates() 

            ###
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  DMT catalog returned for year : {d_year:d}...')
            ###


        elif os.path.isfile(f_loc_iris):
            # Data downloaded using iris FetchData
            
            ###
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Reading IRIS event file : {f_loc_iris:s}')
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  For year : {d_year:d}')
            ###


            try:
                # Read cmt table, sort and re-index.
                nms = ['N', 'DATETIME', 'LAT', 'LON', 'DEP', 'CAT', 'TEMP_CAT', 'CMT_ID', 'MAG', 'LOC']
                data_cat = pd.read_csv(f'{f_loc_iris}', delimiter="|", header=None, names = nms)
                data_cat = data_cat.sort_values(by='DATETIME').reset_index(drop=True)

                # Add ev_ids
                ev_df = pd.DataFrame(columns = ['EV_ID'], index = range(0,len(data_cat)))

                for index, row in data_cat.iterrows(): 
                    eq_time = UTCDateTime(row['DATETIME'])
                    ev_df.loc[index, 'EV_ID'] = f'{eq_time.year:4d}{eq_time.month:02d}{eq_time.day:02d}{eq_time.hour:02d}{eq_time.minute:02d}{eq_time.second:02d}'

                data_cat = pd.merge(data_cat, ev_df, left_index = True, right_index = True)    

                # Remove the "MW,"
                # Function to extract float from the string
                def extract_float(entry):
                    return float(entry.split(',')[1])

                data_cat['MAG'] = data_cat['MAG'].apply(lambda x: extract_float(x))

                ###
                tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  IRIS data catalog returned for year : {d_year:d}...')
                ###

            except:
                data_cat = pd.DataFrame()


        else:
            # Return empty dataframe

            ###
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  No file found at : {f_loc_dmt:s} or {f_loc_iris:s}')
            tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Skipping...')
            ###

            data_cat = pd.DataFrame()

        return data_cat


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def get_match_catalog(self):
        """
        For every line in the data catalog, search for matched event 
        in cmt catalog. Using error values in params (time, lat, lon, 
        depth, magnitude). If match found, output to dataframe

        :param self.data_cat: data catalog dataframe
        :type self.data_cat: pd.df
        :param self.cmt_cat: cmt catalog dataframe
        :type self.cmt_cat: pd.df
        :param match_error_lat: acceptable for lat difference from params file
        :type match_error_lat: float
        :param match_error_lon: acceptable for lon difference from params file
        :type match_error_lon: float
        :param match_error_dep: acceptable for dep difference from params file
        :type match_error_dep: float
        :param match_error_mag: acceptable for mag difference from params file
        :type match_error_mag: float

        :return df_out: matched dataframe with all necessary details
        :type df_out: pd.df
        """
        from toolkit import Toolkit
        tk = Toolkit()

        log_statement = f'* {os.path.basename(__file__).split('.')[0]}, {str(inspect.stack()[0][3])}'
        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Attempting CMT-DMT catalog match...')
        ###

        out_cat = []

        for index, row in self.data_cat.iterrows():

            ind = np.where(np.array(np.abs(row['LAT']-self.cmt_cat.evla) <= self.params['match_error_lat']) &
                        np.array(np.abs(row['LON']-self.cmt_cat.evlo) <= self.params['match_error_lon']) &
                        np.array(np.abs(row['DEP']-self.cmt_cat.evdp) <= self.params['match_error_dep']) &
                        np.array(np.abs(row['MAG']-self.cmt_cat.mag)  <= self.params['match_error_mag']))[0]
            if len(ind) < 1:
                
                ###
                tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  **No Match #1 for data_cat: {index:d}')
                ###

            else:
                # Now check the times.
                ind2 = []
                for m in ind:
                    if np.abs(UTCDateTime(row['DATETIME'])-self.cmt_cat['evtm'][m]) <= self.params['match_error_time']:
                        ind2.append(m)

                if len(ind2) < 1:
                    # Shouldnt really get here.

                    ###
                    tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  **No Match #2 for data_cat: {index:d}')
                    ###

                    continue

                elif len(ind2) > 1:
                    # Take first entry:
                    ind2 = ind2[0]
                    
                    ###
                    tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Taking first entry... {index:d}')
                    tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  {str(self.cmt_cat.loc[ind2]):s}')
                    ###

                else:

                    ###
                    tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Match found data_cat: {index:d}, {str(ind2):s}')
                    ###

                # Add matched event to list
                # Require the outputs as params['match_cat_outcols']
                out_cat.append([str(self.cmt_cat.evtm[ind2].to_numpy()[0]),float(self.cmt_cat.evla[ind2].to_numpy()[0]),
                                float(self.cmt_cat.evlo[ind2].to_numpy()[0]),float(self.cmt_cat.evdp[ind2].to_numpy()[0]),
                                float(self.cmt_cat.mag[ind2].to_numpy()[0]),str(self.cmt_cat.evid[ind2].to_numpy()[0]),
                                str(self.cmt_cat.id[ind2].to_numpy()[0]),str(row['EV_ID']),
                                int(self.cmt_cat['id_fmt_ctm'][ind2].to_numpy()[0]),int(self.cmt_cat['id_fmt_htm'][ind2].to_numpy()[0])])

        df_out = pd.DataFrame(sorted(out_cat),columns = self.params['match_cat_outcols'])

        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Match event catalog returned...')
        ###

        return df_out


    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
    def write_match_catalog(self):
        """
        Takes input matched event catalog dataframe
        Saves matched event catalog as csv, returns nothing

        :param self.match_catalog: input matched_catalog
        :type self.match_catalog: pd.df
        """
        from toolkit import Toolkit
        tk = Toolkit()

        log_statement = f'* {os.path.basename(__file__).split('.')[0]}, {str(inspect.stack()[0][3])}'

        cmt_outfile_name = f'{self.params['cmt_outfile']}_{str(self.params['year'])}.csv'

        ###
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ')
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ,  Saving the following DataFrame to: {str(cmt_outfile_name)}')
        tk.print_log(self.params, self.logfile, f'{log_statement:s}  ')
        tk.print_log(self.params, self.logfile, self.match_catalog)
        ###

        self.match_catalog.to_csv(str(cmt_outfile_name), index = False)
        return

    def execute(self):
        """ 
        Execute all functions in main, using arguments in params

        :param self.params:cdict from main code
        :type self.params: dict
        
        :return self.match_catalog: matched catalog dataframe object
        :type self.match_catalog: pd.df
        """
        # Execute already defined within __init__

        # Get CMT catalog
        self.cmt_cat = self.get_cmt_catalog()

        # Get DMT catalog
        self.data_cat = self.get_data_catalog_year()

        # Match the two catalogs
        self.match_catalog = self.get_match_catalog()

        # Write matched catalog to csv
        self.write_match_catalog()

        return self.match_catalog

def main():
    
    from toolkit import Toolkit
    tk = Toolkit()

    params_in = tk.get_params("params_in.yaml")
    match_catalog = MatchCatalog(params_in)
    evt_id_tab = match_catalog.execute()

if __name__ == '__main__':
    main()
    