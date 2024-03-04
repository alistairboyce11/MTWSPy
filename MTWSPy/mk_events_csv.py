import pandas as pd
import os
from obspy import UTCDateTime
import numpy as np
import toolkit
import inspect

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_log_file(params):
    '''
    Return an open log file in log_loc/'code_start_time'/'filename'/mk_csv_events_df.log
    '''

    lf_loc = params['home'] + '/' + params['log_loc'] + '/'+str(params['code_start_time']) + '/'+os.path.basename(__file__).split('.')[0]
    if not os.path.exists(lf_loc):
        os.makedirs(lf_loc)

    lf_name = lf_loc + '/mk_csv_events_df.log'

    logfile = open(lf_name,'w')
    
    return logfile

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #

def write_params_logfile(params, logfile):
    '''
    write key params to logfile. To be changed for each script
    '''
    justify = 30

    logfile.write(' ')
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

    params_list = ['min_year', 'max_year', 'cmt_infile', 'cmt_outfile', 'match_error_time', 'match_error_lat', 'match_error_lon', 'match_error_dep', 'match_error_mag']
    for k, param in enumerate(params_list):
        logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format(param,' : ',str(params[param]), x = justify) )
    
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_cmt_catalog(params, logfile):
    '''
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
    '''
    log_statement = '* ' + os.path.basename(__file__).split('.')[0] + ', '+str(inspect.stack()[0][3])
    
    ###
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Reading CMT catalog')
    ###
    
    cat = pd.read_csv(params['home']+'/'+params['cmt_infile'])

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
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  CMT catalog returned...')
    ###

    return cat


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_dmt_catalog_year(params, d_year, logfile):
    '''
    Read obspy dmt catalog_table.txt file into pandas dataframe
    With the following columns:

    'N','LAT', 'LON', 'DEP', 'DATETIME', 'MAG', 'EV_ID'

    Return dataframe
    '''
    log_statement = '* '+os.path.basename(__file__).split('.')[0]+', '+str(inspect.stack()[0][3])

    f_loc = params['home']+'/dmt/e'+str(d_year)+'/EVENTS-INFO/catalog_table.txt'

    if os.path.isfile(f_loc):  
        
        ###
        toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Reading DMT file : {f_loc:s}')
        toolkit.print_log(params, logfile, f'{log_statement:s}  ,  For year : {d_year:d}')
        ###

        file = open(f_loc,'r')
        lines = file.readlines()
        file.close()

        out_cat = []
        N = []; LAT = []; LON = []; DEP = []; DATETIME = []; MAG = []; EV_ID = []
        for line in lines:
            if 'Command' in line or 'pyDMT' in line or 'DATETIME' in line or '--' in line or len(line) <= 1:
                # Headers
                continue
            else:
                N.append(int(line.split()[0]))
                LAT.append(float(line.split()[1]))
                LON.append(float(line.split()[2]))
                DEP.append(float(line.split()[3]))
                DATETIME.append(UTCDateTime(str(line.split()[4])))
                MAG.append(float(line.split()[5]))
                EV_ID.append(str(line.split()[7]))

                out_cat.append([int(line.split()[0]), float(line.split()[1]), 
                                float(line.split()[2]), float(line.split()[3]), 
                                UTCDateTime(str(line.split()[4])), float(line.split()[5]),str(line.split()[7])])

        column_names = ['N','LAT', 'LON', 'DEP', 'DATETIME', 'MAG', 'EV_ID']
        dmt_cat = pd.DataFrame(sorted(out_cat),columns = column_names)

    else:
        # Return empty dataframe

        ###
        toolkit.print_log(params, logfile, f'{log_statement:s}  ,  No file found at : {f_loc:s}')
        toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Skipping...')
        ###

        dmt_cat = pd.DataFrame()

    return dmt_cat



# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_dmt_catalog(params,logfile):
    '''
    Determine the number of years required in matched catalog
    Use get_dmt_catalog_year to get dataframe for each year requested from dmt file.
    Sort by EV_ID and return dmt dataframe
    '''

    log_statement = '* ' + os.path.basename(__file__).split('.')[0] + ', '+str(inspect.stack()[0][3])
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Get sorted DMT catalog...')

    if params['min_year'] == params['max_year']:
        d_year = params['min_year']
        dmt_cat = get_dmt_catalog_year(params, d_year, logfile)
    else:
        data_years = range(params['min_year'], params['max_year']+1,1)
        for d, d_year in enumerate(data_years):

            ###
            toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Requesting DMT catalog year {d_year:d}...')
            ###

            out_cat = get_dmt_catalog_year(params, d_year, logfile)
            if d != 0:
                out_cat = pd.concat([dmt_cat,out_cat]).reset_index(drop = True)

            dmt_cat = out_cat

    dmt_cat = dmt_cat.sort_values(by = ['EV_ID'])

    ###
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Concatenated DMT catalog returned...')
    ###

    return dmt_cat

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def get_match_catalog(params, cat, dmt_cat, logfile):
    '''
    For every line in DMT catalog, search for matched event in cmt catalog
    Using error values in params (time, lat, lon, depth, magnitude)
    If match found, output to dataframe
    '''
    log_statement = '* ' + os.path.basename(__file__).split('.')[0] + ', '+str(inspect.stack()[0][3])
    
    ###
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Attempting CMT-DMT catalog match...')
    ###

    out_cat = []

    for index, row in dmt_cat.iterrows():

        ind = np.where(np.array(np.abs(row['LAT']-cat.evla) <= params['match_error_lat']) &
                       np.array(np.abs(row['LON']-cat.evlo) <= params['match_error_lon']) &
                       np.array(np.abs(row['DEP']-cat.evdp) <= params['match_error_dep']) &
                       np.array(np.abs(row['MAG']-cat.mag)  <= params['match_error_mag']))[0]
        if len(ind) < 1:
            
            ###
            toolkit.print_log(params, logfile, f'{log_statement:s}  ,  **No Match #1 for dmt_cat: {index:d}')
            ###

        else:
            # Now check the times.
            ind2 = []
            for m in ind:
                if np.abs(row['DATETIME']-cat['evtm'][m]) <= params['match_error_time']:
                    ind2.append(m)

            if len(ind2) < 1:
                # Shouldnt really get here.

                ###
                toolkit.print_log(params, logfile, f'{log_statement:s}  ,  **No Match #2 for dmt_cat: {index:d}')
                ###

                continue

            elif len(ind2) > 1:
                # Take first entry:
                ind2 = ind2[0]
                
                ###
                toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Taking first entry... {index:d}')
                toolkit.print_log(params, logfile, f'{log_statement:s}  ,  {str(cat.loc[ind2]):s}')
                ###

            else:

                ###
                toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Match found dmt_cat: {index:d}, {str(ind2):s}')
                ###

            # Add matched event to list
            # Require the outputs as params['match_cat_outcols']
            out_cat.append([str(cat.evtm[ind2].to_numpy()[0]),float(cat.evla[ind2].to_numpy()[0]),
                            float(cat.evlo[ind2].to_numpy()[0]),float(cat.evdp[ind2].to_numpy()[0]),
                            float(cat.mag[ind2].to_numpy()[0]),str(cat.evid[ind2].to_numpy()[0]),
                            str(cat.id[ind2].to_numpy()[0]),str(row['EV_ID']),
                            int(cat['id_fmt_ctm'][ind2].to_numpy()[0]),int(cat['id_fmt_htm'][ind2].to_numpy()[0])])

    df_out = pd.DataFrame(sorted(out_cat),columns = params['match_cat_outcols'])

    ###
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Match event catalog returned...')
    ###

    return df_out

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def write_match_catalog(params, df_out,logfile):
    '''
    Input: matched event catalog dataframe
    Saves: matched event catalog as csv
    Returns: Nothing
    '''
    log_statement = '* ' + os.path.basename(__file__).split('.')[0] + ', '+str(inspect.stack()[0][3])

    ###
    toolkit.print_log(params, logfile, f'{log_statement:s}  ')
    toolkit.print_log(params, logfile, f'{log_statement:s}  ,  Saving the following DataFrame to: '+str(params['cmt_outfile']))
    toolkit.print_log(params, logfile, f'{log_statement:s}  ')
    toolkit.print_log(params, logfile, df_out)
    ###

    df_out.to_csv(str(params['cmt_outfile']), index = False)
    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def execute(params_in):
    ''' 
    Execute all functions in main, using arguments in params

    Input: params dict from main code.
    
    Returns: matched catalog dataframe object.
    '''
    logfile = open_log_file(params_in)
    logfile = write_params_logfile(params_in,logfile)

    # Get CMT catalog
    cmt_cat = get_cmt_catalog(params_in,logfile)

    # Get DMT catalog
    dmt_cat = get_dmt_catalog(params_in,logfile)

    # Match the two catalogs
    match_catalog = get_match_catalog(params_in, cmt_cat, dmt_cat, logfile)

    # Write matched catalog to csv
    write_match_catalog(params_in, match_catalog, logfile)

    return match_catalog

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def main():
    '''
    Using values in params_in.yaml
    Match CMT and DMT catalog for given year
    Save output events csv file.
    '''
    # Params go here.
    params_in = toolkit.get_params('params_in.yaml')

    logfile = open_log_file(params_in)
    logfile = write_params_logfile(params_in,logfile)

    # Get CMT catalog
    cmt_cat = get_cmt_catalog(params_in,logfile)

    # Get DMT catalog
    dmt_cat = get_dmt_catalog(params_in,logfile)

    # Match the two catalogs
    match_catalog = get_match_catalog(params_in, cmt_cat, dmt_cat, logfile)  

    # Write matched catalog to csv
    write_match_catalog(params_in, match_catalog, logfile)
    return

if __name__ == '__main__':
    main()

