import os, time, glob
import toolkit
import inspect
import numpy as np

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def find_common_events(params_in):
    '''
    find common twin files for obs and syn
    - list obs twin files
    - list syn twin files
    - find common evid and save to file
    
    Inputs:
    params_in : dict
        For location of twin files based on paths set in .yaml file.

    Outputs:
    Text file : outfile (e.g., EVID.OST-MXT)
        List of matching event names from twin files.

    '''
    # Initiate logfile
    logfile = open_log_file(params_in)
    logfile = write_params_logfile(params_in, logfile)
    log_statement = str(os.path.basename(__file__).split('.')[0])

    # Counting stats
    io_logfile = toolkit.open_io_log_file(params_in)
    num_files_in  = 0
    num_obj_in    = 0
    num_files_out = 0
    num_obj_out   = 0



    ###
    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  ----------')
    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  PROCESSING.....')
    ###

    # Get root folder, channel names, outfile name
    fold_root = params_in['home'] + '/' + params_in['twin_loc'] + '/raw'
    cha_obs, cha_syn = str(params_in['twin_obs_out_loc'][-2:]) + str(params_in['component']), str(params_in['twin_syn_out_loc'][-2:]) + str(params_in['component'])
    outfile = params_in['mtf_outfilename']

    ###
    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  sending output to: {outfile:s}')
    ###

    # Read existing outfile, do not overwrite with new.
    evid_ext = []
    if os.path.exists(outfile):
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  outfile exists (remove it manually if outdated)')
        ###

        with open(outfile, 'r') as file:
            evid_ext = [line.strip() for line in file]

        n_ext = len(evid_ext)
        
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  {n_ext} evids loaded from existing {outfile}')
        ###
    else:
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  creating new outfile....')
        ###


    # Find obs and syn twin files
    fold_obs = fold_root + '/' +  'obs' + '/' +  cha_obs
    fold_syn = fold_root + '/' +  'syn' + '/' +  cha_syn

    files_obs = sorted(glob.glob(fold_obs + '/' + '*.' + cha_obs + '.twin'))
    files_syn = sorted(glob.glob(fold_syn + '/' + '*.' + cha_syn + '.twin'))

    # Find common evid
    # First take full path-filename, remove path, remove cha.twin
    # Find common event names using set.
    filenames_obs, filenames_syn = [x.split('/')[-1] for x in files_obs], [x.split('/')[-1] for x in files_syn]
    evnames_obs, evnames_syn = [x.split('.')[0] for x in filenames_obs], [x.split('.')[0] for x in filenames_syn]
    evid_com = list(set(evnames_obs) & set(evnames_syn))

    n_com = len(evid_com)

    ###
    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  {n_com} common evids found from {len(evnames_obs)} obs and {len(evnames_syn)} syn')
    ###

    # Combining with existing ones
    if evid_ext:
        # Doesnt produce a sorted list.....
        evid_com = list(set(evid_ext) | set(evid_com))
        n_com = len(evid_com)
        
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  {n_com} after joining with existing ones')
        ###
        
        num_files_in += 1
        num_obj_in = np.max([len(evnames_obs), len(evnames_syn)])


    # Sort the list:
    evid_com = sorted(evid_com)

    if n_com > 0:
        # Save to file
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  add the evid list to {outfile}')
        ###

        outfile = open(outfile,'w')
        for i in range(len(evid_com)):
            toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  adding: {evid_com[i]}')
            outfile.write(str(evid_com[i]) + '\n')
        outfile.close()

        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  {n_com} / {np.max([len(evnames_obs), len(evnames_syn)])} possible common evids written...')

        num_obj_out = n_com
        num_files_out += 1

    else:
        # No common ids found
        
        ###
        toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  **no common ids found ....')
        ###

    ###
    toolkit.print_log(params_in, logfile, f'{log_statement:s}  ,  FINISHED.....\n')
    ###


    statement  = toolkit.get_io_statement(log_statement, num_files_in, num_obj_in, num_files_out, num_obj_out)
    ###
    toolkit.print_log(params_in, io_logfile, statement)
    ###

    io_logfile.close()

    logfile.close()
    return


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 
def open_log_file(params_in):
    '''
    Return an open log file in log_loc/'code_start_time'/'filename'/match_twin_files.log
    '''
    # Get location of logfile
    lf_loc = params_in['home'] + '/' + params_in['log_loc'] + '/' + str(params_in['code_start_time']) + '/' + os.path.basename(__file__).split('.')[0]
    if not os.path.exists(lf_loc):
        os.makedirs(lf_loc)

    lf_name = lf_loc + '/match_twin_files.log'

    logfile = open(lf_name,'w')
    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #

def write_params_logfile(params_in, logfile):
    '''
    write key params to logfile. To be changed for each script
    '''
    justify = 30

    logfile.write(' ')
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('matched twin prefix',' : ',str(params_in['mtf_prefix']), x = justify) )
    logfile.write('{0:>{x}s} {1:s} {2:s}\n'.format('matched twin outfile',' : ',str(params_in['mtf_outfilename']), x = justify) )
    logfile.write('----------////               INPUT PARAMETERS                ////----------\n')

    return logfile


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ # 

def main():

    start_time = time.time()

    # Params go here.
    params_in = toolkit.get_params('params_in.yaml')

    # Make the list of common evids and write to file.
    find_common_events(params_in)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ ==  '__main__':
    main()



