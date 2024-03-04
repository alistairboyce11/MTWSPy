#!/usr/bin/env python3
'''

Called using /home/aboyce/bin/obspy_rotate_comps.py $form_event_dir $CHANNEL
Specific exact path to each.

'''
import obspy
from obspy import read
from obspy.core import Stream
from obspy import UTCDateTime
import obspy.signal
import obspy.signal.rotate
import obspy.geodetics.base
import obspy.geodetics
import numpy as np
import sys
import os.path
import glob
import shutil
import os
import time

# For some numpy thing:
# DeprecationWarning: `alltrue` is deprecated as of NumPy 1.25.0, and will be removed in NumPy 2.0. Please use `all` instead.

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# Get unique network and station values
def unique(list1):
    
    # initialize a null list
    unique_list = []
    stat_list = []
 
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x[:-6] not in stat_list:
            unique_list.append(x)
            stat_list.append(x[:-6])
   
    return stat_list, unique_list


def flatten_comprehension(matrix):
    # Flatten a list of lists.
    return [item for row in matrix for item in row]


def loc_chan_comp_priority():
    # This gives the priority of finding files to use for each station in the rotation to ZRT.
    LOCS = '00,10,20'.split(',')
    LH_CHA = 'LH,BH,HH'.split(',')
    BH_CHA = 'BH,HH'.split(',')
    HH_CHA = 'HH'.split(',')
    COMPS = 'Z,E,N,1,2'.split(',')

    LH_LOC_CHA_COMP = {'loc': LOCS, 'cha': LH_CHA, 'comp': COMPS}
    BH_LOC_CHA_COMP = {'loc': LOCS, 'cha': BH_CHA, 'comp': COMPS}
    HH_LOC_CHA_COMP = {'loc': LOCS, 'cha': HH_CHA, 'comp': COMPS}

    output_dict = {'LH': LH_LOC_CHA_COMP, 'BH': BH_LOC_CHA_COMP, 'HH': HH_LOC_CHA_COMP}

    return output_dict


# Collect list of vertical component files to see rotation scripts.
def collect_files(form_event_dir, channel):
    
    file_list = sorted(glob.glob(form_event_dir+'/??????????????_??_*.?H?'))

    stat_list, unique_list = unique(file_list)
    return stat_list, unique_list



def rotate_components(form_event_dir, channel, stat_list, channel_dict):

    # Make a directory for unused files...
    Not_useable = form_event_dir + '/Not_useable/'
    if not os.path.exists(Not_useable):
        os.makedirs(Not_useable)

    # Loop through the unique net-sta combinations.
    for s,stat in enumerate(stat_list):

        print(stat)
        # Make a list of files to test for the rotation given priority in channel_dict
        to_test=[]
        for loc in channel_dict['loc']:
            for cha in channel_dict['cha']:
                # print(loc,cha)

                f_list=[]
                for comp in channel_dict['comp']:
                    # Get all files for given station, with loc, cha, comps we need.
                    f_list.append(glob.glob(str(stat)+str(loc)+'.'+str(cha)+str(comp)))
                f_list=flatten_comprehension(f_list)
                # Dont move on if there is less than 3.

                if len(f_list) >= 3:
                    to_test.append(f_list)

        # Want to continue the loop when each set of files fails (==1)
        # Else (=0) and we are successful and so exit the while
        # Then delete other available files for that net-sta combo.
        # Iterator to start at zero
        i = -1
        FAIL_FLAG = 1
        while FAIL_FLAG == 1:
            if i < len(to_test)-1:
                # print('i',i,'i=i+1',i+1)
                i=i+1
                # Reset fail flag for set of files in to_test
                FAIL_FLAG = 0
            else:
                # Either nothing in to_test or 
                # reached end of options for while loop without success
                break

            # print('Trial files: ',stat,i,to_test[i])

            vertical=to_test[i][0][:-1]+'Z'
            horizontals=to_test[i].copy()
            horizontals.remove(vertical)


            # Read in the components
            if len(horizontals) >= 2: 
                vert_component = read(vertical, format='sac')

                if str(vertical[:-1]) + 'E' in horizontals:
                    east_component = read(str(vertical[:-1]) + 'E', format='sac')
                else:
                    east_component = read(str(vertical[:-1]) + '1', format='sac')

                if str(vertical[:-1]) + 'N' in horizontals:
                    north_component = read(str(vertical[:-1]) + 'N', format='sac')
                else:
                    north_component = read(str(vertical[:-1]) + '2', format='sac')

                # Check trace length and sample rate
                try:
                    if vert_component[0].stats.npts % 10 == 1 & east_component[0].stats.npts % 10 == 1 & north_component[0].stats.npts % 10 == 1:
                        print('Caught XXXXX1 sample points')
                        vert_component[0].data=vert_component[0].data[:-1]
                        east_component[0].data=east_component[0].data[:-1]
                        north_component[0].data=north_component[0].data[:-1]
                except:
                    pass

            #        Checks if sample rate is <10 and equalises
                try:
                    if channel == 'LH':
                        if vert_component[0].stats.sampling_rate < 1.0 or east_component[0].stats.sampling_rate < 1.0 or north_component[0].stats.sampling_rate < 1.0:
                            vert_component.resample(1)
                            east_component.resample(1)
                            north_component.resample(1)
                            print('Some issue with sample rate <1.0 for channel: ', channel)
                    if channel == 'BH':
                        if vert_component[0].stats.sampling_rate < 0.02 or east_component[0].stats.sampling_rate < 0.02 or north_component[0].stats.sampling_rate < 0.02:
                            vert_component.resample(0.02)
                            east_component.resample(0.02)
                            north_component.resample(0.02)
                            print('Some issue with sample rate <0.02 for channel: ', channel)
                    if channel == 'HH':
                        if vert_component[0].stats.sampling_rate < 0.005 or east_component[0].stats.sampling_rate < 0.005 or north_component[0].stats.sampling_rate < 0.005:
                            vert_component.resample(0.005)
                            east_component.resample(0.005)
                            north_component.resample(0.005)
                            print('Some issue with sample rate <0.005 for channel: ', channel)
                except:
                    pass

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
                    vert_component.trim(starttime=start_cut, endtime=end_cut)
                    east_component.trim(starttime=start_cut, endtime=end_cut)
                    north_component.trim(starttime=start_cut, endtime=end_cut)
                    # print('Trimmed to :' + str(vert_component[0].stats['npts']))

                except:
                    print('Trace error in seisN/seisE/seisZ')
                    print('Station: ' + str(vertical[:-1]) + ', Event: ' + str(s) + '\n')
                    print('Issues with starttime/endtime/npts' + '\n')
                    FAIL_FLAG = 1
                    continue

                loc_Z=vert_component[0].stats.location
                seisZ = vert_component.select(channel='*HZ',location=loc_Z)
                seisN_channel  = north_component[0].stats['channel']
                seisE_channel  = east_component[0].stats['channel']
                seisN = north_component.select(channel=seisN_channel,location=loc_Z)
                seisE = east_component.select(channel=seisE_channel,location=loc_Z)

                # print(seisZ[0])
                # print(seisN[0])
                # print(seisE[0])

                ori_seisN_temp=seisN[0].stats['sac']['cmpaz'] # .stats['orientation']
                ori_seisE_temp=seisE[0].stats['sac']['cmpaz'] # .stats['orientation']
                dip_seisE_temp=seisE[0].stats['sac']['cmpinc'] # .stats['dip']
                dip_seisN_temp=seisN[0].stats['sac']['cmpinc'] # .stats['dip']
                ori_seisZ_temp=seisZ[0].stats['sac']['cmpaz'] # .stats['orientation']
                dip_seisZ_temp=seisZ[0].stats['sac']['cmpinc'] # .stats['dip']

        #        make sure the orientations are different and dip is not vertical
                if ori_seisN_temp == ori_seisE_temp:            
                    print('Trace error in seisN/seisE')
                    print('Failed on N/E component')
                    FAIL_FLAG = 1
                    continue
                elif dip_seisE_temp == -90:
                    print('Trace error in seisN/seisE')
                    print('Failed on N/E component')
                    FAIL_FLAG = 1
                    continue
                elif dip_seisN_temp == -90:
                    print('Trace error in seisN/seisE')
                    print('Failed on N/E component')
                    FAIL_FLAG = 1
                    continue
                elif dip_seisZ_temp != -90.0:
                    print('Trace error in seisZ')
                    print('Failed on Z component')
                    FAIL_FLAG = 1
                    continue
                else:
                    FAIL_FLAG=0
      

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

                        FAIL_FLAG = 1
                        continue

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
                    seisN.stats['channel'] = str(channel) + 'N'
                    seisN.data = seisNtmp
                    seisE = seisE[0].copy()
                    seisE.stats['channel'] = str(channel) + 'E'
                    seisE.data = seisEtmp    

                    # should be setting the new orientation and dip to 0,90,-90 etc

                    # rotate components to from North and East to Radial and Transverse
                
                    # print('Rotating to RT using BAZ: ',BAZ)
                    [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(
                        seisN.data, seisE.data, BAZ)
                    seisR = seisN.copy()
                    seisT = seisN.copy()

                    seisR.stats['channel'] = str(channel) + 'R'
                    seisT.stats['channel'] = str(channel) + 'T'
                
                    seisR.data = seisRtmp
                    seisT.data = seisTtmp


                    ########## Sort out SAC headers #################

                    # Radial
                    # seisR.stats['channel'] = str(channel) + 'R'
                    seisR.stats['sac']['cmpaz']=AZ
                    seisR.stats['sac']['az']=AZ
                    seisR.stats['sac']['baz']=BAZ
                    seisR.stats['sac']['gcarc']=DISTDG
                    seisR.stats['sac']['kcmpnm']=str(channel)+'R'

                    # Tangential
                    # seisT.stats['channel'] = str(channel) + 'T'
                    seisT.stats['sac']['cmpaz']=AZ+90.0
                    seisT.stats['sac']['az']=AZ
                    seisT.stats['sac']['baz']=BAZ
                    seisT.stats['sac']['gcarc']=DISTDG
                    seisT.stats['sac']['kcmpnm']=str(channel)+'T'

                    # East
                    # seisE.stats['channel']=str(channel)+'E'
                    seisE.stats['sac']['cmpaz']=90.0
                    seisE.stats['sac']['az']=AZ
                    seisE.stats['sac']['baz']=BAZ
                    seisE.stats['sac']['gcarc']=DISTDG
                    seisE.stats['sac']['kcmpnm']=str(channel)+'E'

                    # North
                    # seisN.stats['channel']=str(channel)+'N'
                    seisN.stats['sac']['cmpaz']=0.0
                    seisN.stats['sac']['az']=AZ
                    seisN.stats['sac']['baz']=BAZ
                    seisN.stats['sac']['gcarc']=DISTDG
                    seisN.stats['sac']['kcmpnm']=str(channel)+'N'

                    # Vertical
                    seisZ[0].stats['sac']['az']=AZ
                    seisZ[0].stats['sac']['baz']=BAZ
                    seisZ[0].stats['sac']['gcarc']=DISTDG


                    filename_E = str(vertical[:-1]) + 'E'
                    filename_N = str(vertical[:-1]) + 'N'
                    filename_V = str(vertical[:-1]) + 'Z'
                    filename_R = str(vertical[:-1]) + 'R'
                    filename_T = str(vertical[:-1]) + 'T'

                    # Write over sac files as we now have everything sorted :)
                    seisE.write(filename_E, 'sac')
                    seisN.write(filename_N, 'sac')
                    seisZ.write(filename_V, 'sac')
                    seisR.write(filename_R, 'sac')
                    seisT.write(filename_T, 'sac')

            else:
                print('Not enough components found for: ' + str(vertical[:-1]))
                # Move horizontals - not sure we will ever be here...
                for file in horizontals:
                    try:
                        # print('Moving :' +str(file)+' -> '+str(Not_useable))
                        shutil.move(file, Not_useable)
                    except:
                        pass
                # Move Vertical - not sure we will ever be here...
                try:
                    # print('Moving :' +str(vertical)+' -> '+str(Not_useable))
                    shutil.move(vertical, Not_useable)
                except:
                    pass



        try:
            # Try to keep all but the ZRT component for each processed station.
            to_delete = sorted(list(set(glob.glob(stat+'*')) - set(vertical[:-1]+'Z') - set(vertical[:-1]+'R') - set(vertical[:-1]+'T')))
            # Keep ZRT components
            to_delete = sorted(glob.glob(stat+'*'))
            to_delete.remove(vertical[:-1]+'Z')
            to_delete.remove(vertical[:-1]+'R')
            to_delete.remove(vertical[:-1]+'T')

        except:
            to_delete = sorted(glob.glob(stat+'*'))

        # print(to_delete)

        if len(to_delete) > 1:
            for file in to_delete:
                try:
                    # print('Moving :' +str(file)+' -> '+str(Not_useable))
                    shutil.move(file, Not_useable)
                except:
                    pass




def main():

    num_args=len(sys.argv)

    if num_args != 3:
        print('Required :       Script.py <form_event_dir> <channel>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     /home/aboyce/d_data_and_docs/dmt/e2008/formatted/20081230194956')
        print('Options [2] :     LH')
        print('Number of arguments (' + str(num_args) +') incorrect... exit')
        exit('exiting....')

    form_event_dir=str(sys.argv[1])
    channel=str(sys.argv[2])
    # print(form_event_dir, channel)

    stat_list, unique_list=collect_files(form_event_dir, channel)
    # print(stat_list)

    channel_dict=loc_chan_comp_priority()[channel]



    rotate_components(form_event_dir, channel, stat_list, channel_dict)


if __name__ == '__main__':
    main()
