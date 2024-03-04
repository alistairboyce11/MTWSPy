#!/usr/bin/env python3

'''

Called using obspyDMT_rem_inst_resp.py $sacfile $respfile $outfile $channel
Specific exact path to each.
'''


from obspy import read
from obspy import read_inventory
import sys

st_file_in=str(sys.argv[1])
resp_file=str(sys.argv[2])
st_file_out=str(sys.argv[3])
channel=str(sys.argv[4])

if channel == "LH":
    # For LH data:
    pre_filt = [0.01, 0.02, 0.1, 0.4]

if channel == "BH":
    # For BH data:
    pre_filt = [0.01, 0.02, 10, 20]

if channel == "HH":
    # For HH data:
    pre_filt = [0.01, 0.02, 40, 60]

trace=read(st_file_in, debug_headers=True)
stxml_file=resp_file
# simxml_file=sim_resp_file

# print('stationXML file: %s' % stxml_file)

inv_st = read_inventory(stxml_file, format="stationxml")

#### OUTPUT TO DISPLACEMENT ###
try:
    trace.remove_response(inventory=inv_st, pre_filt=pre_filt, output="DISP", water_level=None)  
except:
    print('Failed to remove response for: ' + st_file_in)
    print(st_file_out + ' not saved....')

trace.write(st_file_out, format='SAC')