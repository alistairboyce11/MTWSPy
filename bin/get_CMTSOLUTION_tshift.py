from obspy.core import UTCDateTime
import sys
if len(sys.argv) != 3:
    print('Provide correct arguements....')
    print('>> python3 get_CMTSOLUTION_tshift.py <TIME> <TSHIFT>')
    print('exiting...')
    print(len(sys.argv))
    print(sys.argv)
    sys.exit()

time = str(sys.argv[1])
tshift = float(sys.argv[2])
#print(time)
time_shifted = UTCDateTime(time) + tshift
# print(time_shifted.year, time_shifted.month, time_shifted.day, time_shifted.julday, time_shifted.hour, time_shifted.minute, time_shifted.second, time_shifted.microsecond)
print('{0:4d} {1:02d} {2:02d} {3:03d} {4:02d} {5:02d} {6:02d} {7:6d}'.format(time_shifted.year, time_shifted.month, time_shifted.day, time_shifted.julday, time_shifted.hour, time_shifted.minute, time_shifted.second, time_shifted.microsecond))

