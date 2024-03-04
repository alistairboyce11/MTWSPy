from obspy.core import UTCDateTime
import sys
if len(sys.argv) != 2:
    print('Provide depth....')
    print('>> python3 get_julday.py <TIME>')
    print('exiting...')
    sys.exit()

time = str(sys.argv[1])
#print(time)
julday=UTCDateTime(time).julday
print(julday)


