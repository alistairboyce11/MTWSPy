import glob,sys
import numpy as np

def get_sta_files(home, year):

    sta_evt_files = sorted(glob.glob(home+'/e'+str(year)+'/'+str(year)+'????_??????.a/info/station_event'))
    # print(sta_evt_files)
    print(len(sta_evt_files))

    sf=dict()

    for sta_evt in sta_evt_files:
        print(sta_evt)
        file=open(sta_evt,'r')
        lines=file.readlines()
        file.close()
        for line in lines[0:]:
            data=line.split(',')
            # print(data)
            nwk=str(data[0])
            sta=str(data[1])
            loc=str(data[2])
            cha=str(data[3])[0:2]
            lat=np.round(float(data[4]),4)
            lon=np.round(float(data[5]),4)
            elev=np.round(float(data[6]),1)
            dep=np.round(float(data[7]),1)
            # print(nwk, sta, loc, cha, lat, lon, elev, dep)

            # print(nwk in sf)
            if not nwk in sf:
                sf[nwk]=dict()
                sf = dict(sorted(sf.items()))

            if not sta in sf[nwk]:
                sf[nwk][sta]=dict()
                sf[nwk] = dict(sorted(sf[nwk].items()))

            if not cha in sf[nwk][sta]:
                sf[nwk][sta][cha]=dict()
                sf[nwk][sta] = dict(sorted(sf[nwk][sta].items()))

            if not loc in sf[nwk][sta][cha]:
                sf[nwk][sta][cha][loc]=dict()
                sf[nwk][sta][cha] = dict(sorted(sf[nwk][sta][cha].items()))

            sf[nwk][sta][cha][loc]['lat']=lat
            sf[nwk][sta][cha][loc]['lon']=lon
            sf[nwk][sta][cha][loc]['elev']=elev
            sf[nwk][sta][cha][loc]['dep']=dep

    return sf

def write_sta_dict(home, year, sf):
    outfile=home+'/e'+str(year)+'/STATIONS_'+str(year)+'_dmt'
    file=open(outfile,'w')

    for nwk in sf:
        # print(sf[nwk])
        for sta in sf[nwk]:
            for cha in sf[nwk][sta]:
                for loc in sf[nwk][sta][cha]:
                    # print(nwk, sta, loc, cha, sf[nwk][sta][cha][loc]['lat'], sf[nwk][sta][cha][loc]['lon'], sf[nwk][sta][cha][loc]['elev'], sf[nwk][sta][cha][loc]['dep'])
                    # print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta, nwk, str(sf[nwk][sta][cha][loc]['lat']), str(sf[nwk][sta][cha][loc]['lon']), str(sf[nwk][sta][cha][loc]['elev']), str(sf[nwk][sta][cha][loc]['dep'])))
                    file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta, nwk, str(sf[nwk][sta][cha][loc]['lat']), str(sf[nwk][sta][cha][loc]['lon']), str(sf[nwk][sta][cha][loc]['elev']), str(sf[nwk][sta][cha][loc]['dep'])))
    file.close()


def main():

    num_args=len(sys.argv)

    if num_args <= 2:
        print('Required :       Script.py <channel> <year1> <year2> <year3> <etc...>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     LH')
        print('Options [2] :     2008')
        print('Options [3] :     2009')
        print('Number of arguments (' + str(num_args) +') incorrect... exit')
        exit('exiting....')

    else:
        print('Making SPECFEM STATIONS file...')
    listyears = []

    if len(sys.argv) > 2:
        print('Years requested: '+str(sys.argv[2:]))
        for i in range(1,len(sys.argv[2:])+1):
            print('Adding: '+str(sys.argv[i+1]))
            listyears.append(sys.argv[i+1])

    # Set params:
    home='/home/aboyce/d_data_and_docs/dmt'

    # Loop through years
    for year in listyears:
        sta_dict=get_sta_files(home, year)
        
        write_sta_dict(home, year, sta_dict)


if __name__ == '__main__':
    main()