from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import URL_MAPPINGS
from obspy import UTCDateTime
import numpy as np
import sys,os
import subprocess

# for key in sorted(URL_MAPPINGS.keys()):
#     print("{0:<11} {1}".format(key,  URL_MAPPINGS[key])) 

class SpecfemStations:
    """
    Class to handle generation of specfem station files
    """

    def __init__(self):
        pass

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def unique(self, list1):
        """
        Get unique network and station values
        
        :param list1: list of network and station pairs
        :type list1: list
        :param unique_list: unique list of network and station pairs
        :type unique_list: list
        """
        # initialize a null list
        unique_list = []
        # traverse for all elements
        for x in list1:
            # check if exists in unique_list or not
            if x not in unique_list:
                unique_list.append(x)
    
        return unique_list

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
    def get_inventory(self, home, year, clients, network, station, 
                      location, channel):
        """
        Get inventory of stations for given year, data clients, 
        networks, locations, channels etc
        
        :param home: base string where inv will be saved
        :type home: str
        :param year: year of data 
        :type year: str/int
        :param clients: list of data clients/centers, "*"
        :type clients: list
        :param network: network string "*"
        :type network: str
        :param station station string "*"
        :type station: str
        :param location: location string "*"
        :type location: str
        :param channel: channel string "LH"
        :type channel: str
        """

        starttime = UTCDateTime(str(year)+"-01-01T00:00:00")
        endtime = UTCDateTime(str(year)+"-12-31T23:59:59")

        for i, client_name in enumerate(clients):
            print (i, str(client_name))
            client = Client(str(client_name))
            try:
                inv_1 = client.get_stations(network=network, station=station, 
                                            loc=location, channel=channel, 
                                            starttime=starttime, 
                                            endtime=endtime, level='channel')
            except:
                continue

            if i != 0:
                inventory.__add__(inv_1)
            else:
                inventory = inv_1

        return inventory


    def write_inventory(self, home, year, inventory):
        """
        Write inventory to station xml and Specfem format.

        :param home: base string where inv will be saved
        :type home: str
        :param year: year of data 
        :type year: str/int
        :param inventory: obspy station inventory to write to file
        :type inventory: obspy inventory
        """

        # Make the directory for the outputs
        if not os.path.exists(home+'/e'+str(year)):
            os.makedirs(home+'/e'+str(year))

        outfile = home + '/e' + str(year) + '/STATIONS_' + str(year) + '_obspy'

        # Write an output xml file    
        outfile_xml = outfile + '.xml'
        inventory.write(outfile_xml, format="STATIONXML")  

        file=open(outfile,'w')
        nwk_sta_list=[]
        # Loop through networks
        for i in range(len(inventory.networks)):
            net_name = inventory.networks[i]._code
            print(net_name)

            # Loop through stations and get name, lat, lon, elev
            for j in range(len(inventory.networks[i].stations)):
                sta_name = inventory.networks[i].stations[j]._code
                sta_lat = np.round(inventory.networks[i].stations[j]._latitude,4)
                sta_lon = np.round(inventory.networks[i].stations[j]._longitude,4)
                sta_elev = np.round(inventory.networks[i].stations[j]._elevation,1)
                sta_name = inventory.networks[i].stations[j]._code
                
                if str(net_name)+'_'+str(sta_name) in nwk_sta_list:
                    print('Found duplicate: '+net_name, sta_name)
                    continue
                else:
                    nwk_sta_list.append(str(net_name)+'_'+str(sta_name))
                    # print(net_name, sta_name)

                    channel_name = []
                    channel_elev = []
                    channel_depth = []
                    channel_loc = []

                    # Get unique channel info
                    for k in range(len(inventory.networks[i].stations[j].channels)):
                        
                        channel_depth.append(str(np.round(inventory.networks[i].stations[j].channels[k]._depth,1)))
                        channel_elev.append(str(np.round(inventory.networks[i].stations[j].channels[k]._elevation,1)))
                        channel_name.append(str(inventory.networks[i].stations[j].channels[k]._code))
                        channel_loc.append(str(inventory.networks[i].stations[j].channels[k]._location_code))
                    
                    channel_depth = self.unique(channel_depth)
                    channel_elev = self.unique(channel_elev)
                    channel_name = self.unique(channel_name)
                    channel_loc = self.unique(channel_loc)

                    # Check for duplicates
                    if len(channel_depth) == 0 or len(channel_elev) == 0 or len(channel_name) == 0 or len(channel_loc) == 0:
                        print('WE HAVE AN ISSUE WITH: '+net_name+'_'+sta_name)
                        print('Missing info... ')
                        print('channel_name: '+str(channel_name))
                        print('channel_loc: '+str(channel_loc))
                        print('channel_elev: '+str(channel_elev))
                        print('channel_depth: '+str(channel_depth))

                    elif len(channel_depth) > 1:
                        print('WE HAVE A DEPTH ISSUE WITH: '+net_name+'_'+sta_name)
                        print('Depths: '+str(channel_depth[:]))
                        print('Taking '+channel_depth[0])
                        print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
                        file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))

                    elif len(channel_elev) > 1:
                        print('WE HAVE A ELEVATION ISSUE WITH: '+net_name+'_'+sta_name)
                        print('Depths: '+str(channel_elev[:]))
                        print('Taking '+str(channel_elev[0]))
                        print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
                        file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))

                    elif len(channel_name) > 1:
                        print('WE HAVE A CHANNEL ISSUE WITH: '+net_name+'_'+sta_name)
                        print('Depths: '+str(channel_name[:]))
                        print('Taking '+str(channel_name[0]))
                        print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
                        file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))

                    elif len(channel_loc) > 1:
                        print('WE HAVE A LOCATION ISSUE WITH: '+net_name+'_'+sta_name)
                        print('Depths: '+str(channel_loc[:]))
                        print('Taking '+str(channel_loc[0]))
                        print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
                        file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))

                    # print to file
                    else:
                        print("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
                        file.write("{0:<5s}      {1:2s}      {2:>8s}   {3:>9s}   {4:>7s}  {5:>6s}\n".format(sta_name, net_name, str(sta_lat), str(sta_lon), str(channel_elev[0]), str(channel_depth[0])))
        file.close()

        
        # Try to sort the file by the network then the station name.
        try:
            subprocess.run(['sort', '-k2,2', '-k1,1', '-o', outfile, outfile], check=True)
            print(f"File '{outfile}' successfully sorted and overwritten.")
        except subprocess.CalledProcessError as e:
            print(f"Error sorting the file: {e}")

        return




def main():

    specfem_stations = SpecfemStations()
    num_args=len(sys.argv)

    if num_args <= 2:
        print('Required :       Script.py <channel> <year1> <year2> <year3> <etc...>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     HH?,BH?,LH?')
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
    home = os.getcwd()

    # clients = ['AUSPASS','ETH','GEOFON','GEONET','ICGC','INGV','IPGP','IRIS','KNMI','KOERI','LMU','NIEP','NOA','ORFEUS','RESIF','UIB-NORSAR','USP']
    clients = ['IRIS']

    channel = str(sys.argv[1])
    location = '*'
    network = "*"
    station = "*"
 
    # Loop through years
    for year in listyears:
        inv = specfem_stations.get_inventory(home, year, clients, network, station, location, channel)
        
        specfem_stations.write_inventory(home, year, inv)


if __name__ == '__main__':
    main()
