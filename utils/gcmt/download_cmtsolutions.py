import requests
from bs4 import BeautifulSoup
import pandas as pd 
import os,sys

#url="https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=1976&mo=1&day=1&otype=ymd&oyr=1976&omo=1&oday=10&jyr=1976&jday=1&ojyr=1976&ojday=1&nday=3&lmw=0&umw=10&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd=0&uhd=1000&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"

class DownloadCMTSolutions:

    def __init__(self):
        pass

    def get_CMTSOLUTIONS(self, listyear, min_depth, max_depth, min_mag, max_mag):
        """
        Takes a list of years, depth and magnitude range
        Downloads all CMT solutions for each year
        Formats for SPECFEM3D 
   
        listyear :List of years for which to fetch CMT solutions
        listyear : list
        max_depth: Depth constraints on CMTs
        max_depth: float/int
        min_depth: Depth constraints on CMTs
        min_depth: float/int
        min_mag: Magnitude constraints on CMTs
        min_mag: float/int
        max_mag: Magnitude constraints on CMTs
        max_mag: float/int

        """

        for year in listyear:
            j=0
            outdir = './' + str(year) + '_CMT'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if not os.path.exists(outdir+'/HD0'):
                os.makedirs(outdir+'/HD0')
            
            file = open(outdir+"/"+str(year)+"_data.txt","a") # file containing all CMT files that fit criteria

            for i in range(1,13): #Loop over months
            # In the url, we can modify lat max & min, long min & max, depth lhd (min) & uhd (max), Mw max & min, Ms & Mb same for Mw. 
                print('Requesting year ' + str(year) + ', month ' + str(i) + ' for depths: '+str(min_depth)+'-'+str(max_depth)+'km, magnitudes: '+str(min_mag)+'-'+str(max_mag)+'Mw')

                if i == 4 or i == 6 or i == 9 or i == 11:
                    num_days=30
                elif i == 2: # February, careful with leap years.
                    num_days=29
                else:
                    num_days=31
                    
                url="https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr="+str(year)+"&mo="+str(i)+"&day=1&otype=ymd&oyr="+str(year)+"&omo="+str(i)+"&oday="+str(num_days)+"&jyr=1976&jday=1&ojyr=1976&ojday=1&nday=1&lmw="+str(min_mag)+"&umw="+str(max_mag)+"&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd="+str(min_depth)+"&uhd="+str(max_depth)+"&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"

                page = requests.get(url)
                soup = BeautifulSoup(page.content, "html.parser")
                table = soup.select_one('pre:nth-of-type(2)').text.splitlines()
                
                # make a table with the available CMTs corresponding to our criteria per month.
                
                fnum = int(j/13)+1
                CMT = open(outdir+"/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum)),"a")
                CMT_SPECFEM = open(outdir+"/HD0/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum))+"_hd0","a")

                for l, line  in enumerate(table):
                    if line != '':
                        j += 1
                        if j%13 == 0:
                            file.write(line+"\n")
                            CMT.write(line+"\n")
                            CMT_SPECFEM.write(line+"\n")
                            if l != len(table)-2:
                                fnum = int(j/13)+1
                                CMT = open(outdir+"/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum)),"a")
                                CMT_SPECFEM = open(outdir+"/HD0/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum))+"_hd0","a")

                        else:
                            if line[0:5] == ' PDEW':
                                line1 = line.replace(" PDEW"," PDE ")# SPECFEM3D globe cant read PDEW
                                file.write(line1+"\n")
                                CMT.write(line1+"\n")
                                CMT_SPECFEM.write(line1+"\n")
                            elif line[0:5] == ' PDEQ':
                                line1 = line.replace(" PDEQ"," PDE ")# SPECFEM3D globe cant read PDEQ
                                file.write(line1+"\n")
                                CMT.write(line1+"\n")
                                CMT_SPECFEM.write(line1+"\n")

                            elif line[0:5] == ' SWEQ':
                                line1 = line.replace(" SWEQ"," SWE ")# SPECFEM3D globe cant read SWEQ
                                file.write(line1+"\n")
                                CMT.write(line1+"\n")
                                CMT_SPECFEM.write(line1+"\n")

                            # Uncomment the following if half duration = 0 is required (As SPECFEM3D manual)	
                            elif line[0:4] == 'half': 
                                line2 = 'half duration:   0.0000'
                                file.write(line+"\n")
                                CMT.write(line+"\n")
                                CMT_SPECFEM.write(line2+"\n")
                            else:
                                file.write(line+"\n")
                                CMT.write(line+"\n")
                                CMT_SPECFEM.write(line+"\n")

            file.close()

        return

def main():

    download_CMT_solutions = DownloadCMTSolutions()
    num_args=len(sys.argv)

    if num_args <= 1:
        print('Required :       Script.py <year1> <year2> <year3> <etc...>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     2007')
        print('Options [2] :     2008')
        print('Options [3] :     2009')
        print('Number of arguments (' + str(num_args) +') incorrect... exit')
        exit('exiting....')

    listyear = []

    if len(sys.argv) > 1:
        print('Years requested: '+str(sys.argv[1:]))
        for i in range(1,len(sys.argv[1:])+1):
            print('Adding: '+str(sys.argv[i]))
            listyear.append(sys.argv[i])

    min_depth = 0
    max_depth = 700
    min_mag = 5.5
    max_mag = 7.5

    download_CMT_solutions.get_CMTSOLUTIONS(listyear, min_depth, max_depth, min_mag, max_mag)


if __name__  ==  '__main__':
    main()

