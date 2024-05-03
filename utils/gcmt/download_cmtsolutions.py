import requests
from bs4 import BeautifulSoup
import pandas as pd 
import os,sys

# ATTENTION : défaut de ce programme c'est qu'il ne prend pas en compte quand nos critères sont trop larges. En effet sur le site de GlobalCMT
# il va y avoir plusieurs pages de CMT possible et avec mes critères je n'avais qu'une page à chaque fois. Mais il faut faire attention dans 
# certains cas. 

#url="https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=1976&mo=1&day=1&otype=ymd&oyr=1976&omo=1&oday=10&jyr=1976&jday=1&ojyr=1976&ojday=1&nday=3&lmw=0&umw=10&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd=0&uhd=1000&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"
#listannee=[2012,2013,2014,2015,2016,2017,2018,2019,2021,2022]

def get_CMTSOLUTIONS(listannee, min_depth, max_depth, min_mag, max_mag):

    for annee in listannee:
        j=0
        outdir = './' + str(annee) + '_CMT'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(outdir+'/HD0'):
            os.makedirs(outdir+'/HD0')
        
        fichier=open(outdir+"/"+str(annee)+"_data.txt","a")#fichier qui aura tous les CMT correspondant à mes critères. 
        # CMT=open(outdir+"/CMTSOLUTION_run001","a")
        # CMT_SPECFEM=open(outdir+"/HD0/CMTSOLUTION_run001_hd0","a")

        for i in range(1,13):#boucle sur les mois pour pouvoir gérer le nombre de jours dans chacun.
        # Dans url, on peut modifier lat max et min, long min et max, la profondeur lhd (min) et uhd (max), Mw max et min,Ms et Mb pareil que pour Mw. 
            print('Requesting year ' + str(annee) + ', month ' + str(i) + ' for depths: '+str(min_depth)+'-'+str(max_depth)+'km, magnitudes: '+str(min_mag)+'-'+str(max_mag)+'Mw')

            if i==4 or i==6 or i==9 or i==11:
                num_days=30
            elif i==2:#mois de Février, attention aux années bisextiles.
                num_days=29
            else:
                num_days=31
                
            url="https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr="+str(annee)+"&mo="+str(i)+"&day=1&otype=ymd&oyr="+str(annee)+"&omo="+str(i)+"&oday="+str(num_days)+"&jyr=1976&jday=1&ojyr=1976&ojday=1&nday=1&lmw="+str(min_mag)+"&umw="+str(max_mag)+"&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd="+str(min_depth)+"&uhd="+str(max_depth)+"&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"

            page = requests.get(url)
            soup = BeautifulSoup(page.content, "html5lib")
            tableau = soup.select_one('pre:nth-of-type(2)').text.splitlines()
            # En gros ici ça fait un tableau avec les CMT disponible correspondant à nos critères par mois. 
            # Le programme parcourt chaque ligne du tableu et on obtient des CMTsolution qu'on peut utiliser dans specfem
            
            # print(len(tableau))

            fnum=int(j/13)+1
            CMT=open(outdir+"/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum)),"a")
            CMT_SPECFEM=open(outdir+"/HD0/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum))+"_hd0","a")


            for l, ligne  in enumerate(tableau):
                if ligne != '':
                    j+=1
                    if j%13==0:
                        fichier.write(ligne+"\n")
                        CMT.write(ligne+"\n")
                        CMT_SPECFEM.write(ligne+"\n")
                        if l != len(tableau)-2:
                            fnum=int(j/13)+1
                            CMT=open(outdir+"/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum)),"a")
                            CMT_SPECFEM=open(outdir+"/HD0/CMTSOLUTION_run"+str("{0:0=3d}".format(fnum))+"_hd0","a")

                    else:
                        if ligne[0:5]==' PDEW':
                            ligne1=ligne.replace(" PDEW"," PDE ")# SPECFEM3D globe ne lit pas PDEW => ça crée une erreur
                            fichier.write(ligne1+"\n")
                            CMT.write(ligne1+"\n")
                            CMT_SPECFEM.write(ligne1+"\n")
                        elif ligne[0:5]==' PDEQ':
                            ligne1=ligne.replace(" PDEQ"," PDE ")# SPECFEM3D globe ne lit pas PDEQ => ça crée une erreur
                            fichier.write(ligne1+"\n")
                            CMT.write(ligne1+"\n")
                            CMT_SPECFEM.write(ligne1+"\n")

                        elif ligne[0:5]==' SWEQ':
                            ligne1=ligne.replace(" SWEQ"," SWE ")# SPECFEM3D globe ne lit pas SWEQ => ça crée une erreur
                            fichier.write(ligne1+"\n")
                            CMT.write(ligne1+"\n")
                            CMT_SPECFEM.write(ligne1+"\n")

                        # à décommenter si on a besoin d'un half duration =0 dans SPECFEM3D (ce qui est recommandé dans le manuel)	
                        elif ligne[0:4]=='half': 
                            ligne2='half duration:   0.0000'
                            fichier.write(ligne+"\n")
                            CMT.write(ligne+"\n")
                            CMT_SPECFEM.write(ligne2+"\n")
                        else:
                            fichier.write(ligne+"\n")
                            CMT.write(ligne+"\n")
                            CMT_SPECFEM.write(ligne+"\n")

        fichier.close()


def main():

    num_args=len(sys.argv)

    if num_args <= 1:
        print('Required :       Script.py <year1> <year2> <year3> <etc...>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     2007')
        print('Options [2] :     2008')
        print('Options [3] :     2009')
        print('Number of arguments (' + str(num_args) +') incorrect... exit')
        exit('exiting....')

    listannee = []

    if len(sys.argv) > 1:
        print('Years requested: '+str(sys.argv[1:]))
        for i in range(1,len(sys.argv[1:])+1):
            print('Adding: '+str(sys.argv[i]))
            listannee.append(sys.argv[i])


    min_depth = 0
    max_depth = 700
    min_mag = 5.5
    max_mag = 7.5

    get_CMTSOLUTIONS(listannee, min_depth, max_depth, min_mag, max_mag)


if __name__ == '__main__':
    main()

