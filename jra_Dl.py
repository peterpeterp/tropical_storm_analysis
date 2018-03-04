#! /usr/bin/env python
#
# python script to download selected files from rda.ucar.edu
# after you save the file, don't forget to make it executable
#   i.e. - "chmod 755 <name_of_script>"
#
import sys
import os
import urllib2
import cookielib
from subprocess import Popen
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage



# if (len(sys.argv) != 2):
#   print "usage: "+sys.argv[0]+" [-q] password_on_RDA_webserver"
#   print "-q suppresses the progress message for each file that is downloaded"
#   sys.exit(1)
# #
passwd_idx=1
verbose=True
# if (len(sys.argv) == 3 and sys.argv[1] == "-q"):
#   passwd_idx=2
#   verbose=False
#
cj=cookielib.MozillaCookieJar()
opener=urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
#
# check for existing cookies file and authenticate if necessary
do_authentication=False
if (os.path.isfile("auth.rda.ucar.edu")):
  cj.load("auth.rda.ucar.edu",False,True)
  for cookie in cj:
    if (cookie.name == "sess" and cookie.is_expired()):
      do_authentication=True
else:
  do_authentication=True
if (do_authentication):
  login=opener.open("https://rda.ucar.edu/cgi-bin/login","email=peter.pfleiderer@climateanalytics.org&password=b3u!L1ronhave&action=login")



#
# save the authentication cookies for future downloads
# NOTE! - cookies are saved for future sessions because overly-frequent authentication to our server can cause your data access to be blocked
  cj.clear_session_cookies()
  cj.save("auth.rda.ucar.edu",True,True)
#
# download the data file(s)
print('still running')
os.chdir('/p/projects/tumble/carls/shared_folder/reanalysis/JRA55/')
sys.path.append('/p/projects/tumble/carls/shared_folder/reanalysis/')

from add_hybrid_coordinates import add_hybrid

try:
    years=[str(sys.argv[1])]
except:
    years=[str(yr) for yr in range(1979,2018)]

print(years)

decade_info=["***060100_***061018","***061100_***062018","***062100_***063018","***070100_***071018","***071100_***072018","***072100_***073118","***080100_***081018","***081100_***082018","***082100_***083118","***090100_***091018","***091100_***092018","***092100_***093018","***100100_***101018","***101100_***102018","***102100_***103118","***110100_***111018","***111100_***112018","***112100_***113018"]

month_info=["***060100_***063018","***070100_***073118","***080100_***083118","***090100_***093018","***100100_***103118","***110100_***113018"]

year_info="***010100_***123118"

mslp_loc="anl_surf/***/anl_surf.001_pres.reg_tl319."
u10_loc="anl_surf/***/anl_surf.033_ugrd.reg_tl319."
v10_loc="anl_surf/***/anl_surf.034_vgrd.reg_tl319."
t_loc="anl_mdl/***/anl_mdl.011_tmp.reg_tl319."
u_loc="anl_mdl/***/anl_mdl.033_ugrd.reg_tl319."
v_loc="anl_mdl/***/anl_mdl.034_vgrd.reg_tl319."


for year in years:
    print(year)
    print('***** surface *****')
    for name,loc in zip(['Psurf','U10','V10'],[mslp_loc,u10_loc,v10_loc]):
        # try if yearly files
        try:
            file=(loc+year_info).replace('***',year)
            idx=file.rfind("/")
            if (idx > 0):
                ofile=file[idx+1:]
            else:
                ofile=file
            sys.stdout.write("downloading "+ofile+"...")
            sys.stdout.flush()
            infile=opener.open("http://rda.ucar.edu/data/ds628.0/"+file)
            outfile=open(ofile,"wb")
            outfile.write(infile.read())
            outfile.close()
            Popen('cdo -f nc copy '+ofile+' '+ofile+'.nc',shell=True).wait()
            Popen('cdo -sellonlatbox,-110,-10,0,50 -selmon,6/11 '+ofile+'.nc atl_'+year+'_'+name+'.nc',shell=True).wait()
            Popen('rm '+ofile+' '+ofile+'.nc',shell=True).wait()
        # if not must be monthly
        except:
            for mon in month_info:
                file=(loc+mon).replace('***',year)
                idx=file.rfind("/")
                if (idx > 0):
                    ofile=file[idx+1:]
                else:
                    ofile=file
                sys.stdout.write("downloading "+ofile+"...")
                sys.stdout.flush()
                infile=opener.open("http://rda.ucar.edu/data/ds628.0/"+file)
                outfile=open(ofile,"wb")
                outfile.write(infile.read())
                outfile.close()
                Popen('cdo -f nc copy '+ofile+' '+ofile+'.nc',shell=True).wait()
                Popen('cdo -sellonlatbox,-110,-10,0,50 '+ofile+'.nc atl_'+year+'_'+name+'_'+mon+'.nc',shell=True).wait()
                Popen('rm '+ofile+' '+ofile+'.nc',shell=True).wait()

            Popen('cdo mergetime  atl_'+year+'_'+name+'* atl_'+year+'_'+name+'.nc',shell=True).wait()
            Popen('rm atl_'+year+'_'+name+'_*',shell=True).wait()

    # temperature
    print('***** upper *****')
    for name,loc,pressures,var_grib,var_128 in zip(['T','U','V'],[t_loc,u_loc,v_loc],[[85000,50000,20000],[85000],[85000]],['var11','var33','var34'],['130','131','132']):
        for dec_i,dec in enumerate(decade_info):
            file=(loc+dec).replace('***',year)
            idx=file.rfind("/")
            if (idx > 0):
                ofile=file[idx+1:]
            else:
                ofile=file
            sys.stdout.write("downloading "+ofile+"...")
            sys.stdout.flush()
            infile=opener.open("http://rda.ucar.edu/data/ds628.0/"+file)
            outfile=open(ofile,"wb")
            outfile.write(infile.read())
            outfile.close()
            Popen('cdo -f nc copy '+ofile+' '+ofile+'.nc',shell=True).wait()
            Popen('cdo -sellonlatbox,-110,-10,0,50 '+ofile+'.nc tmp_'+year+'_'+str(dec_i)+'_'+name+'.nc',shell=True).wait()
            add_hybrid('tmp_'+year+'_'+str(dec_i)+'_'+name+'.nc',var_grib,'hy_'+year+'_'+str(dec_i)+'_'+name+'.nc')
            Popen('cdo merge hy_'+year+'_'+str(dec_i)+'_'+name+'.nc oro.nc atl_'+year+'_Psurf.nc oro_'+year+'_'+str(dec_i)+'_'+name+'.nc',shell=True).wait()
            Popen('cdo ml2pl,'+','.join([str(prr) for prr in pressures])+' -settabnum,128 -chcode,1,134,6,129,11,130,33,131,34,132 oro_'+year+'_'+str(dec_i)+'_'+name+'.nc atl_'+year+'_'+str(dec_i)+'_'+name+'.nc',shell=True).wait()

            Popen('rm '+ofile+' '+ofile+'.nc tmp_'+year+'_'+str(dec_i)+'_'+name+'.nc hy_'+year+'_'+str(dec_i)+'_'+name+'.nc oro_'+year+'_'+str(dec_i)+'_'+name+'.nc',shell=True).wait()
            sys.stdout.write("done.\n")


        Popen('cdo mergetime  atl_'+year+'*_'+name+'.nc atl_'+year+'_'+name+'.nc',shell=True).wait()
        Popen('rm atl_'+year+'_*_'+name+'.nc',shell=True).wait()
