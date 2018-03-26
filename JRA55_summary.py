from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage


sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False

if os.path.isfile('detection/JRA55/JRA55_obs_track_info.nc')==False or overwrite:
    # check for duplicates
    obs_tracks={}
    for identifier in [str(yr) for yr in range(1979,2017)]:
        print(obs_tracks)
        tmp=da.read_nc('detection/JRA55/'+str(identifier)+'_JRA55/obs_track_info.nc')['obs_track_info']
        for storm in tmp.storm:
            obs_tracks[storm]=tmp[storm,:,:]

    obs_tracks=da.Dataset({'obs_tracks':obs_tracks})
    obs_tracks.write_nc('detection/JRA55/JRA55_obs_track_info.nc',mode='w')

for style in ['contours','knutson2007']:
    if os.path.isfile('detection/JRA55/JRA55_all_tracks_'+style+'.nc')==False or overwrite:
        # check for duplicates
        all_tracks={}
        for identifier in [str(yr) for yr in range(1979,2017)]:
            print(all_tracks)
            tmp=da.read_nc('detection/JRA55/'+str(identifier)+'_JRA55/track_info_'+style+'.nc')['track_info_'+style]
            for storm in tmp.storm:
                all_tracks[storm]=tmp[storm,:,:]

        all_tracks=da.Dataset({'all_tracks':all_tracks})
        all_tracks.write_nc('detection/JRA55/JRA55_all_tracks_'+style+'.nc',mode='w')
