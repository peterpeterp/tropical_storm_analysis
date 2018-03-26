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

for style in ['contours','knutson2007']:
    if os.path.isfile('detection/CAR25/CAR25_all_tracks_'+style+'.nc')==False or overwrite:
        # check for duplicates
        all_tracks={}
        for identifier in [str(yr) for yr in range(1979,2017)]:
            tmp=da.read_nc('detection/CAR25/'+str(identifier)+'_CAR25/track_info_'+style+'.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    all_tracks[id_]=track

        all_tracks=da.Dataset({'all_tracks':all_tracks})
        all_tracks.write_nc('detection/CAR25/CAR25_all_tracks_'+style+'.nc',mode='w')
