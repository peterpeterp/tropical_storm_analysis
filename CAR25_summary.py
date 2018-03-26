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
    data_path='data/CAR25/'
    local=True
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_717/region/'
    local=False
    
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
        for identifier in [ff.split('_')[-3] for ff in glob.glob(data_path+'/item16222_6hrly_inst/item16222_6hrly_inst*')]:
            tmp=da.read_nc('detection/CAR25/'+str(identifier)+'_CAR25/track_info_'+style+'.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    all_tracks[id_]=track

        all_tracks=da.Dataset({'all_tracks':all_tracks})
        all_tracks.write_nc('detection/CAR25/CAR25_all_tracks_'+style+'.nc',mode='w')
