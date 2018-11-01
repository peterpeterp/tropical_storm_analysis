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

for base_path in ['/Users/peterpfleiderer/Projects/tropical_cyclones','tropical_cyclones','/p/projects/tumble/carls/shared_folder/tropical_cyclones']:
    try:
        os.chdir(base_path)
    except:
        pass

sys.path.append('tropical_storm_analysis')
sys.path.append('tropical_storm_analysis/feature_tracking')

import TC_support ;  TC_support = reload(TC_support)
import TC_detection;    TC_detection = reload(TC_detection)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False

identifiers = [nn.split('_')[-3] for nn in glob.glob('../data/WAH/batch_755/region/item16222_6hrly_inst/*')]

for style in ['contours']:
    if os.path.isfile('detection/ATL/ATL_all_tracks_'+style+'.nc')==False or overwrite:
        # check for duplicates
        all_tracks={}
        for identifier in identifiers:
            tmp=da.read_nc('detection/ATL/'+str(identifier)+'/track_info_'+style+'.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    all_tracks[id_]=track

        all_tracks=da.Dataset({'all_tracks':all_tracks})
        all_tracks.write_nc('detection/ATL/ATL_all_tracks_'+style+'.nc',mode='w')
