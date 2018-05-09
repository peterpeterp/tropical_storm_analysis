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
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    local=True
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    local=False

sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
args = parser.parse_args()

for style in ['belanger']:
    if os.path.isfile('aew_detection/ERAint/ERAint_all_tracks_'+style+'.nc')==False or args.overwrite:
        # check for duplicates
        all_tracks={}
        for file_name in glob.glob('aew_detection/ERAint/*/track_info_belanger.nc'):
            tmp=da.read_nc(file_name)
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    all_tracks[id_]=track

        all_tracks=da.Dataset({'all_tracks':all_tracks})
        all_tracks.write_nc('aew_detection/ERAint/ERAint_all_tracks_'+style+'.nc',mode='w')
