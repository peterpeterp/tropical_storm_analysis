from __future__ import print_function
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

from cdo import *   # python version

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    cdo = Cdo()
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *


if os.path.isfile('detection/ERA5_obs_track_info.nc')==False:
    # check for duplicates
    obs_tracks={}
    for identifier in [str(yr) for yr in range(2010,2017)]:
        print(obs_tracks)
        obs_tracks.update(da.read_nc('detection/'+str(identifier)+'_ERA5/obs_track_info.nc'))

    da.Dataset({'obs_tracks':obs_tracks}).write_nc('detection/ERA5_obs_track_info.nc',mode='w')

else:
    obs_tracks=da.read_nc('detection/ERA5_obs_track_info.nc')['obs_tracks']['obs_track_info']


categories={}
for cat in range(6):
    pos=np.where(obs_tracks==cat)
    categories[cat]=obs_tracks.ix[pos[0],pos[1],:]
