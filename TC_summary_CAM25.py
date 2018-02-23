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

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAM25/'
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *


try:
    identifiers=[sys.argv[1]]
except:
    identifieres=[ff.split('_')[-3] for ff in glob.glob(data_path+'/item3225_daily_mean/item3225_daily*')]

if os.path.isfile('detection/CAM25_summary.nc')==False:
    found_tracks={}
    for identifier in identifieres:
        found_tracks.update(da.read_nc('detection/'+str(identifier)+'_CAM25/track_info.nc'))

    summary=da.array(np.zeros([len(found_tracks),2]),axes=[found_tracks.keys(),['category','duration']],dims=['track','stat'])
    for id_,track in found_tracks.items():
        summary[id_,'category']=np.nanmax(track[:,'cat'])
        summary[id_,'duration']=len(np.where((track[:,'cat']>0) & (np.isfinite(track[:,'cat'])))[0])

    da.Dataset({'summary':summary}).write_nc('detection/CAM25_summary.nc',mode='w')

else:
    summary=da.read_nc('detection/CAM25_summary.nc')['summary']
