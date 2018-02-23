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

if os.path.isfile('detection/CAM25_all_tracks.nc')==False:
    found_tracks={}
    for identifier in identifieres:
        found_tracks.update(da.read_nc('detection/'+str(identifier)+'_CAM25/track_info.nc'))


    da.Dataset(found_tracks).write_nc('detection/CAM25_all_tracks.nc',mode='w')
else:
    summary=da.read_nc('detection/CAM25_all_tracks.nc')['summary']


# plt.close('all')
# fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(7,4))
# ax=axes[0]
# for i in range(1,6):
#     ax.bar([i-0.5],len(np.where(summary[:,'category'].values==i)[0])/float(summary.shape[0]),width=1)
# ax.set_xlim(0.5,5.5)
# ax.set_ylabel('TCs per season')
# ax.set_xlabel('TC Category')
#
# ax=axes[1]
# tmp=np.histogram(summary[:,'duration'].values)
# ax.hist(summary[:,'duration'].values,10,normed=True)
# ax.set_ylabel('TCs density')
# ax.set_xlabel('TC duration')
#
# plt.suptitle('Tropical cyclones detected in CAM25')
# plt.tight_layout(rect=(0,0,1,0.95))
# plt.savefig('detection/CAM25_summary.png')
