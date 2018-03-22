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

obs_tracks=da.read_nc('detection/JRA55/JRA55_obs_track_info.nc')['obs_tracks']

categories={}
for cat in range(6):
    pos=np.where(obs_tracks==cat)
    categories[cat]=obs_tracks.ix[pos[0],pos[1],:]


tc_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
fig,axes=plt.subplots(nrows=2,ncols=2)
for ax,vari in zip(axes.flatten(),['MSLP','VO','Wind10','T500']):
    for cat in range(5):
        tmp=categories[cat][:,:,vari].values
        tmp=tmp[np.isfinite(tmp)]
        hist=np.histogram(tmp,bins=20,density=True)
        ax.plot(hist[1][:-1],hist[0],color=tc_colors[cat])
    ax.set_xlabel(vari)

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_hists.png')
