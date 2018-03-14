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
    data_path='data/CAM25/'
    cdo = Cdo()
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *

all_tracks=da.read_nc('detection/CAM25/CAM25_all_tracks.nc')['all_tracks']

# intensity relations
plt.close('all')
fig,axes=plt.subplots(nrows=1,ncols=3,figsize=(9,2.5))
tmp=[]
for id_ in all_tracks.ID:
    track=all_tracks[id_]
    track=track[np.isfinite(track[:,'t']),:]
    track=track[track[:,'tc_cond']==3,:]
    tmp+=[axes[0].plot(range(track.shape[0]),track[:,'MSLP'],c='k',alpha=0.1)]
    tmp+=[axes[1].plot(range(track.shape[0]),track[:,'Wind10'],c='k',alpha=0.1)]
    tmp+=[axes[2].plot(track[:,'Wind10'],track[:,'MSLP'],'.k',linestyle='',alpha=0.05)]

y=np.arange(91000,101000,10)
axes[2].plot(3.5*((101000-y)/100.)**0.64,y,c='b')
axes[2].plot(2.3*((101000-y)/100.)**0.76,y,c='r')

axes[0].set_ylim((91000,101000))
axes[0].set_xlim((0,14))

axes[1].set_ylim((5,55))
axes[1].set_xlim((0,14))

axes[2].set_ylim((91000,101000))
axes[2].set_xlim((5,55))

plt.tight_layout(rect=(0,0,1,0.95))

plt.savefig('detection/CAM25/CAM25_intensities.png')


summary={'cat':[],'duration':[]}
for id_ in all_tracks.ID:
    summary['cat'].append(np.nanmax(all_tracks[id_,:,'cat']))
    summary['duration'].append(len(np.where((all_tracks[id_,:,'cat']>0) & (np.isfinite(all_tracks[id_,:,'cat'])))[0]))

runs=set([id_.split('_')[0] for id_ in all_tracks.ID])
n_runs=len(runs)

# show summary plots
plt.close('all')
fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(7,4))
ax=axes[0]
for i in range(1,6):
    ax.bar([i-0.5],len(np.where(np.array(summary['cat'])==i)[0])/float(len(all_tracks.ID)),width=1)
ax.set_xlim(0.5,5.5)
ax.set_ylabel('probability density')
ax.set_yscale('log')
ax.set_xlabel('TC Category')

ax=axes[1]
tmp=np.histogram(summary['duration'])
ax.hist(summary['duration'],10,normed=True)
ax.set_ylabel('probability density')
ax.set_xlabel('TC duration')

plt.suptitle(str(len(all_tracks.ID))+' Tropical cyclones detected in CAM25')
plt.tight_layout(rect=(0,0,1,0.95))
plt.savefig('detection/CAM25/CAM25_summary.png')

# show tracks on map
nc = Dataset(data_path+'/item16222_daily_mean/item16222_daily_mean_o6uk_2017-06_2017-10.nc')
lats = nc.variables['global_latitude0'][:]
lons = nc.variables['global_longitude0'][:]
rlats = nc.variables['latitude0'][:]
rlons = nc.variables['longitude0'][:]
rlons, rlats = np.meshgrid(rlons, rlats)
data = nc.variables['item16222_daily_mean'][:,0,:,:].squeeze()
data = np.ma.masked_values(data,-999.)
rotpole = nc.variables['rotated_pole0']

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
o_lon_p = rotpole.grid_north_pole_longitude
o_lat_p = rotpole.grid_north_pole_latitude

plt.close('all')
tmp,txt=[],[]
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
m= Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
           llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
           llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
x,y = m(lons,lats)
m.drawmapboundary(fill_color='1.')
m.drawmapboundary(fill_color='darkblue')
m.fillcontinents(color='darkgreen',lake_color='darkblue')
m.drawcoastlines(linewidth=0.3)
m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
ax.invert_yaxis()
ax.set_title('CAM25 '+str(n_runs)+' runs')

for id_ in all_tracks.ID:
    track=all_tracks[id_]
    track=track[np.isfinite(track[:,'t']),:]
    x_=[int(xx) for xx in track[:,'x']]
    y_=[int(yy) for yy in track[:,'y']]
    x,y = m(lons[y_,x_],lats[y_,x_])
    tmp+=tc_plot(m,x,y,track[:,'MSLP'].values,color_type=tc_pressure_color,alpha=0.4)


plt.tight_layout()
# for storm in tmp:
#     storm[0].set_alpha(0.1)
#     #tmp.set_linewidth(0.001)

plt.savefig('detection/CAM25/CAM25_summary_map.png')
