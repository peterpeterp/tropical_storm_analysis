from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy

from scipy import fftpack

from windspharm.standard import VectorWind
from skimage.feature import peak_local_max

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
import TC_support ;  TC_support = reload(TC_support)
import tc_detection;    tc_detection = reload(tc_detection)

from scipy.signal import butter, lfilter


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--u_wind','-u',help='u wind',required=True)
parser.add_argument('--v_wind','-v',help='v wind',required=True)

parser.add_argument('--out_file','-o',help='',default='test.png')

parser.add_argument('--x1',help='',type=int,required=False)
parser.add_argument('--x2',help='',type=int,required=False)
parser.add_argument('--y1',help='',type=int,required=False)
parser.add_argument('--y2',help='',type=int,required=False)
parser.add_argument('--level',help='',type=int,required=False)
parser.add_argument('--time_steps','-ts',help='time steps',nargs='+',type=int,required=False)
args = parser.parse_args()

print(args)

# Read u and v wind components from file.
u = da.read_nc(args.u_wind)['var33'].ix[args.time_steps,args.level,:,:].squeeze()
v = da.read_nc(args.v_wind)['var34'].ix[args.time_steps,args.level,:,:].squeeze()

u.lon[u.lon>180]-=360
v.lon[v.lon>180]-=360
lon_rolling=len(u.lon)-np.where(u.lon<0)[0][0]
u.lon=np.roll(u.lon,lon_rolling,axis=-1)
v.lon=np.roll(v.lon,lon_rolling,axis=-1)
u.values=np.roll(u,lon_rolling,axis=-1)
v.values=np.roll(v,lon_rolling,axis=-1)

u=u[args.y1:args.y2,args.x1:args.x2]
v=v[args.y1:args.y2,args.x1:args.x2]
lons,lats=np.meshgrid(u.lon,u.lat)


#VO=TC_support.rel_vort(u.values[:,:],v.values[:,:],u.lat,u.lon)

u[np.isnan(u)]=0
v[np.isnan(v)]=0

xx,yy = np.meshgrid(u.lon,u.lat)
dx,dy = np.meshgrid(u.lon.copy()*0+np.mean(np.diff(u.lon,1)),u.lat.copy()*0+np.mean(np.diff(u.lat,1)))
dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
dy*=6371000*2*np.pi/360.

W = (u**2+v**2)**0.5

du_dx = (u-np.roll(u,1,axis=1))/dx
du_dy = (u-np.roll(u,1,axis=0))/dy

dv_dx = (v-np.roll(v,1,axis=1))/dx
dv_dy = (v-np.roll(v,1,axis=0))/dy

dW_dx = (W-np.roll(W,1,axis=1))/dx
dW_dy = (W-np.roll(W,1,axis=0))/dy

vo=dv_dx-du_dy

shear_vort=v/W*dW_dx-u/W*dW_dy
curv_vort=vo-shear_vort

asp=abs(args.x1-args.x2)*2/(abs(args.y1-args.y2)*3)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(5*asp,5),subplot_kw={'projection': plate_carree})
for ax in axes.flatten():
    #ax.axis('off')
    ax.set_global()
    ax.coastlines(edgecolor='magenta')
    ax.set_extent([np.min(u.lon),np.max(u.lon),np.min(u.lat),np.max(u.lat)],crs=plate_carree)

axes[0,0].pcolormesh(lons,lats,u,vmin=-30,vmax=30); axes[0,0].annotate('U', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[1,0].pcolormesh(lons,lats,v,vmin=-30,vmax=30); axes[1,0].annotate('V', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[2,0].pcolormesh(lons,lats,W,vmin=-30,vmax=30); axes[2,0].annotate('W', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')

axes[0,1].pcolormesh(lons,lats,vo,vmin=-0.00001,vmax=0.00001); axes[0,1].annotate('rel. vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[1,1].pcolormesh(lons,lats,shear_vort,vmin=-0.00001,vmax=0.00001); axes[1,1].annotate('shear vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[2,1].pcolormesh(lons,lats,curv_vort,vmin=-0.00001,vmax=0.00001); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')

plt.tight_layout()
plt.savefig(args.out_file)
