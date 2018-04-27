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


# Read u and v wind components from file.
u = da.read_nc('easterly_waves/U_2017.nc')['var33'].ix[:,0,:,:].squeeze()
#u=u.ix[:,:,100:160,range(100)+range(500,640)]

v = da.read_nc('easterly_waves/V_2017.nc')['var34'].ix[:,0,:,:].squeeze()
#v=v.ix[:,:,100:160,range(100)+range(500,640)]

u.lon[u.lon>180]-=360
v.lon[v.lon>180]-=360
lon_rolling=len(u.lon)-np.where(u.lon<0)[0][0]
u.lon=np.roll(u.lon,lon_rolling,axis=-1)
v.lon=np.roll(v.lon,lon_rolling,axis=-1)
u.values=np.roll(u,lon_rolling,axis=-1)
v.values=np.roll(v,lon_rolling,axis=-1)

u=u[:,30:0,-100:50]
v=v[:,30:0,-100:50]
lons,lats=np.meshgrid(u.lon,u.lat)


VO=TC_support.rel_vort(u.values[:,:,:],v.values[:,:,:],u.lat,u.lon)


u[np.isnan(u)]=0
v[np.isnan(v)]=0


def rel_vort(U,V,y,x):
    """adapted from wrftools: github.com/keltonhalbert/wrftools"""
    xx,yy = np.meshgrid(x,y)
    xx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
    yy*=6371000*2*np.pi/360.
    dy=np.gradient(yy)
    dx=np.gradient(xx)

    vort=U.copy()*np.nan
    for i in range(U.shape[0]):
        du = np.gradient( U[i,:,:] )
        dv = np.gradient( V[i,:,:] )
        vort[i,:,:]= dv[-1]/dx[-1] - du[-2]/dy[-2]
    return vort


xx,yy = np.meshgrid(u.lon,u.lat)
dx,dy = np.meshgrid(u.lon.copy()*0+np.mean(np.diff(u.lon,1)),u.lat.copy()*0+np.mean(np.diff(u.lat,1)))
dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
dy*=6371000*2*np.pi/360.


t=100
u,v=u.ix[t,:,:],v.ix[t,:,:]

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

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes = plt.subplots(nrows=3,ncols=5,figsize=(10,2),subplot_kw={'projection': plate_carree})
for ax in axes.flatten():
    #ax.axis('off')
    ax.set_global()
    ax.coastlines(edgecolor='magenta')
    ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
    ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
    ax.set_extent([np.min(u.lon),np.max(u.lon),np.min(u.lat),np.max(u.lat)],crs=plate_carree)

axes[0,0].pcolormesh(lons,lats,u,vmin=-30,vmax=30)
axes[1,0].pcolormesh(lons,lats,du_dy,vmin=-0.00001,vmax=0.00001)
axes[2,0].pcolormesh(lons,lats,du_dx,vmin=-0.00001,vmax=0.00001)

axes[0,1].pcolormesh(lons,lats,v,vmin=-30,vmax=30)
axes[1,1].pcolormesh(lons,lats,dv_dy,vmin=-0.00001,vmax=0.00001)
axes[2,1].pcolormesh(lons,lats,dv_dx,vmin=-0.00001,vmax=0.00001)

axes[0,2].pcolormesh(lons,lats,W,vmin=-30,vmax=30)
axes[1,2].pcolormesh(lons,lats,dW_dy,vmin=-0.00001,vmax=0.00001)
axes[2,2].pcolormesh(lons,lats,dW_dx,vmin=-0.00001,vmax=0.00001)

axes[0,3].pcolormesh(lons,lats,vo,vmin=-0.00001,vmax=0.00001)
axes[1,3].pcolormesh(lons,lats,dv_dx-du_dy,vmin=-0.00001,vmax=0.00001)
axes[2,3].pcolormesh(lons,lats,dW_dx.copy()*0,vmin=-0.00001,vmax=0.00001)

axes[0,4].pcolormesh(lons,lats,shear_vort,vmin=-0.00001,vmax=0.00001)
axes[1,4].pcolormesh(lons,lats,curv_vort,vmin=-0.00001,vmax=0.00001)
axes[2,4].pcolormesh(lons,lats,dW_dx.copy()*0,vmin=-0.00001,vmax=0.00001)

plt.tight_layout()
plt.savefig('test.png')
