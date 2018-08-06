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

try:
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/JRA55/'
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/JRA55/'

import aew_detection;    aew_detection = reload(aew_detection)

def toYearFraction(date):
    date=datetime(year=date.year,month=date.month,day=date.day,hour=date.hour)
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed/yearDuration

    return date.year + fraction

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
u_wind='easterly_waves/u_2008_sample.nc'
v_wind='easterly_waves/v_2008_sample.nc'
time_steps=12
level=1
x1,x2,y1,y2=-100,0,40,0
out_file='test'

# Read u and v wind components from file.
u = da.read_nc(u_wind)['var33'].ix[:,level,:,:].squeeze()
v = da.read_nc(v_wind)['var34'].ix[:,level,:,:].squeeze()

u.lon[u.lon>180]-=360
v.lon[v.lon>180]-=360
lon_rolling=len(u.lon)-np.where(u.lon<0)[0][0]
u.lon=np.roll(u.lon,lon_rolling,axis=-1)
v.lon=np.roll(v.lon,lon_rolling,axis=-1)
u.values=np.roll(u,lon_rolling,axis=-1)
v.values=np.roll(v,lon_rolling,axis=-1)

u=u[:,y1:y2,x1:x2]
v=v[:,y1:y2,x1:x2]
lons,lats=np.meshgrid(u.lon,u.lat)

u[np.isnan(u)]=0
v[np.isnan(v)]=0


nc=da.read_nc(u_wind)
time_=nc.time
dates=[num2date(t,units = nc.axes['time'].units,calendar = nc.axes['time'].calendar) for t in time_]

aews=aew_detection.aew_tracks(U=u,V=v,lats=lats,lons=lons,time_=time_,dates=dates,identifier='2008',working_dir='easterly_waves/2008/')
plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
ax = plt.axes(projection=plate_carree)
ax.set_global()
ax.coastlines()
ax.set_xlim(np.min(lons),np.max(lons))
ax.set_ylim(np.min(lats),np.max(lats))
aews.init_map(ax=ax,transform=plate_carree)

aews.detect(overwrite=True,thr_curv_vort=1*10**(-5))
aews.combine_tracks(overwrite=True)
aews.plot_detect_summary()
aews.plot_season()


plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes=plt.subplots(nrows=2,ncols=1,figsize=(10,8),subplot_kw={'projection': plate_carree})
for ax in axes:
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))

aews.plot_surrounding(axes=axes)
