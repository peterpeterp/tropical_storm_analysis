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
from shapely.geometry import Polygon

try:
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='easterly_waves/'
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/ERA_6hourly/u_v/'

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


curv_vort = da.read_nc(data_path+'ano_curv_vort_850-700-600_2008.nc')['curv_vort'][:,70000,35:1,-89:3].squeeze().values
nc = da.read_nc(data_path+'atl_u_v_850-700-600_2008.nc')
u=nc['U'][:,70000,35:1,-89:3].squeeze()
v=nc['V'][:,70000,35:1,-89:3].squeeze().values
lon,lat=u.lon,u.lat
u=u.values
lon[lon>180]-=360
lon[lon>180]-=360
lon_rolling=len(lon)-np.where(lon<0)[0][0]
lon=np.roll(lon,lon_rolling,axis=-1)
lon=np.roll(lon,lon_rolling,axis=-1)
u=np.roll(u,lon_rolling,axis=-1)
v=np.roll(v,lon_rolling,axis=-1)
lons,lats=np.meshgrid(lon,lat)
time_=nc.time
dates=[num2date(t,units = nc.axes['time'].units,calendar = nc.axes['time'].calendar) for t in time_]
xx,yy = np.meshgrid(lon,lat)
dx,dy = np.meshgrid(lon.copy()*0+np.mean(np.diff(lon,1)),lat.copy()*0+np.mean(np.diff(lat,1)))
dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
dy*=6371000*2*np.pi/360.

W = (u**2+v**2)**0.5
du_dx = (u-np.roll(u,1,axis=-1))/dx
du_dy = (u-np.roll(u,1,axis=-2))/dy
dv_dx = (v-np.roll(v,1,axis=-1))/dx
dv_dy = (v-np.roll(v,1,axis=-2))/dy
dW_dx = (W-np.roll(W,1,axis=-1))/dx
dW_dy = (W-np.roll(W,1,axis=-2))/dy
vo=dv_dx-du_dy
dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=-1))/dx
dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=-2))/dy
curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

aews=aew_detection.aew_tracks(U=u,V=v,vo=vo,curv_vort=curv_vort,curv_vort_advect=curv_vort_advect,lats=lats,lons=lons,smoothing_factor=3,coarsening_factor=2,time_=time_,dates=dates,identifier='2008',working_dir='aew_detection/ERAint/2008/')
plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
ax = plt.axes(projection=plate_carree)
ax.set_global()
ax.coastlines()
ax.set_xlim(np.min(lons),np.max(lons))
ax.set_ylim(np.min(lats),np.max(lats))
aews.init_map(ax=ax,transform=plate_carree)

aews.detect(overwrite=False,thr_curv_vort=np.nanpercentile(curv_vort,75))
aews.combine_tracks(overwrite=True)
aews.plot_detect_summary()
aews.plot_season()

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(16,3),subplot_kw={'projection': plate_carree})
for ax in axes:
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))

aews.plot_surrounding(axes=axes,time_steps=range(300,303),thr_curv_vort=np.nanpercentile(curv_vort,75))
