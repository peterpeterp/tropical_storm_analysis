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

import cartopy.crs as ccrs

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import * ; reload(sys.modules['TC_support'])

data_path='data/CAM25/'
identifier='o6uk'
nc=da.read_nc(data_path+'item16222_daily_mean/item16222_daily_mean_'+identifier+'_2017-06_2017-10.nc')

lats = nc['global_latitude0'].values
lons = nc['global_longitude0'].values
lons[lons>180]-=360

o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
lon_0 = normalize180(o_lon_p-180.)

mslp=nc['item16222_daily_mean'].values[10,0,:,:]

rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p,central_rotated_longitude=lon_0)
plate_carree = ccrs.PlateCarree()
plt.figure(figsize=(6, 3))
ax = plt.axes(projection=rot_pole)
ax.set_global()
ax.coastlines()
ax.set_extent([np.min(lons),np.max(lons),np.min(lats),np.max(lats)],crs=plate_carree)
#ax.set_ylim((np.min(lats),np.max(lats)))
#ax.set_xlim((np.min(lons),np.max(lons)))
ax.plot(-74,20,'*r',transform=plate_carree)
ax.contourf(lons, lats, mslp, transform=plate_carree)  # didn't use transform, uh oh!

xs, ys, zs = rot_pole.transform_points(plate_carree,
                                 np.array([np.min(lons),np.max(lons)]),
                                 np.array([np.min(lats),np.max(lats)])).T
#rot_pole.transform_points(plate_carree,
#ax.set_xlim(xs)
#ax.set_ylim(ys)
plt.tight_layout()
plt.savefig('plots/cartopy_test.png')




lon = np.linspace(-80, 80, 25)
lat = np.linspace(30, 70, 25)
lon2d, lat2d = np.meshgrid(lon, lat)

data = np.cos(np.deg2rad(lat2d) * 4) + np.sin(np.deg2rad(lon2d) * 4)
data_crs = ccrs.PlateCarree()
# A rotated pole projection again...
projection = ccrs.RotatedPole(pole_longitude=-177.5, pole_latitude=37.5)
plt.figure(figsize=(6, 3))
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()

# ...but now using the transform argument
ax.contourf(lon, lat, data, transform=data_crs)
plt.savefig('plots/cartopy_test.png')
