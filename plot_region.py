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
from shapely.geometry.polygon import Polygon

import cartopy.crs as ccrs

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/region/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import * ; reload(sys.modules['TC_support'])

# TC=da.read_nc('../data/Allstorms.ibtracs_all.v03r10.nc')
# tc_sel=TC.ix[np.where(TC['basin'][:,0]==0)[0]]
# tc_lat=tc_sel['lat_for_mapping']
# tc_lon=tc_sel['lon_for_mapping']
# tc_lon[tc_lon<0]+=360
# tc_wind=tc_sel['source_wind']
# tc_wind[np.isnan(tc_wind)]=-999

nc=da.read_nc('a0p1ga.pej6mar.nc')
lats = nc['global_latitude0'].values
lons = nc['global_longitude0'].values
grid_lats = nc['latitude0'].values
grid_lons = nc['longitude0'].values
lons[lons>180]-=360
o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
lon_0 = normalize180(o_lon_p-180.)
T=nc['field16'].values[0,0,:,:]

cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}



rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
plate_carree = ccrs.PlateCarree()
globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)

plt.close('all')
plt.figure(figsize=(6, 3))
ax = plt.axes(projection=globe)
ax.set_global()
ax.coastlines()
#ax.set_extent([np.min(grid_lons),np.max(grid_lons),np.min(grid_lats),np.max(grid_lats)],crs=rot_pole)
ax.contourf(grid_lons, grid_lats, T, transform=rot_pole)
ax.gridlines()

plt.tight_layout()
plt.savefig('region.png')

plt.close('all')
plt.figure(figsize=(6, 3))
ax = plt.axes(projection=rot_pole)
ax.set_global()
ax.coastlines()
#ax.set_extent([np.min(grid_lons),np.max(grid_lons),np.min(grid_lats),np.max(grid_lats)],crs=rot_pole)
ax.contourf(grid_lons, grid_lats, T, transform=rot_pole)
ax.gridlines()

plt.tight_layout()
plt.savefig('region_rot_pole.png')


plt.close('all')
plt.figure(figsize=(6, 3))
ax = plt.axes(projection=rot_pole)
ax.set_global()
ax.coastlines()
ax.set_extent([np.min(grid_lons)-10,np.max(grid_lons)+10,np.min(grid_lats)-10,np.max(grid_lats)+10],crs=rot_pole)
border=T.copy()*0+1
border[:,[0,-2]]=0
border[[0,-2],:]=0
ax.pcolormesh(grid_lons, grid_lats, border, cmap="bone",transform=rot_pole)
ax.gridlines()
tmp={}
for storm in tc_sel.storm:
    tmp[storm]=ax.plot(tc_lon[storm,:],tc_lat[storm,:],color=cat_colors[tc_cat(tc_wind[storm,:,0].max(),'wind')],alpha=0.2,linewidth=1,transform=plate_carree)

plt.tight_layout()
plt.savefig('region_tcs.png')
