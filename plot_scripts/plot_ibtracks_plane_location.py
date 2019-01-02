from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import matplotlib.ticker as mticker
import cartopy

import cartopy.crs as ccrs

os.chdir('/Users/peterpfleiderer/Projects/tropical_cyclones/')
sys.path.append('/Users/peterpfleiderer/Projects/tropical_cyclones/TC_scripts')
import TC_support ;  TC_support = reload(TC_support)

if 'TC' not in globals():
    TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
    tc_sel=TC.ix[np.where(TC['basin'][:,0]==0)[0]]
    tc_sel=tc_sel.ix[tc_sel['season']>=1979]
    tc_sel=tc_sel.ix[tc_sel['season']>=1979]
    tc_lat=tc_sel['lat_for_mapping']
    tc_lon=tc_sel['lon_for_mapping']
    tc_lon[tc_lon<0]+=360
    tc_wind=tc_sel['source_wind']
    tc_wind[np.isnan(tc_wind)]=-999

    tc_sel_cat=np.array(TC_support.tc_cat(np.nanmin(tc_sel['source_pres'],axis=(1,2)),'pressure'))
    tc_sel=tc_sel.ix[np.where(tc_sel_cat>0)]

    nc=da.read_nc('data/CAR25/item16222_6hrly_inst/item16222_6hrly_inst_p014_2017-06_2017-10.nc')
    lats = nc['global_latitude1'].values
    lons = nc['global_longitude1'].values
    grid_lats = nc['latitude1'].values
    grid_lons = nc['longitude1'].values
    lons[lons>180]-=360
    o_lon_p = nc['rotated_pole1'].attrs['grid_north_pole_longitude']
    o_lat_p = nc['rotated_pole1'].attrs['grid_north_pole_latitude']
    lon_0 = TC_support.normalize180(o_lon_p-180.)

    rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
    globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)

cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}

# sparial stuff
plt.close('all')
plt.figure(figsize=(10,5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()
gl=ax.gridlines(color='#BF3EFF',linewidth=1)
gl.ylocator = mticker.FixedLocator(np.arange(-20,60,10))
gl.xlocator = mticker.FixedLocator(np.arange(-120,60,10))
for yy in np.arange(-10,40,10):   ax.text(5,yy,str(yy),color='#BF3EFF',transform=ccrs.PlateCarree())
for xx in np.arange(-100,0,10):   ax.text(xx,-8,str(xx),color='#BF3EFF',transform=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
ax.set_xlim(-120,50)
ax.set_ylim(-15,50)

#########
# Houston
#########

reg_atl=Polygon([(-100,25),(-90,25),(-90,35),(-100,35),(-100,25)])
ax.add_geometries([reg_atl], ccrs.PlateCarree(), color='red',alpha=1,facecolor='none')

ax.lines=[]
for storm in tc_sel.storm:
    plot = False
    for ti in range(tc_lon[storm,:].shape[0]):
        xx = tc_lon[storm,ti]
        if xx>180:
            xx-=360
        if reg_atl.contains(Point(xx,tc_lat[storm,ti])):
            plot = True

    if plot:
        ax.plot(tc_lon[storm,:],tc_lat[storm,:],color=cat_colors[TC_support.tc_cat(np.nanmin(tc_sel['source_pres'][storm,:,:]),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_tracks_houston.png',dpi=300)

ax.lines=[]
for storm in tc_sel.storm:
    plot = False
    for ti in range(tc_lon[storm,:].shape[0]):
        xx = tc_lon[storm,ti]
        if xx>180:
            xx-=360
        if reg_atl.contains(Point(xx,tc_lat[storm,ti])):
            plot = True

    if plot:
        ax.plot(tc_lon[storm,0],tc_lat[storm,0],color=cat_colors[TC_support.tc_cat(np.nanmin(tc_sel['source_pres'][storm,:,:]),'pressure')],alpha=0.7,marker='o',transform=ccrs.PlateCarree())

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_genesis_houston.png',dpi=300)


# sparial stuff
plt.close('all')
plt.figure(figsize=(10,5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()
gl=ax.gridlines(color='#BF3EFF',linewidth=1)
gl.ylocator = mticker.FixedLocator(np.arange(-20,60,10))
gl.xlocator = mticker.FixedLocator(np.arange(-120,60,10))
for yy in np.arange(-10,40,10):   ax.text(5,yy,str(yy),color='#BF3EFF',transform=ccrs.PlateCarree())
for xx in np.arange(-100,0,10):   ax.text(xx,-8,str(xx),color='#BF3EFF',transform=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
ax.set_xlim(-120,50)
ax.set_ylim(-15,50)

#########
# North Carolina
#########

reg_atl=Polygon([(-79,30),(-79,40),(-75,40),(-75,30),(-79,30)])
ax.add_geometries([reg_atl], ccrs.PlateCarree(), color='red',alpha=1,facecolor='none')

ax.lines=[]
for storm in tc_sel.storm:
    plot = False
    for ti in range(tc_lon[storm,:].shape[0]):
        xx = tc_lon[storm,ti]
        if xx>180:
            xx-=360
        if reg_atl.contains(Point(xx,tc_lat[storm,ti])):
            plot = True

    if plot:
        ax.plot(tc_lon[storm,:],tc_lat[storm,:],color=cat_colors[TC_support.tc_cat(np.nanmin(tc_sel['source_pres'][storm,:,:]),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_tracks_ncarolina.png',dpi=300)

ax.lines=[]
for storm in tc_sel.storm:
    plot = False
    for ti in range(tc_lon[storm,:].shape[0]):
        xx = tc_lon[storm,ti]
        if xx>180:
            xx-=360
        if reg_atl.contains(Point(xx,tc_lat[storm,ti])):
            plot = True

    if plot:
        ax.plot(tc_lon[storm,0],tc_lat[storm,0],color=cat_colors[TC_support.tc_cat(np.nanmin(tc_sel['source_pres'][storm,:,:]),'pressure')],alpha=0.7,marker='o',transform=ccrs.PlateCarree())

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_genesis_ncarolina.png',dpi=300)
