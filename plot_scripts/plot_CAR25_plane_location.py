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
sys.path.append('/Users/peterpfleiderer/Projects/tropical_cyclones/tc_detection')
import TC_support ;  TC_support = reload(TC_support)


nc=da.read_nc('data/CAR25/item16222_6hrly_inst/item16222_6hrly_inst_p014_2017-06_2017-10.nc')
grid_lats = nc['latitude1'].values
grid_lons = nc['longitude1'].values
o_lon_p = nc['rotated_pole1'].attrs['grid_north_pole_longitude']
o_lat_p = nc['rotated_pole1'].attrs['grid_north_pole_latitude']
lon_0 = TC_support.normalize180(o_lon_p-180.)

rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)


lats = nc['global_latitude1'].values
lons = nc['global_longitude1'].values
lons[lons>180]-=360

obs_tracks=da.read_nc('detection/CAR25/CAR25_all_tracks_knutson2007.nc')['all_tracks']
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
reg=Polygon([(grid_lons[0],grid_lats[0]),(grid_lons[-1],grid_lats[0]),(grid_lons[-1],grid_lats[-1]),(grid_lons[0],grid_lats[-1]),(grid_lons[0],grid_lats[0])])
ax.add_geometries([reg], rot_pole, color='lightgreen',alpha=1,facecolor='none')

############
# Houston
############

reg_atl=Polygon([(-100,25),(-90,25),(-90,35),(-100,35),(-100,25)])
ax.add_geometries([reg_atl], ccrs.PlateCarree(), color='red',alpha=1,facecolor='none')

for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    plot = False
    for ti,tt in enumerate(track[:,'t']):
        if reg_atl.contains(Point(lons[int(track[tt,'y']),int(track[tt,'x'])],lats[int(track[tt,'y']),int(track[tt,'x'])])):
            plot = True

    if plot:
        ax.plot(lons[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_tracks_houston.png',dpi=300)



ax.lines=[]
for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    plot = False
    for ti,tt in enumerate(track[:,'t']):
        if reg_atl.contains(Point(lons[int(track[tt,'y']),int(track[tt,'x'])],lats[int(track[tt,'y']),int(track[tt,'x'])])):
            plot = True

    if plot:
        ax.plot(lons[int(track.ix[0,1]),int(track.ix[0,2])],lats[int(track.ix[0,1]),int(track.ix[0,2])],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,marker='o',transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_genesis_houston.png',dpi=300)


############
# North Carolina
############
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
reg=Polygon([(grid_lons[0],grid_lats[0]),(grid_lons[-1],grid_lats[0]),(grid_lons[-1],grid_lats[-1]),(grid_lons[0],grid_lats[-1]),(grid_lons[0],grid_lats[0])])
ax.add_geometries([reg], rot_pole, color='lightgreen',alpha=1,facecolor='none')

reg_atl=Polygon([(-79,30),(-79,40),(-75,40),(-75,30),(-79,30)])
ax.lines=[]
ax.add_geometries([reg_atl], ccrs.PlateCarree(), color='red',alpha=1,facecolor='none')
for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    plot = False
    for ti,tt in enumerate(track[:,'t']):
        if reg_atl.contains(Point(lons[int(track[tt,'y']),int(track[tt,'x'])],lats[int(track[tt,'y']),int(track[tt,'x'])])):
            plot = True

    if plot:
        ax.plot(lons[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_tracks_ncarolina.png',dpi=300)



ax.lines=[]
for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    plot = False
    for ti,tt in enumerate(track[:,'t']):
        if reg_atl.contains(Point(lons[int(track[tt,'y']),int(track[tt,'x'])],lats[int(track[tt,'y']),int(track[tt,'x'])])):
            plot = True

    if plot:
        ax.plot(lons[int(track.ix[0,1]),int(track.ix[0,2])],lats[int(track.ix[0,1]),int(track.ix[0,2])],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,marker='o',transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_genesis_ncarolina.png',dpi=300)



#


#
