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
import matplotlib.ticker as mticker
import cartopy

import cartopy.crs as ccrs

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
import TC_support ;  TC_support = reload(TC_support)


nc=da.read_nc('data/CAR25/item16222_6hrly_inst/item16222_6hrly_inst_p014_2017-06_2017-10.nc')
grid_lats = nc['latitude1'].values
grid_lons = nc['longitude1'].values
o_lon_p = nc['rotated_pole1'].attrs['grid_north_pole_longitude']
o_lat_p = nc['rotated_pole1'].attrs['grid_north_pole_latitude']
lon_0 = TC_support.normalize180(o_lon_p-180.)

rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
plate_carree = ccrs.PlateCarree()
globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)


lats = nc['global_latitude1'].values
lons = nc['global_longitude1'].values
lons[lons>180]-=360

obs_tracks=da.read_nc('detection/CAR25/CAR25_all_tracks_knutson2007.nc')['all_tracks']
cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}


# sparial stuff
plt.close('all')
plt.figure(figsize=(10,5))
ax = plt.axes(projection=rot_pole)
ax.set_global()
ax.coastlines()
gl=ax.gridlines(color='#8470ff',linewidth=1)
gl.ylocator = mticker.FixedLocator(np.arange(-20,60,10))
gl.xlocator = mticker.FixedLocator(np.arange(-120,20,10))
for yy in np.arange(-10,40,10):   ax.text(-15,yy,str(yy),color='#8470ff',transform=plate_carree)
for xx in np.arange(-100,0,10):   ax.text(xx,8,str(xx),color='#8470ff',transform=plate_carree)
ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
ax.set_xlim(np.min(grid_lons)-5,np.max(grid_lons)+20)
ax.set_ylim(np.min(grid_lats)-5,np.max(grid_lats)+5)
reg=Polygon([(grid_lons[0],grid_lats[0]),(grid_lons[-1],grid_lats[0]),(grid_lons[-1],grid_lats[-1]),(grid_lons[0],grid_lats[-1]),(grid_lons[0],grid_lats[0])])
ax.add_geometries([reg], rot_pole, color='lightgreen',alpha=1,facecolor='none')

for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    ax.plot(lons[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int),np.array(track[:,'x'],int)],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,linewidth=2,transform=plate_carree)

plt.title('CAR25')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_tracks.png',dpi=300)

ax.lines=[]
for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    ax.plot(lons[int(track.ix[0,1]),int(track.ix[0,2])],lats[int(track.ix[0,1]),int(track.ix[0,2])],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,marker='o',transform=plate_carree)

plt.title('CAR25')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_genesis.png',dpi=300)


# statistics
identifiers=list(set([id_.split('_')[0] for id_ in obs_tracks.unnamed]))
hurrs_in_seas=[]
tcs_in_season=[]
for year in identifiers:
    track_ids=[id_ for id_ in obs_tracks.unnamed if id_.split('_')[0]==year]
    clim_cat=np.array(TC_support.tc_cat(np.nanmin(obs_tracks[track_ids,:,'MSLP'],axis=1),'pressure'))
    hurrs_in_seas.append(np.sum(clim_cat>0))
    tcs_in_season.append(np.sum(clim_cat>-1))

plt.close('all')
plt.figure(figsize=(5,3))
clim_cat=np.array(TC_support.tc_cat(np.nanmin(obs_tracks[:,:,'MSLP'],axis=1),'pressure'))
for cat in range(1,6):
    plt.bar(cat,np.sum(clim_cat==cat)/float(len(identifiers)),color=cat_colors[cat])
plt.annotate('total: '+str(int(round(np.mean(hurrs_in_seas))))+' $\pm$ '+str(int(round(np.std(hurrs_in_seas)))),xy=(0.75,0.85),xycoords='axes fraction')
plt.ylabel('hurricanes per season')
plt.xlabel('hurricane category')
plt.title('CAR25')
plt.tight_layout()
plt.savefig('plots/CAR25/CAR25_cat_hist.png',dpi=300)





#
