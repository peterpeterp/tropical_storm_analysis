from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
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

cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}

obs_tracks=da.read_nc('aew_detection/ERAint/ERAint_all_tracks_'+'belanger'+'.nc')['all_tracks']
nc=da.read_nc('easterly_waves/ano_curv_vort_850-700-600_2008.nc')
lons,lats=nc.lon,nc.lat
lons[lons>180]-=360

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
    ax.plot(lons[np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int)],color='orange s',alpha=0.7,linewidth=2,transform=plate_carree)

plt.title('ERAint AEW')
plt.tight_layout()
plt.savefig('plots/ERAint/ERAint_AEW_tracks.png',dpi=300)

ax.lines=[]
for storm in obs_tracks.unnamed:
    track=obs_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    ax.plot(lons[int(track.ix[0,2])],lats[int(track.ix[0,1])],color='orange',alpha=0.7,marker='o',transform=plate_carree)

plt.title('ERAint EAW genesis')
plt.tight_layout()
plt.savefig('plots/ERAint_AEW_genesis.png',dpi=300)








#
