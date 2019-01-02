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

for base_path in ['/Users/peterpfleiderer/Projects/tropical_cyclones','tropical_cyclones','/p/projects/tumble/carls/shared_folder/tropical_cyclones']:
    try:
        os.chdir(base_path)
    except:
        pass

sys.path.append('tropical_storm_analysis')
sys.path.append('tropical_storm_analysis/feature_tracking')

import TC_support ;  TC_support = reload(TC_support)
import TC_detection;    TC_detection = reload(TC_detection)

cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}

found_tracks=da.read_nc('detection/ATL/ATL_all_tracks_contours.nc')['all_tracks']
nc=da.read_nc('../data/WAH/batch_755/region/item16222_6hrly_inst/item16222_6hrly_inst_p034_2017-06_2017-10.nc')
lons,lats=nc.longitude1,nc.latitude1
lons -= 25
lons[lons>180]-=360

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
reg=Polygon([(max(lons),max(lats)),(min(lons),max(lats)),(min(lons),min(lats)),(max(lons),min(lats)),(max(lons),max(lats))])
ax.add_geometries([reg], ccrs.PlateCarree(), color='lightgreen',alpha=1,facecolor='none')

for storm in found_tracks.unnamed:
    track=found_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    ax.plot(lons[np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int)],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/ATL/ATL_tracks.png',dpi=300)

ax.lines=[]
for storm in found_tracks.unnamed:
    track=found_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    ax.plot(lons[int(track.ix[0,2])],lats[int(track.ix[0,1])],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,marker='o',transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/ATL/ATL_genesis.png',dpi=300)


# statistics
identifiers=list(set([id_.split('_')[0] for id_ in found_tracks.unnamed]))
hurrs_in_seas=[]
tcs_in_season=[]
for year in identifiers:
    track_ids=[id_ for id_ in found_tracks.unnamed if id_.split('_')[0]==year]
    clim_cat=np.array(TC_support.tc_cat(np.nanmin(found_tracks[track_ids,:,'MSLP'],axis=1),'pressure'))
    hurrs_in_seas.append(np.sum(clim_cat>0))
    tcs_in_season.append(np.sum(clim_cat>-1))

sns.reset_orig()

plt.close('all')
plt.figure(figsize=(5,3))
clim_cat=np.array(TC_support.tc_cat(np.nanmin(found_tracks[:,:,'MSLP'],axis=1),'pressure'))
for cat in range(1,6):
    plt.bar(cat,np.sum(clim_cat==cat)/float(len(identifiers)),color=cat_colors[cat],label=cat)
plt.annotate('total: '+str(int(round(np.mean(hurrs_in_seas))))+' $\pm$ '+str(int(round(np.std(hurrs_in_seas)))),xy=(0.75,0.85),xycoords='axes fraction')
plt.ylabel('hurricanes per season')
plt.xlabel('hurricane category')
plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/ATL/ATL_cat_hist.png',dpi=300)


# for one locations
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
reg=Polygon([(max(lons),max(lats)),(min(lons),max(lats)),(min(lons),min(lats)),(max(lons),min(lats)),(max(lons),max(lats))])
ax.add_geometries([reg], ccrs.PlateCarree(), color='lightgreen',alpha=1,facecolor='none')


reg_atl=Polygon([(-67,20),(-10,20),(-10,5),(-67,5),(-67,20)])
ax.add_geometries([reg_atl], ccrs.PlateCarree(), color='red',alpha=1,facecolor='none')

for storm in found_tracks.unnamed:
    track=found_tracks[storm,:,:]
    track=track[np.isfinite(track[:,'t']),:]
    plot = False
    for ti,tt in enumerate(track[:,'t']):
        print(reg_atl.contains(Point(track[tt,'x'],track[tt,'y'])))
        if reg_atl.contains(Point(track[tt,'x']-25,track[tt,'y'])):
            plot = True

    if plot:
        ax.plot(lons[np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int)],color=cat_colors[TC_support.tc_cat(track[:,'MSLP'].min(),'pressure')],alpha=0.7,linewidth=2,transform=ccrs.PlateCarree())

plt.title('ATL')
plt.tight_layout()
plt.savefig('plots/ATL/ATL_tracks_houston.png',dpi=300)


#
