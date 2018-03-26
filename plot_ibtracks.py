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

# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
# tc_sel=TC.ix[np.where(TC['basin'][:,0]==0)[0]]
# tc_sel=tc_sel.ix[tc_sel['season']>=1979]
# tc_lat=tc_sel['lat_for_mapping']
# tc_lon=tc_sel['lon_for_mapping']
# tc_lon[tc_lon<0]+=360
# tc_wind=tc_sel['source_wind']
# tc_wind[np.isnan(tc_wind)]=-999
#
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
plate_carree = ccrs.PlateCarree()
globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)

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

for storm in tc_sel.storm:
    ax.plot(tc_lon[storm,:],tc_lat[storm,:],color=cat_colors[TC_support.tc_cat(tc_wind[storm,:,0].max(),'wind')],alpha=0.3,linewidth=2,transform=plate_carree)

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_tracks.png',dpi=300)

ax.lines=[]
for storm in tc_sel.storm:
    ax.plot(tc_lon[storm,0],tc_lat[storm,0],color=cat_colors[TC_support.tc_cat(tc_wind[storm,:,0].max(),'wind')],alpha=0.3,marker='o',transform=plate_carree)

plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_genesis.png',dpi=300)


# statistics
hurrs_in_seas=[]
for year in range(1979,2018):
    tc_clim=tc_sel.ix[tc_sel['season']==year]
    print(tc_clim['source_pres'].values)
    clim_cat=np.array(TC_support.tc_cat(np.nanmin(tc_clim['source_pres'],axis=(1,2)),'pressure'))
    hurrs_in_seas.append(np.sum(clim_cat>0))

tc_clim=tc_sel.ix[tc_sel['season']>=1979]
plt.close('all')
plt.figure(figsize=(5,3))
clim_cat=np.array(TC_support.tc_cat(np.nanmin(tc_clim['source_pres'],axis=(1,2)),'pressure'))
for cat in range(1,6):
    plt.bar(cat,np.sum(clim_cat==cat)/float(2017-1979),color=cat_colors[cat])
plt.text(4.2,2.5,'total: '+str(int(round(np.mean(hurrs_in_seas))))+' $\pm$ '+str(int(round(np.std(hurrs_in_seas)))))
plt.ylabel('hurricanes per season')
plt.xlabel('hurricane category')
plt.title('ibtracks')
plt.tight_layout()
plt.savefig('plots/ibtracks/ibtracks_cat_hist.png',dpi=300)





#
