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
mpl.use('Agg')

sns.set_palette(sns.color_palette("plasma"))

try:
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='easterly_waves/'
    local=True
    years=[2008]
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/ERA_6hourly/u_v/'
    local=False
    years=range(1979,2017)

nc = da.read_nc('aew_detection/ERAint/2008/info.nc')
lat=nc['lats'].values
lon=nc['lons'].values
lat=lat[:,0]; lat=lat[np.isfinite(lat)]
lon=lon[0,:]; lon=lon[np.isfinite(lon)]
lons,lats=np.meshgrid(lon,lat)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes=plt.subplots(nrows=4,ncols=1,figsize=(8,6),subplot_kw={'projection': plate_carree},gridspec_kw = {'height_ratios':[3,1,3,1]})
for ax in axes[[0,2]]:
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lon),np.max(lon))
    ax.set_ylim(np.min(lat),np.max(lat))

occurence=lats.copy()*0.0
genesis=lats.copy()*0.0
for style in ['belanger']:
    if os.path.isfile('aew_detection/ERAint/ERAint_all_tracks_AEW_'+style+'.nc')==False or overwrite:
        # check for duplicates
        for identifier in [str(yr) for yr in years]:
            tmp=da.read_nc('aew_detection/ERAint/'+str(identifier)+'/track_info_'+style+'.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    track=track[np.isfinite(track[:,'t']),:]
                    genesis[int(track.ix[0,1]),int(track.ix[0,2])]+=1
                    for point in track.values.tolist():
                        occurence[int(point[1]),int(point[2])]+=1

occurence/=float(len(years))
cs=axes[0].contourf(lons,lats,occurence)
axes[1].outline_patch.set_edgecolor('white')
cbar_ax=fig.add_axes([0,0.6,1,0.15])
cbar_ax.axis('off')
cb=fig.colorbar(cs,orientation='horizontal',label='seasonal occurence frequency',ax=cbar_ax)

genesis/=float(len(years))
cs=axes[2].contourf(lons,lats,genesis)
axes[3].outline_patch.set_edgecolor('white')
cbar_ax=fig.add_axes([0,0.1,1,0.15])
cbar_ax.axis('off')
cb=fig.colorbar(cs,orientation='horizontal',label='seasonal genesis frequency',ax=cbar_ax)


plt.tight_layout()
plt.savefig('plots/ERAint/ERAint_AEW_tracks.png',dpi=300)
