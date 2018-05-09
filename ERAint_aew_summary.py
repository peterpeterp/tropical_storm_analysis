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
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/ERA_6hourly/u_v/'
    local=False

nc = da.read_nc('aew_detection/ERAint/2008/info.nc')
lat=nc['lats'].values
lon=nc['lons'].values
lat=lat[:,0]; lat=lat[np.isfinite(lat)]
lon=lon[0,:]; lon=lon[np.isfinite(lon)]
lons,lats=np.meshgrid(lon,lat)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes=plt.subplots(nrows=3,ncols=1,figsize=(10,10),subplot_kw={'projection': plate_carree})
for ax in axes:
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lon),np.max(lon))
    ax.set_ylim(np.min(lat),np.max(lat))

occurence=lats.copy()*0.0
for style in ['belanger']:
    if os.path.isfile('aew_detection/ERAint/ERAint_all_tracks_AEW_'+style+'.nc')==False or overwrite:
        # check for duplicates
        for identifier in [str(yr) for yr in range(1979,2017)]:
            try:
                tmp=da.read_nc('aew_detection/ERAint/'+str(identifier)+'/track_info_'+style+'.nc')
                for id_,track in tmp.items():
                    if id_ not in ['z','time']:
                        track=track[np.isfinite(track[:,'t']),:]
                        axes[1].plot(lon[np.array(track.ix[0,2],int)],lat[np.array(track.ix[0,1],int)],'.m',alpha=0.5,linestyle='',transform=plate_carree)
                        axes[2].plot(lon[np.array(track[:,'x'],int)],lat[np.array(track[:,'y'],int)],color='orange',alpha=0.5,linewidth=2,transform=plate_carree)
                        for point in track.values.tolist():
                            occurence[int(point[1]),int(point[2])]+=1
            except:
                pass
axes[0].contourf(lons,lats,occurence)
plt.tight_layout()
plt.savefig('plots/ERAint/ERAint_AEW_tracks.png',dpi=300)
