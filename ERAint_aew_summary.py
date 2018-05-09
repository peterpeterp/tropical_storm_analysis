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

nc=da.read_nc(data_path+'ano_curv_vort_850-700-600_2008.nc')
lons,lats=nc.lon,nc.lat
lons[lons>180]-=360

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes=plt.subplots(nrows=2,ncols=1,figsize=(10,10),subplot_kw={'projection': plate_carree})
for ax in axes:
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))


for style in ['belanger']:
    if os.path.isfile('aew_detection/ERAint/ERAint_all_tracks_AEW_'+style+'.nc')==False or overwrite:
        # check for duplicates
        for identifier in [str(yr) for yr in range(1979,2017)]:
            tmp=da.read_nc('aew_detection/ERAint/'+str(identifier)+'/track_info_'+style+'.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    track=track[np.isfinite(track[:,'t']),:]
                    axes[0].plot(lons[np.array(track[:,'x'],int)],lats[np.array(track[:,'y'],int)],color='orange',alpha=0.5,linewidth=2,transform=plate_carree)
                    axes[1].plot(lons[np.array(track[0,'x'],int)],lats[np.array(track[0,'y'],int)],color='orange',alpha=0.5,linestyle='',transform=plate_carree)


plt.tight_layout()
plt.savefig('plots/ERAint/ERAint_AEW_tracks.png',dpi=300)
