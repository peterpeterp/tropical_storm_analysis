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

sns.set_palette(sns.color_palette("plasma"))

try:
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    local=True
    years=[2008]
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/tropical_cyclones/tropical_storm_analysis')
    os.chdir('/p/projects/tumble/carls/shared_folder/tropical_cyclones/')
    local=False
    years=range(1979,2017)
    data_path='/p/projects/tumble/carls/shared_folder/data/reanalysis/ERAint/'
    plt.switch_backend('agg')

yrFr = np.array([xx/365./4. +152./365.  for xx in np.arange(0,732,1)])
day = np.array([xx/4 for xx in np.arange(0,732,1)])
month = day.copy() * 0
month[day<30] = 6
month[(day>=30) & (day<61)] = 7
month[(day>=61) & (day<92)] = 8
month[(day>=92) & (day<122)] = 9
month[(day>=122) & (day<153)] = 10
month[day>=153] = 11

cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue","green","yellow","red","brown","violet","white"])

for style in ['belanger']: #'dieng'
    print(style)
    nc = da.read_nc('aew_detection/ERAint/2008/2008_'+style+'_grid.nc')
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
        ax.set_xlim(-70,30)
        ax.set_ylim(0,30)

    occurence=lats.copy()*0.0
    genesis=lats.copy()*0.0

    _time = np.zeros([100,1]) * np.nan
    _month = np.zeros([100,1]) * np.nan
    _lon = np.zeros([100,1]) * np.nan
    _lat = np.zeros([100,1]) * np.nan
    _dummy = np.zeros([100,1]) * np.nan
    _ids = []

    if os.path.isfile('aew_detection/ERAint/ERAint_all_tracks_AEW_'+style+'.nc')==False or True:
        # check for duplicates
        for identifier in [str(yr) for yr in years]:
            tmp=da.read_nc('aew_detection/ERAint/'+str(identifier)+'/'+str(identifier)+'_'+style+'_track_info.nc')
            for id_,track in tmp.items():
                if id_ not in ['z','time']:
                    _ids.append(id_)
                    track=track[np.isfinite(track[:,'t']),:]
                    genesis[int(track.ix[0,1]),int(track.ix[0,2])]+=1
                    for point in track.values.tolist():
                        occurence[int(point[1]),int(point[2])]+=1

                    _lon = np.vstack((_lon,_dummy))
                    _lon[:track.shape[0],-1] = lons[np.array(track.ix[:,1],np.int),np.array(track.ix[:,2],np.int)]
                    _lat = np.vstack((_lat,_dummy))
                    _lat[:track.shape[0],-1] = lats[np.array(track.ix[:,1],np.int),np.array(track.ix[:,2],np.int)]

                    _time = np.vstack((_time,_dummy))
                    _time[:track.shape[0],-1] = yrFr[np.array(track[:,'t'],np.int)]+float(identifier)
                    _month = np.vstack((_month,_dummy))
                    _month[:track.shape[0],-1] = month[np.array(track[:,'t'],np.int)]

        ds = da.Dataset({
            'lon':da.DimArray(_lon,axes=[np.arange(0,100,1,np.short),np.array(_ids)],dims=['step','ID']),
            'lat':da.DimArray(_lat,axes=[np.arange(0,100,1,np.short),np.array(_ids)],dims=['step','ID']),
            'time':da.DimArray(_time,axes=[np.arange(0,100,1,np.short),np.array(_ids)],dims=['step','ID']),
            'month':da.DimArray(_month,axes=[np.arange(0,100,1,np.short),np.array(_ids)],dims=['step','ID']),
        })
        ds.write_nc('aew_detection/ERAint/ERAint_all_tracks_AEW_'+style+'_clean.nc')

    occurence/=float(len(years))
    cs=axes[0].contourf(lons,lats,occurence,20,cmap=cmap)
    axes[1].outline_patch.set_edgecolor('white')
    cbar_ax=fig.add_axes([0,0.6,1,0.15])
    cbar_ax.axis('off')
    cb=fig.colorbar(cs,orientation='horizontal',label='seasonal occurence frequency',ax=cbar_ax)

    genesis/=float(len(years))
    cs=axes[2].contourf(lons,lats,genesis,20,cmap=cmap)
    axes[3].outline_patch.set_edgecolor('white')
    cbar_ax=fig.add_axes([0,0.1,1,0.15])
    cbar_ax.axis('off')
    cb=fig.colorbar(cs,orientation='horizontal',label='seasonal genesis frequency',ax=cbar_ax)


    plt.tight_layout()
    plt.savefig('plots/ERAint/ERAint_AEW_tracks_'+style+'.png',dpi=300)
