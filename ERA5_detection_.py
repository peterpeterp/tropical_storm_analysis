from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max
import cartopy.crs as ccrs
import cartopy
sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAM25/'
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')

from TC_support import * ; reload(sys.modules['TC_support'])
from tc_detection import * ; reload(sys.modules['tc_detection'])

def date_(t):
    yr=int(t/10000.)
    mn=int((t-yr*10000)/100.)
    day=int((t-yr*10000-mn*100))
    hr=int((t-yr*10000-mn*100-day)*24)
    return(datetime(yr,mn,day,hr,0))

# start = time.time()
# identifier='2016'
# print('*** started run '+identifier+' ***')
# os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
# data_path='data/ERA5/'
# wind_nc=da.read_nc(data_path+'atl_'+identifier+'_10mWind.nc')
# Wind10=wind_nc['ws'].values
# time_=wind_nc.time
# MSLP=da.read_nc(data_path+'atl_'+identifier+'_surface.nc')['MSL'].values/100.
# T=da.read_nc(data_path+'atl_'+identifier+'_T.nc')['T'].values
# SST=da.read_nc(data_path+'atl_'+identifier+'_surface.nc')['SSTK'].values-273.15
# dates=[date_(t) for t in time_]
#
# U=da.read_nc(data_path+'atl_'+identifier+'_UV.nc')['U']
# V=da.read_nc(data_path+'atl_'+identifier+'_UV.nc')['V']
# VO=rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.lat,U.lon)
#
# lons,lats=np.meshgrid(wind_nc.lon,wind_nc.lat)
# lons[lons>180]-=360
#
# TC_season=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['season']
# TC_basin=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['basin']
# tc_sel=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc').ix[np.where((TC_season==int(identifier)) & (TC_basin[:,0]==0))[0]]
#
elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
ax = plt.axes(projection=plate_carree)
ax.set_global()
ax.coastlines()
ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
ax.set_extent([np.min(lons),np.max(lons),np.min(lats),np.max(lats)],crs=plate_carree)

working_dir='detection/ERA5/'+str(identifier)+'_ERA5_/'
found_tcs=tc_tracks(Wind10=Wind10,MSLP=MSLP,MSLP_smoothed=ndimage.gaussian_filter(MSLP,sigma=(0,3,3)),SST=SST,VO=VO,T=T[:,1,:,:],lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
found_tcs.init_map(ax=ax,transform=plate_carree)
found_tcs.init_obs_tcs(tc_sel)
elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)


# thresholds method
found_tcs.detect_knutson2007(overwrite=True)
found_tcs.plot_detect_summary(thr_wind=0)
found_tcs.combine_tracks(overwrite=True,thr_wind=0,search_radius=6,total_steps=5,warm_steps=3,consecutive_warm_strong_steps=0,plot=False)
found_tcs.plot_season()

# contours method
found_tcs.detect_contours(overwrite=False,p_radius=27,neighborhood_size=3,warm_core_size=3,cores_distance=1)
found_tcs.plot_detect_summary(thr_wind=17.5)
found_tcs.combine_tracks(overwrite=True,thr_wind=17.5,search_radius=6,total_steps=12,warm_steps=8,consecutive_warm_strong_steps=6)
found_tcs.plot_season()
elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)

# plt.close('all')
# fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,5),subplot_kw={'projection': plate_carree})
# axes=axes.flatten()
# for ax in axes:
#     ax.set_global()
#     ax.coastlines(edgecolor='magenta')
#     #ax.stock_img()
#     ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
#     ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
#     ax.set_extent([np.min(lons),np.max(lons),np.min(lats),np.max(lats)],crs=plate_carree)
#
# found_tcs.plot_surrounding(axes=axes,time_steps=range(600,640))#; convert -delay 50 track_surrounding/{94..127}* TC.gif
