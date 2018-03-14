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
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/JRA55/'
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/JRA55/'

import TC_support ;  TC_support = reload(TC_support)
import tc_detection;    tc_detection = reload(tc_detection)
# import tc_detection;    if is_changed(tc_detection): foo = reload(tc_detection)
# import TC_support;    if is_changed(TC_support): foo = reload(TC_support)



import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--verbosity",'-v', help="increase output verbosity",action="store_true")
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
parser.add_argument('--year','-y',help='year',required=False)
parser.add_argument('--surrounding','-s',help='time_steps for which the surroundings are plotted',nargs='+',required=False, type=int)
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False

if args.year is not None:
    identifiers=[args.year]
else:
    identifiers=[str(yr) for yr in range(1979,2018)]
    identifiers=['2010']

print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')

    U10=da.read_nc(data_path+'atl_'+identifier+'_U10.nc')['var33'].values.squeeze()
    V10=da.read_nc(data_path+'atl_'+identifier+'_V10.nc')['var34'].values.squeeze()
    T=da.read_nc(data_path+'atl_'+identifier+'_T_ana.nc')['var11'].values[:,1:3,:,:].mean(axis=1)

    Wind10=(U10**2+V10**2)**0.5
    nc=da.read_nc(data_path+'atl_'+identifier+'_MSLP.nc')
    MSLP=nc['psl'].values/100.
    time_=nc.time
    dates=[num2date(t,units = nc.axes['time'].units,calendar = nc.axes['time'].calendar) for t in time_]

    U=da.read_nc(data_path+'atl_'+identifier+'_U.nc')['var33']
    V=da.read_nc(data_path+'atl_'+identifier+'_V.nc')['var34']
    VO=TC_support.rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.lat,U.lon)
    VO[np.isnan(VO)]=-999

    lons,lats=np.meshgrid(nc.lon,nc.lat)
    lons[lons>180]-=360

    TC_season=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['season']
    TC_basin=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['basin']
    tc_sel=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc').ix[np.where((TC_season==int(identifier)) & (TC_basin[:,0]==0))[0]]

    plt.close('all')
    plate_carree = ccrs.PlateCarree()
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    ax = plt.axes(projection=plate_carree)
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
    ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))

    working_dir='detection/JRA55/'+str(identifier)+'_JRA55/'
    found_tcs=tc_detection.tc_tracks(Wind10=Wind10,MSLP=MSLP,MSLP_smoothed=ndimage.gaussian_filter(MSLP,sigma=(0,2,2)),SST=None,VO=VO,T=T,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
    found_tcs.init_map(ax=ax,transform=plate_carree)
    found_tcs.init_obs_tcs(tc_sel)
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)

    # contours method
    found_tcs.detect_contours(overwrite=True,p_radius=27,dis_mslp_min=3,warm_core_size=3,dis_cores=1)
    found_tcs.plot_detect_summary(thr_wind=10)
    found_tcs.combine_tracks(overwrite=True,thr_wind=17.5,search_radius=6,total_steps=12,warm_steps=8,consecutive_warm_strong_steps=4,plot=False)
    found_tcs.plot_season()
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)

    # thresholds method
    found_tcs.detect_knutson2007(overwrite=True,thr_vort=3.5*10**(-5),dis_vort_max=4,dis_cores=2,thr_MSLP_inc=2,dis_MSLP_inc=5,thr_T_drop=0.8,dis_T_drop=5,tc_size=7)
    found_tcs.plot_detect_summary(thr_wind=15)
    found_tcs.combine_tracks(overwrite=True,thr_wind=15,search_radius=6,total_steps=8,strong_steps=8,warm_steps=8,consecutive_warm_strong_steps=0,plot=False)
    found_tcs.plot_season()

    found_tcs.obs_track_info(overwrite=True,core_radius=3,full_radius=7)

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
