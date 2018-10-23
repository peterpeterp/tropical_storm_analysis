from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max
import cartopy.crs as ccrs
import cartopy
sns.set_palette(sns.color_palette("plasma"))

try:
    sys.path.append('/Users/peterpfleiderer/Projects/tropical_cyclones/TC_scripts')
    sys.path.append('/Users/peterpfleiderer/Projects/tropical_cyclones/TC_scripts/feature_tracking')
    os.chdir('/Users/peterpfleiderer/Projects/tropical_cyclones/')
    data_path='data/MEX/'
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/TC_scripts')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')

import TC_support ;  TC_support = reload(TC_support)
import TC_detection;    TC_detection = reload(TC_detection)

identifiers = ['a0ks']
print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')

    data_path = 'data/MEX/'+identifier+'/'


    U10=da.read_nc(data_path+'item3225_daily_mean_'+identifier+'_2004-05_2005-10.nc')['item3225_daily_mean'].values.squeeze()
    V10=da.read_nc(data_path+'item3226_daily_mean_'+identifier+'_2004-05_2005-10.nc')['item3226_daily_mean'].values.squeeze()
    Wind10=(U10**2+V10**2)**0.5

    MSLP=da.read_nc(data_path+'item16222_daily_mean_'+identifier+'_2004-05_2005-10.nc')['item16222_daily_mean'].squeeze().values[:,:-1,:]/100.

    nc=da.read_nc(data_path+'item15201_daily_mean_'+identifier+'_2004-05_2005-10.nc')
    time_=nc.time0
    dates=num2date(time_,units = nc['time0'].units)

    U=da.read_nc(data_path+'item15201_daily_mean_'+identifier+'_2004-05_2005-10.nc')['item15201_daily_mean']
    V=da.read_nc(data_path+'item15202_daily_mean_'+identifier+'_2004-05_2005-10.nc')['item15202_daily_mean']
    VO_850=TC_support.rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.latitude0,U.longitude0)
    VO_500=TC_support.rel_vort(U.values[:,1,:,:],V.values[:,1,:,:],U.latitude0,U.longitude0)
    VO_av=(VO_850.squeeze() + VO_500.squeeze()) * 0.5

    land_mask=np.array(da.read_nc('data/MEX/land_mask.nc')['lsm'].values.squeeze()[:-1]*-1,np.int)

    # prepare map
    lats = nc['global_latitude0'].values
    lons = nc['global_longitude0'].values
    grid_lats = nc['latitude0'].values
    grid_lons = nc['longitude0'].values
    lons[lons>180]-=360
    o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
    o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
    lon_0 = TC_support.normalize180(o_lon_p-180.)

    rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
    plate_carree = ccrs.PlateCarree()

    working_dir='detection/MEX/'+str(identifier)+'_MEX/'



    import TC_detection;    TC_detection = reload(TC_detection)
    found_tcs=TC_detection.tc_tracks(working_dir=working_dir,
                                 lats=lats,
                                 lons=lons,
                                 time_=time_,
                                 dates=dates,
                                 identifier=identifier)
    found_tcs.plot_settings(ccrs.PlateCarree(),
                            extent=[np.min(lons),np.max(lons),np.min(lats),np.max(lats)])



    input_data={'VO':ndimage.gaussian_filter(VO_av,sigma=(0,4,4)),'MSLP':MSLP,'Wind10':Wind10,'land_mask':land_mask}
    found_tcs.add_data(input_data)
    found_tcs.detect_daily(overwrite=True)

    found_tcs.combine_tracks(overwrite=True,
                             thr_wind=15.0,search_radius=6,strong_steps=2,lat_formation_cutoff=30,
                             total_steps=2,warm_steps=0,consecutive_warm_strong_steps=0,plot=True)
    ax=found_tcs.plot_all_tracks()

    found_tcs.plot_surrounding()

    #
