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
    os.chdir('/Users/peterpfleiderer/Projects/')
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/tropical_cyclones/TC_detection/TC_scripts')
    sys.path.append('/p/projects/tumble/carls/shared_folder/tropical_cyclones/TC_detection/TC_scripts/feature_tracking')
    os.chdir('/p/projects/tumble/carls/shared_folder/')

import TC_support ;  TC_support = reload(TC_support)
import TC_detection;    TC_detection = reload(TC_detection)

identifiers = [nn.split('_')[-3] for nn in glob.glob('data/WAH/batch_755/region/item16222_6hrly_inst/*')]
print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')

    data_path = 'data/WAH/batch_755/region/'

    nc=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')
    time_=nc.time0
    dates=num2date(time_,units = nc['time0'].units)

    MSLP=da.read_nc(data_path+'item16222_6hrly_inst/item16222_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16222_6hrly_inst'].squeeze().values[:,:-1,:]/100.
    T=da.read_nc(data_path+'item16203_6hrly_inst/item16203_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16203_6hrly_inst'].values[:,1:3,:-1,:].mean(axis=1).squeeze()

    U=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15201_6hrly_inst']
    V=da.read_nc(data_path+'item15202_6hrly_inst/item15202_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15202_6hrly_inst']
    VO_850=TC_support.rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.latitude0,U.longitude0)
    VO_500=TC_support.rel_vort(U.values[:,1,:,:],V.values[:,1,:,:],U.latitude0,U.longitude0)
    VO_av=(VO_850.squeeze() + VO_500.squeeze()) * 0.5

    Wind10=np.array(np.sqrt(U.values[:,0,:,:]**2+V.values[:,0,:,:]**2))

    #land_mask=np.array(da.read_nc('data/MEX/land_mask.nc')['lsm'].values.squeeze()[:-1]*-1,np.int)

    # prepare map
    lats = nc['latitude0'].values
    lons = nc['longitude0'].values -25
    lons[lons>180]-=360
    lons,lats = np.meshgrid(lons,lats)
    plate_carree = ccrs.PlateCarree()

    working_dir='tropical_cyclones/detection/ATL/'+str(identifier)+'/'

    found_tcs=TC_detection.tc_tracks(working_dir=working_dir,
                                 lats=lats,
                                 lons=lons,
                                 time_=time_,
                                 dates=dates,
                                 identifier=identifier)
    found_tcs.plot_settings(ccrs.PlateCarree(),
                            extent=[np.min(lons),np.max(lons),np.min(lats),np.max(lats)])



    input_data={'VO':ndimage.gaussian_filter(VO_av,sigma=(0,4,4)),'MSLP':MSLP,'MSLP_smoothed':ndimage.gaussian_filter(MSLP,sigma=(0,3,3)),'T':T,'Wind10':Wind10}
    found_tcs.add_data(input_data)
    # contours method
    detected=found_tcs.detect_contours(overwrite=True,p_radius=27,dis_mslp_min=3,warm_core_size=3,dis_cores=1)
    ax = found_tcs.plot_detect_summary(thr_wind=3)


    found_tcs.combine_tracks(overwrite=True,
                             thr_wind=5,search_radius=6,
                             total_steps=12,warm_steps=8,consecutive_warm_strong_steps=4,plot=False)
    ax=found_tcs.plot_all_tracks()

    found_tcs.plot_surrounding()

    #
