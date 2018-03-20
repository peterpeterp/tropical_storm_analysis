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

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAR25/'
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')

import TC_support ;  TC_support = reload(TC_support)
import tc_detection;    tc_detection = reload(tc_detection)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--verbosity",'-v', help="increase output verbosity",action="store_true")
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
parser.add_argument('--portion','-p',help='tenth of the available files to treat',required=False)
parser.add_argument('--surrounding','-s',help='time_steps for which the surroundings are plotted',nargs='+',required=False, type=int)
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False

identifiers=[ff.split('_')[-3] for ff in glob.glob(data_path+'/item3225_daily_mean/item3225_daily*')]
portion=int(len(identifiers)/10)

if args.portion is not None:
    if (int(sys.argv[1])+1)*portion>=len(identifiers):
        identifiers=identifiers[int(sys.argv[1])*portion:len(identifiers)-1]
    else:
        identifiers=identifiers[int(sys.argv[1])*portion:(int(sys.argv[1])+1)*portion]


# identifiers=['p014']

print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')
    MSLP=da.read_nc(data_path+'item16222_6hrly_inst/item16222_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16222_6hrly_inst'].values[:,0,1:,:]
    nc=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')
    U=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15201_6hrly_inst'].ix[:,0,:,:]
    V=da.read_nc(data_path+'item15202_6hrly_inst/item15202_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15202_6hrly_inst'].ix[:,0,:,:]
    VO=rel_vort(U.values[:,:,:],V.values[:,:,:],U.latitude0,U.longitude0)
    Wind10=np.array(np.sqrt(U**2+V**2))
    T=da.read_nc(data_path+'item16203_6hrly_inst/item16203_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16203_6hrly_inst'].ix[:,1:3,:,:].values.mean(axis=1)
    time_=nc.time0
    dates=[num2date(t,units = nc.axes['time0'].units,calendar = nc.axes['time0'].calendar) for t in time_]

    # prepare map
    lats = nc['global_latitude0'].values
    lons = nc['global_longitude0'].values
    grid_lats = nc['latitude0'].values
    grid_lons = nc['longitude0'].values
    lons[lons>180]-=360
    o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
    o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
    lon_0 = normalize180(o_lon_p-180.)

    rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
    plate_carree = ccrs.PlateCarree()
    globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)

    plt.close('all')
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=plate_carree)
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
    ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))

    working_dir='detection/CAR25/'+str(identifier)+'_CAR25/'
    elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
    found_tcs=tc_detection.tc_tracks(Wind10=Wind10,MSLP=MSLP,MSLP_smoothed=None,SST=None,VO=VO,T=T,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
    found_tcs.init_map(ax=ax,transform=rot_pole)
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
