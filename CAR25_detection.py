from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAR25/'
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_717/region/'
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

identifiers=[ff.split('_')[-3] for ff in glob.glob(data_path+'/item16222_6hrly_inst/item16222_6hrly_inst*')]
portion=int(len(identifiers)/10)

if args.portion is not None:
    if (int(sys.argv[1])+1)*portion>=len(identifiers):
        identifiers=identifiers[int(sys.argv[1])*portion:len(identifiers)-1]
    else:
        identifiers=identifiers[int(sys.argv[1])*portion:(int(sys.argv[1])+1)*portion]


# identifiers=['p014']
# overwrite=False

print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')
    MSLP=da.read_nc(data_path+'item16222_6hrly_inst/item16222_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16222_6hrly_inst'].values[:,0,1:,:]/100.
    nc=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')
    U=da.read_nc(data_path+'item15201_6hrly_inst/item15201_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15201_6hrly_inst'].ix[:,0,:,:]
    V=da.read_nc(data_path+'item15202_6hrly_inst/item15202_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item15202_6hrly_inst'].ix[:,0,:,:]
    VO=TC_support.rel_vort(U.values[:,:,:],V.values[:,:,:],U.latitude0,U.longitude0)
    Wind10=np.array(np.sqrt(U**2+V**2))
    T=da.read_nc(data_path+'item16203_6hrly_inst/item16203_6hrly_inst_'+identifier+'_2017-06_2017-10.nc')['item16203_6hrly_inst'].ix[:,1:3,1:,:].values.mean(axis=1)
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
    lon_0 = TC_support.normalize180(o_lon_p-180.)

    rot_pole = ccrs.RotatedPole(pole_longitude=o_lon_p, pole_latitude=o_lat_p)
    plate_carree = ccrs.PlateCarree()
    globe= ccrs.Orthographic(central_longitude=-60.0, central_latitude=20.0, globe=None)

    plt.close('all')
    plt.figure(figsize=(10,5))
    ax = plt.axes(projection=rot_pole)
    ax.set_global()
    ax.coastlines()
    gl=ax.gridlines(color='lightblue',linewidth=1)
    gl.ylocator = mticker.FixedLocator(np.arange(-10,60,10))
    gl.xlocator = mticker.FixedLocator(np.arange(-110,0,10))
    for yy in np.arange(0,40,10):   ax.text(-35,yy,str(yy),color='lightblue',transform=plate_carree)
    for xx in np.arange(-90,-20,10):   ax.text(xx,8,str(xx),color='lightblue',transform=plate_carree)
    ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
    ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
    ax.set_xlim(np.min(grid_lons),np.max(grid_lons))
    ax.set_ylim(np.min(grid_lats),np.max(grid_lats))

    # #ax.pcolormesh(lons,lats,MSLP[0,:,:],transform=plate_carree)
    # ax.plot(lons[50,50],lats[50,50],'og',transform=plate_carree)
    # plt.savefig('test.png')

    working_dir='detection/CAR25/'+str(identifier)+'_CAR25/'
    elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
    found_tcs=tc_detection.tc_tracks(Wind10=Wind10,MSLP=MSLP,MSLP_smoothed=None,SST=None,VO=VO,T=T,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
    found_tcs.init_map(ax=ax,transform=plate_carree)
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)

    # contours method
    found_tcs.detect_contours(overwrite=overwrite,p_radius=27,dis_mslp_min=3,warm_core_size=3,dis_cores=1)
    found_tcs.plot_detect_summary(thr_wind=15)
    found_tcs.combine_tracks(overwrite=overwrite,thr_wind=17.5,search_radius=6,total_steps=12,warm_steps=8,consecutive_warm_strong_steps=4,lat_formation_cutoff=30,plot=False)
    found_tcs.plot_season()
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)

    # thresholds method
    found_tcs.detect_knutson2007(overwrite=overwrite,thr_vort=3.5*10**(-5),dis_vort_max=4,dis_cores=2,thr_MSLP_inc=2,dis_MSLP_inc=5,thr_T_drop=0.8,dis_T_drop=5,tc_size=7)
    found_tcs.plot_detect_summary(thr_wind=15)
    found_tcs.combine_tracks(overwrite=overwrite,thr_wind=15,search_radius=6,total_steps=8,strong_steps=8,warm_steps=8,consecutive_warm_strong_steps=0,lat_formation_cutoff=30,plot=False)
    found_tcs.plot_season()


    # plt.close('all')
    # fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(10,5),subplot_kw={'projection': rot_pole})
    # axes=axes.flatten()
    # for ax in axes:
    #     ax.set_global()
    #     ax.coastlines(color='magenta')
    #     ax.set_xlim(np.min(grid_lons),np.max(grid_lons))
    #     ax.set_ylim(np.min(grid_lats),np.max(grid_lats))
    #
    # found_tcs.plot_surrounding(axes=axes,time_steps=range(500,540))#; convert -delay 50 track_surrounding/{94..127}* TC.gif
