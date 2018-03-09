from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage

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

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--verbosity",'-v', help="increase output verbosity",action="store_true")
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
parser.add_argument('--portion','-p',help='tenth of the available files to treat',required=False)
parser.add_argument('--surrounding', ,help='time_steps for which the surroundings are plotted',nargs='+',required=False, type=int)
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

print(identifiers)

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')
    MSLP=da.read_nc(data_path+'item16222_daily_mean/item16222_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item16222_daily_mean'].values[:,0,1:,:]
    nc=da.read_nc(data_path+'item3225_daily_mean/item3225_daily_mean_'+identifier+'_2017-06_2017-10.nc')
    U=da.read_nc(data_path+'item3225_daily_mean/item3225_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item3225_daily_mean'].ix[:,0,:,:]
    V=da.read_nc(data_path+'item3226_daily_mean/item3226_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item3226_daily_mean'].ix[:,0,:,:]
    VO=rel_vort(U.values[:,:,:],V.values[:,:,:],U.latitude0,U.longitude0)
    Wind10=np.array(np.sqrt(U**2+V**2))

    time_=nc.time0
    dates=[num2date(t,units = nc.axes['time0'].units,calendar = nc.axes['time0'].calendar) for t in time_]

    # prepare map
    lats = nc['global_latitude0'].values
    lons = nc['global_longitude0'].values

    o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
    o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
    lon_0 = normalize180(o_lon_p-180.)

    # init map that is used for plotting
    plt.close('all')
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))
    m= Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
               llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
               llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
    m.drawmapboundary(fill_color='1.')
    m.drawmapboundary(fill_color='darkblue')
    m.fillcontinents(color='darkgreen',lake_color='darkblue')
    m.drawcoastlines(linewidth=0.3)
    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
    ax.invert_yaxis()
    plot_lon,plot_lat = m(lons,lats)

    working_dir='detection/CAM25/'+str(identifier)+'_CAM25/'
    elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
    found_tcs=tc_tracks(Wind10=Wind10,MSLP=MSLP,SST=None,VO=VO,T=None,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
    found_tcs.init_map(m=m,ax=ax,plot_lat=plot_lat,plot_lon=plot_lon)

    found_tcs.set_thresholds(thr_wind=15,thr_vort=5*10**(-5),thr_mslp=101500,thr_ta=0,thr_sst=26.5,win1=7,win2=12,win_step=20,neighborhood_size=8,min_time_steps=2)
    found_tcs.detect(overwrite=overwrite)
    found_tcs.combine_tracks(overwrite=overwrite)
    #found_tcs.gather_info_track(overwrite=False)
    #track_info,track=found_tcs.plot_track_evolution()
    #found_tcs.plot_season()

    if args.surrounding is not None:
        plt.close('all')
        fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,5))
        axes=axes.flatten()
        maps=[]
        for ax in axes:
            mm=Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
                       llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
                       llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
            mm.drawcoastlines(linewidth=0.7,color='m')
            mm.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
            mm.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
            ax.invert_yaxis()
            maps.append(mm)
        found_tcs.plot_surrounding(maps=maps,axes=axes,time_steps=args.surrounding)#; convert -delay 50 track_surrounding/{94..127}* TC.gif
        elapsed = time.time() - start;  print('Done with plotting %.3f seconds.' % elapsed)
    print('memory in use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10.**6))
