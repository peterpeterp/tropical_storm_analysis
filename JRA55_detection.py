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
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/JRA55/'
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='../reanalysis/JRA55/'

import TC_support; TC_support=reload(TC_support)
import tc_detection; tc_detection=reload(tc_detection)

try:
    identifiers=[sys.argv[1]]
except:
    identifiers=[ff.split('_')[1] for ff in glob.glob(data_path+'atl_*_MSLP.nc')]


for identifier in ['2010']:
    start = time.time()
    print('*** started run '+identifier+' ***')

    U10=da.read_nc(data_path+'atl_'+identifier+'_U10.nc')['var33'].values.squeeze()
    V10=da.read_nc(data_path+'atl_'+identifier+'_V10.nc')['var34'].values.squeeze()
    T=da.read_nc(data_path+'atl_'+identifier+'_T.nc')['var11'].values

    Wind10=(U10**2+V10**2)**0.5
    nc=da.read_nc(data_path+'atl_'+identifier+'_MSLP.nc')
    MSLP=nc['psl'].values
    time_=nc.time
    dates=[num2date(t,units = nc.axes['time'].units,calendar = nc.axes['time'].calendar) for t in time_]

    U=da.read_nc(data_path+'atl_'+identifier+'_U.nc')['var33']
    V=da.read_nc(data_path+'atl_'+identifier+'_V.nc')['var34']
    VO=TC_support.rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.lat,U.lon)

    lons,lats=np.meshgrid(nc.lon,nc.lat)
    lons[lons>180]-=360

    plt.close('all')
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    m = Basemap(ax=ax,llcrnrlon=np.min(lons),urcrnrlon=np.max(lons),llcrnrlat=np.min(lats),urcrnrlat=np.max(lats),resolution="l",projection='cyl')
    m.drawmapboundary(fill_color='1.')
    m.drawmapboundary(fill_color='darkblue')
    m.fillcontinents(color='darkgreen',lake_color='darkblue')
    m.drawcoastlines(linewidth=0.3)
    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)

    plot_lon,plot_lat=m(lons,lats)

    TC_season=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['season']
    TC_basin=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['basin']
    tc_sel=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc').ix[np.where((TC_season==int(identifier)) & (TC_basin[:,0]==0))[0]]

    working_dir='detection/JRA55/'+str(identifier)+'_JRA55/'
    elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
    found_tcs=tc_detection.tc_tracks(Wind10=Wind10,MSLP=MSLP,SST=None,VO=VO,T=T,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
    found_tcs.init_map(m=m,ax=ax,plot_lat=plot_lat,plot_lon=plot_lon)
    found_tcs.init_obs_tcs(tc_sel)
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)
    found_tcs.set_thresholds(thr_wind=15,thr_vort=5*10**(-5),thr_mslp=101500,thr_ta=0,thr_sst=26.5,win1=7,win2=12,win_step=10,neighborhood_size=8,min_time_steps=6)
    found_tcs.detect(overwrite=False)
    found_tcs.combine_tracks(overwrite=True)
    found_tcs.obs_track_info(overwrite=False)
    #found_tcs.gather_info_track(overwrite=False)
    #found_tcs.plot_track_evolution()
    found_tcs.plot_season()
    elapsed = time.time() - start;  print('Done with plotting %.3f seconds.' % elapsed)
    print('memory in use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10.**6))
    #
    # plt.close('all')
    # fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,5))
    # axes=axes.flatten()
    # maps=[]
    # for ax in axes:
    #     mm = Basemap(ax=ax,llcrnrlon=np.min(lons),urcrnrlon=np.max(lons),llcrnrlat=np.min(lats),urcrnrlat=np.max(lats),resolution="l",projection='cyl')
    #     mm.drawcoastlines(linewidth=0.7,color='m')
    #     mm.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
    #     mm.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
    #     maps.append(mm)
    #
    # found_tcs.plot_surrounding(maps=maps,axes=axes,time_steps=range(470,500))#; convert -delay 50 track_surrounding/{94..127}* TC.gif
