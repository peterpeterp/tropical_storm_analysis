from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy
from shapely.geometry import Polygon
mpl.use('Agg')

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

import aew_detection;    aew_detection = reload(aew_detection)

def toYearFraction(date):
    date=datetime(year=date.year,month=date.month,day=date.day,hour=date.hour)
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed/yearDuration

    return date.year + fraction

if local==False:
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
        #identifiers=['2010']

    print(identifiers)

if local:
    identifiers=['2008']

for identifier in identifiers:
    curv_vort = da.read_nc(data_path+'ano_curv_vort_850-700-600_'+identifier+'.nc')['curv_vort'][:,70000,35:1,-89:30].squeeze().values
    nc = da.read_nc(data_path+'atl_u_v_850-700-600_'+identifier+'.nc')
    u=nc['U'][:,70000,35:1,-89:30].squeeze()
    v=nc['V'][:,70000,35:1,-89:30].squeeze().values
    lon,lat=u.lon,u.lat
    u=u.values
    lon[lon>180]-=360
    lon_rolling=len(lon)-np.where(lon<0)[0][0]
    lon=np.roll(lon,lon_rolling,axis=-1)
    lon=np.roll(lon,lon_rolling,axis=-1)
    u=np.roll(u,lon_rolling,axis=-1)
    v=np.roll(v,lon_rolling,axis=-1)
    lons,lats=np.meshgrid(lon,lat)
    time_=nc.time
    dates=[num2date(t,units = nc.axes['time'].units,calendar = nc.axes['time'].calendar) for t in time_]
    xx,yy = np.meshgrid(lon,lat)
    dx,dy = np.meshgrid(lon.copy()*0+np.mean(np.diff(lon,1)),lat.copy()*0+np.mean(np.diff(lat,1)))
    dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
    dy*=6371000*2*np.pi/360.

    W = (u**2+v**2)**0.5
    du_dx = (u-np.roll(u,1,axis=-1))/dx
    du_dy = (u-np.roll(u,1,axis=-2))/dy
    dv_dx = (v-np.roll(v,1,axis=-1))/dx
    dv_dy = (v-np.roll(v,1,axis=-2))/dy
    dW_dx = (W-np.roll(W,1,axis=-1))/dx
    dW_dy = (W-np.roll(W,1,axis=-2))/dy
    vo=dv_dx-du_dy
    dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=-1))/dx
    dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=-2))/dy
    curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

    aews=aew_detection.aew_tracks(U=u,V=v,vo=vo,curv_vort=curv_vort,curv_vort_advect=curv_vort_advect,lats=lats,lons=lons,smoothing_factor=3,coarsening_factor=2,time_=time_,dates=dates,identifier=identifier,working_dir='aew_detection/ERAint/'+identifier+'/')
    plt.close('all')
    plate_carree = ccrs.PlateCarree()
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    ax = plt.axes(projection=plate_carree)
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))
    aews.init_map(ax=ax,transform=plate_carree)

    aews.detect(overwrite=True,thr_curv_vort=3.2*10**(-6)
    aews.combine_tracks(overwrite=True)
    aews.plot_detect_summary()
    aews.plot_season()

    if args.surrounding is not None:
        plt.close('all')
        plate_carree = ccrs.PlateCarree()
        fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(16,3),subplot_kw={'projection': plate_carree})
        for ax in axes:
            ax.set_global()
            ax.coastlines()
            ax.set_xlim(np.min(lons),np.max(lons))
            ax.set_ylim(np.min(lats),np.max(lats))

        aews.plot_surrounding(axes=axes,time_steps=args.surrounding,thr_curv_vort=3.2*10**(-6)
