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

try:
    sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/ERAint/'
    local=True
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/reanalysis/ERAint/'
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

    surrounding=args.surrounding

    print(identifiers)

if local:
    identifiers=['2008']
    surrounding=range(280,310)

for identifier in identifiers:
    RH = da.read_nc(data_path+'atmos/atl_atmos_850-700-600_'+identifier+'.nc')['R'][:,70000,:,:].squeeze().values
    nc = da.read_nc(data_path+'u_v/atl_u_v_850-700-600_'+identifier+'.nc')
    u=nc['U'][:,70000,:,:].squeeze()
    v=nc['V'][:,70000,:,:].squeeze().values
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

    curv_vort = da.read_nc(data_path+'u_v/ano_curv_vort_850-700-600_'+identifier+'.nc')['curv_vort'][:,70000,:,:].squeeze().values

    dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=-1))/dx
    dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=-2))/dy
    curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

    plt.close('all')
    plate_carree = ccrs.PlateCarree()
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    ax = plt.axes(projection=plate_carree)
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(-70,30)
    ax.set_ylim(0,30)

    # dieng
    dieng=aew_detection.aew_tracks(identifier=identifier+'_dieng',working_dir='aew_detection/ERAint/'+identifier+'/')
    dieng.prepare_data(lats=lats,lons=lons,time_=time_,dates=dates,smoothing_factor=2,coarsening_factor=2)
    dieng.add_fields(VO=vo,RH=RH)
    dieng.init_map(ax=ax,transform=plate_carree)

    dieng.detect_dieng(overwrite=True,dis_VO_max=8,contour_radius=25,min_number_cells=3,thr_VO=1*10**(-5),thr_RH=50)
    dieng.plot_detect_summary()
    dieng.combine_tracks(overwrite=True)
    dieng.plot_season()

    # belanger
    belanger=aew_detection.aew_tracks(identifier=identifier+'_belanger',working_dir='aew_detection/ERAint/'+identifier+'/')
    belanger.prepare_data(lats=lats,lons=lons,time_=time_,dates=dates,smoothing_factor=3,coarsening_factor=2)
    belanger.add_fields(U=u,V=v,CURV_VORT=curv_vort,CURV_VORT_ADVECT=curv_vort_advect,VO=vo)
    belanger.init_map(ax=ax,transform=plate_carree)

    belanger.detect_belanger(overwrite=True,thr_curv_vort=3.2*10**(-6))
    belanger.plot_detect_summary()
    belanger.combine_tracks(overwrite=True)
    belanger.plot_season()


    if surrounding is not None:
        for t in surrounding:
            plt.close('all')
            plate_carree = ccrs.PlateCarree()
            fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(16,6),subplot_kw={'projection': plate_carree})
            for ax in axes.flatten():
                ax.set_global()
                ax.coastlines()
                ax.set_xlim(-70,30)
                ax.set_ylim(0,30)

            dieng.add_surounding_dieng(axes=axes[0,:],t=t,thr_VO=1*10**(-5))
            belanger.add_surounding_belanger(axes=axes[1,:],t=t,thr_curv_vort=3.2*10**(-6))
            plt.suptitle(str(aews._dates[t]))
            plt.tight_layout()
            plt.savefig(aews._working_dir+'track_surrounding/combined_'+str(t)+'.png', bbox_inches = 'tight')
