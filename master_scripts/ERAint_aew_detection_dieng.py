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
    data_path='easterly_waves/'
    local=True
except:
    sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/reanalysis/ERAint/u_v/'
    local=False

import aew_detection_dieng;    aew_detection_dieng = reload(aew_detection_dieng)

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
    surrounding=range(290,400)

for identifier in identifiers:
    RH = da.read_nc(data_path+'atl_atmos_850-700-600_'+identifier+'.nc')['R'][:,70000,:,:].squeeze().values
    nc = da.read_nc(data_path+'atl_u_v_850-700-600_'+identifier+'.nc')
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

    aews=aew_detection_dieng.aew_tracks(VO=vo,RH=RH,lats=lats,lons=lons,time_=time_,dates=dates,smoothing_factor=2,coarsening_factor=2,identifier=identifier,working_dir='aew_detection/ERAint/'+identifier+'/')
    plt.close('all')
    plate_carree = ccrs.PlateCarree()
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    ax = plt.axes(projection=plate_carree)
    ax.set_global()
    ax.coastlines()
    ax.set_xlim(np.min(lons),np.max(lons))
    ax.set_ylim(np.min(lats),np.max(lats))
    aews.init_map(ax=ax,transform=plate_carree)

    #aews.detect_dieng_groups(overwrite=True,dis_VO_max=8,max_extend=25,min_number_cells=6,thr_VO=1*10**(-5),thr_RH=50)
    aews.detect_dieng(overwrite=True,dis_VO_max=8,contour_radius=25,min_number_cells=3,thr_VO=1*10**(-5),thr_RH=50)
    aews.combine_tracks(overwrite=True,search_radius=8)
    aews.plot_detect_summary()
    aews.plot_season()

    if surrounding is not None:
        plt.close('all')
        plate_carree = ccrs.PlateCarree()
        fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,3),subplot_kw={'projection': plate_carree})
        for ax in [axes]:
            ax.set_global()
            ax.coastlines()
            ax.set_xlim(np.min(lons),np.max(lons))
            ax.set_ylim(np.min(lats),np.max(lats))

        aews.plot_surrounding(axes=[axes],time_steps=surrounding,thr_VO=1*10**(-5))
