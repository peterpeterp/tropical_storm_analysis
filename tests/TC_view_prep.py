import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np

def toYearFraction(date):
    date=datetime(year=date.year,month=date.month,day=date.day,hour=date.hour)
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def yearFractiontoDate(yr_fr):
    year=int(yr_fr)
    fr=yr_fr-year
    startOfEpoch = datetime(year=1970, month=1, day=1)
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    seconds=fr*yearDuration+(startOfThisYear - startOfEpoch).total_seconds()
    dt=time.localtime(seconds)
    return datetime(year=dt.tm_year, month=dt.tm_mon, day=dt.tm_mday, hour=dt.tm_hour)

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

year='1979'
sel_timesteps=range(244,488)

nc=Dataset('detection/ERA_interim/atl_interim_ocean_'+year+'.nc')
yr_fr=np.array([toYearFraction(num2date(ttt,units = nc.variables['time'].units,calendar = nc.variables['time'].calendar)) for ttt in nc.variables['time'][:]])

ERA_atmos=da.read_nc('detection/ERA_interim/atl_interim_atmos_'+year+'.nc')
ERA_ocean=da.read_nc('detection/ERA_interim/atl_interim_ocean_'+year+'.nc')

lat_train=ERA_ocean.lat
lon_train=ERA_ocean.lon
lon_train[lon_train>180]-=360

if sel_timesteps is not None:
    sel_time=ERA_atmos.time[sel_timesteps]
    time_train=yr_fr[sel_timesteps]

train={}
for var in ERA_atmos.keys():
    train[var]=da.DimArray(axes=[time_train,lat_train,lon_train,ERA_atmos.lev],dims=['time','lat','lon','lev'])
    tmp=ERA_atmos[var][sel_time,:,:,:].values
    tmp=np.rollaxis(tmp,1,4)
    train[var][:,:,:,:]=tmp

for var in ERA_ocean.keys():
    train[var]=da.DimArray(axes=[time_train,lat_train,lon_train],dims=['time','lat','lon'])
    tmp=ERA_ocean[var][sel_time,:,:].values
    train[var][:,:,:]=tmp

TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
tc_sel=TC.ix[np.where((TC['season']==1979) & (TC['basin'][:,0]==0))[0]]
t=tc_sel['source_time'].values
lat=tc_sel['lat_for_mapping'].values
lon=tc_sel['lon_for_mapping'].values
lon[lon>180]-=360

train['TC']=da.DimArray(axes=[time_train,lat_train,lon_train],dims=['time','lat','lon'])
train['TC'].values=0

for storm in range(t.shape[0]):
    for ts in range(t.shape[1]):
        if np.isfinite(t[storm,ts]):
            tmp_time=toYearFraction(num2date(t[storm,ts],units = 'days since 1858-11-17 00:00:00'))
            if tmp_time in time_train:
                tmp_lat=lat[storm,ts]
                tmp_lon=lon[storm,ts]
                if np.min(abs(lat_train-tmp_lat))<0.75 and np.min(abs(lon_train-tmp_lon))<0.75:
                    approx_lat=lat_train[np.argmin(abs(lat_train-tmp_lat))]
                    approx_lon=lon_train[np.argmin(abs(lon_train-tmp_lon))]
                    train['TC'][tmp_time,approx_lat,approx_lon]=storm
                else:
                    print tmp_lat,lat_train[np.argmin(abs(lat_train-tmp_lat))],tmp_lon,lon_train[np.argmin(abs(lon_train-tmp_lon))]


ds=da.Dataset(train)
ds.write_nc('detection/tc_env_'+year+'.nc')
