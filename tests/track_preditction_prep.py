import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np

def plot_boxplot(ax,distr,x,width,color):
    box=np.nanpercentile(distr,[25,75])
    lim=np.nanpercentile(distr,[0,100])
    ax.fill_between([x-width,x+width],[box[0],box[0]],[box[1],box[1]],color=color,alpha=0.3)
    ax.plot([x,x],[box[0],lim[0]],color=color,linewidth=0.3)
    ax.plot([x,x],[box[1],lim[1]],color=color,linewidth=0.3)
    median=np.nanpercentile(distr,50)
    ax.plot([x-width,x+width],[median,median],color=color)
    ax.plot([x],[np.nanmean(distr)],marker='*',color=color)

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

def load_season(year,sel_timesteps=None):
    nc=Dataset('detection/ERA_interim/atl_interim_ocean_'+year+'.nc')
    yr_fr=np.array([toYearFraction(num2date(ttt,units = nc.variables['time'].units,calendar = nc.variables['time'].calendar)) for ttt in nc.variables['time'][:]])

    ERA_atmos=da.read_nc('detection/ERA_interim/atl_interim_atmos_'+year+'.nc')
    ERA_ocean=da.read_nc('detection/ERA_interim/atl_interim_ocean_'+year+'.nc')

    lat_train=ERA_ocean.lat+0.375
    lon_train=ERA_ocean.lon
    lon_train[lon_train>180]-=360
    lon_train+= -0.375


    if sel_timesteps is not None:
        sel_time=ERA_atmos.time[sel_timesteps]
        time_train=yr_fr[sel_timesteps]

    train={}
    for var in ERA_atmos.keys():
        train[var]=da.DimArray(axes=[time_train,lat_train,lon_train,ERA_atmos.lev,['left_top','left_bottom','right_bottom','right_top']],dims=['time','lat','lon','lev','place'])
        tmp=ERA_atmos[var][sel_time,:,:,:].values
        tmp=np.rollaxis(tmp,1,4)

        train[var][:,:,:,:,'right_bottom']=tmp
        train[var][:,:,:,:,'right_top']=np.roll(tmp,1,axis=1)
        train[var][:,:,:,:,'left_bottom']=np.roll(tmp,1,axis=2)
        train[var][:,:,:,:,'left_top']=np.roll(np.roll(tmp,1,axis=2),1,axis=1)

    for var in ERA_ocean.keys():
        train[var]=da.DimArray(axes=[time_train,lat_train,lon_train,['left_top','left_bottom','right_bottom','right_top']],dims=['time','lat','lon','place'])
        tmp=ERA_ocean[var][sel_time,:,:].values

        train[var][:,:,:,'right_bottom']=tmp
        train[var][:,:,:,'right_top']=np.roll(tmp,1,axis=1)
        train[var][:,:,:,'left_bottom']=np.roll(tmp,1,axis=2)
        train[var][:,:,:,'left_top']=np.roll(np.roll(tmp,1,axis=2),1,axis=1)


    TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
    tc_sel=TC.ix[np.where((TC['season']==1979) & (TC['basin'][:,0]==0))[0]]
    t=tc_sel['source_time'].values
    lat=tc_sel['lat_for_mapping'].values
    lon=tc_sel['lon_for_mapping'].values
    lon[lon>180]-=360

    train['TC']=da.DimArray(axes=[time_train,lat_train,lon_train],dims=['time','lat','lon'])
    train['TC'].values=0

    for i in range(t.shape[0]):
        for j in range(t.shape[1]):
            if np.isfinite(t[i,j]):
                tmp_time=toYearFraction(num2date(t[i,j],units = 'days since 1858-11-17 00:00:00',calendar = 'julian'))
                if tmp_time in time_train:
                    tmp_lat=lat[i,j]
                    tmp_lon=lon[i,j]
                    if np.min(abs(lat_train-tmp_lat))<0.75 and np.min(abs(lon_train-tmp_lon))<0.75:
                        approx_lat=lat_train[np.argmin(abs(lat_train-tmp_lat))]
                        approx_lon=lon_train[np.argmin(abs(lon_train-tmp_lon))]
                        train['TC'][tmp_time,approx_lat,approx_lon]=i
                    else:
                        print tmp_lat,lat_train[np.argmin(abs(lat_train-tmp_lat))],tmp_lon,lon_train[np.argmin(abs(lon_train-tmp_lon))]


    ds=da.Dataset(train)
    ds.write_nc('detection/training_'+year+'.nc')

    # X=np.hstack((train_atmos.values.reshape((240*67*133,112)),train_ocean.values.reshape((240*67*133,4))))
    # X[np.isnan(X)]=-99
    # Y=tc_.values.reshape((240*67*133))
    #
    # np.save('detection/training_'+year+'_X',X)
    # np.save('detection/training_'+year+'_Y',Y)

load_season('1979',sel_timesteps=range(240,480))



train=da.read_nc('detection/training_1979.nc')
y=train['TC'].values.reshape((240*67*133))
indices_tc=np.where(y!=0)[0]
indices_no=np.where(y==0)[0]
indices=np.concatenate((indices_tc,indices_no[np.random.randint(len(indices_no), size=100000)]))

# check
for t,y,x in zip(id_[0],id_[1],id_[2]):
    print yearFractiontoDate(train.time[t]),train.lat[y],train.lon[x]

import matplotlib.pyplot as plt
import seaborn as sn

# plot all variables
plt.close('all')
fig,axes=plt.subplots(nrows=6,ncols=5,figsize=(12,12))
axes=axes.flatten()
count=0
for var in ['T','R','U','V']:
    for lvl in train.lev:
        ax=axes[count]
        tmp=train[var][:,:,:,lvl,:].values.reshape((240*67*133,4))
        for i in range(4):
            plot_boxplot(ax,tmp[indices_no,i],i/4.,0.05,'b')
            plot_boxplot(ax,tmp[indices_tc,i],0.1+i/4.,0.05,'r')
        ax.set_title(var+' '+str(lvl))
        count+=1

ax=axes[count]
var='SSTK'
tmp=train[var][:,:,:,:].values.reshape((240*67*133,4))
for i in range(4):
    plot_boxplot(ax,tmp[indices_no,i],i/4.,0.05,'b')
    plot_boxplot(ax,tmp[indices_tc,i],0.1+i/4.,0.05,'r')
ax.set_title(var+' '+str(lvl))
count+=1
axes[count].axis('off')
for ax in axes: ax.get_xaxis().set_visible(False)
fig.tight_layout()
plt.savefig('plots/all_vars.png',dpi=300)



#
# plt.close('all')
# fig,axes=plt.subplots(nrows=3,ncols=1,figsize=(10,6))
# count=0
# ax=axes[0]
#
# tmp=train['SSTK'][:,:,:,:].values.reshape((240*67*133,4))
# tmp-=np.nanmean(tmp)
# for i in range(4):
#     plot_boxplot(ax,tmp[indices_no,i],count+i/4.,0.05,'b')
#     plot_boxplot(ax,tmp[indices_tc,i],count+0.1+i/4.,0.05,'r')
# count+=1
#
# tmp=train['T'][:,:,:,100000,:].values.reshape((240*67*133,4))-train['T'][:,:,:,10000,:].values.reshape((240*67*133,4))
# tmp-=np.nanmean(tmp)
# for i in range(4):
#     plot_boxplot(ax,tmp[indices_no,i],count+i/4.,0.05,'b')
#     plot_boxplot(ax,tmp[indices_tc,i],count+0.1+i/4.,0.05,'r')
# count+=1
#
# tmp=train['R'][:,:,:,85000,:].values.reshape((240*67*133,4))
# tmp-=np.nanmean(tmp)
# for i in range(4):
#     plot_boxplot(ax,tmp[indices_no,i],count+i/4.,0.05,'b')
#     plot_boxplot(ax,tmp[indices_tc,i],count+0.1+i/4.,0.05,'r')
# count+=1
#
# tmp=train['U'][:,:,:,100000,:].values.reshape((240*67*133,4))-train['U'][:,:,:,10000,:].values.reshape((240*67*133,4))
# tmp-=np.nanmean(tmp)
# for i in range(4):
#     plot_boxplot(ax,tmp[indices_no,i],count+i/4.,0.05,'b')
#     plot_boxplot(ax,tmp[indices_tc,i],count+0.1+i/4.,0.05,'r')
# count+=1
#
# tmp=train['V'][:,:,:,100000,:].values.reshape((240*67*133,4))-train['V'][:,:,:,10000,:].values.reshape((240*67*133,4))
# tmp-=np.nanmean(tmp)
# for i in range(4):
#     plot_boxplot(ax,tmp[indices_no,i],count+i/4.,0.05,'b')
#     plot_boxplot(ax,tmp[indices_tc,i],count+0.1+i/4.,0.05,'r')
# count+=1
#
# plt.tight_layout()
# plt.savefig('plots/test.png',dpi=300)




#sad
