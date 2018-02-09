import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
import seaborn as sn
import webcolors
import numpy as np
import dimarray as da
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
from support_functions import *

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
tc_sel=TC.ix[np.where((TC['season']==2008) & (TC['basin'][:,0]==0))[0]]
t=tc_sel['source_time'].values

tc_time=t.copy()*np.nan
for i in range(t.shape[0]):
    for j in range(t.shape[1]):
        if np.isfinite(t[i,j]):
            tc_time[i,j]=toYearFraction(num2date(t[i,j],units = 'days since 1858-11-17 00:00:00'))

lat_tc=tc_sel['lat_for_mapping'].values
lon_tc=tc_sel['lon_for_mapping'].values
lon_tc[lon_tc<0]+=360


year=2008
os.system('mkdir plots/detected/view_dayX_'+str(year))
vort_nc=da.read_nc('data/ERA_dayX/atl_atmos_'+str(year)+'.nc')
vort=vort_nc['VO']
temp=vort_nc['T']
wind=da.read_nc('data/ERA_dayX/atl_10mWind_'+str(year)+'.nc')['ws']
mslp=da.read_nc('data/ERA_dayX/atl_MSL_'+str(year)+'.nc')['MSL']
sst=da.read_nc('data/ERA_dayX/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
time,lat,lon=wind.time,wind.lat,wind.lon
dates=[num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar) for t in time]
year_fraction=[toYearFraction(dd) for dd in dates]


for i in range(len(wind.time)):


    plt.close()
    fig,axes = plt.subplots(nrows=2,ncols=2)
    axes=axes.flatten()

    ax=axes[0]
    ax.set_title('10m wind speed')
    data=wind.values[i,:,:]
    y_w,x_w = local_extrem(data,'max',threshold=12,neighborhood_size=3)
    im=ax.imshow(data,vmin=0,vmax=13)
    im.set_cmap('bone')
    ax.autoscale(False)
    ax.axis('off')
    ax.scatter(x_w,y_w,color='r',s=50,facecolors='none')
    for y,x in zip(y_w,x_w):
        # y_=lat[y]
        # x_=lon[x]
        y,x=np.argmin(mslp.ix[i,max([0,y-3]):min([y+3,len(lat)-1]),max([0,x-3]):min([x+3,len(lon)-1])])
        y=np.where(lat==y)[0][0]
        x=np.where(lon==x)[0][0]
        if np.max(wind.values[i,max([0,y-3]):min([y+3,len(lat)-1]),max([0,x-3]):min([x+3,len(lon)-1])])>16:
            ax.plot([x],[y], 'm^')
        if np.max(sst.values[i,max([0,y-3]):min([y+3,len(lat)-1]),max([0,x-3]):min([x+3,len(lon)-1])])>6:
            ax.plot([x],[y], 'c*')

    ax=axes[1]
    ax.set_title('rel. Vorticity')
    data=vort.values[i,0,:,:]
    y_v,x_v = local_extrem(data,'max',threshold=2.5*10**(-5))
    im=ax.imshow(data,vmin=-9.5*10**(-5),vmax=0.0002)
    im.set_cmap('bone')
    ax.autoscale(False)
    ax.axis('off')
    ax.scatter(x_v,y_v,color='r',s=50,facecolors='none')

    ax=axes[2]
    ax.set_title('detected (MSLP)')
    data=mslp.values[i,:,:]
    y_p,x_p = local_extrem(data,'min',threshold=101500,neighborhood_size=5)
    im=ax.imshow(data,vmin=100360,vmax=103000)
    im.set_cmap('bone')
    ax.autoscale(False)
    ax.axis('off')
    ax.scatter(x_p,y_p,color='r',s=50,facecolors='none')

    ax=axes[3]
    # aplot showing max wind, min mslp and max vort next to tc loc
    im=ax.imshow(wind.values[i,:,:]*0,vmin=0,vmax=13)
    im.set_cmap('bone')
    ax.autoscale(False)
    ax.axis('off')
    # add ibtracks info
    obs_tc=np.where(abs(tc_time-year_fraction[i])<0.002)
    if len(obs_tc[0])>0:
        for storm in set(obs_tc[0]):
            first_day=obs_tc[1][obs_tc[0]==storm][0]
            # if np.isfinite(tc_sel['source_wind'].ix[storm,first_day,0]):
            x_,y_=np.argmin(abs(lon-lon_tc[storm,first_day])),np.argmin(abs(lat-lat_tc[storm,first_day]))
            ax.text(x_,y_,'wind: '+str(round(np.max(wind[wind.time[i],lat[y_]+3:lat[y_]-3,lon[x_]-3:lon[x_]+3]),01))+'\nmslp: '+str(round(np.min(mslp[mslp.time[i],lat[y_]+3:lat[y_]-3,lon[x_]-3:lon[x_]+3]),06))+'\nvort: '+str(np.max(vort[vort.time[i],85000,lat[y_]+3:lat[y_]-3,lon[x_]-3:lon[x_]+3])),color='white',va='top',fontsize=7)

    # add ibtracks info
    obs_tc=np.where(abs(tc_time-year_fraction[i])<0.002)
    if len(obs_tc[0])>0:
        for oo in range(len(obs_tc[0])):
            if np.isfinite(tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]):
                for ax in axes:
                    ax.plot(np.argmin(abs(lon-lon_tc[obs_tc[0][oo],obs_tc[1][oo]])),np.argmin(abs(lat-lat_tc[obs_tc[0][oo],obs_tc[1][oo]])),color=get_tc_color(tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]),marker='.')


    plt.suptitle(str(dates[i]))
    plt.tight_layout()
    plt.savefig('plots/detected/view_dayX_'+str(year)+'/'+str(i)+'.png', bbox_inches = 'tight')

#
# import imageio
# images = []
# for filename in sorted(glob.glob('plots/detected/view_'+str(year)+'/*')):
#     images.append(imageio.imread(filename))
# imageio.mimsave('plots/detected/daily_'+str(year)+'.gif', images)
