import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import seaborn as sn
import webcolors
import matplotlib
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def latlon_to_meshgrid(lat,lon):
    x,y=lon.copy(),lat.copy()
    #x-=np.diff(x,1)[0]/2.
    #y-=np.diff(y,1)[0]/2.
    x=np.append(x,[x[-1]+np.diff(x,1)[0]])
    y=np.append(y,[y[-1]+np.diff(y,1)[0]])
    return np.meshgrid(x,y)

def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

def unique(l):
    return list(uniq(sorted(l, reverse=True)))

os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
from support_functions import *

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

year=2016
vort_nc=da.read_nc('data/detection/atl_atmos_'+str(year)+'.nc')
vort=vort_nc['VO']
wind=da.read_nc('data/detection/atl_10mWind_'+str(year)+'.nc')['ws']
mslp=da.read_nc('data/detection/atl_MSL_'+str(year)+'.nc')['MSL']
sst=da.read_nc('data/detection/ERA_interim/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
time,lat,lon=wind.time,wind.lat,wind.lon

tcs={}
all_pos=np.zeros([1,3])
id_=0
for t,t_i in zip(time,range(len(time))):

    start_pos=np.where(wind[t,:,:]>13.5)
    y,x=0,0
    for y_i,x_i in zip(lat[start_pos[0]],lon[start_pos[1]]):
        y_,x_=np.argmax(wind[t,y_i+5:y_i-5,x_i-5:x_i+5])
        if y_!=y and x_!=x:
            y,x=y_,x_
            y_vort,x_vort=np.argmax(vort[t,100000,y+2.5:y-2.5,x-2.5:x+2.5])
            y_mslp,x_mslp=np.argmin(mslp[t,y+2.5:y-2.5,x-2.5:x+2.5])
            if sst[t,y_mslp,x_mslp]>26.5:
                # TC found

                if [t,y_mslp,x_mslp] not in all_pos.tolist():
                    #print [t,y_mslp,x_mslp]
                    #print all_pos.tolist()

                    tk={'lat':[y_mslp],'lon':[x_mslp],'time':[t]}
                    i=0
                    for tt in time[t_i+1:]:
                        y_0,x_0=tk['lat'][-1],tk['lon'][-1]
                        if y_0<lat[4] and y_0>lat[-4] and x_0<lon[-4] and x_0>lon[4]:
                            y_1,x_1=np.argmax(vort[tt,100000,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                            if i>=2:
                                y_2=2*tk['lat'][i-1]-tk['lat'][i-2]
                                x_2=2*tk['lon'][i-1]-tk['lon'][i-2]
                                if mslp[tt,y_1,x_1]<mslp[tt,y_2,x_2]:
                                    y_2,x_2=y_1,x_1
                            else:
                                y_2,x_2=y_1,x_1
                            y_3,x_3=np.argmin(mslp[tt,y_2+2.5:y_2-2.5,x_2-2.5:x_2+2.5])
                            if np.nanmax(wind[tt,y_3+2.5:y_3-2.5,x_3-2.5:x_3+2.5])>13.9:
                                tk['lat'].append(y_3)
                                tk['lon'].append(x_3)
                                tk['time'].append(tt)
                                i+=1
                            else:
                                pos_ar=np.array([tk['time'],tk['lat'],tk['lon']]).T
                                if i>2 and len(unique(all_pos.tolist()+pos_ar.tolist()))-all_pos.shape[0]>3:
                                    plt.close()
                                    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
                                    m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
                                    m.drawmapboundary(fill_color='1.')
                                    # show coastlines and borders
                                    m.drawcoastlines(linewidth=0.3)
                                    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
                                    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
                                    m.plot(tk['lon'],tk['lat'],'^k')
                                    ax.set_title(str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)))
                                    plt.tight_layout()
                                    plt.savefig('plots/detected/track_'+str(year)+'_'+str(t_i)+'_'+str(i)+'.png')

                                    tcs[str(id_)]=da.DimArray(pos_ar,axes=[tk['time'],['time','lat','lon']],dims=['time','pos'])
                                    all_pos=np.concatenate((all_pos,pos_ar))
                                    id_+=1
                                break
                        else:
                            pos_ar=np.array([tk['time'],tk['lat'],tk['lon']]).T
                            if i>2 and len(unique(all_pos.tolist()+pos_ar.tolist()))-all_pos.shape[0]>3:
                                plt.close()
                                fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
                                m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
                                m.drawmapboundary(fill_color='1.')
                                # show coastlines and borders
                                m.drawcoastlines(linewidth=0.3)
                                m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
                                m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
                                m.plot(tk['lon'],tk['lat'],'^k')
                                ax.set_title(str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)))
                                plt.tight_layout()
                                plt.savefig('plots/detected/track_'+str(year)+'_'+str(t_i)+'_'+str(i)+'.png')

                                tcs[str(id_)]=da.DimArray(pos_ar,axes=[tk['time'],['time','lat','lon']],dims=['time','pos'])
                                all_pos=np.concatenate((all_pos,pos_ar))
                                id_+=1
                            break



plt.close()
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
m.drawmapboundary(fill_color='1.')
# show coastlines and borders
m.drawcoastlines(linewidth=0.3)
m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
for i in range(len(tcs)):
    m.plot(tcs[str(i)].values[:,2],tcs[str(i)].values[:,1],color='k',marker='^')
ax.set_title(str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)))
plt.tight_layout()
plt.savefig('plots/detected/track_'+str(year)+'_found_tracks.png')



#asdasd
