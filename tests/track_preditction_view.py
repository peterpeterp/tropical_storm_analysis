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

def tc_plot(m, x, y, z, treshs=[82,95,112,136], colors=['#ffffcc','#ffe775','#ffc148','#ff8f20','#ff6060']):
    for i in range(len(x)-1):
        if np.isfinite(x[i+1]):
            if z[i]<=82: color='#ffffcc'
            if z[i]>82: color='#ffe775'
            if z[i]>95: color='#ffc148'
            if z[i]>112: color='#ff8f20'
            if z[i]>136: color='#ff6060'
            m.plot(x[i:i+2],y[i:i+2],color=color)
    return 0


os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
import plot_map as plot_map; reload(plot_map)
from plot_map import col_conv

def plot_boxplot(ax,distr,x,width,color):
    box=np.nanpercentile(distr,[25,75])
    lim=np.nanpercentile(distr,[0,100])
    ax.fill_between([x-width,x+width],[box[0],box[0]],[box[1],box[1]],color=color,alpha=0.3)
    ax.plot([x,x],[box[0],lim[0]],color=color,linewidth=0.3)
    ax.plot([x,x],[box[1],lim[1]],color=color,linewidth=0.3)
    median=np.nanpercentile(distr,50)
    ax.plot([x-width,x+width],[median,median],color=color)
    ax.plot([x],[np.nanmean(distr)],marker='*',color=color)

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

def latlon_to_meshgrid(lat,lon):
    x,y=lon.copy(),lat.copy()
    #x-=np.diff(x,1)[0]/2.
    #y-=np.diff(y,1)[0]/2.
    x=np.append(x,[x[-1]+np.diff(x,1)[0]])
    y=np.append(y,[y[-1]+np.diff(y,1)[0]])
    return np.meshgrid(x,y)


TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
tc_sel=TC.ix[np.where((TC['season']==1979) & (TC['basin'][:,0]==0))[0]]

train=da.read_nc('detection/tc_env_'+year+'.nc')
tc=train['TC']
pos='right_bottom'
lat_tc,lon_tc=train.lat,train.lon
lat_map,lon_map=train.lat,train.lon

for storm in set(tc[tc!=0]):

    storm=9

    i_ = np.where(tc==storm)
    print storm

    for i in range(len(i_[0])):
        t=train.time[i_[0][i]]
        yy=lat_tc[i_[1][i]]
        xx=lon_tc[i_[2][i]]

        plt.close()
        fig,axes=plt.subplots(nrows=3,ncols=2,figsize=(12,8))

        for ax,var,lvl in zip(axes.flatten(),['T','T','SSTK','R','U','V'],[100000,10000,'',100000,100000,100000]):
            if lvl!='':
                to_plot=train[var][t,:,:,lvl].values
            if lvl=='':
                to_plot=train[var][t,:,:].values
            m = Basemap(ax=ax,llcrnrlon=min(lon_map),urcrnrlon=max(lon_map),llcrnrlat=min(lat_map),urcrnrlat=max(lat_map),resolution="l",projection='cyl')
            m.drawmapboundary(fill_color='1.')
            # show coastlines and borders
            m.drawcoastlines(linewidth=0.3)
            m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
            m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)

            color_range=np.percentile(to_plot[np.isfinite(to_plot)],[0,100])
            xi,yi=latlon_to_meshgrid(lat_map,lon_map)
            im = m.pcolormesh(xi,yi,to_plot,cmap=plt.cm.jet,vmin=color_range[0],vmax=color_range[1],latlon=True)
            cb = m.colorbar(im,'bottom', size="5%", pad="2%")

            m.plot([xx],[yy],'*m')
            tc_plot(m,tc_sel['lon_for_mapping'].ix[storm,:].values,tc_sel['lat_for_mapping'].ix[storm,:].values,tc_sel['source_wind'].ix[storm,:,0].values)
            ax.set_title(var+' '+str(lvl))

            axins = zoomed_inset_axes(ax, 5, loc=1)
            plt.xticks(visible=False)
            plt.yticks(visible=False)
            map2 = Basemap(ax=axins,llcrnrlon=xx-3,urcrnrlon=xx+3,llcrnrlat=yy-3,urcrnrlat=yy+3,resolution="l",projection='cyl')
            lon_i=range(np.where(lon_tc==xx)[0]-10,np.where(lon_tc==xx)[0]+10)
            lat_i=range(np.where(lat_tc==yy)[0]-10,np.where(lat_tc==yy)[0]+10)
            xi,yi=latlon_to_meshgrid(lat_map[lat_i],lon_map[lon_i])
            to_plot_=to_plot[lat_i[0]:lat_i[-1]+1,lon_i[0]:lon_i[-1]+1]
            #color_range=[np.min(to_plot[np.isfinite(to_plot_)]),np.max(to_plot[np.isfinite(to_plot_)])]
            im2 = map2.pcolormesh(xi,yi,to_plot_,cmap=plt.cm.jet,vmin=color_range[0],vmax=color_range[1],latlon=True)
            map2.drawmapboundary(fill_color='1.')
            map2.drawcoastlines(linewidth=0.3)
            tc_plot(map2,tc_sel['lon_for_mapping'].ix[storm,:].values,tc_sel['lat_for_mapping'].ix[storm,:].values,tc_sel['source_wind'].ix[storm,:,0].values)
            map2.plot([xx],[yy],'*m')
            mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

        plt.tight_layout()
        plt.suptitle(str(yearFractiontoDate(t)))
        plt.savefig('plots/tc_view/'+str(storm)+'_'+str(i)+'.png')

    asdas


#asdas
