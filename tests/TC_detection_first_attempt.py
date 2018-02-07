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

for t,t_i in zip(time,range(len(time))):
    start_pos=np.where(wind[t,:,:]>17.5)
    y,x=0,0
    for i in range(len(start_pos[0])):
        y_i,x_i=lat[start_pos[0][i]],lon[start_pos[1][i]]
        y_,x_=np.argmax(wind[t,y_i+5:y_i-5,x_i-5:x_i+5])
        if y_!=y and x_!=x:
            y,x=y_,x_

            min_vort=np.argmax(vort[t,100000,y+2.5:y-2.5,x-2.5:x+2.5])
            min_mslp=np.argmin(mslp[t,y+2.5:y-2.5,x-2.5:x+2.5])
            if sst[t,min_vort[0],min_vort[1]]>26.5:
                print t_i,min_vort,min_mslp,wind[t,y,x],str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)),sst[t,min_vort[0],min_vort[1]]
                plt.close()
                fig,axes=plt.subplots(nrows=1,ncols=3,figsize=(12,3))
                for ax,to_plot,title in zip(axes.flatten(),[wind[t,:,:],vort[t,100000,:,:],mslp[t,:,:]],['10m wind speed','relative vorticity','mean sea level pressure']):
                    m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
                    m.drawmapboundary(fill_color='1.')
                    # show coastlines and borders
                    m.drawcoastlines(linewidth=0.3)
                    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
                    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
                    color_range=np.percentile(to_plot[np.isfinite(to_plot)],[0,100])
                    xi,yi=latlon_to_meshgrid(lat,lon)
                    im = m.pcolormesh(xi,yi,to_plot,cmap=plt.cm.jet,vmin=color_range[0],vmax=color_range[1],latlon=True)
                    cb = m.colorbar(im,'bottom', size="5%", pad="2%")
                    m.plot([min_vort[1]],[min_vort[0]],'^k')
                    m.plot([min_mslp[1]],[min_mslp[0]],'vw')
                    ax.set_title(title)
                    #tc_plot(m,tc_sel['lon_for_mapping'].ix[storm,:].values,tc_sel['lat_for_mapping'].ix[storm,:].values,tc_sel['source_wind'].ix[storm,:,0].values)

                plt.tight_layout()
                plt.suptitle(str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)))
                plt.savefig('plots/detected/'+str(year)+'_'+str(t_i)+'_'+str(i)+'.png')
