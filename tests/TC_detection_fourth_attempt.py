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


class tc_tracks(object):
    def __init__(self,wind,vort,mslp,sst,dates,year):
        self.year=year
        self.wind=wind
        self.vort=vort
        self.mslp=mslp
        self.sst=sst
        self.time=wind.time
        self.lat=wind.lat
        self.lon=wind.lon
        self.dates=dates
        # initialize outputs
        self.all_pos=np.zeros([1,3])
        self.id_=0
        self.tcs={}

    def save_track(self,tk,t_i,i_f,i_b,y_ini,x_ini):
        pos_ar=np.array([tk['time'],tk['lat'],tk['lon']]).T
        #if True:
        if i_f+i_b>4 and len(unique(self.all_pos.tolist()+pos_ar.tolist()))-self.all_pos.shape[0]>3:
            self.tcs[str(self.id_)]=da.DimArray(pos_ar,axes=[tk['time'],['time','lat','lon']],dims=['time','pos'])
            self.all_pos=np.concatenate((self.all_pos,pos_ar))
            self.plot_track(tk,t_i,i_f,i_b,y_ini,x_ini)
            self.id_+=1
        else:
            pass
            #print 'droped because i='+str(i)+' or overlap='+str(len(unique(self.all_pos.tolist()+pos_ar.tolist()))-self.all_pos.shape[0])

    def plot_track(self,tk,t_i,i_f,i_b,y_ini,x_ini):
        plt.close()
        fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
        m = Basemap(ax=ax,llcrnrlon=min(self.lon),urcrnrlon=max(self.lon),llcrnrlat=min(self.lat),urcrnrlat=max(self.lat),resolution="l",projection='cyl')
        m.drawmapboundary(fill_color='1.')
        # show coastlines and borders
        m.drawcoastlines(linewidth=0.3)
        m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
        m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
        m.plot(tk['lon'],tk['lat'],'-k')
        m.plot(tk['lon'][0:i_b],tk['lat'][0:i_b],'^b')
        m.plot(tk['lon'][i_b:],tk['lat'][i_b:],'^r')
        m.plot([x_ini],[y_ini],'*m')
        ax.set_title(str(dates[t_i]))
        plt.tight_layout()
        plt.savefig('plots/detected/track_'+str(self.year)+'_'+str(t_i)+'_'+str(self.id_)+'_'+str(i_b+i_f)+'_vort.png')

    def plot_season(self):
        plt.close()
        fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
        m = Basemap(ax=ax,llcrnrlon=min(self.lon),urcrnrlon=max(self.lon),llcrnrlat=min(self.lat),urcrnrlat=max(self.lat),resolution="l",projection='cyl')
        m.drawmapboundary(fill_color='1.')
        # show coastlines and borders
        m.drawcoastlines(linewidth=0.3)
        m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
        m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
        for i in range(len(self.tcs)):
            m.plot(self.tcs[str(i)].values[:,2],self.tcs[str(i)].values[:,1],marker='^')
        plt.tight_layout()
        plt.savefig('plots/detected/track_'+str(self.year)+'_found_tracks_vort.png')

    def complete_track(self,t_i,t,y_ini,x_ini):
        tk={'lat':[y_ini],'lon':[x_ini],'time':[t]}
        i_f=0
        # forewards
        for tt in self.time[t_i+1:]:
            y_0,x_0=tk['lat'][-1],tk['lon'][-1]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmax(self.vort[tt,50000,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if i_f>=2:
                    y_2=2*tk['lat'][i_f-1]-tk['lat'][i_f-2]
                    x_2=2*tk['lon'][i_f-1]-tk['lon'][i_f-2]
                    if self.mslp[tt,y_1,x_1]<self.mslp[tt,y_2,x_2]:
                        y_2,x_2=y_1,x_1
                else:
                    y_2,x_2=y_1,x_1
                y_3,x_3=np.argmin(self.mslp[tt,y_2+1.25:y_2-1.25,x_2-1.25:x_2+1.25])
                if np.nanmax(self.wind[tt,y_3+1.25:y_3-1.25,x_3-1.25:x_3+1.25])>10:
                    tk['lat']+=[y_3]
                    tk['lon']+=[x_3]
                    tk['time']+=[tt]
                    i_f+=1
                else:
                    break
            else:
                break

        i_b=0
        # backwards
        for tt in self.time[t_i-1::-1]:
            y_0,x_0=tk['lat'][0],tk['lon'][0]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmax(self.vort[tt,50000,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if i_b+i_f>=2:
                    y_2=2*tk['lat'][0]-tk['lat'][1]
                    x_2=2*tk['lon'][0]-tk['lon'][1]
                    if self.mslp[tt,y_1,x_1]<self.mslp[tt,y_2,x_2]:
                        y_2,x_2=y_1,x_1
                else:
                    y_2,x_2=y_1,x_1
                y_3,x_3=np.argmin(self.mslp[tt,y_2+1.25:y_2-1.25,x_2-1.25:x_2+1.25])
                if np.nanmax(self.wind[tt,y_3+1.25:y_3-1.25,x_3-1.25:x_3+1.25])>10:
                    tk['lat']=[y_3]+tk['lat']
                    tk['lon']=[x_3]+tk['lon']
                    tk['time']=[tt]+tk['time']
                    i_b+=1
                else:
                    break
            else:
                break

        self.save_track(tk,t_i,i_f,i_b,y_ini,x_ini)

    def complete_track_vort(self,t_i,t,y_ini,x_ini):
        tk={'lat':[y_ini],'lon':[x_ini],'time':[t]}
        i_f=0
        # forewards
        for tt in self.time[t_i+1:]:
            y_0,x_0=tk['lat'][-1],tk['lon'][-1]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmax(self.vort[tt,85000,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if np.nanmax(self.wind[tt,y_1+1.25:y_1-1.25,x_1-1.25:x_1+1.25])>10:
                    tk['lat']+=[y_3]
                    tk['lon']+=[x_3]
                    tk['time']+=[tt]
                    i_f+=1
                else:
                    break
            else:
                break

        i_b=0
        # backwards
        for tt in self.time[t_i-1::-1]:
            y_0,x_0=tk['lat'][0],tk['lon'][0]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmax(self.vort[tt,85000,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if np.nanmax(self.wind[tt,y_1+1.25:y_1-1.25,x_1-1.25:x_1+1.25])>10:
                    tk['lat']+=[y_3]
                    tk['lon']+=[x_3]
                    tk['time']+=[tt]
                    i_b+=1
                else:
                    break
            else:
                break

        self.save_track(tk,t_i,i_f,i_b,y_ini,x_ini)

    def detect_TC(self,):
        # go through time
        for t,t_i in zip(self.time,range(len(self.time))):
            # detect strong wind
            start_pos=np.where(self.wind[t,:,:]>16)
            y,x=0,0
            for y_i,x_i in zip(self.lat[start_pos[0]],self.lon[start_pos[1]]):
                # detect max wind
                y_,x_=np.argmax(self.wind[t,y_i+5:y_i-5,x_i-5:x_i+5])
                if y_!=y and x_!=x:
                    y,x=y_,x_
                    # detect max vort and min mslp
                    #y_vort,x_vort=np.argmax(self.vort[t,100000,y+2.5:y-2.5,x-2.5:x+2.5])
                    y_mslp,x_mslp=np.argmin(self.mslp[t,y+2.5:y-2.5,x-2.5:x+2.5])
                    # check if over warm water
                    if self.sst[t,y_mslp,x_mslp]>24:
                        # check if TC already detected
                        if [t,y_mslp,x_mslp] not in self.all_pos.tolist():
                            # complete the track of the cyclone
                            self.complete_track(t_i,t,y_mslp,x_mslp)
                    else:
                        pass
                        #print str(t_i)+' '+str(y)+' '+str(x)+' droped because of '+str(self.sst[t,y_mslp,x_mslp])

    def detect_TC_thresh(self,):
        # go through time
        for t,t_i in zip(self.time,range(len(self.time))):
            # detect strong wind
            start_pos=np.where((self.wind[t,:,:]>16) && (self.mslp[t,]) )
            y,x=0,0
            for y_i,x_i in zip(self.lat[start_pos[0]],self.lon[start_pos[1]]):
                # detect max wind
                y_,x_=np.argmax(self.wind[t,y_i+5:y_i-5,x_i-5:x_i+5])
                if y_!=y and x_!=x:
                    y,x=y_,x_
                    # detect max vort and min mslp
                    #y_vort,x_vort=np.argmax(self.vort[t,100000,y+2.5:y-2.5,x-2.5:x+2.5])
                    y_mslp,x_mslp=np.argmin(self.mslp[t,y+2.5:y-2.5,x-2.5:x+2.5])
                    # check if over warm water
                    if self.sst[t,y_mslp,x_mslp]>24:
                        # check if TC already detected
                        if [t,y_mslp,x_mslp] not in self.all_pos.tolist():
                            # complete the track of the cyclone
                            self.complete_track(t_i,t,y_mslp,x_mslp)
                    else:
                        pass
                        #print str(t_i)+' '+str(y)+' '+str(x)+' droped because of '+str(self.sst[t,y_mslp,x_mslp])

    def detect_TC_vort(self,):
        # go through time
        for t,t_i in zip(self.time,range(len(self.time))):
            # detect strong wind
            start_pos=np.where(self.vort[t,85000,:,:]>2.5*10**(-5))
            y,x=0,0
            for y_i,x_i in zip(self.lat[start_pos[0]],self.lon[start_pos[1]]):
                # detect max wind
                y_,x_=np.argmax(self.vort[t,85000,y_i+5:y_i-5,x_i-5:x_i+5])
                if y_!=y and x_!=x:
                    y,x=y_,x_
                    if np.max(self.wind[t,y+2.5:y-2.5,x-2.5:x+2.5])>15 and self.sst[t,y,x]>26.5:
                        # check if TC already detected
                        if [t,y,x] not in self.all_pos.tolist():
                            # complete the track of the cyclone
                            self.complete_track(t_i,t,y,x)
                    else:
                        pass
                        #print str(t_i)+' '+str(y)+' '+str(x)+' droped because of '+str(self.sst[t,y_mslp,x_mslp])

    def detect_TC__(self,):
        # go through time
        for t,t_i in zip(self.time[500:600],range(len(self.time[500:600]))):
        #for t,t_i in zip(self.time,range(len(self.time))):
            # detect strong wind
            tmp_wind=self.wind[t,:,:].copy()
            while True:
                if np.nanmax(tmp_wind)>16:
                    print np.nanmax(tmp_wind)
                    y,x=np.argmax(tmp_wind)
                    tmp_wind[y+2.5:y-2.5,x-2.5:x+2.5]=0
                    # detect max vort and min mslp
                    y_mslp,x_mslp=np.argmin(self.mslp[t,y+2.5:y-2.5,x-2.5:x+2.5])
                    # check if over warm water
                    if self.sst[t,y_mslp,x_mslp]>24:
                        # check if TC already detected
                        if [t,y_mslp,x_mslp] not in self.all_pos.tolist():
                            # complete the track of the cyclone
                            print t_i,t,y_mslp,x_mslp
                            self.complete_track(t_i,t,y_mslp,x_mslp)
                    else:
                        pass
                        #print str(t_i)+' '+str(y)+' '+str(x)+' droped because of '+str(self.sst[t,y_mslp,x_mslp])
                else:
                    break



tc_2016=tc_tracks(wind,vort,mslp,sst,dates,year)
tc_2016.detect_TC()
tc_2016.plot_season()







year=2017
vort_nc=da.read_nc('data/detection/atl_atmos_'+str(year)+'.nc')
vort=vort_nc['VO']
wind=da.read_nc('data/detection/atl_10mWind_'+str(year)+'.nc')['ws']
mslp=da.read_nc('data/detection/atl_MSL_'+str(year)+'.nc')['MSL']
sst=da.read_nc('data/detection/ERA_interim/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
time,lat,lon=wind.time,wind.lat,wind.lon
dates=[num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar) for t in time]



# plt.close()
# fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
# m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
# m.drawmapboundary(fill_color='1.')
# # show coastlines and borders
# m.drawcoastlines(linewidth=0.3)
# m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
# m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
# for i in range(len(tcs)):
#     m.plot(tcs[str(i)].values[:,2],tcs[str(i)].values[:,1],color='k',marker='^')
# ax.set_title(str(num2date(t,units = vort_nc.axes[3].units,calendar = vort_nc.axes[3].calendar)))
# plt.tight_layout()
# plt.savefig('plots/detected/track_'+str(year)+'_found_tracks.png')



#asdasd
