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
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
from support_functions import *

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')


class tc_tracks(object):
    def __init__(self,wind,mslp,sst,dates,year,tc_sel,working_dir):
        self.year=year
        self.working_dir=working_dir
        # input fields
        self.wind=wind
        #self.vort=vort
        self.mslp=mslp
        self.sst=sst
        self.time=wind.time
        self.lat=wind.lat
        self.lon=wind.lon
        self.dates=dates
        self.yr_frac=[toYearFraction(dd) for dd in dates]
        # initialize outputs
        self.all_pos=np.zeros([1,3])
        self.id_=0
        self.tcs={}
        # ibtracks
        self.tc_sel=tc_sel
        tmp_time=tc_sel['source_time'].values
        self.tc_time=tmp_time.copy()*np.nan
        for i in range(tmp_time.shape[0]):
            for j in range(tmp_time.shape[1]):
                if np.isfinite(tmp_time[i,j]):
                    self.tc_time[i,j]=toYearFraction(num2date(tmp_time[i,j],units = 'days since 1858-11-17 00:00:00'))
        self.tc_lat=tc_sel['lat_for_mapping'].values
        self.tc_lon=tc_sel['lon_for_mapping'].values
        self.tc_lon[self.tc_lon<0]+=360
        self.tc_intens=tc_sel['source_wind'].values
        # init map that is used for plotting
        fig,self.ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
        self.m = Basemap(ax=self.ax,llcrnrlon=min(self.lon),urcrnrlon=max(self.lon),llcrnrlat=min(self.lat),urcrnrlat=max(self.lat),resolution="l",projection='cyl')
        self.m.drawmapboundary(fill_color='1.')
        self.m.drawmapboundary(fill_color='#99ffff')
        self.m.fillcontinents(color='#cc9966',lake_color='#99ffff')
        self.m.drawcoastlines(linewidth=0.3)
        self.m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
        self.m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
        for storm in range(len(self.tc_sel.storm)):
            self.m.plot(self.tc_lon[storm,:],self.tc_lat[storm,:],color='gray')

    def save_track(self,tk,t_i,i_f,i_b,y_ini,x_ini,plot):
        pos_ar=np.array([tk['time'],tk['lat'],tk['lon']]).T
        #if True:
        if i_f+i_b>10 and len(unique(self.all_pos.tolist()+pos_ar.tolist()))-self.all_pos.shape[0]>3:
            self.tcs[str(self.id_)]=da.DimArray(pos_ar,axes=[tk['time'],['time','lat','lon']],dims=['time','pos'])
            self.all_pos=np.concatenate((self.all_pos,pos_ar))
            if plot:    self.plot_track(tk,t_i,i_f,i_b,y_ini,x_ini)
            self.id_+=1
        else:
            pass
            #print 'droped because i='+str(i)+' or overlap='+str(len(unique(self.all_pos.tolist()+pos_ar.tolist()))-self.all_pos.shape[0])

    def plot_track(self,tk,t_i,i_f,i_b,y_ini,x_ini):
        tmp,text=[],[]
        tmp.append(self.m.plot(tk['lon'],tk['lat'],'-k'))
        tmp.append(self.m.plot(tk['lon'][0:i_b],tk['lat'][0:i_b],'.b'))
        tmp.append(self.m.plot(tk['lon'][i_b:],tk['lat'][i_b:],'.r'))
        tmp.append(self.m.plot([x_ini],[y_ini],'*m'))
        self.ax.set_title(str(dates[t_i]))


        storms=[np.where(self.tc_time==zz)[0][0] for zz in self.yr_frac[max(0,t_i-i_b-5):min(t_i+i_f+5,len(self.yr_frac)-1)] if zz in self.tc_time]
        for storm in set(storms):
            tmp+=tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens[storm,:,0])
            last_pos=np.where(np.isfinite(self.tc_lon[storm,:]))[0][-1]
            text.append(self.ax.text(self.tc_lon[storm,last_pos],self.tc_lat[storm,last_pos],''.join(self.tc_sel['name'].ix[storm,:])))

        plt.tight_layout()
        plt.savefig(self.working_dir+'track_'+str(self.year)+'_'+str(t_i)+'_'+str(self.id_)+'_'+str(i_b+i_f)+'.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in text:
            element.remove()

    def plot_season(self,out_name=None):
        if out_name is None:
            out_name=self.working_dir+'track_'+str(self.year)+'_found_tracks_vort.png'

        for i in range(len(self.tcs)):
            self.m.plot(self.tcs[str(i)].values[:,2],self.tcs[str(i)].values[:,1],marker='.')

        for storm in range(len(self.tc_sel.storm)):
            tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens[storm,:,0])
        plt.tight_layout()
        plt.savefig(out_name)

    def complete_track(self,t_i,t,y_ini,x_ini,plot):
        tk={'lat':[y_ini],'lon':[x_ini],'time':[t]}
        i_f=0
        # search forewards
        for tt in self.time[t_i+1:]:
            y_0,x_0=tk['lat'][-1],tk['lon'][-1]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmin(self.mslp[tt,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if np.nanmax(self.wind[tt,y_1+2.5:y_1-2.5,x_1-2.5:x_1+2.5])>10:
                    tk['lat']+=[y_1]
                    tk['lon']+=[x_1]
                    tk['time']+=[tt]
                    i_f+=1
                else:
                    break
            else:
                break

        # search backwards
        i_b=0
        for tt in self.time[t_i-1::-1]:
            y_0,x_0=tk['lat'][0],tk['lon'][0]
            if y_0<self.lat[4] and y_0>self.lat[-4] and x_0<self.lon[-4] and x_0>self.lon[4]:
                y_1,x_1=np.argmin(self.mslp[tt,y_0+2.5:y_0-2.5,x_0-2.5:x_0+2.5])
                if np.nanmax(self.wind[tt,y_1+2.5:y_1-2.5,x_1-2.5:x_1+2.5])>8:
                    tk['lat']=[y_1]+tk['lat']
                    tk['lon']=[x_1]+tk['lon']
                    tk['time']=[tt]+tk['time']
                    i_b+=1
                else:
                    break
            else:
                break

        self.save_track(tk,t_i,i_f,i_b,y_ini,x_ini,plot)

    def detect_TC(self,wind_thresh=15,sst_thresh=26,detection_window=None,plot=True):

        if detection_window is None:
            detection_window=self.time

        # go through time
        for t,t_i in zip(detection_window,range(len(detection_window))):
            # detect strong wind
            y,x=local_extrem(self.wind[t,:,:],'max',threshold=wind_thresh)
            for y_,x_ in zip(self.lat[y],self.lon[x]):
                # find the eye of the cyclone
                y_c,x_c=np.argmin(self.mslp[t,y_+2.5:y_-2.5,x_-2.5:x_+2.5])
                # check if over warm water (or some other conditions)
                if self.sst[t,y_c,x_c]>sst_thresh:
                    # check if TC already detected
                    if [t,y_c,x_c] not in self.all_pos.tolist():
                        # complete the track of the cyclone
                        self.complete_track(t_i,t,y_c,x_c,plot)


# # read ibtracks
# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')

found_tracks={}
for year in [2008,2016]:
    start = time.time()
    # read ERA interim
    wind_nc=da.read_nc('data/detection/atl_10mWind_'+str(year)+'.nc')
    wind=wind_nc['ws']
    mslp=da.read_nc('data/detection/atl_MSL_'+str(year)+'.nc')['MSL']
    sst=da.read_nc('data/detection/ERA_interim/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
    dates=[num2date(t,units = wind_nc.axes[0].units,calendar = wind_nc.axes[0].calendar) for t in wind.time]

    tc_sel=TC.ix[np.where((TC['season']==year) & (TC['basin'][:,0]==0))[0]]
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)

    os.system('mkdir plots/detected/'+str(year))
    found_tracks[year]=tc_tracks(wind,mslp,sst,dates,year,tc_sel,working_dir='plots/detected/'+str(year)+'/')
    found_tracks[year].detect_TC(wind_thresh=15,sst_thresh=26,plot=True)
    found_tracks[year].plot_season()
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)
