import os,sys,glob,time,collections,gc,calendar,weakref
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

from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches

os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
from support_functions import *

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')





class tc_tracks(object):
    def __init__(self,wind,mslp,sst,dates,year,tc_sel,working_dir,time_steps=None):
        self.year=year
        self.working_dir=working_dir
        # input fields
        if time_steps is None:
            time_steps=range(len(wind.time))
        self.time_i=time_steps
        self.wind=wind
        self.vort=vort
        self.mslp=mslp
        self.sst=sst
        self.time=wind.time
        self.lat=wind.lat
        self.lon=wind.lon
        self.dates=dates
        self.yr_frac=[toYearFraction(dd) for dd in self.dates]
        # initialize outputs
        self.all_pos=np.zeros([1,3])
        self.id_=0
        self.tcs={}
        self.detection_t_i=[]
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
        self.tc_intens=tc_sel['source_wind']
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

    def local_max(self,data,threshold,neighborhood_size):
        y,x=np.where(data>threshold)
        ex_x,ex_y=[],[]
        for y_,x_ in zip(y,x):
            if data[y_,x_]==np.max(data[max([0,y_-neighborhood_size]):min([y_+neighborhood_size,len(self.lat)-1]),max([0,x_-neighborhood_size]):min([x_+neighborhood_size,len(self.lon)-1])]):
                ex_y.append(y_)
                ex_x.append(x_)
        return ex_y,ex_x

    def local_min(self,data,threshold,neighborhood_size):
        y,x=np.where(data<threshold)
        ex_x,ex_y=[],[]
        for y_,x_ in zip(y,x):
            if data[y_,x_]==np.min(data[max([0,y_-neighborhood_size]):min([y_+neighborhood_size,len(lat)-1]),max([0,x_-neighborhood_size]):min([x_+neighborhood_size,len(lon)-1])]):
                ex_y.append(y_)
                ex_x.append(x_)
        return ex_y,ex_x

    def save_track(self,tk,t_i,i_f,i_b,y_ini,x_ini,plot):
        pos_ar=np.array([tk['time'],tk['lat'],tk['lon']]).T
        #if i_f+i_b>1 and len(unique(self.all_pos.tolist()+pos_ar.tolist()))-self.all_pos.shape[0]>4:
        if i_f+i_b>2:
            if sum([bb in self.all_pos.tolist() for bb in pos_ar.tolist()])<2:

                self.tcs[str(self.id_)]=da.DimArray(pos_ar,axes=[tk['time'],['time','lat','lon']],dims=['time','pos'])
                self.all_pos=np.concatenate((self.all_pos,pos_ar))
                if plot:    self.plot_track(tk,t_i,i_f,i_b,y_ini,x_ini)
                self.detection_t_i.append(t_i)
                self.id_+=1

    def plot_track(self,tk,t_i,i_f,i_b,y_ini,x_ini):
        tmp,text=[],[]
        tmp.append(self.m.plot(tk['lon'],tk['lat'],'-k'))
        tmp.append(self.m.plot(tk['lon'][0:i_b],tk['lat'][0:i_b],'.b'))
        tmp.append(self.m.plot(tk['lon'][i_b:],tk['lat'][i_b:],'.r'))
        tmp.append(self.m.plot([x_ini],[y_ini],'*m'))
        self.ax.set_title(str(dates[t_i]))


        storms=[np.where(self.tc_time==zz)[0][0] for zz in self.yr_frac[max(0,t_i-i_b-5):min(t_i+i_f+5,len(self.yr_frac)-1)] if zz in self.tc_time]
        print t_i,self.dates[t_i]
        print set(storms)
        storms=np.where(abs(self.tc_time-self.yr_frac[t_i])<0.002)[0]
        print set(storms)
        for storm in set(storms):
            tmp+=tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.ix[storm,:,0])
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
            out_name=self.working_dir+'track_'+str(self.year)+'_found_tracks.png'

        for i in range(len(self.tcs)):
            self.m.plot(self.tcs[str(i)].values[:,2],self.tcs[str(i)].values[:,1],marker='.')

        for storm in range(len(self.tc_sel.storm)):
            tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.values[storm,:,0])

        tmp=np.array(self.detected)

        self.m.scatter(tmp[:,2],tmp[:,1],marker='.')

        plt.tight_layout()
        plt.savefig(out_name)

    def plot_surrounding(self,win1=7,win2=4,thr_wind=13,thr_sst=26):
        for t,t_i in zip(self.time[self.time_i],self.time_i):
            plt.close()
            fig,axes = plt.subplots(nrows=2,ncols=2)
            axes=axes.flatten()

            ax=axes[0]
            ax.set_title('10m wind speed')
            data=self.wind.values[t_i,:,:]
            y_w,x_w = self.local_max(data,threshold=wind_thresh,neighborhood_size=3)
            im=ax.imshow(data,vmin=0,vmax=13)
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_w,y_w,color='r',s=50,facecolors='none')

            ax=axes[1]
            ax.set_title('rel. Vorticity')
            data=self.vort.values[t_i,0,:,:]
            y_v,x_v = local_extrem(data,'max',threshold=2*10**(-5),neighborhood_size=15)
            im=ax.imshow(data,vmin=-9.5*10**(-5),vmax=0.0002)
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_v,y_v,color='r',s=50,facecolors='none')

            ax=axes[2]
            # aplot showing max wind, min mslp and max vort next to tc loc
            im=ax.imshow(self.wind.values[t_i,:,:]*0,vmin=0,vmax=13)
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            # add ibtracks info
            obs_tc=np.where(abs(self.tc_time-self.yr_frac[t_i])<0.002)
            if len(obs_tc[0])>0:
                for storm in set(obs_tc[0]):
                    first_day=obs_tc[1][obs_tc[0]==storm][0]
                    # if np.isfinite(tc_sel['source_wind'].ix[storm,first_day,0]):
                    x_,y_=np.argmin(abs(self.lon-self.tc_lon[storm,first_day])),np.argmin(abs(self.lat-self.tc_lat[storm,first_day]))
                    ax.text(x_,y_,'wind: '+str(round(np.max(self.wind[t,self.lat[y_]+3:self.lat[y_]-3,self.lon[x_]-3:self.lon[x_]+3]),01))+'\nmslp: '+str(round(np.min(mslp[t,self.lat[y_]+3:self.lat[y_]-3,self.lon[x_]-3:self.lon[x_]+3]),06))+'\nvort: '+str(np.max(vort[t,85000,self.lat[y_]+3:self.lat[y_]-3,self.lon[x_]-3:self.lon[x_]+3])),color='white',va='top',fontsize=7)

            # add ibtracks info
            obs_tc=np.where(abs(self.tc_time-self.yr_frac[t_i])<0.002)
            if len(obs_tc[0])>0:
                for oo in range(len(obs_tc[0])):
                    if np.isfinite(self.tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]):
                        for ax in axes:
                            ax.plot(np.argmin(abs(self.lon-self.tc_lon[obs_tc[0][oo],obs_tc[1][oo]])),np.argmin(abs(self.lat-self.tc_lat[obs_tc[0][oo],obs_tc[1][oo]])),color=get_tc_color(self.tc_intens.ix[obs_tc[0][oo],obs_tc[1][oo],0]),marker='.')

            ax=axes[3]
            ax.set_title('detected (MSLP)')
            data=self.mslp.values[t_i,:,:]
            y_p,x_p = local_extrem(data,'min',threshold=101500,neighborhood_size=10)
            im=ax.imshow(data,vmin=100360,vmax=103000)
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_p,y_p,marker='o',color='r',s=50,facecolors='none')
            ax.scatter(x_v,y_v,marker='*',color='m')
            ax.scatter(x_w,y_w,marker='*',color='c')


            for y,x in zip(y_p,x_p):
                for y_,x_ in zip(y_v,x_v):
                    #print x,y,x_,y_,x-win1,x+win1
                    if x_<x+win2 and x_>x-win2 and y_<y+win2 and y_>y-win2:
                        ax.plot([x-win1*1.1,x-win1*1.1,x+win1*1.1,x+win1*1.1,x-win1*1.1],[y-win1*1.1,y+win1*1.1,y+win1*1.1,y-win1*1.1,y-win1*1.1],'-m')
                        if np.max(self.wind.ix[t_i,max(0,y-win1):min(y+win1,len(self.lat)-1),max(x-win1,0):min(x+win1,len(self.lon)-1)])>thr_wind:
                            ax.plot([x-win1*1.2,x-win1*1.2,x+win1*1.2,x+win1*1.2,x-win1*1.2],[y-win1*1.2,y+win1*1.2,y+win1*1.2,y-win1*1.2,y-win1*1.2],'-c')
                        if self.sst.ix[t_i,y,x]>thr_sst or np.isnan(self.sst.ix[t_i,y,x]):
                            ax.plot([x-win1,x-win1,x+win1,x+win1,x-win1],[y-win1,y+win1,y+win1,y-win1,y-win1],'--y')

            plt.suptitle(str(dates[t_i]))
            plt.tight_layout()
            plt.savefig(self.working_dir+str(t_i)+'.png', bbox_inches = 'tight')

    def detect(self,win1=7,win2=4,thr_wind=13,thr_sst=26):
        self.detected=[]
        self.detected_i=[]
        for t,t_i in zip(self.time[self.time_i],self.time_i):
            y_p,x_p = local_extrem(self.mslp[t,:,:].values,'min',threshold=101500,neighborhood_size=10)
            for y,x in zip(y_p,x_p):
                for y_,x_ in zip(y_v,x_v):
                    if x_<x+win1 and x_>x-win1 and y_<y+win1 and y_>y-win1:
                        if np.max(self.wind.ix[t_i,max(0,y-win1):min(y+win1,len(self.lat)-1),max(x-win1,0):min(x+win1,len(self.lon)-1)])>thr_wind:
                            if self.sst.ix[t_i,y,x]>thr_sst or np.isnan(self.sst.ix[t_i,y,x]):
                                self.detected_i.append((t_i,y,x))
                                self.detected.append((t,self.lat[y],self.lon[x]))


        '''implement temperature condition!!!!
        asdasd
        '''

# # read ibtracks
# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')

found_tracks={}
for year in [2008]:
    start = time.time()
    # read ERA interim
    wind_nc=da.read_nc('data/ERA_dayX/atl_10mWind_'+str(year)+'.nc')
    wind=wind_nc['ws']
    mslp=da.read_nc('data/ERA_dayX/atl_MSL_'+str(year)+'.nc')['MSL']
    sst=da.read_nc('data/ERA_dayX/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
    dates=[num2date(t,units = wind_nc.axes['time'].units,calendar = wind_nc.axes['time'].calendar) for t in wind.time]

    tc_sel=TC.ix[np.where((TC['season']==year) & (TC['basin'][:,0]==0))[0]]
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)

    working_dir='plots/detected/'+str(year)+'_dayX/'
    os.system('mkdir '+working_dir)
    found_tracks[year]=tc_tracks(wind,mslp,sst,dates,year,tc_sel,working_dir=working_dir,time_steps=range(90,120))
    #found_tracks[year].detect_TC(wind_thresh=15,sst_thresh=26,plot=True)
    found_tracks[year].detect()
    found_tracks[year].plot_season()
    found_tracks[year].plot_surrounding(win2=5)
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)
