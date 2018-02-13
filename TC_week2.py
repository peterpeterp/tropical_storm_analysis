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

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches

os.chdir('/Users/peterpfleiderer/Documents/Projects/Scripts/allgemeine_scripte')
from support_functions import *

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')



class tc_tracks(object):
    def __init__(self,wind,mslp,sst,vort,ta,dates,year,tc_sel,working_dir,time_steps=None):
        self.year=year
        self.working_dir=working_dir
        # input fields
        if time_steps is None:
            time_steps=range(len(wind.time))
        self.time_i=time_steps
        self.wind=wind
        self.vort=vort
        self.ta=ta
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

    def set_thresholds(self,thr_wind,thr_sst,thr_vort,thr_mslp,win1,win2,win3,neighborhood_size):
        self.thr_wind=thr_wind
        self.thr_sst=thr_sst
        self.thr_vort=thr_vort
        self.thr_mslp=thr_mslp
        self.win1=win1
        self.win2=win2
        self.win3=win3
        self.neighborhood_size=neighborhood_size

    def time_ind(self,t):
        return np.where(self.time==t)[0][0]

    def get_surrounding(self,y,x,window,latlon=False):
        if latlon:
            y=np.where(self.lat==y)[0][0]
            x=np.where(self.lon==x)[0][0]
        y_=range(max(0,y-window),min(len(self.lat)-1,y+window))
        x_=range(max(0,x-window),min(len(self.lon)-1,x+window))
        return y_,x_

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
            if data[y_,x_]==np.min(data[max([0,y_-neighborhood_size]):min([y_+neighborhood_size,len(self.lat)-1]),max([0,x_-neighborhood_size]):min([x_+neighborhood_size,len(self.lon)-1])]):
                ex_y.append(y_)
                ex_x.append(x_)
        return ex_y,ex_x

    def plot_surrounding(self,win1=7,win2=2,thr_wind=13,thr_sst=26):
        for t,t_i in zip(self.time[self.time_i],self.time_i):
            plt.close()
            fig,axes = plt.subplots(nrows=2,ncols=2)
            axes=axes.flatten()

            ax=axes[0]
            ax.set_title('10m wind speed')
            data=self.wind.values[t_i,:,:]

            y_w,x_w = self.local_max(data,threshold=self.thr_wind,neighborhood_size=self.neighborhood_size)
            im=ax.imshow(data,vmin=0,vmax=13,interpolation='none')
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_w,y_w,color='r',s=50,facecolors='none')

            ax=axes[1]
            ax.set_title('rel. Vorticity')
            data=self.vort.values[t_i,0,:,:]
            y_v,x_v = self.local_max(data,threshold=self.thr_vort,neighborhood_size=self.neighborhood_size)
            im=ax.imshow(data,vmin=-9.5*10**(-5),vmax=0.0002,interpolation='none')
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_v,y_v,color='r',s=50,facecolors='none')

            ax=axes[2]
            # aplot showing max wind, min mslp and max vort next to tc loc
            im=ax.imshow(self.ta.values[t_i,1,:,:],vmin=255,vmax=275,interpolation='none')
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
            y_p,x_p = self.local_min(data,threshold=self.thr_mslp,neighborhood_size=self.neighborhood_size)
            im=ax.imshow(data,vmin=100360,vmax=103000,interpolation='none')
            im.set_cmap('bone')
            ax.autoscale(False)
            ax.axis('off')
            ax.scatter(x_p,y_p,marker='o',color='r',s=50,facecolors='none')
            ax.scatter(x_v,y_v,marker='*',color='m')
            ax.scatter(x_w,y_w,marker='*',color='c')

            p_mins=[(y,x) for y,x in zip(y_p,x_p)]

            # i vort max
            y_v,x_v = self.local_max(self.vort.values[t_i,0,:,:],threshold=self.thr_vort,neighborhood_size=self.neighborhood_size)
            for y,x in zip(y_v,x_v):
                # ax.plot([x-self.win2,x-self.win2,x+self.win2,x+self.win2,x-self.win2],[y-self.win2,y+self.win2,y+self.win2,y-self.win2,y-self.win2],'-b')
                surrounding=[]
                for i in range(max(0,y-self.win1),min(len(self.lat)-1,y+self.win1)):
                    for j in range(max(0,x-self.win1),min(len(self.lon)-1,x+self.win1)):
                        surrounding.append((i,j))

                y_win1=range(max(0,y-self.win1),min(len(self.lat)-1,y+self.win1))
                x_win1=range(max(0,x-self.win1),min(len(self.lon)-1,x+self.win1))
                y_win2=range(max(0,y-self.win2),min(len(self.lat)-1,y+self.win2))
                x_win2=range(max(0,x-self.win2),min(len(self.lon)-1,x+self.win2))
                # ii pressure min
                if sum([(yy,xx) in p_mins for yy,xx in surrounding]):
                    ax.plot([x-self.win2*1.1,x-self.win2*1.1,x+self.win2*1.1,x+self.win2*1.1,x-self.win2*1.1],[y-self.win2*1.1,y+self.win2*1.1,y+self.win2*1.1,y-self.win2*1.1,y-self.win2*1.1],'-g')

                    # iii wind speed
                    if self.wind.ix[t_i,y_win2,x_win2].max()>14:
                        ax.plot([x-self.win2*1.2,x-self.win2*1.2,x+self.win2*1.2,x+self.win2*1.2,x-self.win2*1.2],[y-self.win2*1.2,y+self.win2*1.2,y+self.win2*1.2,y-self.win2*1.2,y-self.win2*1.2],'-c')

                        # iv warm core
                        if self.wind.ix[t_i,y_win1,x_win1].max()-self.wind.ix[t_i,y_win2,x_win2].mean()>2:
                            ax.plot([x-self.win2*0.9,x-self.win2*0.9,x+self.win2*0.9,x+self.win2*0.9,x-self.win2*0.9],[y-self.win2*0.9,y+self.win2*0.9,y+self.win2*0.9,y-self.win2*0.9,y-self.win2*0.9],'-r')


            plt.suptitle(str(dates[t_i]))
            plt.tight_layout()
            plt.savefig(self.working_dir+str(t_i)+'.png', bbox_inches = 'tight')

    def save_track(self,track,t_i,plot=True):
        if track.shape[0]>1:
            if self.sst[track[0,0],self.lat[np.argmin(abs(self.lat-track[0,1]))],self.lon[np.argmin(abs(self.lon-track[0,2]))]]>self.thr_sst:
                self.tcs[str(self.id_)]=da.DimArray(track,axes=[track[:,0],['time','lat','lon','t_i']],dims=['time','pos'])
                if plot:    self.plot_track(track,t_i)
                self.id_+=1

    def plot_track(self,track,t_i):
        print t_i
        tmp,text=[],[]
        points=np.array(self.detected[:])
        tmp.append(self.m.plot(points[:,2],points[:,1],'.g'))
        tmp.append(self.m.plot(track[:,2],track[:,1],'-k'))
        tmp.append(self.m.plot(track[0,2],track[0,1],'*b'))
        self.ax.set_title(str(dates[t_i]))

        storms=[np.where(self.tc_time==zz)[0][0] for zz in self.yr_frac[max(0,t_i):min(t_i+track.shape[0]+5,len(self.yr_frac)-1)] if zz in self.tc_time]
        storms=np.where(abs(self.tc_time-self.yr_frac[t_i])<0.002)[0]
        for storm in set(storms):
            tmp+=tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.ix[storm,:,0])
            last_pos=np.where(np.isfinite(self.tc_lon[storm,:]))[0][-1]
            text.append(self.ax.text(self.tc_lon[storm,last_pos],self.tc_lat[storm,last_pos],''.join(self.tc_sel['name'].ix[storm,:])))

        plt.tight_layout()
        plt.savefig(self.working_dir+'track_'+str(self.year)+'_'+str(t_i)+'_'+str(self.id_)+'_.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in text:
            element.remove()

    def plot_season(self,out_name=None):
        if out_name is None:
            out_name=self.working_dir+'track_'+str(self.year)+'_found_tracks.png'

        for storm in range(len(self.tc_sel.storm)):
            tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.values[storm,:,0])

        for i in range(len(self.tcs)):
            self.m.plot(self.tcs[str(i)].values[:,2],self.tcs[str(i)].values[:,1],marker='.')

        tmp=np.array(self.detected)

        self.m.scatter(tmp[:,2],tmp[:,1],marker='.')

        plt.tight_layout()
        plt.savefig(out_name)

    def extend_tracks(self):
        points=np.array(self.detected[:])
        for id_ in self.tcs.keys():
            tmp=self.tcs[id_].values.copy()
            t_i = tmp[0,3]
            while t_i>0:
                t_i-=1
                y_p,x_p = self.local_min(self.mslp.values[t_i,:,:],threshold=self.thr_mslp,neighborhood_size=self.neighborhood_size)
                p_mins=[(y,x) for y,x in zip(y_p,x_p)]

                y_0,x_0=np.where(self.lat==tmp[0,1])[0][0],np.where(self.lon==tmp[0,2])[0][0]
                window3 = Polygon([(y_0-self.win3, x_0-self.win3), (y_0-self.win3, x_0+self.win3), (y_0+self.win3, x_0+self.win3), (y_0+self.win3, x_0-self.win3)])
                # ii pressure min
                if sum([window3.contains(Point(yy,xx)) for yy,xx in p_mins]):
                    y_win3,x_win3=self.get_surrounding(y_0,x_0,self.win3)
                    y_1,x_1=np.argmin(self.mslp.ix[t_i,y_win3,x_win3])
                    y_win2,x_win2=self.get_surrounding(y_1,x_1,self.win2,latlon=True)
                    if self.wind.ix[t_i,y_win2,x_win2].max()>10:
                        tmp=np.vstack((np.array([self.time[t_i],y_1,x_1,t_i]),tmp))

                else:
                    break

            t_i = tmp[-1,3]
            while t_i<len(self.time)-1:
                t_i+=1
                y_p,x_p = self.local_min(self.mslp.values[t_i,:,:],threshold=self.thr_mslp,neighborhood_size=self.neighborhood_size)
                p_mins=[(y,x) for y,x in zip(y_p,x_p)]

                y_0,x_0=np.where(self.lat==tmp[-1,1])[0][0],np.where(self.lon==tmp[-1,2])[0][0]
                window3 = Polygon([(y_0-self.win3, x_0-self.win3), (y_0-self.win3, x_0+self.win3), (y_0+self.win3, x_0+self.win3), (y_0+self.win3, x_0-self.win3)])
                # ii pressure min
                if sum([window3.contains(Point(yy,xx)) for yy,xx in p_mins]):
                    print sum([window3.contains(Point(yy,xx)) for yy,xx in p_mins])
                    y_win3,x_win3=self.get_surrounding(y_0,x_0,self.win3)
                    y_1,x_1=np.argmin(self.mslp.ix[t_i,y_win3,x_win3])
                    y_win2,x_win2=self.get_surrounding(y_1,x_1,self.win2,latlon=True)
                    if self.wind.ix[t_i,y_win2,x_win2].max()>10:
                        tmp=np.vstack((tmp,np.array([self.time[t_i],y_1,x_1,t_i])))

                else:
                    break

            self.tcs[id_]=da.DimArray(tmp,axes=[tmp[:,0],['time','lat','lon','t_i']],dims=['time','pos'])

                    # if tmp.shape[0]>1:
                    #     y_00,x_00=np.where(self.lat==tmp[-2,1])[0][0],np.where(self.lon==tmp[-2,2])[0][0]
                    #     y_1=2*y_0-y_00
                    #     x_1=2*x_0-x_00
                    #     y_win1,x_win1=self.get_surrounding(y_1,x_1,self.win1)
                        # window1 = Polygon([(y_1-self.win1, x_1-self.win1), (y_1-self.win1, x_1+self.win1), (y_1+self.win1, x_1+self.win1), (y_1+self.win1, x_1-self.win1)])
                        # for tr_i in tcs:
                        #     tr=tcs[tr_i].ix[0,:].values
                        #     print tr,t_i,y_1,x_1
                        #     if tr[3]==t_i and window1.contains(Point(np.where(self.lat==tr[1])[0][0],np.where(self.lon==tr[2])[0][0])):
                        #         tmp=np.vstack((tmp,tcs[tr_i].values))
                        #         print 'hey'

                        # if len(y_win1)>0 and len(x_win1)>0:
                        #     if self.wind.ix[t_i,y_win1,x_win1].max()>10:
                        #         y_win1,x_win1=self.get_surrounding(y_1,x_1,self.win1)
                        #         if len(y_win1)>0 and len(x_win1)>0:
                        #             y_2,x_2=np.argmin(self.mslp.ix[t_i,y_win1,x_win1])
                        #             tmp=np.vstack((tmp,np.array([self.time[t_i],y_2,x_2,t_i])))



    def combine_tracks(self):
        points=np.array(self.detected[:])
        used_pos=[]
        a=0
        for p in points.tolist():
            a+=1
            #if point not already recorded:
            step=0
            #tmp=[]
            #tmp.append(self.m.plot(p[2],p[1],'^k'))
            if p not in used_pos:
                used_pos.append(p)
                track=[]
                all_locs={step:[p]}

                for p_1 in points[points[:,0]==p[0],:].tolist():
                    if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<3:
                        all_locs[step].append(p_1)
                        used_pos.append(p_1)
                track.append(np.array(all_locs[step]).mean(axis=0).tolist()+[self.time_ind(p_1[0])])

                while True:
                    step+=1
                    all_locs[step]=[]
                    p=track[-1]
                    #tmp.append(self.m.plot(p[2],p[1],'^k'))
                    tmp_points=points[points[:,0]==p[0]+24,:]
                    #tmp.append(self.m.plot(tmp_points[:,2],tmp_points[:,1],'*g'))
                    for p_1 in points[points[:,0]==p[0]+24,:].tolist():
                        if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<9:
                            if p_1 not in all_locs[step]:
                                all_locs[step].append(p_1)
                                used_pos.append(p_1)
                                #tmp.append(self.m.plot(p_1[2],p_1[1],'*m'))

                    #plt.savefig(self.working_dir+'__'+str(np.where(self.time==p[0])[0][0])+'_'+str(a)+'_'+str(step)+'_.png')
                    # clean map
                    #for element in tmp:
                    #    if len(element)>0:
                    #        l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

                    if len(all_locs[step])!=0:
                        track.append(np.array(all_locs[step]).mean(axis=0).tolist()+[self.time_ind(p_1[0])])
                    if len(all_locs[step])==0:
                        break

                self.save_track(np.array(track),np.where(self.time==track[0][0])[0][0])

    def detect(self,plot=True):
        self.detected=[]
        self.detected_i=[]
        for t,t_i in zip(self.time[self.time_i],self.time_i):
            y_p,x_p = self.local_min(self.mslp.values[t_i,:,:],threshold=self.thr_mslp,neighborhood_size=self.neighborhood_size)
            p_mins=[(y,x) for y,x in zip(y_p,x_p)]
            # i vort max
            y_v,x_v = self.local_max(self.vort.values[t_i,0,:,:],threshold=self.thr_vort,neighborhood_size=self.neighborhood_size)
            for y,x in zip(y_v,x_v):
                surrounding=[]
                for i in range(max(0,y-self.win1),min(len(self.lat)-1,y+self.win1)):
                    for j in range(max(0,x-self.win1),min(len(self.lon)-1,x+self.win1)):
                        surrounding.append((i,j))


                y_win1=range(max(0,y-self.win1),min(len(self.lat)-1,y+self.win1))
                x_win1=range(max(0,x-self.win1),min(len(self.lon)-1,x+self.win1))
                y_win2=range(max(0,y-self.win2),min(len(self.lat)-1,y+self.win2))
                x_win2=range(max(0,x-self.win2),min(len(self.lon)-1,x+self.win2))
                # ii pressure min
                if sum([(yy,xx) in p_mins for yy,xx in surrounding]):

                    # iii wind speed
                    if self.wind.ix[t_i,y_win2,x_win2].max()>self.thr_wind:

                        # iv warm core
                        if self.wind.ix[t_i,y_win1,x_win1].max()-self.wind.ix[t_i,y_win2,x_win2].mean()>0:

                            self.detected_i.append((t_i,y,x))
                            self.detected.append((t,self.lat[y],self.lon[x]))




# # read ibtracks
# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')

found_tracks={}
for year in [2008]:
    start = time.time()
    # read ERA interim
    wind_nc=da.read_nc('data/ERA_dayX/atl_10mWind_'+str(year)+'.nc')
    wind=wind_nc['ws']
    mslp=da.read_nc('data/ERA_dayX/atl_MSL_'+str(year)+'.nc')['MSL']
    vort=da.read_nc('data/ERA_dayX/atl_atmos_'+str(year)+'.nc')['VO']
    ta=da.read_nc('data/ERA_dayX/atl_atmos_'+str(year)+'.nc')['T']
    sst=da.read_nc('data/ERA_dayX/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
    dates=[num2date(t,units = wind_nc.axes['time'].units,calendar = wind_nc.axes['time'].calendar) for t in wind.time]

    tc_sel=TC.ix[np.where((TC['season']==year) & (TC['basin'][:,0]==0))[0]]
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)

    working_dir='plots/detected/'+str(year)+'_dayX/'
    os.system('mkdir '+working_dir)
    found_tracks[year]=tc_tracks(wind,mslp,sst,vort,ta,dates,year,tc_sel,working_dir=working_dir,time_steps=range(90,100))
    found_tracks[year].set_thresholds(thr_wind=14,thr_vort=1*10**(-4),thr_mslp=101500,thr_sst=26,win1=3,win2=5,win3=7,neighborhood_size=4)
    found_tracks[year].detect()
    found_tracks[year].combine_tracks()
    found_tracks[year].extend_tracks()
    found_tracks[year].plot_season()
    #found_tracks[year].plot_surrounding()
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)
