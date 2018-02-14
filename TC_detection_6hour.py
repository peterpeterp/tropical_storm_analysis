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
        self.ta=ta.mean(axis=1)
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

    def set_thresholds(self,thr_wind,thr_sst,thr_vort,thr_mslp,thr_ta,win1,win2,win_step,neighborhood_size):
        self.thr_wind=thr_wind
        self.thr_sst=thr_sst
        self.thr_vort=thr_vort
        self.thr_mslp=thr_mslp
        self.thr_ta=thr_ta
        self.win1=win1
        self.win2=win2
        self.win_step=win_step
        self.neighborhood_size=neighborhood_size

    def time_ind(self,t):
        return np.where(self.time==t)[0][0]

    def plot_on_map(self,x,y,**kwargs):
        if isinstance(x,np.ndarray) or isinstance(x,list) or isinstance(x,da.core.dimarraycls.DimArray):
            x=self.lon[[int(xx) for xx in x]]
            y=self.lat[[int(yy) for yy in y]]
        elif isinstance(x,np.float64) or isinstance(x,int):
            x=self.lon[int(x)]
            y=self.lat[int(y)]
        return self.m.plot(x,y,**kwargs)

    def get_surrounding(self,y,x,window):
        y_=range(max(0,y-window),min(len(self.lat),y+window+1))
        x_=range(max(0,x-window),min(len(self.lon),x+window+1))
        return y_,x_

    def points_in_box(self,box):
        points=[]
        for i in range(box[0],box[1]):
            for j in range(box[2],box[3]):
                points.append((i,j))
        return points

    def get_box(self,y,x,window):
        y_min=max(0,y-window)
        y_max=min(len(self.lat),y+window+1)
        x_min=max(0,x-window)
        x_max=min(len(self.lon),x+window+1)
        return (y_min,y_max,x_min,x_max)

    def local_max(self,data,threshold,neighborhood_size):
        y,x=np.where(data>threshold)
        ex_x,ex_y=[],[]
        for y_,x_ in zip(y,x):
            box=self.get_box(y_,x_,neighborhood_size)
            if data[y_,x_]==np.max(data[box[0]:box[1],box[2]:box[3]]):
                ex_y.append(y_)
                ex_x.append(x_)
        return ex_y,ex_x

    def local_min(self,data,threshold,neighborhood_size):
        y,x=np.where(data<threshold)
        ex_x,ex_y=[],[]
        for y_,x_ in zip(y,x):
            box=self.get_box(y_,x_,neighborhood_size)
            if data[y_,x_]==np.min(data[box[0]:box[1],box[2]:box[3]]):
                ex_y.append(y_)
                ex_x.append(x_)
        return ex_y,ex_x

    def plot_surrounding(self):
        os.system('mkdir '+self.working_dir+'surrounding')
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
            im=ax.imshow(self.ta.values[t_i,:,:],vmin=255,vmax=275,interpolation='none')
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
            obs_tc=np.where(self.tc_time==self.yr_frac[t_i])
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


            if t_i in self.detect.keys():
                for point in self.detect[t_i]:
                    if 'vort' in self.detect[t_i][point]:
                        ax.plot(point[1]+2,point[0],'^w')
                    if 'mslp' in self.detect[t_i][point]:
                        ax.plot(point[1]+4,point[0],'*g')
                    if 'wind' in self.detect[t_i][point]:
                        ax.plot(point[1]+6,point[0],'vm')
                    if 'ta' in self.detect[t_i][point]:
                        ax.plot(point[1]+8,point[0],'*r')

            plt.suptitle(str(dates[t_i]))
            plt.tight_layout()
            plt.savefig(self.working_dir+'surrounding/'+str(t_i)+'.png', bbox_inches = 'tight')

    def plot_track(self,track):
        t=int(track.ix[0,0])
        tmp,text=[],[]
        #points=np.array(self.detected[:])
        #tmp.append(self.m.plot(points[:,2],points[:,1],'.g'))
        tmp.append(self.plot_on_map(track[:,'x'],track[:,'y'],c='k'))
        tmp.append(self.plot_on_map(track[track[:,'tc_cond']==4,:].ix[0,2],track[track[:,'tc_cond']==4,:].ix[0,1],marker='*',c='b'))
        self.ax.set_title(str(dates[t]))

        storms=[np.where(self.tc_time==zz)[0][0] for zz in self.yr_frac[max(0,t):min(t+track.shape[0]+5,len(self.yr_frac)-1)] if zz in self.tc_time]
        storms=np.where(abs(self.tc_time-self.yr_frac[t])<0.002)[0]
        for storm in set(storms):
            tmp+=tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.ix[storm,:,0])
            last_pos=np.where(np.isfinite(self.tc_lon[storm,:]))[0][-1]
            text.append(self.ax.text(self.tc_lon[storm,last_pos],self.tc_lat[storm,last_pos],''.join(self.tc_sel['name'].ix[storm,:])))

        plt.tight_layout()
        plt.savefig(self.working_dir+'track_'+str(self.year)+'_'+str(t)+'_'+str(self.id_)+'_.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in text:
            element.remove()

    def plot_season(self,out_name=None):
        tmp=[]
        if out_name is None:
            out_name=self.working_dir+'track_'+str(self.year)+'_found_tracks.png'

        for storm in range(len(self.tc_sel.storm)):
            tmp+=tc_plot(self.m,self.tc_lon[storm,:],self.tc_lat[storm,:],self.tc_intens.values[storm,:,0])

        for i in range(len(self.tcs)):
            tmp.append(self.plot_on_map(self.tcs[i][:,'x'],self.tcs[i][:,'y'],marker='.'))

        #tmp=np.array(self.detected)
        #self.m.scatter(tmp[:,2],tmp[:,1],marker='.')

        plt.tight_layout()
        plt.savefig(out_name)

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def combine_tracks(self,plot=True):
        detected=self.detect[self.detect[:,'tc_cond']==4].values
        mslp_found=self.detect[self.detect[:,'mslp']==1].values
        used_pos=[]
        for p in detected.tolist():
            if p[7] and p[8]:
                if p not in used_pos:
                    used_pos.append(p)
                    track=[p]

                    # go foreward
                    while True:
                        p=track[-1]
                        end=True
                        candidates=[]
                        for p_1 in detected[detected[:,0]==p[0]+1,:].tolist():
                            if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self.win_step:
                                candidates.append(p_1)
                                end=False

                        if end:
                            # check for points with fewer conditions fullfilled
                            for p_1 in mslp_found[mslp_found[:,0]==p[0]+1,:].tolist():
                                if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self.win_step:
                                    candidates.append(p_1)
                                    end=False
                            if end:
                                break
                            else:
                                track.append(candidates[0])
                                used_pos.append(candidates[0])
                        else:
                            track.append(candidates[0])
                            used_pos.append(candidates[0])

                    # search origin of storm
                    while p[0]>0:
                        p=track[0]
                        end=True
                        candidates=[]
                        for p_1 in mslp_found[mslp_found[:,0]==p[0]-1,:].tolist():
                            if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self.win_step:
                                candidates.append(p_1)
                                end=False

                        if end:
                            break
                        else:
                            track=[candidates[0]]+track

                    track=da.DimArray(np.array(track),axes=[np.array(track)[:,0],['t','y','x','vort','mslp','wind','ta','sst','tropical','tc_cond']],dims=['time','z'])
                    if track[track[:,'tc_cond']==4].shape[0]>6 or track.shape[0]>10:
                        self.tcs[self.id_]=track
                        if plot:    self.plot_track(track)
                        self.id_+=1

    def detect(self):
        detect=np.array([[np.nan]*10])
        for t in self.time_i:
            y_p,x_p = self.local_min(self.mslp.values[t,:,:],threshold=self.thr_mslp,neighborhood_size=self.neighborhood_size)
            p_mins=[(y,x) for y,x in zip(y_p,x_p)]
            # i vort max
            y_v,x_v = self.local_max(self.vort.values[t,0,:,:],threshold=self.thr_vort,neighborhood_size=self.neighborhood_size)
            for y,x in zip(y_v,x_v):
                box_1=self.get_box(y,x,self.win1)
                box_2=self.get_box(y,x,self.win2)
                points_1=self.points_in_box(box_1)
                tmp=[t,y,x,1,0,0,0,0,0,0]
                # ii pressure min
                if sum([(yy,xx) in p_mins for yy,xx in points_1]):
                    tmp[4]=1
                # iii wind speed
                if self.wind.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()>self.thr_wind:
                    tmp[5]=1
                # iv warm core
                if self.ta.ix[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()-self.ta.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].mean()>self.thr_ta:
                    tmp[6]=1
                # v warm sea
                if self.sst.ix[t,y,x]>self.thr_sst:
                    tmp[7]=1
                # vi tropical
                if self.lat[y]<=30:
                    tmp[8]=1
                tmp[9]=sum(tmp[3:7])
                detect=np.concatenate((detect,np.array([tmp])))

        self.detect=da.DimArray(detect[1:,:],axes=[range(detect.shape[0]-1),['t','y','x','vort','mslp','wind','ta','sst','tropical','tc_cond']],dims=['ID','z'])


# # read ibtracks
# TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')

found_tracks={}
for year in range(2016,2018):
    start = time.time()
    # read ERA interim
    wind_nc=da.read_nc('data/ERA_6hourly/atl_10mWind_'+str(year)+'.nc')
    wind=wind_nc['ws']
    mslp=da.read_nc('data/ERA_6hourly/atl_MSL_'+str(year)+'.nc')['MSL']
    vort=da.read_nc('data/ERA_6hourly/atl_atmos_'+str(year)+'.nc')['VO']
    ta=da.read_nc('data/ERA_6hourly/atl_atmos_'+str(year)+'.nc')['T']
    sst=da.read_nc('data/ERA_6hourly/atl_interim_ocean_'+str(year)+'.nc')['SSTK']-273.15
    dates=[num2date(t,units = wind_nc.axes['time'].units,calendar = wind_nc.axes['time'].calendar) for t in wind.time]
    tc_sel=TC.ix[np.where((TC['season']==year) & (TC['basin'][:,0]==0))[0]]
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)

    working_dir='plots/detected/'+str(year)+'_6hourly/'
    os.system('mkdir '+working_dir)
    found_tracks[year]=tc_tracks(wind,mslp,sst,vort,ta,dates,year,tc_sel,working_dir=working_dir)#,time_steps=range(360,400)

    found_tracks[year].set_thresholds(thr_wind=14,thr_vort=1*10**(-4),thr_mslp=101500,thr_ta=0,thr_sst=26.5,win1=3,win2=5,win_step=6,neighborhood_size=4)
    found_tracks[year].detect(); saved_detect=found_tracks[year].detect
    #found_tracks[year].detect=saved_detect
    found_tracks[year].combine_tracks()
    found_tracks[year].plot_season()
    #found_tracks[year].plot_surrounding()
    elapsed = time.time() - start
    print('Elapsed %.3f seconds.' % elapsed)
