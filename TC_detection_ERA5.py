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


try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import *

class tc_tracks(object):
    def __init__(self,VO,Wind10,MSLP,SST,T,dates,year,tc_sel,plot_dir,time_steps=None):
        self._year=year
        self._plot_dir=plot_dir
        # input fields
        self._lat=U.lat
        self._lon=U.lon

        if time_steps is None:
            time_steps=range(len(U.time))
        self._time_i=time_steps
        self._time=U.time
        self._dates=dates
        self._yr_frac=[toYearFraction(dd) for dd in self._dates]

        self._VO=VO
        self._Wind10=Wind10
        self._T=T.mean(axis=1)
        self._MSLP=MSLP
        self._SST=SST


        # initialize outputs
        self._id=0
        self._tcs={}
        # ibtracks
        self._tc_sel=tc_sel
        tmp_time=tc_sel['source_time'].values
        self._tc_time=tmp_time.copy()*np.nan
        for i in range(tmp_time.shape[0]):
            for j in range(tmp_time.shape[1]):
                if np.isfinite(tmp_time[i,j]):
                    self._tc_time[i,j]=toYearFraction(num2date(tmp_time[i,j],units = 'days since 1858-11-17 00:00:00'))
        self._tc_lat=tc_sel['lat_for_mapping'].values
        self._tc_lon=tc_sel['lon_for_mapping'].values
        self._tc_lon[self._tc_lon<0]+=360
        self._tc_intens=tc_sel['source_wind']
        # init map that is used for plotting
        fig,self._ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
        self._m = Basemap(ax=self._ax,llcrnrlon=min(self._lon),urcrnrlon=max(self._lon),llcrnrlat=min(self._lat),urcrnrlat=max(self._lat),resolution="l",projection='cyl')
        self._m.drawmapboundary(fill_color='1.')
        self._m.drawmapboundary(fill_color='#99ffff')
        self._m.fillcontinents(color='#cc9966',lake_color='#99ffff')
        self._m.drawcoastlines(linewidth=0.3)
        self._m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
        self._m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
        for storm in range(len(self._tc_sel.storm)):
            self._m.plot(self._tc_lon[storm,:],self._tc_lat[storm,:],color='gray')

    def set_thresholds(self,thr_wind,thr_sst,thr_vort,thr_mslp,thr_ta,win1,win2,win_step,neighborhood_size):
        self._thr_wind=thr_wind
        self._thr_sst=thr_sst
        self._thr_vort=thr_vort
        self._thr_mslp=thr_mslp
        self._thr_T=thr_ta
        self._win1=win1
        self._win2=win2
        self._win_step=win_step
        self._neighborhood_size=neighborhood_size

    def plot_on_map(self,x,y,**kwargs):
        if isinstance(x,np.ndarray) or isinstance(x,list) or isinstance(x,da.core.dimarraycls.DimArray):
            x=self._lon[[int(xx) for xx in x]]
            y=self._lat[[int(yy) for yy in y]]
        elif isinstance(x,np.float64) or isinstance(x,int):
            x=self._lon[int(x)]
            y=self._lat[int(y)]
        return self._m.plot(x,y,**kwargs)

    def points_in_box(self,box):
        points=[]
        for i in range(box[0],box[1]):
            for j in range(box[2],box[3]):
                points.append((i,j))
        return points

    def get_box(self,y,x,window):
        y_min=max(0,y-window)
        y_max=min(len(self._lat),y+window+1)
        x_min=max(0,x-window)
        x_max=min(len(self._lon),x+window+1)
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
        os.system('mkdir '+self._plot_dir+'surrounding')
        for t in self._time_i:
            plt.close()
            fig,axes = plt.subplots(nrows=2,ncols=2)
            axes=axes.flatten()

            ax=axes[0]; ax.set_title('rel. Vorticity')
            y_v,x_v = self.local_max(self._VO.values[t,:,:],threshold=self._thr_vort,neighborhood_size=self._neighborhood_size)
            ax.plot(x_v,y_v,'or')

            im=ax.imshow(self._VO.values[t,:,:],vmin=-9.5*10**(-5),vmax=0.0002,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[1]; ax.set_title('mean sea level pressure')
            y_p,x_p = self.local_min(self._MSLP.values[t,:,:],threshold=self._thr_mslp,neighborhood_size=self._neighborhood_size)
            ax.plot(x_p,y_p,'or')

            im=ax.imshow(self._MSLP.values[t,:,:],vmin=100360,vmax=103000,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[2]; ax.set_title('10m wind speed')
            im=ax.imshow(self._Wind10.values[t,:,:],vmin=0,vmax=13,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            # add ibtracks info
            obs_tc=np.where(abs(self._tc_time-self._yr_frac[t])<0.001)
            if len(obs_tc[0])>0:
                for oo in range(len(obs_tc[0])):
                    if np.isfinite(self._tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]):
                        for ax in axes:
                            ax.plot(np.argmin(abs(self._lon-self._tc_lon[obs_tc[0][oo],obs_tc[1][oo]])),np.argmin(abs(self._lat-self._tc_lat[obs_tc[0][oo],obs_tc[1][oo]])),color=get_tc_color(self._tc_intens.ix[obs_tc[0][oo],obs_tc[1][oo],0]),marker='.')


            for point in self._detect[self._detect[:,'t']==t].values.tolist():
                if point[3]:
                    axes[0].plot(point[2]+2,point[1],'*b')
                if point[4]:
                    axes[1].plot(point[2]+2,point[1],'*g')
                if point[5]:
                    axes[1].plot(point[2]+6,point[1],'*m')

            ax=axes[3]; ax.set_title('ibtracks')
            # aplot showing max wind, min mslp and max vort next to tc loc
            im=ax.imshow(self._T.values[t,:,:],vmin=260,vmax=275,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')
            # add ibtracks info
            obs_tc=np.where(abs(self._tc_time-self._yr_frac[t])<0.001)
            if len(obs_tc[0])>0:
                for storm in set(obs_tc[0]):
                    first_day=obs_tc[1][obs_tc[0]==storm][0]
                    # if np.isfinite(tc_sel['source_wind'].ix[storm,first_day,0]):
                    x_,y_=np.argmin(abs(self._lon-self._tc_lon[storm,first_day])),np.argmin(abs(self._lat-self._tc_lat[storm,first_day]))
                    box_1=self.get_box(y_,x_,self._win1)
                    box_2=self.get_box(y_,x_,self._win2)
                    ax.text(x_,y_,'wind: '+str(round(self._Wind10.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max(),01))+'\nmslp: '+str(round(self._MSLP.ix[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].min(),06))+'\nvort: '+str(self._VO.ix[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()),color='white',va='top',fontsize=7)

            plt.suptitle(str(dates[t]))
            plt.tight_layout()
            plt.savefig(self._plot_dir+'surrounding/'+str(t)+'.png', bbox_inches = 'tight')

    def plot_track(self,track):
        t=int(track.ix[0,0])
        tmp,text=[],[]
        #points=np.array(self._detected[:])
        #tmp.append(self._m.plot(points[:,2],points[:,1],'.g'))
        tmp.append(self.plot_on_map(track[:,'x'],track[:,'y'],c='k'))
        tmp.append(self.plot_on_map(track[track[:,'tc_cond']==3,:].ix[0,2],track[track[:,'tc_cond']==3,:].ix[0,1],marker='*',c='b'))
        self._ax.set_title(str(dates[t]))

        storms=[np.where(self._tc_time==zz)[0][0] for zz in self._yr_frac[max(0,t):min(t+track.shape[0]+5,len(self._yr_frac)-1)] if zz in self._tc_time]
        storms=np.where(abs(self._tc_time-self._yr_frac[t])<0.002)[0]
        for storm in set(storms):
            tmp+=tc_plot(self._m,self._tc_lon[storm,:],self._tc_lat[storm,:],self._tc_intens.ix[storm,:,0])
            last_pos=np.where(np.isfinite(self._tc_lon[storm,:]))[0][-1]
            text.append(self._ax.text(self._tc_lon[storm,last_pos],self._tc_lat[storm,last_pos],''.join(self._tc_sel['name'].ix[storm,:])))

        plt.tight_layout()
        plt.savefig(self._plot_dir+'track_'+str(self._year)+'_'+str(t)+'_'+str(self._id)+'_.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in text:
            element.remove()

    def plot_season(self,out_name=None):
        tmp=[]
        if out_name is None:
            out_name=self._plot_dir+'track_'+str(self._year)+'_found_tracks.png'

        for storm in range(len(self._tc_sel.storm)):
            tmp+=tc_plot(self._m,self._tc_lon[storm,:],self._tc_lat[storm,:],self._tc_intens.values[storm,:,0])

        for i in range(len(self._tcs)):
            tmp.append(self.plot_on_map(self._tcs[i][:,'x'],self._tcs[i][:,'y'],marker='.'))

        plt.tight_layout()
        plt.savefig(out_name)

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def combine_tracks(self,plot=True):
        detected=self._detect[self._detect[:,'tc_cond']==3].values
        mslp_found=self._detect[self._detect[:,'mslp']==1].values
        used_pos=[]
        for p in detected[detected[:,-2].argsort()[::-1],:].tolist():
            if p[7] and p[8]:
                if p not in used_pos:
                    track=[p]

                    # go foreward
                    while True:
                        p=track[-1]
                        end=True
                        candidates=[]
                        for p_1 in detected[detected[:,0]==p[0]+1,:].tolist():
                            if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._win_step:
                                candidates.append(p_1)
                                end=False

                        if end:
                            # check for points with fewer conditions fullfilled
                            for p_1 in mslp_found[mslp_found[:,0]==p[0]+1,:].tolist():
                                if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._win_step:
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
                            if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._win_step:
                                candidates.append(p_1)
                                end=False

                        if end:
                            break
                        else:
                            track=[candidates[0]]+track

                    if sum([pp in used_pos for pp in track])/float(len(track))>0.2:
                        used_pos+=track
                        track=da.DimArray(np.array(track),axes=[np.array(track)[:,0],['t','y','x','cd_mslp','cd_wind','cd_ta','cd_sst','cd_tropical','tc_cond','vort','mslp','wind']],dims=['time','z'])
                        if track[track[:,'tc_cond']==3].shape[0]>6 or track.shape[0]>10:
                            self._tcs[self._id]=track
                            if plot:    self.plot_track(track)
                            self._id+=1

    def detect(self):
        detect=np.array([[np.nan]*12])
        for t in self._time_i:
            print t
            y_p,x_p = self.local_min(self._MSLP.values[t,:,:],threshold=self._thr_mslp,neighborhood_size=self._neighborhood_size)
            p_mins=[(y,x) for y,x in zip(y_p,x_p)]
            # i vort max
            y_v,x_v = self.local_max(self._VO.values[t,:,:],threshold=self._thr_vort,neighborhood_size=self._neighborhood_size)
            for y,x in zip(y_v,x_v):
                box_1=self.get_box(y,x,self._win1)
                box_2=self.get_box(y,x,self._win2)
                points_1=self.points_in_box(box_1)
                # ii pressure min
                if sum([(yy,xx) in p_mins for yy,xx in points_1]):
                    tmp=[t,y,x,1,0,0,0,0,0,self._VO.ix[t,y,x],self._MSLP.ix[t,y,x],self._Wind10.ix[t,y,x]]
                    # iii wind speed
                    if self._Wind10.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()>self._thr_wind:
                        tmp[4]=1
                    # iv warm core
                    if self._T.ix[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()-self._T.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].mean()>self._thr_T:
                        tmp[5]=1
                    # v warm sea
                    if self._SST.ix[t,y,x]>self._thr_sst:   # or np.isnan(self._SST.ix[t,y,x])
                        tmp[6]=1
                    # vi tropical
                    if self._lat[y]<=30:
                        tmp[7]=1
                    tmp[8]=sum(tmp[3:6])
                    detect=np.concatenate((detect,np.array([tmp])))

        self._detect=da.DimArray(detect[1:,:],axes=[range(detect.shape[0]-1),['t','y','x','cd_mslp','cd_wind','cd_ta','cd_sst','cd_tropical','tc_cond','vort','mslp','wind']],dims=['ID','z'])


def date_(t):
    yr=int(t/10000.)
    mn=int((t-yr*10000)/100.)
    day=int((t-yr*10000-mn*100))
    hr=int((t-yr*10000-mn*100-day))*24
    return(datetime(yr,mn,day,hr,0))

found_tracks={}
for year in range(2016,2017):
    start = time.time()
    # read ERA interim
    # wind_nc=da.read_nc('data/ERA5/atl_'+str(year)+'_10mWind.nc')
    # Wind10=wind_nc['ws']
    # MSLP=da.read_nc('data/ERA5/atl_'+str(year)+'_surface.nc')['MSL']
    # T=da.read_nc('data/ERA5/atl_'+str(year)+'_T.nc')['T']
    # SST=da.read_nc('data/ERA5/atl_'+str(year)+'_surface.nc')['SSTK']-273.15
    # dates=[date_(t) for t in wind.time]
    #
    # U=da.read_nc('data/ERA5/atl_'+str(year)+'_UV.nc')['U']
    # V=da.read_nc('data/ERA5/atl_'+str(year)+'_UV.nc')['V']
    # VO=da.DimArray(rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.lon,U.lat),axes=[U.time,U.lat,U.lon],dims=['time','lat','lon'])
    #
    # # read ibtracks
    # TC_season=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['season']
    # TC_basin=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['basin']
    # tc_sel=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc').ix[np.where((TC_season==year) & (TC_basin[:,0]==0))[0]]
    # elapsed = time.time() - start;  print('Elapsed %.3f seconds.' % elapsed)

    plot_dir='plots/detection/'+str(year)+'_ERA5/'
    os.system('mkdir '+plot_dir)
    found_tracks[year]=tc_tracks(Wind10=Wind10,MSLP=MSLP,SST=SST,VO=VO,T=T,dates=dates,year=year,tc_sel=tc_sel,plot_dir=plot_dir)#,time_steps=range(470,520))
    found_tracks[year].set_thresholds(thr_wind=17,thr_vort=1*10**(-5),thr_mslp=101500,thr_ta=0,thr_sst=26.5,win1=6,win2=10,win_step=10,neighborhood_size=8)
    found_tracks[year].detect(); saved_detect=found_tracks[year]._detect
    found_tracks[year]._detect=saved_detect
    #found_tracks[year].plot_surrounding()

    start = time.time()
    found_tracks[year].combine_tracks()
    found_tracks[year].plot_season('plots/detection/season_'+str(year)+'.png')
    elapsed = time.time() - start;  print('Elapsed %.3f seconds.' % elapsed)
