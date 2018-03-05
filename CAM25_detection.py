from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAM25/'
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *

class tc_tracks(object):
    def __init__(self,VO,Wind10,MSLP,SST,T,nc,identifier,working_dir,time_steps=None):
        self._identifier=identifier
        self._working_dir=working_dir
        if os.path.isdir(working_dir)==False:
            os.system('mkdir '+working_dir)
        if os.path.isdir(working_dir+'/track_surrounding')==False:
            os.system('mkdir '+working_dir+'/track_surrounding')
        if os.path.isdir(working_dir+'/track_path')==False:
            os.system('mkdir '+working_dir+'/track_path')
        if os.path.isdir(working_dir+'/track_evolution')==False:
            os.system('mkdir '+working_dir+'/track_evolution')

        # input fields
        self._lats=nc['global_latitude0']
        self._lons=nc['global_longitude0']
        self._lons[self._lons>180]-=360

        self._VO=VO
        self._Wind10=Wind10
        self._MSLP=MSLP

        self._time=nc.time0
        if time_steps is None:
            time_steps=range(len(self._time))
        self._time_i=time_steps
        self._dates=[num2date(t,units = nc.axes['time0'].units,calendar = nc.axes['time0'].calendar) for t in self._time]
        self._yr_frac=np.array([toYearFraction(dd) for dd in self._dates])

        if T is not None:
            self._T=T.mean(axis=1)
        else:
            self._T=None
        self._SST=SST

        # initialize outputs
        self._id=0
        self._tcs={}

        # tc cat dict
        self._cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
        self._cat_names={0:'tropical storm',1:'Category 1',2:'Category 2',3:'Category 3',4:'Category 4',5:'Category 5'}

    def prepare_map(self,nc):
        ''' adapted from https://github.com/matplotlib/basemap/blob/master/examples/test_rotpole.py'''

        def normalize180(lon):
            """Normalize lon to range [180, 180)"""
            lower = -180.; upper = 180.
            if lon > upper or lon == lower:
                lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
            if lon < lower or lon == upper:
                lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
            return lower if lon == upper else lon

        rlats = nc.latitude0[:]
        rlons = nc.longitude0[:]
        rlons, rlats = np.meshgrid(rlons, rlats)

        o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
        o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
        lon_0 = normalize180(o_lon_p-180.)

        # init map that is used for plotting
        plt.close('all')
        fig,self._ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))
        self._m= Basemap(ax=self._ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
                   llcrnrlat = self._lats.ix[0,0], urcrnrlat = self._lats.ix[-1,-1],\
                   llcrnrlon = self._lons.ix[0,0], urcrnrlon = self._lons.ix[-1,-1],resolution='c')
        self._plot_lon,self._plot_lat = self._m(self._lons.values,self._lats.values)
        self._m.drawmapboundary(fill_color='1.')
        self._m.drawmapboundary(fill_color='darkblue')
        self._m.fillcontinents(color='darkgreen',lake_color='darkblue')
        self._m.drawcoastlines(linewidth=0.3)
        self._m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
        self._m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
        self._ax.invert_yaxis()

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

    def tc_cat(self,z,method='pressure'):
        if method=='wind':
            if z<=64: cat= 0
            if z>64: cat= 1
            if z>82: cat= 2
            if z>95: cat= 3
            if z>112: cat= 4
            if z>136: cat= 5
            if np.isnan(z): cat= 0
            return cat
        if method=='pressure':
            if z>=1020000: cat= 0
            if z<1020000: cat= 1
            if z<98000: cat= 2
            if z<96500: cat= 3
            if z<94500: cat= 4
            if z<92000: cat= 5
            if np.isnan(z): cat= 0
            return cat

    def plot_on_map(self,m,x_in,y_in,z=None,tc_cat_method='pressure',**kwargs):
        if isinstance(x_in,np.ndarray) or isinstance(x_in,list) or isinstance(x_in,da.core.dimarraycls.DimArray):
            x=self._plot_lon[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
            y=self._plot_lat[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
        elif isinstance(x_in,np.float64) or isinstance(x_in,int) or isinstance(x_in,float):
            x=self._plot_lon[int(y_in),int(x_in)]
            y=self._plot_lat[int(y_in),int(x_in)]
        if z is not None:
            tmp=[]
            for i in range(len(x)-1):
                if np.isfinite(x[i+1]):
                    tmp.append(m.plot(x[i:i+2],y[i:i+2],color=self._cat_colors[self.tc_cat(z[i])],**kwargs))
            return tmp
        else:
            return m.plot(x,y,**kwargs)

    def points_in_box(self,box):
        points=[]
        for i in range(box[0],box[1]):
            for j in range(box[2],box[3]):
                points.append((i,j))
        return points

    def get_box(self,y,x,window):
        y_min=int(max(0,y-window))
        y_max=int(min(self._lats.shape[0],y+window+1))
        x_min=int(max(0,x-window))
        x_max=int(min(self._lats.shape[1],x+window+1))
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

    def plot_surrounding(self,time_steps=None):
        if time_steps is None:
            time_steps=self._time_i

        def normalize180(lon):
            """Normalize lon to range [180, 180)"""
            lower = -180.; upper = 180.
            if lon > upper or lon == lower:
                lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
            if lon < lower or lon == upper:
                lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
            return lower if lon == upper else lon

        rlats = nc.latitude0[:]
        rlons = nc.longitude0[:]
        rlons, rlats = np.meshgrid(rlons, rlats)

        o_lon_p = nc['rotated_pole0'].attrs['grid_north_pole_longitude']
        o_lat_p = nc['rotated_pole0'].attrs['grid_north_pole_latitude']
        lon_0 = normalize180(o_lon_p-180.)


        plt.close('all')
        fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,5))
        axes=axes.flatten()
        maps=[]
        for ax in axes:
            mm=Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
                       llcrnrlat = self._lats.ix[0,0], urcrnrlat = self._lats.ix[-1,-1],\
                       llcrnrlon = self._lons.ix[0,0], urcrnrlon = self._lons.ix[-1,-1],resolution='c')
            mm.drawcoastlines(linewidth=0.7,color='m')
            mm.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
            mm.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
            ax.invert_yaxis()
            maps.append(mm)
        #plt.tight_layout()

        for t in time_steps:
            tmp,txt=[],[]
            ax=axes[0]; ax.set_title('rel. Vorticity'); mm=maps[0]
            im=mm.imshow(self._VO[t,:,:],vmin=-9.5*10**(-5),vmax=0.0002,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')
            y_v,x_v = self.local_max(self._VO[t,:,:],threshold=self._thr_vort,neighborhood_size=self._neighborhood_size)
            tmp.append(self.plot_on_map(mm,x_v,y_v,c='r',marker='*',linestyle=''))

            ax=axes[1]; ax.set_title('mean sea level pressure'); mm=maps[1]
            im=mm.imshow(self._MSLP[t,:,:],vmin=100360,vmax=103000,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[2]; ax.set_title('10m wind speed'); mm=maps[2]
            im=mm.imshow(self._Wind10[t,:,:],vmin=0,vmax=15,interpolation='none')
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            for point in self._detected[self._detected[:,'t']==t].values.tolist():
                box_1=self.get_box(point[1],point[2],self._win1)
                box_2=self.get_box(point[1],point[2],self._win2)
                if point[3]:
                    tmp.append(self.plot_on_map(maps[1],point[2],point[1],c='b',marker='*'))
                    stats='wind: '+str(round(self._Wind10[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max(),01))+'\nmslp: '+str(round(self._MSLP[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].min()/100.,01))
                    txt.append(axes[3].text(self._plot_lon[point[1],point[2]],self._plot_lat[point[1],point[2]],stats,color='red',va='bottom',ha='right',fontsize=7))
                if point[4]:
                    tmp.append(self.plot_on_map(maps[2],point[2],point[1],c='g',marker='*'))

            ax=axes[3]; ax.set_title('10m wind [m/s] and mslp [mbar]')

            plt.suptitle(str(self._dates[t]))
            plt.savefig(self._working_dir+'track_surrounding/'+str(t)+'.png', bbox_inches = 'tight')

            # clean map
            for ax in axes:
                for imm in ax.images:
                    ax.images.remove(imm)
            for element in tmp:
                l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
            for element in txt:
                element.remove()

    def plot_track_path(self,track):
        t=int(track.ix[0,0])
        tmp,txt=[],[]
        #points=np.array(self._detecteded[:])
        #tmp.append(self._m.plot(points[:,2],points[:,1],'.g'))
        tmp.append(self.plot_on_map(self._m,track[:,'x'],track[:,'y'],c='orange'))
        tmp.append(self.plot_on_map(self._m,track[track[:,'tc_cond']==3,:].ix[:,2],track[track[:,'tc_cond']==3,:].ix[:,1],marker='*',c='m',linestyle=''))
        self._ax.set_title(str(self._dates[t]))

        plt.tight_layout()
        plt.savefig(self._working_dir+'track_path/'+str(self._identifier)+'_'+str(t)+'_'+str(self._id)+'_.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in txt:
            element.remove()

    def plot_season(self,out_name=None):
        tmp=[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_tracks.png'

        self._ax.set_title('season '+self._identifier)#

        summary={0:[],1:[],2:[],3:[],4:[],5:[]}
        for id_,track in self._tcs.items():
            track=np.array(track[np.isfinite(track[:,'t']),:],dtype=np.int)
            tmp.append(self.plot_on_map(self._m,track[0,2],track[0,1],linestyle='',marker='o',c='r'))
            tmp.append(self.plot_on_map(self._m,track[:,2],track[:,1],linestyle='-'))
            tmp+=self.plot_on_map(self._m,track[:,2],track[:,1],z=self._MSLP[track[:,0],track[:,1],track[:,2]],linestyle='',marker='*')
            print(self._MSLP[track[:,0],track[:,1],track[:,2]].min())
            summary[self.tc_cat(self._MSLP[track[:,0],track[:,1],track[:,2]].min())].append(id_)

        summary.pop(0)
        print (summary)
        txt=[]
        for cat,y in zip(summary.keys(),[0.99,0.95,0.91,0.87,0.83]):
            txt.append(self._ax.text(0.005,y,self._cat_names[cat]+': '+''.join(['X']*len(summary[cat])),transform=self._ax.transAxes,color=self._cat_colors[cat],va='top',ha='left',fontsize=12))
        plt.tight_layout()
        plt.savefig(out_name)

        # clean map
        for element in txt:
            element.remove()
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def plot_track_evolution(self):
        for id_,track in self._tcs.items():
            track=track[np.isfinite(track[:,'t']),:]
            track_info=self._track_info[id_]
            track_info=track_info[np.isfinite(track_info[:,'Wind10',12,12]),:,:,:]
            print(track_info.shape,track.shape)

            # find historic storm
            max_wind=track_info[:,'Wind10',:,:].values.max(axis=(-1,-2))
            p=track.values[np.nanargmax(max_wind),:]
            tc_found=np.where(abs(self._tc_time-self._yr_frac[int(p[0])])<0.0004)
            storm=tc_found[0][np.argmin((self._tc_lat[tc_found]-self._lat[p[1]])**2+(self._tc_lon[tc_found]-self._lon[p[2]])**2)]

            plt.close('all')
            fig,axes = plt.subplots(nrows=2,ncols=2)
            axes=axes.flatten()
            asa=np.array(track[:,'t'],'int')
            axes[0].plot(self._yr_frac[asa],max_wind)
            axes[0].plot(self._tc_time[storm,:],self._tc_intens.ix[storm,:,0])

            plt.savefig(self._working_dir+'track_evolution/'+str(self._identifier)+'_'+str(id_)+'_.png')

    def gather_info_track(self,overwrite=False):
        out_file=self._working_dir+'surrounding_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._track_info=da.read_nc(out_file)
            return self._track_info

        track_info={}
        for id_,track in self._tcs.items():
            track=track[np.isfinite(track[:,'t']),:]
            info=np.zeros([track.shape[0],5,self._win2*2+1,self._win2*2+1])*np.nan
            for i,p in enumerate(track.values.tolist()):
                box_2=[int(bb) for bb in self.get_box(p[1],p[2],self._win2)]
                info[i,0,abs(p[1]-box_2[0]-12):box_2[1]-p[1]+12,abs(p[2]-box_2[2]-12):box_2[3]-p[2]+12]=self._VO[int(p[0]),box_2[0]:box_2[1],box_2[2]:box_2[3]]
                info[i,1,abs(p[1]-box_2[0]-12):box_2[1]-p[1]+12,abs(p[2]-box_2[2]-12):box_2[3]-p[2]+12]=self._MSLP[int(p[0]),box_2[0]:box_2[1],box_2[2]:box_2[3]]
                info[i,2,abs(p[1]-box_2[0]-12):box_2[1]-p[1]+12,abs(p[2]-box_2[2]-12):box_2[3]-p[2]+12]=self._Wind10[int(p[0]),box_2[0]:box_2[1],box_2[2]:box_2[3]]
                if self._T is not None:
                    info[i,3,abs(p[1]-box_2[0]-12):box_2[1]-p[1]+12,abs(p[2]-box_2[2]-12):box_2[3]-p[2]+12]=self._T.ix[int(p[0]),box_2[0]:box_2[1],box_2[2]:box_2[3]]
                if self._SST is not None:
                    info[i,4,abs(p[1]-box_2[0]-12):box_2[1]-p[1]+12,abs(p[2]-box_2[2]-12):box_2[3]-p[2]+12]=self._SST.ix[int(p[0]),box_2[0]:box_2[1],box_2[2]:box_2[3]]
            track_info[str(id_)]=da.DimArray(info,axes=[range(len(track.time)),['VO','MSLP','Wind10','T','SST'],range(self._win2*2+1),range(self._win2*2+1)],dims=['time_id','variable','lat','lon'])
        self._track_info=da.Dataset(track_info)
        self._track_info.write_nc(out_file,mode='w')

    def combine_tracks(self,plot=True,overwrite=False):
        out_file=self._working_dir+'track_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._tcs=da.read_nc(out_file)
            return self._tcs

        detected=self._detected[self._detected[:,'tc_cond']==3].values
        mslp_found=self._detected[self._detected[:,'cd_mslp']==1].values
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

                        # add info
                        info=np.zeros([len(track),4])*np.nan
                        for i,p in enumerate(track):
                            info[i,0]=self._VO[int(p[0]),int(p[1]),int(p[2])]
                            info[i,1]=self._MSLP[int(p[0]),int(p[1]),int(p[2])]
                            info[i,2]=self._Wind10[int(p[0]),int(p[1]),int(p[2])]
                            info[i,3]=self.tc_cat(self._MSLP[int(p[0]),int(p[1]),int(p[2])])

                        track=np.hstack((np.array(track),info))
                        track=da.DimArray(track,axes=[np.array(track)[:,0],['t','y','x','cd_mslp','cd_wind','cd_ta','cd_sst','cd_tropical','tc_cond','VO','MSLP','Wind10','cat']],dims=['time','z'])
                        print(track)
                        if track[track[:,'tc_cond']==3].shape[0]>6 or track.shape[0]>10:
                            self._tcs[self._identifier+'_'+str(self._id)]=track
                            if plot:    self.plot_track_path(track)
                            self._id+=1

        self._tcs=da.Dataset(self._tcs)
        self._tcs.write_nc(out_file,mode='w')
        return self._tcs

    def detect(self,overwrite=False):
        out_file=self._working_dir+'detected_positions.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected

        detect=np.array([[np.nan]*9])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*7]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()
            # i vort max
            y_v,x_v = self.local_max(self._VO[t,:,:],threshold=self._thr_vort,neighborhood_size=self._neighborhood_size)
            for y,x in zip(y_v,x_v):
                box_1=self.get_box(y,x,self._win1)
                tmp=self._MSLP[t,box_1[0]:box_1[1],box_1[2]:box_1[3]]
                y,x=np.where(tmp==tmp.min()); y,x=box_1[0]+y[0],box_1[2]+x[0]
                box_2=self.get_box(y,x,self._win2)
                # ii relative pressure min
                if self._MSLP[t,y,x]==self._MSLP[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].min():
                    box_1=self.get_box(y,x,self._win1)
                    box_2=self.get_box(y,x,self._win2)
                    tmp=[t,y,x,1,0,0,0,0,0]
                    # iii wind speed
                    if self._Wind10[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()>self._thr_wind:
                        tmp[4]=1
                    # iv warm core
                    if self._T is None:
                        tmp[5]=1
                    elif self._T.ix[t,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()-self._T.ix[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].mean()>self._thr_T:
                        tmp[5]=1
                    # v warm sea
                    if self._SST is None:
                        tmp[6]=1
                    elif self._SST.ix[t,y,x]>self._thr_sst:   # or np.isnan(self._SST.ix[t,y,x])
                        tmp[6]=1
                    # vi tropical
                    if self._lats.ix[y,x]<=30:
                        tmp[7]=1
                    tmp[8]=sum(tmp[3:6])
                    detect=np.concatenate((detect,np.array([tmp])))

        self._detected=da.DimArray(np.array(detect[1:,:],dtype=np.int),axes=[range(detect.shape[0]-1),['t','y','x','cd_mslp','cd_wind','cd_ta','cd_sst','cd_tropical','tc_cond']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('done')
        return self._detected


try:
    identifiers=[sys.argv[1]]
except:
    identifiers=[ff.split('_')[-3] for ff in glob.glob(data_path+'/item3225_daily_mean/item3225_daily*')]

for identifier in identifiers:
    start = time.time()
    print('*** started run '+identifier+' ***')
    MSLP=ndimage.gaussian_filter(da.read_nc(data_path+'item16222_daily_mean/item16222_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item16222_daily_mean'].ix[:,0,1:,:],sigma=(0,2,2))
    nc=da.read_nc(data_path+'item3225_daily_mean/item3225_daily_mean_'+identifier+'_2017-06_2017-10.nc')
    U=da.read_nc(data_path+'item3225_daily_mean/item3225_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item3225_daily_mean'].ix[:,0,:,:]
    V=da.read_nc(data_path+'item3226_daily_mean/item3226_daily_mean_'+identifier+'_2017-06_2017-10.nc')['item3226_daily_mean'].ix[:,0,:,:]
    VO=ndimage.gaussian_filter(rel_vort(U.values[:,:,:],V.values[:,:,:],U.latitude0,U.longitude0),sigma=(0,1,1))
    Wind10=np.array(np.sqrt(U**2+V**2))


    working_dir='detection/CAM25/'+str(identifier)+'_CAM25/'
    elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
    found_tcs=tc_tracks(Wind10=Wind10,MSLP=MSLP,SST=None,VO=VO,T=None,nc=nc,identifier=identifier,working_dir=working_dir)#,time_steps=range(470,520))
    found_tcs.prepare_map(nc)
    elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)
    found_tcs.set_thresholds(thr_wind=15,thr_vort=5*10**(-5),thr_mslp=101500,thr_ta=0,thr_sst=26.5,win1=7,win2=12,win_step=20,neighborhood_size=8)
    found_tcs.detect(overwrite=False)
    found_tcs.combine_tracks(overwrite=False)
    #found_tcs.gather_info_track(overwrite=False)
    #track_info,track=found_tcs.plot_track_evolution()
    found_tcs.plot_season()
    found_tcs.plot_surrounding(range(100,130))#; convert -delay 50 track_surrounding/{94..127}* TC.gif
    elapsed = time.time() - start;  print('Done with plotting %.3f seconds.' % elapsed)
    print('memory in use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10.**6))
