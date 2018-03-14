from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max

sns.set_palette(sns.color_palette("plasma"))

sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import * ; reload(sys.modules['TC_support'])

'''
Distances in km not in grid cells - conversion within init
Murakami's detection scheme
'''


class tc_tracks(object):
    def __init__(self,VO,Wind10,MSLP,MSLP_smoothed,SST,T,lats,lons,time_,dates,identifier,working_dir,time_steps=None):
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
        self._lats=lats
        self._lons=lons

        self._VO=VO
        self._Wind10=Wind10
        self._MSLP=MSLP
        if MSLP_smoothed is not None:
            self._MSLP_smoothed=MSLP_smoothed
        else:
            self._MSLP_smoothed=MSLP

        self._time=time_
        if time_steps is None:
            time_steps=range(len(self._time))
        self._time_i=time_steps
        self._dates=dates
        self._yr_frac=np.array([toYearFraction(dd) for dd in self._dates])

        self._T=T
        self._SST=SST

        # initialize outputs
        self._id=0
        self._tcs={}

        # tc cat dict
        self._obs_tc=False
        self._cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
        self._cat_names={0:'tropical storm',1:'Category 1',2:'Category 2',3:'Category 3',4:'Category 4',5:'Category 5'}

    def init_map(self,m,ax,plot_lat,plot_lon):
        self._m=m
        self._ax=ax
        self._plot_lat=plot_lat
        self._plot_lon=plot_lon

        # for storm in range(len(self._tc_sel.storm)):
        #     self._m.plot(self._tc_lon[storm,:],self._tc_lat[storm,:],color='gray')

    def init_obs_tcs(self,tc_sel):
        self._tc_sel=tc_sel
        tmp_time=tc_sel['source_time'].values
        self._tc_time=tmp_time.copy()*np.nan
        for i in range(tmp_time.shape[0]):
            for j in range(tmp_time.shape[1]):
                if np.isfinite(tmp_time[i,j]):
                    self._tc_time[i,j]=toYearFraction(num2date(tmp_time[i,j],units = 'days since 1858-11-17 00:00:00'))
        self._tc_lat=tc_sel['lat_for_mapping'].values
        self._tc_lon=tc_sel['lon_for_mapping'].values
        self._tc_lon[self._tc_lon>180]-=360
        self._tc_intens=np.nanmean(tc_sel['source_wind'],axis=-1)
        self._obs_tc=True

    def set_thresholds(self,thr_wind,thr_mslp,p_radius,neighborhood_size,warm_core_size,cores_distance,search_radius,min_time_steps):
        self._thr_wind=thr_wind
        self._thr_mslp=thr_mslp
        self._min_time_steps=min_time_steps

        self._p_radius=self.degree_to_step(p_radius)
        self._neighborhood_size=self.degree_to_step(neighborhood_size)
        self._warm_core_size=self.degree_to_step(warm_core_size)**2
        self._cores_distance=self.degree_to_step(cores_distance)
        self._search_radius=self.degree_to_step(search_radius)

    def degree_to_step(self,degree):
        y_step=abs(np.diff(self._lats[:,0],1).mean())
        x_step=abs(np.diff(self._lons[0,:],1).mean())
        return round(degree/(y_step+x_step)*2)

    def tc_cat(self,z,method='pressure'):
        def cat__(zz):
            if method=='wind':
                if zz<=64: cat= 0
                if zz>64: cat= 1
                if zz>82: cat= 2
                if zz>95: cat= 3
                if zz>112: cat= 4
                if zz>136: cat= 5
                if np.isnan(zz): cat= 0
                return cat
            if method=='pressure':
                if zz>=1020: cat= 0
                if zz<1020: cat= 1
                if zz<980: cat= 2
                if zz<965: cat= 3
                if zz<945: cat= 4
                if zz<920: cat= 5
                if np.isnan(zz): cat= 0
                return cat
        if isinstance(z,np.ndarray) or isinstance(z,list) or isinstance(z,da.core.dimarraycls.DimArray):
            return [cat__(zz) for zz in z]
        else:
            return cat__(z)

    def plot_on_map(self,m,x_in,y_in,cat=None,latlon=False,**kwargs):
        if latlon:
            x,y=x_in,y_in
        if latlon==False:
            if isinstance(x_in,np.ndarray) or isinstance(x_in,list) or isinstance(x_in,da.core.dimarraycls.DimArray):
                x=self._plot_lon[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
                y=self._plot_lat[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
            elif isinstance(x_in,np.float64) or isinstance(x_in,int) or isinstance(x_in,float):
                x=self._plot_lon[int(y_in),int(x_in)]
                y=self._plot_lat[int(y_in),int(x_in)]
        if cat is not None:
            tmp=[]
            for i in range(len(x)-1):
                if np.isfinite(x[i+1]):
                    tmp.append(m.plot(x[i:i+2],y[i:i+2],color=self._cat_colors[cat[i]],**kwargs))
            return tmp
        else:
            return m.plot(x,y,**kwargs)

    def find_closed_contours(self,field,y,x,step=2,search_radius=30,n_contours=None,method='min'):
        '''
        finds closed contours around center
        center is supposed to be a minimum
        '''
        y_min=int(max(0,y-search_radius))
        y_max=int(min(field.shape[0],y+search_radius+1))
        x_min=int(max(0,x-search_radius))
        x_max=int(min(field.shape[1],x+search_radius+1))

        if method=='min':
            im=-field.copy()
            threshold=-field[y,x]
        else:
            im=field.copy()
            threshold=field[y,x]

        ny,nx=im.shape
        running=True
        ncont=0
        cont=np.zeros((ny, nx))
        while n_contours is None or ncont<n_contours:
            threshold-=step
            th, im_floodfill = cv2.threshold(im, threshold, 255, cv2.THRESH_BINARY_INV);
            mask = np.zeros((ny+2, nx+2), np.uint8)
            cv2.floodFill(im_floodfill, mask, (x,y), 1);
            y_,x_=np.where(mask[1:-1,1:-1]==1)
            if y_min in y_ or x_min in x_ or y_max-1 in y_ or x_max-1 in x_:
                cont[cont==0]=np.nan
                return cont,ncont
            cont=np.array(mask[1:-1,1:-1],np.float)
            ncont+=1

        cont[cont==0]=np.nan
        return cont,ncont

    def get_box(self,y,x,window):
        y_min=int(max(0,y-window))
        y_max=int(min(self._lats.shape[0],y+window+1))
        x_min=int(max(0,x-window))
        x_max=int(min(self._lats.shape[1],x+window+1))
        return (y_min,y_max,x_min,x_max)

    def plot_surrounding(self,maps,axes,time_steps=None):
        if time_steps is None:
            time_steps=self._time_i

        #plt.tight_layout()
        for t in time_steps:
            tmp,txt=[],[]
            ax=axes[0]; ax.set_title('mean sea level pressure'); mm=maps[0]
            im=mm.pcolormesh(self._plot_lon,self._plot_lat,self._MSLP[t,:,:])
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[1]; ax.set_title('temperature'); mm=maps[1]
            im=mm.pcolormesh(self._plot_lon,self._plot_lat,self._T[t,1,:,:])
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[2]; ax.set_title('10m wind speed'); mm=maps[2]
            im=mm.pcolormesh(self._plot_lon,self._plot_lat,self._Wind10[t,:,:],vmin=0,vmax=15)
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            for point in self._detected[self._detected[:,'t']==t].values.tolist():
                if point[3]:
                    tmp.append(self.plot_on_map(maps[0],point[2],point[1],c='b',marker='.'))
                    stats='wind: '+str(round(point[6],01))+'\nmslp: '+str(round(point[5],01))
                    txt.append(axes[3].text(self._plot_lon[point[1],point[2]],self._plot_lat[point[1],point[2]],stats,color='red',va='bottom',ha='right',fontsize=7))
                if point[4]:
                    tmp.append(self.plot_on_map(maps[1],point[2],point[1],c='g',marker='*'))

            ax=axes[3]; ax.set_title('10m wind [m/s] and mslp [mbar]')

            if self._obs_tc:
                obs_tc=np.where(abs(self._tc_time-self._yr_frac[t])<0.0004)
                if len(obs_tc[0])>0:
                    for oo in range(len(obs_tc[0])):
                        if np.isfinite(self._tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]):
                            tmp.append(maps[3].plot(self._tc_lon[obs_tc[0][oo],obs_tc[1][oo]],self._tc_lat[obs_tc[0][oo],obs_tc[1][oo]],color=self._cat_colors[self.tc_cat(self._tc_intens[obs_tc[0][oo],obs_tc[1][oo]],method='wind')],marker='.'))


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
        t=int(track.ix[np.nanargmin(track[:,'MSLP'].values),0])
        tmp,txt=[],[]

        if self._obs_tc:
            storms=np.where(abs(self._tc_time-self._yr_frac[t])<0.002)[0]
            for storm in set(storms):
                tmp+=self.plot_on_map(self._m,self._tc_lon[storm,:],self._tc_lat[storm,:],cat=self.tc_cat(self._tc_intens[storm,:],method='wind'),latlon=True)
                last_pos=np.where(np.isfinite(self._tc_lon[storm,:]))[0][-1]
                txt.append(self._ax.text(self._tc_lon[storm,last_pos],self._tc_lat[storm,last_pos],''.join(self._tc_sel['name'].ix[storm,:])))

        #tmp.append(self.plot_on_map(self._m,track[:,'x'],track[:,'y'],c='orange'))
        tmp+=self.plot_on_map(self._m,track[:,'x'],track[:,'y'],cat=self.tc_cat(track[:,'MSLP'].values),marker='.',linestyle='')
        self._ax.set_title(str(self._dates[t]))

        #tmp+=self.plot_on_map(self._m,self._detected[:,'x'],self._detected[:,'y'],marker='.',linestyle='',color='m')

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
            track=track[np.isfinite(track[:,'t']),:]
            tmp.append(self.plot_on_map(self._m,track.ix[0,2],track.ix[0,1],linestyle='',marker='o',c='r'))
            tmp.append(self.plot_on_map(self._m,track[:,'x'],track[:,'y'],linestyle='-',linewidth=0.5,c='r'))
            tmp+=self.plot_on_map(self._m,track[:,'x'],track[:,'y'],cat=self.tc_cat(track[:,'MSLP'].values),marker='.',linestyle='')
            summary[max(self.tc_cat(track[:,'MSLP'].values))].append(id_)

        if self._obs_tc:
            for storm in range(len(self._tc_sel.storm)):
                tmp+=self.plot_on_map(self._m,self._tc_lon[storm,:],self._tc_lat[storm,:],cat=self.tc_cat(self._tc_intens[storm,:],method='wind'),latlon=True)


        summary.pop(0)
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

    def obs_track_info(self,overwrite=False):
        out_file=self._working_dir+'obs_track_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._obs_track_info=da.read_nc(out_file)
            return self._obs_track_info


        obs_summary=np.zeros([len(self._tc_sel.storm),200,7])*np.nan
        for i,storm in enumerate(self._tc_sel.storm):
            tmp_t=self._tc_time[i,:]
            last_val=len(np.where(np.isfinite(tmp_t))[0])
            obs_summary[i,0:last_val,0]=[self.tc_cat(z,method='wind') for z in np.nanmean(self._tc_sel['source_wind'].values[i,0:last_val,:],axis=-1)]

            for t in range(last_val):
                t_=np.where(abs(self._yr_frac-self._tc_time[i,t])<0.0004)[0]
                if len(t_)!=0:
                    t_=t_[0]
                    y,x=np.argmin(abs(self._lats[:,0]-self._tc_lat[i,t])),np.argmin(abs(self._lons[0,:]-self._tc_lon[i,t]))
                    box_1=[int(bb) for bb in self.get_box(y,x,self._win1)]
                    box_2=[int(bb) for bb in self.get_box(y,x,self._win2)]
                    obs_summary[i,t,1]=self._VO[t_,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()
                    obs_summary[i,t,2]=self._MSLP[t_,box_1[0]:box_1[1],box_1[2]:box_1[3]].min()
                    obs_summary[i,t,3]=self._Wind10[t_,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()
                    obs_summary[i,t,4]=self._T[t_,0,y,x]
                    obs_summary[i,t,5]=self._T[t_,1,y,x]
                    if self._SST is not None:
                        obs_summary[i,t,6]=self._SST[t_,y,x]

        obs_summary=obs_summary[:,np.isfinite(np.nanmean(obs_summary,axis=(0,-1))),:]
        self._obs_track_info=da.DimArray(obs_summary,axes=[self._tc_sel.storm,range(obs_summary.shape[1]),['cat','VO','MSLP','Wind10','T850','T500','SST']],dims=['storm','time','variable'])

        da.Dataset({'obs_track_info':self._obs_track_info}).write_nc(out_file)

        # print summary
        sys.stdout.write('Category:\t0\t\t1\t\t2\t\t3\t\t4\t\t5'); sys.stdout.flush()
        for vari,name in zip(range(2,7),['VO','MSLP','Wind10','T850','T500','SST']):
            sys.stdout.write('\n'+name+'\t\t'); sys.stdout.flush()
            for cat in range(6):
                pos=np.where(obs_summary==cat)
                sys.stdout.write(str(np.nanmean(obs_summary[pos[0],pos[1],vari]))+'\t'); sys.stdout.flush()

        return self._obs_track_info

    def combine_tracks(self,plot=True,overwrite=False):
        out_file=self._working_dir+'track_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._tcs=da.read_nc(out_file)
            return self._tcs

        warm_core=self._detected.values
        pressure_low=self._detected[self._detected[:,'pressure_low']==1].values
        used_pos=[]
        for p in warm_core[warm_core[:,-1].argsort()[::-1],:].tolist():
            if p not in used_pos:
                track=[p]

                # go foreward
                while True:
                    p=track[-1]
                    end=True
                    candidates=[]
                    for p_1 in warm_core[warm_core[:,0]==p[0]+1,:].tolist():
                        if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._search_radius:
                            candidates.append(p_1)
                            end=False

                    if end:
                        # check for points with fewer conditions fullfilled
                        for p_1 in pressure_low[pressure_low[:,0]==p[0]+1,:].tolist():
                            if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._search_radius:
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
                    for p_1 in pressure_low[pressure_low[:,0]==p[0]-1,:].tolist():
                        if ((p[1]-p_1[1])**2+(p[2]-p_1[2])**2)**0.5<self._search_radius:
                            candidates.append(p_1)
                            end=False

                    if end:
                        break
                    else:
                        track=[candidates[0]]+track

                if sum([pp in used_pos for pp in track])/float(len(track))>0.2:
                    used_pos+=track

                    track=da.DimArray(track,axes=[np.array(track)[:,0],['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['time','z'])
                    if track.shape[0]>=12:
                        if track[track[:,'warm_core']==1].shape[0]>=8:
                            if track[track[:,'Wind10']>=self._thr_wind].shape[0]>=6:
                                if self._lats[track[track[:,'warm_core']==1].ix[0,1],0]<40:
                                    self._tcs[self._identifier+'_'+str(self._id)]=track
                                    if plot:    self.plot_track_path(track)
                                    self._id+=1

        self._tcs=da.Dataset(self._tcs)
        self._tcs.write_nc(out_file,mode='w')
        return self._tcs

    def plot_detect_summary(self,out_name=None):
        tmp=[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_positions.png'

        self._ax.set_title('season '+self._identifier)#

        detect=self._detected
        tmp.append(self.plot_on_map(self._m,detect[:,'x'],detect[:,'y'],linestyle='',marker='o',c='c'))
        warm_core=detect[detect[:,'warm_core']==1]
        tmp.append(self.plot_on_map(self._m,warm_core[:,'x'],warm_core[:,'y'],linestyle='',marker='*',c='r'))
        strong_wind=detect[detect[:,'Wind10']>=self._thr_wind]
        tmp.append(self.plot_on_map(self._m,strong_wind[:,'x'],strong_wind[:,'y'],linestyle='',marker='^',c='m'))

        plt.tight_layout()
        plt.savefig(out_name)
        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def detect(self,overwrite=False):
        out_file=self._working_dir+'detected_positions.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected

        detect=np.array([[np.nan]*7])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()
            coords=peak_local_max(-self._MSLP_smoothed[t,:,:], min_distance=int(self._neighborhood_size))
            for y_p,x_p in zip(coords[:,0],coords[:,1]):
                tc_area,ncont=self.find_closed_contours(self._MSLP_smoothed[t,:,:],y_p,x_p,step=1,search_radius=self._p_radius,method='min')
                if ncont>-999:
                    tmp=[t,y_p,x_p,1,0,0,0]
                    # have to check boundary issues here
                    box=self.get_box(y_p,x_p,self._cores_distance)
                    y_,x_=np.where(self._T[t,1,box[0]:box[1],box[2]:box[3]]==self._T[t,1,box[0]:box[1],box[2]:box[3]].max())
                    y_t,x_t=box[0]+y_[0],box[2]+x_[0]
                    warm_core_area,ncont=self.find_closed_contours(self._T[t,1,:,:],y_t,x_t,step=1,search_radius=self._p_radius,n_contours=1,method='max')
                    yy,xx=np.where(warm_core_area==1)
                    if len(np.where(warm_core_area==1)[0])<self._warm_core_size and ncont==1:
                        tmp[4]=1

                    tmp[5]=self._MSLP[t,y_p,x_p]
                    tmp[6]=np.nanmax(tc_area*self._Wind10[t,:,:])
                    detect=np.concatenate((detect,np.array([tmp])))

        self._detected=da.DimArray(np.array(detect[1:,:],dtype=np.int),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('done')
        return self._detected



def date_(t):
    yr=int(t/10000.)
    mn=int((t-yr*10000)/100.)
    day=int((t-yr*10000-mn*100))
    hr=int((t-yr*10000-mn*100-day)*24)
    return(datetime(yr,mn,day,hr,0))

start = time.time()
identifier='2016'
print('*** started run '+identifier+' ***')
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
data_path='data/ERA5/'
# wind_nc=da.read_nc(data_path+'atl_'+identifier+'_10mWind.nc')
# Wind10=wind_nc['ws'].values
# time_=wind_nc.time
# MSLP=da.read_nc(data_path+'atl_'+identifier+'_surface.nc')['MSL'].values/100.
# T=da.read_nc(data_path+'atl_'+identifier+'_T.nc')['T'].values
# SST=da.read_nc(data_path+'atl_'+identifier+'_surface.nc')['SSTK'].values-273.15
# dates=[date_(t) for t in time_]
#
# U=da.read_nc(data_path+'atl_'+identifier+'_UV.nc')['U']
# V=da.read_nc(data_path+'atl_'+identifier+'_UV.nc')['V']
# VO=rel_vort(U.values[:,0,:,:],V.values[:,0,:,:],U.lat,U.lon)
#
# lons,lats=np.meshgrid(wind_nc.lon,wind_nc.lat)
# lons[lons>180]-=360
#
#
# TC_season=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['season']
# TC_basin=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')['basin']
# tc_sel=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc').ix[np.where((TC_season==int(identifier)) & (TC_basin[:,0]==0))[0]]
#

plt.close('all')
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
m = Basemap(ax=ax,llcrnrlon=np.min(lons),urcrnrlon=np.max(lons),llcrnrlat=np.min(lats),urcrnrlat=np.max(lats),resolution="l",projection='cyl')
m.drawmapboundary(fill_color='1.')
m.drawmapboundary(fill_color='darkblue')
m.fillcontinents(color='darkgreen',lake_color='darkblue')
m.drawcoastlines(linewidth=0.3)
m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)

plot_lon,plot_lat=m(lons,lats)

working_dir='detection/ERA5/'+str(identifier)+'_ERA5_/'
elapsed = time.time() - start;  print('Data loaded %.3f seconds.' % elapsed)
found_tcs=tc_tracks(Wind10=Wind10,MSLP=MSLP,MSLP_smoothed=ndimage.gaussian_filter(MSLP,sigma=(0,3,3)),SST=SST,VO=None,T=T,lats=lats,lons=lons,time_=time_,dates=dates,identifier=identifier,working_dir=working_dir)
found_tcs.init_map(m=m,ax=ax,plot_lat=plot_lat,plot_lon=plot_lon)
found_tcs.init_obs_tcs(tc_sel)
elapsed = time.time() - start;  print('Done with preparations %.3f seconds.' % elapsed)
found_tcs.set_thresholds(thr_wind=17.5,thr_mslp=101500,p_radius=27,neighborhood_size=2,warm_core_size=3,cores_distance=1,search_radius=9,min_time_steps=10)
found_tcs.detect(overwrite=False)
found_tcs.plot_detect_summary()
found_tcs.combine_tracks(overwrite=True)
found_tcs.plot_season()

# plt.close('all')
# fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,5))
# axes=axes.flatten()
# maps=[]
# for ax in axes:
#     mm = Basemap(ax=ax,llcrnrlon=np.min(lons),urcrnrlon=np.max(lons),llcrnrlat=np.min(lats),urcrnrlat=np.max(lats),resolution="l",projection='cyl')
#     mm.drawcoastlines(linewidth=0.7,color='m')
#     mm.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
#     mm.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
#     maps.append(mm)
#
# found_tcs.plot_surrounding(maps=maps,axes=axes,time_steps=range(700,715))#; convert -delay 50 track_surrounding/{94..127}* TC.gif
