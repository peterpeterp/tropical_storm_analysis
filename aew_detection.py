from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max
import cartopy.crs as ccrs
import cartopy
from itertools import combinations
from itertools import permutations
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
mpl.use('Agg')

sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import * ; reload(sys.modules['TC_support'])

def coarsener(data,coarseness = 2):
    if len(data.shape)==2:
        data=data.reshape((1,data.shape[0],data.shape[1]))
    while data.shape[-2] % coarseness!=0:
        data=data[:,:-1,:]
    while data.shape[-1] % coarseness != 0:
        data=data[:,:,:-1]
    tmp = data.reshape((data.shape[0],data.shape[-2]//coarseness,coarseness,data.shape[-1]//coarseness,coarseness))
    return(np.nanmean(tmp, axis=(-3,-1)).squeeze())

def smoother(data,kernel_size=3):
    kernel = np.ones([kernel_size,kernel_size])/float(kernel_size*kernel_size)
    tmp=ndimage.convolve(data,kernel)
    tmp=ndimage.convolve(tmp,kernel)
    return(tmp)

class aew_tracks(object):
    def __init__(self,identifier,working_dir):
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


    def prepare_data(self,lats,lons,time_,dates,smoothing_factor=1,coarsening_factor=1,land_mask=None,time_steps=None):
        self._time=time_
        if time_steps is None:
            time_steps=range(len(self._time))
        self._time_i=time_steps
        self._dates=dates
        self._yr_frac=np.array([toYearFraction(dd) for dd in self._dates])

        self._coarsening_factor=coarsening_factor
        self._smoothing_factor=smoothing_factor

        # input fields
        self._lats_fine=lats
        self._lons_fine=lons
        self._lats=smoother(coarsener(lats,coarsening_factor),smoothing_factor)
        self._lons=smoother(coarsener(lons,coarsening_factor),smoothing_factor)
        self._lat=self._lats[:,0]
        self._lon=self._lons[0,:]

        if land_mask is not None:
            self._land_mask=land_mask
        else:
            self._land_mask=np.ones([len(self._time),lats.shape[0],lats.shape[1]])

        info=da.Dataset({'lats':self._lats,'lons':self._lons})
        info.write_nc(self._working_dir+str(self._identifier)+'_grid.nc',mode='w')


    def add_fields(self,U=None,V=None,VO=None,RH=None,CURV_VORT=None,CURV_VORT_ADVECT=None):
        if VO is not None:
            self._VO_fine=VO
            self._VO=coarsener(VO,self._coarsening_factor)
            for t in range(len(self._time)):
                self._VO[t,:,:]=smoother(self._VO[t,:,:],self._smoothing_factor)

        if RH is not None:
            self._RH_fine=RH
            self._RH=coarsener(RH,self._coarsening_factor)
            for t in range(len(self._time)):
                self._RH[t,:,:]=smoother(self._RH[t,:,:],self._smoothing_factor)

        if U is not None:
            self._U_fine=U
            self._U=coarsener(U,self._coarsening_factor)
            for t in range(len(self._time)):
                self._U[t,:,:]=smoother(self._U[t,:,:],self._smoothing_factor)

        if V is not None:
            self._V_fine=V
            self._V=coarsener(V,self._coarsening_factor)
            for t in range(len(self._time)):
                self._V[t,:,:]=smoother(self._V[t,:,:],self._smoothing_factor)

        if CURV_VORT is not None:
            self._CURV_VORT_fine=CURV_VORT
            self._CURV_VORT=coarsener(CURV_VORT,self._coarsening_factor)
            for t in range(len(self._time)):
                self._CURV_VORT[t,:,:]=smoother(self._CURV_VORT[t,:,:],self._smoothing_factor)

        if CURV_VORT_ADVECT is not None:
            self._CURV_VORT_ADVECT_fine=CURV_VORT_ADVECT
            self._CURV_VORT_ADVECT=coarsener(CURV_VORT_ADVECT,self._coarsening_factor)
            for t in range(len(self._time)):
                self._CURV_VORT_ADVECT[t,:,:]=smoother(self._CURV_VORT_ADVECT[t,:,:],self._smoothing_factor)




    def init_map(self,ax,transform):
        self._ax=ax
        self._transform=transform

        # for storm in range(len(self._tc_sel.storm)):
        #     self._m.plot(self._tc_lon[storm,:],self._tc_lat[storm,:],color='gray')

    def degree_to_step(self,degree):
        y_step=abs(np.diff(self._lats[:,0],1).mean())
        x_step=abs(np.diff(self._lons[0,:],1).mean())
        return round(degree/(y_step+x_step)*2)

    def step_to_distance(self,y1,y2,x1,x2):
        y_step=self._lats[int(y1),int(x1)]-self._lats[int(y2),int(x2)]*6371000*2*np.pi/360.
        x_step=self._lons[int(y1),int(x1)]-self._lons[int(y2),int(x2)]*np.cos(np.radians(self._lats[int(y1),int(x1)]))*6371000*2*np.pi/360.
        return(y_step,x_step)

    def point_to_latlon(self,point):
        return([self._lat[point[0]],self._lon[point[1]]])

    # plotting
    def plot_on_map(self,ax,x_in,y_in,latlon=False,**kwargs):
        if latlon:
            x,y=x_in,y_in
        if latlon==False:
            if isinstance(x_in,np.ndarray) or isinstance(x_in,list) or isinstance(x_in,da.core.dimarraycls.DimArray):
                x=self._lons[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
                y=self._lats[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
            elif isinstance(x_in,np.float64) or isinstance(x_in,int) or isinstance(x_in,float):
                x=self._lons[int(y_in),int(x_in)]
                y=self._lats[int(y_in),int(x_in)]
        return ax.plot(x,y,transform=self._transform,**kwargs)

    def plot_track_path(self,track):
        t=int(track.ix[0,0])
        tmp,txt=[],[]

        tmp.append(self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],c='orange'))
        self._ax.set_title(str(self._dates[t]))

        #tmp+=self.plot_on_map(self._m,self._detected[:,'x'],self._detected[:,'y'],marker='.',linestyle='',color='m')

        plt.savefig(self._working_dir+'track_path/'+str(self._identifier)+'_'+str(t)+'_'+str(self._id)+'.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in txt:
            element.remove()

    def plot_season(self,out_name=None,start_point=True,facecolor='none'):
        tmp,txt=[],[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_tracks.png'

        self._ax.set_title('season '+self._identifier)#

        for id_,track in self._aews.items():
            track=track[np.isfinite(track[:,'t']),:]
            if start_point:
                tmp.append(self.plot_on_map(self._ax,track.ix[0,2],track.ix[0,1],linestyle='',marker='o',c='r'))
                txt.append(self._ax.text(self._lons[int(track.ix[0,1]),int(track.ix[0,2])]-1,self._lats[int(track.ix[0,1]),int(track.ix[0,2])],str(track.ix[0,0]),color='red',va='bottom',ha='right',fontsize=7,transform=self._transform))
            tmp.append(self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],linestyle='-',linewidth=0.5,c='r'))


        plt.tight_layout()
        plt.savefig(out_name,facecolor=facecolor)

        # clean map
        for element in txt:
            element.remove()
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def plot_detect_summary(self,out_name=None):
        tmp=[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_positions.png'

        self._ax.set_title('season '+self._identifier)#

        detect=self._detected.copy()
        tmp.append(self.plot_on_map(self._ax,detect[:,'x'],detect[:,'y'],linestyle='',marker='o',c='g'))

        plt.tight_layout()
        plt.savefig(out_name)
        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def add_surounding_dieng(self,axes,t,thr_VO=1*10**(-5)):
        ax=axes[0]
        ax.contourf(self._lons,self._lats,self._VO[t,:,:],np.arange(1,8,0.5)*10**(-5),cmap=plt.cm.YlOrRd,transform=self._transform)
        ax.contourf(self._lons,self._lats,self._RH[t,:,:], levels=[20, 35, 50, 100],colors='none', hatches=['//////','----','....'   ,'++++++',None],transform=self._transform)

        if hasattr(self,'_detected'):
            for point in self._detected[self._detected[:,'t']==t].values.tolist():
                self.plot_on_map(ax,int(point[2]),int(point[1]),c='c',marker='*',linestyle='')
                #yy,xx=self.find_group(field=self._VO[t,:,:],y=int(point[1]),x=int(point[2]),thresh=thr_VO)
                #self.plot_on_map(ax,xx,yy,c='g',marker='.',linestyle='')

        if hasattr(self,'_aews'):
            for wave in self._aews.values():
                wave=wave[np.isfinite(wave[:,'t']),:]
                if t in wave[:,'t'].values:
                    self.plot_on_map(ax,wave[:,'x'].values.astype(int),wave[:,'y'].values.astype(int),c='m')

    def add_surounding_belanger(self,axes,t,thr_curv_vort=1*10**(-5),thr_u=2.5):
        ax=axes[0]
        ax.contourf(self._lons,self._lats,self._CURV_VORT[t,:,:],np.arange(-3,3,0.5)*10**(-5),cmap=plt.cm.bwr,transform=self._transform)
        ax.contour(self._lons_fine,self._lats_fine,self._VO_fine[t,:,:],np.arange(0,10,2)*10**(-5),colors='k',transform=self._transform)
        ax.contour(self._lons,self._lats,self._CURV_VORT_ADVECT[t,:,:],[-1,0,1],colors='g',transform=self._transform)
        asign = np.sign(self._CURV_VORT_ADVECT[t,:,:])
        signchange_y = ((np.roll(asign, 1,axis=0) - asign) != 0).astype(int)
        signchange_x = ((np.roll(asign, 1,axis=1) - asign) != 0).astype(int)
        y,x=np.where((signchange_y+signchange_x>0))
        #self.plot_on_map(ax,x,y,linestyle='',color='g',marker='.')

        y,x=np.where((signchange_y+signchange_x>0) & (self._U[t,:,:]<thr_u) & (self._CURV_VORT[t,:,:]>thr_curv_vort))
        self.plot_on_map(ax,x,y,linestyle='',color='m',marker='*')

        ax=axes[1]
        if hasattr(self,'_aews'):
            for wave in self._aews.values():
                wave=wave[np.isfinite(wave[:,'t']),:]
                if t in wave[:,'t'].values:
                    self.plot_on_map(ax,wave[:,'x'].values.astype(int),wave[:,'y'].values.astype(int),c='m')

        for id_,track in self._aews.items():
            track=track[np.isfinite(track[:,'t']),:]
            #self.plot_on_map(ax,track[:,'x'],track[:,'y'],linestyle='-',linewidth=1,c='g')

        ax.contourf(self._lons,self._lats,self._CURV_VORT[t,:,:],np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd,transform=self._transform)

        for point in self._detected[self._detected[:,'t']==t].values.tolist():
            self.plot_on_map(ax,[int(point[i]) for i in [4,2,6]],[int(point[i]) for i in [3,1,5]],c='b',marker='.')
            try:
                self.plot_on_map(ax,[int(point[i]) for i in [15,16,17,18,19]],[int(point[i]) for i in [10,11,12,13,14]],c='g',marker='.')
            except:
                pass

    def plot_surrounding(self,axes,style,time_steps=None,thr_VO=1*10**(-5),thr_curv_vort=1*10**(-5),thr_u=2.5):
        if time_steps is None:
            time_steps=self._time_i

        #plt.tight_layout()
        for t in time_steps:
            if style=='dieng':
                self.add_surounding_dieng(axes,t,thr_VO=thr_VO)

            if style=='belanger':
                self.add_surounding_belanger(axes,t,thr_curv_vort=1*10**(-5),thr_u=2.5)

            plt.suptitle(str(self._dates[t]))
            plt.savefig(self._working_dir+'track_surrounding/'+str(self._identifier)+'_'+str(t)+'.png', bbox_inches = 'tight')

            # clean map
            for ax in axes:
                ax.clear()
                ax.set_global()
                ax.coastlines()
                ax.set_xlim(np.min(self._lons),np.max(self._lons))
                ax.set_ylim(np.min(self._lats),np.max(self._lats))

    # combine detected positions
    def combine_tracks(self,plot=True,search_radius=8,overwrite=False):
        out_file=self._working_dir+str(self._identifier)+'_track_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
            os.system('rm '+self._working_dir+'track_path/'+str(self._identifier)+'_*_*.png')

        elif overwrite==False and os.path.isfile(out_file):
            self._aews=da.read_nc(out_file)
            return self._aews

        def unit_vector(vector):
            return vector / np.linalg.norm(vector)

        def angle_between(v1, v2):
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        def v_len(v):
            return sum([zz**2 for zz in v])**0.5


        # convert distances from degrees into grid-cells
        search_radius=self.degree_to_step(search_radius)

        found_id=0
        found_tracks={}

        postions=self._detected.copy().values
        used_pos=[]
        print('combining tracks\n10------50-------100')
        for p,progress in zip(postions.tolist(),np.array([['-']+['']*(len(postions.tolist())/20+1)]*20).flatten()[0:len(postions.tolist())]):
            sys.stdout.write(progress); sys.stdout.flush()
            if p[0:3] not in used_pos:
                track=[p]

                running=True
                #go backwards
                while True:
                    p=track[0]
                    candidates={}
                    for p_1 in postions[postions[:,0]==p[0]-1,:].tolist():
                        dist=((p_1[1]-p[1])**2+(p_1[2]-p[2])**2)**0.5
                        if dist<search_radius:
                            candidates[dist]=p_1
                            end=False
                    if len(candidates.keys())>0:
                        track=[candidates[np.min(candidates.keys())]]+track
                    else:
                        break

                #go forewards
                while True:
                    p=track[-1]
                    candidates={}
                    for p_1 in postions[postions[:,0]==p[0]+1,:].tolist():
                        dist=((p_1[1]-p[1])**2+(p_1[2]-p[2])**2)**0.5
                        if dist<search_radius:
                            candidates[dist]=p_1
                            end=False
                    if len(candidates.keys())>0:
                        track=track+[candidates[np.min(candidates.keys())]]
                    else:
                        break

                track=da.DimArray(track,axes=[np.array(track)[:,0],self._detected.z],dims=['time','z'])
                start_of_track=track[self._lons[np.array(track[:,'y'],int),np.array(track[:,'x'],int)]>-40]
                if len(start_of_track.time)>=2 and len(track.time)>=8:
                    # propagation speed in starting domain
                    c=np.array([self.step_to_distance(start_of_track[i-1,'y'],start_of_track[i,'y'],start_of_track[i-1,'x'],start_of_track[i,'x'])[0]/(6.*60.*60) for i in start_of_track.time[1:] if i-1 in start_of_track.time])
                    if np.percentile(c,50)<-2:
                        if np.mean(self._lats[np.array(start_of_track[:,'y'],int),np.array(start_of_track[:,'x'],int)])<35:
                            if abs(self._lons[int(track.ix[-1,1]),int(track.ix[-1,2])]-self._lons[int(track.ix[0,1]),int(track.ix[0,2])])>6:
                                keep=True
                                for id__,track__ in found_tracks.items():
                                    if sum([pp in track__.values.tolist() for pp in track.values.tolist()])/float(len(track.values.tolist()))!=0:
                                        if track[:,'members'].sum()>track__[:,'members'].sum():
                                            found_tracks.pop(id__)
                                            break
                                        else:
                                            keep=False
                                            break
                                if keep:
                                    found_tracks[found_id]=track
                                    found_id+=1


        self._aews={}
        self._id=0
        for track in found_tracks.values():
            track=da.DimArray(track,axes=[np.array(track)[:,0],self._detected.z],dims=['time','z'])
            self._aews[self._identifier+'_'+str(self._id)]=track
            if plot:    self.plot_track_path(track)
            self._id+=1

        self._aews=da.Dataset(self._aews)
        self._aews.write_nc(out_file,mode='w')
        print('\ndone')
        return self._aews

    def find_group(self,field,y,x,thresh):
        group=[[y,x]]
        for p in group:
            yy_,xx_=p[0],p[1]
            candidates=np.array([[yy_-1,xx_],[yy_,xx_-1],[yy_+1,xx_],[yy_,xx_+1]])
            for yyy_,xxx_ in zip(candidates[:,0],candidates[:,1]):
                if 0<=yyy_<field.shape[0] and 0<=xxx_<field.shape[1] and [yyy_,xxx_] not in group:
                    if field[yyy_,xxx_]>thresh:
                        group.append([yyy_,xxx_])
        return(np.array(group)[:,0],np.array(group)[:,1])

    # detect positions
    def detect_dieng(self,overwrite=False,dis_VO_max=8,contour_radius=25,min_number_cells=6,thr_VO=1*10**(-5),thr_RH=50):
        out_file=self._working_dir+str(self._identifier)+'_detected_positions.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected

        # convert distances from degrees into grid-cells
        dis_VO_max=self.degree_to_step(dis_VO_max)
        contour_radius=self.degree_to_step(contour_radius)

        detect=np.array([[np.nan]*6])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()
            coords=peak_local_max(self._VO[t,:,:], min_distance=int(dis_VO_max))
            #print(coords)
            for y_,x_ in zip(coords[:,0],coords[:,1]):
                if self._VO[t,y_,x_]>thr_VO:
                    yy,xx=self.find_group(field=self._VO[t,:,:],y=y_,x=x_,thresh=thr_VO)
                    if len(yy)>=min_number_cells:
                        if self._RH[t,y_,x_]>=thr_RH:
                            #y_,x_ = sum(yy) / len(yy), sum(xx) / len(yy)
                            tmp=[t,y_,x_,self._VO[t,y_,x_],self._RH[t,y_,x_],len(yy)]
                            detect=np.concatenate((detect,np.array([tmp])))

        self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','VO','RH','members']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('\ndone')
        return self._detected



    # detect positions
    def detect_belanger(self,overwrite=False,max_extend=25,thr_u=2.5,thr_curv_vort=1*10**(-5)):
        out_file=self._working_dir+str(self._identifier)+'_detected_positions.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected


        def group_points(points):
            used_pos=[]
            groups=[]
            i=0
            while len(used_pos)<len(points):
                p=points[i]
                i+=1
                if p not in used_pos:
                    used_pos.append(p)
                    group=[p]
                    for p in group:
                        yy,xx=p[0],p[1]
                        candidates=[[yy+ystep,xx+xstep] for ystep,xstep in permutations(np.arange(-1,2,1),2)]
                        for pp in candidates:
                            if pp in points and pp not in used_pos:
                                group.append(pp)
                                used_pos.append(pp)

                    groups.append(group)

            return(groups)

        def group_extend(group):
            if len(group)<2:
                return(0,[[np.nan,np.nan],[np.nan,np.nan]])
            max_square_distance = 0
            for pair in combinations(group,2):
                dist=np.sqrt((self._lat[pair[0][0]]-self._lat[pair[1][0]])**2+(self._lon[pair[0][1]]-self._lon[pair[1][1]])**2)
                if dist > max_square_distance:
                    max_square_distance = dist
                    max_pair = pair
            return(max_square_distance,max_pair)

        # convert distances from degrees into grid-cells
        max_extend=self.degree_to_step(max_extend)

        detect=np.array([[np.nan]*20])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()

            asign = np.sign(self._CURV_VORT_ADVECT[t,:,:])
            signchange_y = ((np.roll(asign, 1,axis=0) - asign) != 0).astype(int)
            signchange_x = ((np.roll(asign, 1,axis=1) - asign) != 0).astype(int)
            y,x=np.where((signchange_y+signchange_x>0) & (self._U[t,:,:]<thr_u) & (self._CURV_VORT[t,:,:]>thr_curv_vort))
            x=x[(5<self._lat[y]) & (self._lat[y]<35)]
            y=y[(5<self._lat[y]) & (self._lat[y]<35)]

            points=np.vstack((y,x)).T.tolist()
            groups=group_points(points)

            while True:
                done=True
                for group in groups:
                    done=True
                    if group_extend(group)[0]>max_extend:
                        done=False
                        groups.remove(group)
                        thresh=5*10**(-6)
                        while True:
                            thresh*=1.5
                            for pp in group:
                                if self._CURV_VORT[t,pp[0],pp[1]]<thresh:
                                    group.remove(pp)
                            sub_groups=group_points(group)
                            stop=True
                            for sub in sub_groups:
                                if group_extend(sub)[0]>max_extend:
                                    group=sub
                                    stop=False
                                else:
                                    groups.append(sub)
                            if stop:
                                break
                if done:
                    break

            for group in groups:
                if len(group)>3:
                    dist,pair=group_extend(group)
                    x = [p[0] for p in group]
                    y = [p[1] for p in group]
                    centr = [sum(x) / len(group), sum(y) / len(group)]
                    max_vort=np.max([self._CURV_VORT[t,int(pp[0]),int(pp[1])] for pp in group])
                    #tmp=[t,np.median(np.array(group)[:,0]),np.median(np.array(group)[:,1])]+pair[0]+pair[1]
                    tmp=[t]+centr+pair[0]+pair[1]+[len(group),max_vort,self._U[t,centr[0],centr[1]]]
                    # save a simplified polygon
                    try:
                        hull = ConvexHull(group)
                        y,x=Polygon(np.array(group)[hull.vertices,:]).simplify(1).exterior.xy
                        detect=np.concatenate((detect,np.array([tmp+list(y)+list(x)])))
                    except:
                        detect=np.concatenate((detect,np.array([tmp+[np.nan]*10])))

        self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','y_ext1','x_ext1','y_ext2','x_ext2','members','max_CURV_VORT','u_centroid','y1','y2','y3','y4','y5','x1','x2','x3','x4','x5']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('\ndone')
        return self._detected


    #
    # def group_points(self,points):
    #     used_pos=[]
    #     groups=[]
    #     i=0
    #     while len(used_pos)<len(points):
    #         p=points[i]
    #         i+=1
    #         if p not in used_pos:
    #             used_pos.append(p)
    #             group=[p]
    #             for p in group:
    #                 yy,xx=p[0],p[1]
    #                 candidates=[[yy+ystep,xx+xstep] for ystep,xstep in permutations(np.arange(-2,3,1),2)]
    #                 for pp in candidates:
    #                     if pp in points and pp not in used_pos:
    #                         group.append(pp)
    #                         used_pos.append(pp)
    #
    #             groups.append(group)
    #     return(groups)
    #
    # def group_extend(self,group):
    #     if len(group)<2:
    #         return(0,[[np.nan,np.nan],[np.nan,np.nan]])
    #     max_square_distance = 0
    #     for pair in combinations(group,2):
    #         dist=np.sqrt((self._lat[pair[0][0]]-self._lat[pair[1][0]])**2+(self._lon[pair[0][1]]-self._lon[pair[1][1]])**2)
    #         if dist > max_square_distance:
    #             max_square_distance = dist
    #             max_pair = pair
    #     return(max_square_distance,max_pair)
    #
    # def reduce_groups(self,groups,max_extend,field,thresh,thresh_factor):
    #     while True:
    #         for group in groups:
    #             done=True
    #             if self.group_extend(group)[0]>max_extend:
    #                 done=False
    #                 groups.remove(group)
    #                 thresh__=thresh
    #                 while True:
    #                     thresh__*=thresh_factor
    #                     for pp in group:
    #                         if field[pp[0],pp[1]]<thresh__:
    #                             group.remove(pp)
    #                     sub_groups=self.group_points(group)
    #                     stop=True
    #                     for sub in sub_groups:
    #                         if self.group_extend(sub)[0]>max_extend:
    #                             group=sub
    #                             stop=False
    #                         else:
    #                             groups.append(sub)
    #                     if stop:
    #                         break
    #         if done:
    #             break
    #     return(groups)
    #
    # # detect positions
    # def detect_dieng_groups(self,overwrite=False,dis_VO_max=8,max_extend=25,min_number_cells=6,thr_VO=1*10**(-5),thr_RH=50):
    #     self._add_name='dieng_groups'
    #     out_file=self._working_dir+'detected_positions_'+self._add_name+'.nc'
    #     if overwrite and os.path.isfile(out_file):
    #         os.system('rm '+out_file)
    #     elif overwrite==False and os.path.isfile(out_file):
    #         self._detected=da.read_nc(out_file)['detected']
    #         return self._detected
    #
    #
    #     # convert distances from degrees into grid-cells
    #     dis_VO_max=self.degree_to_step(dis_VO_max)
    #     max_extend=self.degree_to_step(max_extend)
    #
    #     detect=np.array([[np.nan]*6])
    #     print('detecting\n10------50-------100')
    #     for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
    #         sys.stdout.write(progress); sys.stdout.flush()
    #
    #         y,x=np.where( (self._VO[t,:,:]>thr_VO) & (self._RH[t,:,:]>thr_RH))
    #         points=np.vstack((y,x)).T.tolist()
    #         start_time=time.time()
    #         groups=self.group_points(points)
    #         # print(time.time()-start_time,len(groups))
    #         if groups is not None:
    #             groups=self.reduce_groups(groups,max_extend=max_extend,field=self._VO[t,:,:],thresh=thr_VO,thresh_factor=1.1)
    #             # print(time.time()-start_time)
    #             for group in groups:
    #                 yy = [p[0] for p in group]
    #                 xx = [p[1] for p in group]
    #                 y_,x_ = sum(yy) / len(group), sum(xx) / len(group)
    #
    #                 if len(group)>=min_number_cells:
    #                     tmp=[t,y_,x_,self._VO[t,y_,x_],self._RH[t,y_,x_],len(group)]
    #                     detect=np.concatenate((detect,np.array([tmp])))
    #
    #     self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','VO','RH','cells_above_thr']],dims=['ID','z'])
    #     da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
    #     print('done')
    #     return self._detected
