from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max
import cartopy.crs as ccrs
import cartopy
from itertools import combinations

sns.set_palette(sns.color_palette("plasma"))



class aew_tracks(object):
    def __init__(self,U,V,lats,lons,time_,dates,identifier,working_dir,land_mask=None,time_steps=None):
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
        self._lat=lats[:,0]
        self._lon=lons[0,:]
        self._time=time_
        if time_steps is None:
            time_steps=range(len(self._time))
        self._time_i=time_steps
        self._dates=dates

        #self._yr_frac=np.array([toYearFraction(dd) for dd in self._dates])

        U[np.isnan(U)]=0
        V[np.isnan(V)]=0
        self._u=U.values
        self._v=V.values

        xx,yy = np.meshgrid(U.lon,U.lat)
        dx,dy = np.meshgrid(U.lon.copy()*0+np.mean(np.diff(U.lon,1)),U.lat.copy()*0+np.mean(np.diff(U.lat,1)))
        dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
        dy*=6371000*2*np.pi/360.

        self._vo,self._curv_vort,self._curv_vort_advect=self._u.copy(),self._u.copy(),self._u.copy()

        for t in self._time_i:
            u,v=self._u[t,:,:],self._v[t,:,:]
            W = (u**2+v**2)**0.5
            du_dx = (u-np.roll(u,1,axis=1))/dx
            du_dy = (u-np.roll(u,1,axis=0))/dy
            dv_dx = (v-np.roll(v,1,axis=1))/dx
            dv_dy = (v-np.roll(v,1,axis=0))/dy
            dW_dx = (W-np.roll(W,1,axis=1))/dx
            dW_dy = (W-np.roll(W,1,axis=0))/dy

            vo=dv_dx-du_dy
            self._vo[t,:,:]=vo

            shear_vort=v/W*dW_dx-u/W*dW_dy
            curv_vort=vo-shear_vort

            kernel_size=5
            kernel = np.ones([kernel_size,kernel_size])/float(kernel_size*kernel_size)

            curv_vort_coarse=ndimage.filters.gaussian_filter(curv_vort,3)
            curv_vort_smoo=ndimage.convolve(curv_vort_coarse,kernel)
            curv_vort_smoo=ndimage.convolve(curv_vort_smoo,kernel)
            self._curv_vort[t,:,:]=curv_vort_smoo

            dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=1))/dx
            dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=0))/dy
            curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

            curv_vort_advect_coarse=ndimage.filters.gaussian_filter(curv_vort_advect,3)
            curv_vort_advect_smoo=ndimage.convolve(curv_vort_advect_coarse,kernel)
            curv_vort_advect_smoo=ndimage.convolve(curv_vort_advect_smoo,kernel)
            self._curv_vort_advect[t,:,:]=curv_vort_advect_smoo

        if land_mask is not None:
            self._land_mask=land_mask
        else:
            self._land_mask=self._u[0,:,:].copy()*0+1

    def init_map(self,ax,transform):
        self._ax=ax
        self._transform=transform

        # for storm in range(len(self._tc_sel.storm)):
        #     self._m.plot(self._tc_lon[storm,:],self._tc_lat[storm,:],color='gray')

    def degree_to_step(self,degree):
        y_step=abs(np.diff(self._lats[:,0],1).mean())
        x_step=abs(np.diff(self._lons[0,:],1).mean())
        return round(degree/(y_step+x_step)*2)

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

        plt.savefig(self._working_dir+'track_path/'+str(self._identifier)+'_'+self._add_name+'_'+str(t)+'_'+str(self._id)+'.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in txt:
            element.remove()

    def plot_season(self,out_name=None,start_point=True,facecolor='none'):
        tmp,txt=[],[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_tracks_'+self._add_name+'.png'

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
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_positions_'+self._add_name+'.png'

        self._ax.set_title('season '+self._identifier)#

        detect=self._detected.copy()
        tmp.append(self.plot_on_map(self._ax,detect[:,'x'],detect[:,'y'],linestyle='',marker='o',c='g'))

        plt.tight_layout()
        plt.savefig(out_name)
        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def plot_surrounding(self,axes,time_steps=None,thr_u=2.5,thr_curv_vort=1*10**(-5)):
        if time_steps is None:
            time_steps=self._time_i

        #plt.tight_layout()
        for t in time_steps:
            ax=axes[0]

            ax.contourf(self._lons,self._lats,self._curv_vort[t,:,:],np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd,transform=self._transform)
            ax.contour(self._lons,self._lats,self._curv_vort_advect[t,:,:],[-1,0,1],color='k',transform=self._transform)
            asign = np.sign(self._curv_vort_advect[t,:,:])
            signchange_y = ((np.roll(asign, 1,axis=0) - asign) != 0).astype(int)
            signchange_x = ((np.roll(asign, 1,axis=1) - asign) != 0).astype(int)
            y,x=np.where((signchange_y+signchange_x>0))
            self.plot_on_map(ax,x,y,linestyle='',color='g',marker='.')

            y,x=np.where((signchange_y+signchange_x>0) & (self._u[t,:,:]<thr_u) & (self._curv_vort[t,:,:]>thr_curv_vort))
            self.plot_on_map(ax,x,y,linestyle='',color='m',marker='*')

            ax=axes[1]
            for id_,track in self._aews.items():
                track=track[np.isfinite(track[:,'t']),:]
                self.plot_on_map(ax,track[:,'x'],track[:,'y'],linestyle='-',linewidth=1,c='g')

            ax.contourf(self._lons,self._lats,self._curv_vort[t,:,:],np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd,transform=self._transform)

            for point in self._detected[self._detected[:,'t']==t].values.tolist():
                self.plot_on_map(ax,[int(point[i]) for i in [4,2,6]],[int(point[i]) for i in [3,1,5]],c='b',marker='.')

            plt.suptitle(str(self._dates[t]))
            plt.savefig(self._working_dir+'track_surrounding/'+self._add_name+'_'+str(t)+'.png', bbox_inches = 'tight')

            # clean map
            for ax in axes:
                ax.clear()
                ax.set_global()
                ax.coastlines()
                ax.set_xlim(np.min(self._lon),np.max(self._lon))
                ax.set_ylim(np.min(self._lat),np.max(self._lat))

    # combine detected positions
    def combine_tracks(self,plot=True,search_radius=6,overwrite=False):
        out_file=self._working_dir+'track_info_'+self._add_name+'.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
            os.system('rm '+self._working_dir+'track_path/*_'+self._add_name+'_*_*.png')

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

        def consecutive_sequence(zz):
            i,su,out=0,0,[]
            while i <len(zz):
                if zz[i]:
                    su+=1
                else:
                    out.append((i-su,su))
                    su=0
                i+=1
            out.append((i-su,su))
            return np.array(out)

        # convert distances from degrees into grid-cells
        search_radius=self.degree_to_step(search_radius)

        self._id=0
        self._aews={}

        postions=self._detected.copy().values
        used_pos=[]
        for p in postions.tolist():
            if p not in used_pos:
                track=[p]

                running=True
                #go backwards
                while True:
                    p=track[0]
                    if len(track)==1:
                        prev_step=np.array((0,0,0))
                    else:
                        prev_step=np.array((p[1]-track[1][1],p[2]-track[1][2],0))
                    candidates={}
                    for p_1 in postions[postions[:,0]==p[0]-1,:].tolist():
                        new_step=np.array((p_1[1]-p[1],p_1[2]-p[2],0))
                        if v_len(new_step-prev_step)<search_radius:
                            candidates[v_len(new_step-prev_step)]=p_1
                            end=False
                    if len(candidates.keys())>0:
                        track=[candidates[np.min(candidates.keys())]]+track
                    else:
                        break

                #go forewards
                while True:
                    p=track[-1]
                    if len(track)==1:
                        prev_step=np.array((0,0,0))
                    else:
                        prev_step=np.array((p[1]-track[-2][1],p[2]-track[-2][2],0))
                    candidates={}
                    for p_1 in postions[postions[:,0]==p[0]+1,:].tolist():
                        new_step=np.array((p_1[1]-p[1],p_1[2]-p[2],0))
                        if v_len(new_step-prev_step)<search_radius:
                            candidates[v_len(new_step-prev_step)]=p_1
                            end=False
                    if len(candidates.keys())>0:
                        track=track+[candidates[np.min(candidates.keys())]]
                    else:
                        break

                # # delete land steps in the beginning of track
                # while len(track)>0:
                #     if self._land_mask[int(track[0][1]),int(track[0][2])]:
                #         track=track[1:]
                #     else:
                #         break

                if len(track)>=8:
                    if sum([pp in used_pos for pp in track])/float(len(track))<0.3:
                        if -25<np.percentile([self._u[int(pp[0]),int(pp[1]),int(pp[2])] for pp in track],66)<-2:
                            used_pos+=track
                            track=da.DimArray(track,axes=[np.array(track)[:,0],['t','y','x','y_ext1','x_ext1','y_ext2','x_ext2','members','max_curv_vort','u_centroid']],dims=['time','z'])
                            self._aews[self._identifier+'_'+str(self._id)]=track
                            if plot:    self.plot_track_path(track)
                            self._id+=1

        self._aews=da.Dataset(self._aews)
        self._aews.write_nc(out_file,mode='w')
        return self._aews

    # detect positions
    def detect(self,overwrite=False,max_extend=25,thr_u=2.5,thr_curv_vort=1*10**(-5)):
        self._add_name='belanger'
        out_file=self._working_dir+'detected_positions_'+self._add_name+'.nc'
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
                        candidates=[[yy-1,xx-1],[yy-1,xx],[yy-1,xx+1],[yy,xx+1],[yy+1,xx+1],[yy+1,xx],[yy+1,xx-1],[yy,xx-1]]
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

        detect=np.array([[np.nan]*10])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()

            asign = np.sign(self._curv_vort_advect[t,:,:])
            signchange_y = ((np.roll(asign, 1,axis=0) - asign) != 0).astype(int)
            signchange_x = ((np.roll(asign, 1,axis=1) - asign) != 0).astype(int)
            y,x=np.where((signchange_y+signchange_x>0) & (self._u[t,:,:]<thr_u) & (self._curv_vort[t,:,:]>thr_curv_vort))
            x=x[(5<self._lat[y]) & (self._lat[y]<35)]
            y=y[(5<self._lat[y]) & (self._lat[y]<35)]

            points=np.vstack((y,x)).T.tolist()
            groups=group_points(points)

            while True:
                for group in groups:
                    done=True
                    if group_extend(group)[0]>max_extend:
                        done=False
                        groups.remove(group)
                        thresh=5*10**(-6)
                        while True:
                            thresh*=1.5
                            for pp in group:
                                if self._curv_vort[t,pp[0],pp[1]]<thresh:
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
                if len(group)>1:
                    dist,pair=group_extend(group)
                    x = [p[0] for p in group]
                    y = [p[1] for p in group]
                    centroid = [sum(x) / len(group), sum(y) / len(group)]
                    max_vort=np.max([self._curv_vort[t,int(pp[0]),int(pp[1])] for pp in group])
                    #tmp=[t,np.median(np.array(group)[:,0]),np.median(np.array(group)[:,1])]+pair[0]+pair[1]
                    tmp=[t]+centroid+pair[0]+pair[1]+[len(group),max_vort,self._u[t,centroid[0],centroid[1]]]
                    detect=np.concatenate((detect,np.array([tmp])))

        self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','y_ext1','x_ext1','y_ext2','x_ext2','members','max_curv_vort','u_centroid']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('done')
        return self._detected


    # # detect positions
    # def detect_old(self,overwrite=False,loc_max_dist=5,thr_u=2.5,thr_curv_vort=1*10**(-5)):
    #     self._add_name='contours'
    #     out_file=self._working_dir+'detected_positions_'+self._add_name+'.nc'
    #     if overwrite and os.path.isfile(out_file):
    #         os.system('rm '+out_file)
    #     elif overwrite==False and os.path.isfile(out_file):
    #         self._detected=da.read_nc(out_file)['detected']
    #         return self._detected
    #
    #
    #     # convert distances from degrees into grid-cells
    #     loc_max_dist=self.degree_to_step(loc_max_dist)
    #
    #     detect=np.array([[np.nan]*3])
    #     print('detecting\n10------50-------100')
    #     for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
    #         sys.stdout.write(progress); sys.stdout.flush()
    #         coords=peak_local_max(self._curv_vort[t,:,:], min_distance=int(loc_max_dist))
    #         for y,x in coords:
    #             if 5<self._lats[y,x]<35:
    #                 if self._curv_vort[t,y,x]>thr_curv_vort:
    #                     if self._u[t,y,x]<thr_u:
    #                         if np.sign(self._curv_vort_advect[t,y,x])!=np.sign(self._curv_vort_advect[t,y,x+1]):
    #                             tmp=[t,y,x]
    #                             detect=np.concatenate((detect,np.array([tmp])))
    #                         if np.sign(self._curv_vort_advect[t,y,x-1])!=np.sign(self._curv_vort_advect[t,y,x]):
    #                             tmp=[t,y,x]
    #                             detect=np.concatenate((detect,np.array([tmp])))
    #
    #
    #     self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x']],dims=['ID','z'])
    #     da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
    #     print('done')
    #     return self._detected
