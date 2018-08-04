# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

from __future__ import print_function
import os,sys,glob,time,weakref,datetime

from netCDF4 import Dataset,num2date
import numpy as np
import dimarray as da

import scipy.ndimage as ndimage
from skimage.feature import peak_local_max
import cv2

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
import cartopy
sns.set_palette(sns.color_palette("plasma"))

sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
from TC_support import * ; reload(sys.modules['TC_support'])

def unit_vector(vector):
    '''
    get the unit vector from a vector
    '''
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    '''
    get angle between two vectors
    '''
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

class tc_tracks(object):
    """
    This class can be used to detect tropical cyclones. Depending on the method used to detect storm snapshots different input data is required. Currently two storm snapshot detecion methods are implemented: :meth:`detect_contours` or :meth:`detect_knutson2007`. When adding new methods you should stick to the output format of these methods.
    """
    def __init__(self,VO,Wind10,MSLP,MSLP_smoothed,SST,T,T_diff,land_mask,lats,lons,time_,dates,identifier,working_dir,time_steps=None):
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

        if land_mask is not None:
            self._land_mask=land_mask
        else:
            self._land_mask=self._MSLP[0,:,:].copy()*0+1

        self._time=time_
        if time_steps is None:
            time_steps=range(len(self._time))
        self._time_i=time_steps
        self._dates=dates
        self._yr_frac=np.array([toYearFraction(dd) for dd in self._dates])

        self._T=T
        self._T_diff=T_diff
        self._SST=SST

        # tc cat dict
        self._obs_tc=False
        self._cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
        self._cat_names={0:'tropical storm',1:'Category 1',2:'Category 2',3:'Category 3',4:'Category 4',5:'Category 5'}

    def degree_to_step(self,degree):
        """
        Converts distances in degrees into steps on the lat-lon grid

        Parameters
        ----------
            degree: float [deg]
                distance to be converted

        Returns
        -------
            number of steps on the lon-lat grid corresponding to the distance
        """
        y_step=abs(np.diff(self._lats[:,0],1).mean())
        x_step=abs(np.diff(self._lons[0,:],1).mean())
        return round(degree/(y_step+x_step)*2)

    def tc_cat(self,z,method='pressure'):
        """
        Returns the category of a storm snapshot based on 10m wind speed or MSLP

        Parameters
        ----------
            z: float [hPa or m/s]
                value tested against TC-categors thresholds
            method: str ['pressure' or 'wind'], default='pressure'
                specifies if `z` is a wind-speed or a pressure value

        Returns
        -------
            cat: int
                category of the TC form 0-5
        """
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

    # treating ibtracks
    def init_obs_tcs(self,tc_sel):
        """
        This functions transforms ibtracks input into the format used by plot functions

        Parameters
        ----------
            tc_sel: dimarray
                ibtracks loaded as dimarray
        """
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

    def obs_track_info(self,core_radius=3,full_radius=7,overwrite=False):
        """
        Collects and sves info of input fields around observed tracks

        Parameters
        ----------
            core_radius: int, default=3
                radius of grid-cells around storm center which are going to be analyzed for core specific variables
            full_radius: int, default=7
                radius of grid-cells around storm center which are going to be analyzed
            overwrite:  bool, default=False
                if False, existing outputfiles of the function are loaded, if True they are overwritten

        Returns
        -------
            self._obs_track_info: dimarray
                array containg info at time steps and positions of observed storms
        """
        out_file=self._working_dir+'obs_track_info.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._obs_track_info=da.read_nc(out_file)
            return self._obs_track_info


        # convert distances from degrees into grid-cells
        core_radius=int(self.degree_to_step(core_radius))
        full_radius=int(self.degree_to_step(full_radius))

        obs_summary=np.zeros([len(self._tc_sel.storm),200,10])*np.nan
        for i,storm in enumerate(self._tc_sel.storm):
            tmp_t=self._tc_time[i,:]
            last_val=len(np.where(np.isfinite(tmp_t))[0])
            obs_summary[i,0:last_val,0]=[self.tc_cat(z,method='wind') for z in np.nanmean(self._tc_sel['source_wind'].values[i,0:last_val,:],axis=-1)]
            obs_summary[i,0:last_val,1]=np.nanmean(self._tc_sel['source_wind'].values[i,0:last_val,:],axis=-1)
            obs_summary[i,0:last_val,2]=np.nanmean(self._tc_sel['source_pres'].values[i,0:last_val,:],axis=-1)

            for t in range(last_val):
                t_=np.where(abs(self._yr_frac-self._tc_time[i,t])<0.0004)[0]
                if len(t_)!=0:
                    t_=t_[0]
                    y,x=np.argmin(abs(self._lats[:,0]-self._tc_lat[i,t])),np.argmin(abs(self._lons[0,:]-self._tc_lon[i,t]))
                    y_core,x_core=self.area_around(y,x,core_radius)
                    y_full,x_full=self.area_around(y,x,full_radius)
                    obs_summary[i,t,3]=self._VO[t_,y_core,x_core].max()
                    obs_summary[i,t,4]=self._MSLP[t_,y_core,x_core].min()
                    obs_summary[i,t,5]=self._Wind10[t_,y_full,x_full].max()
                    obs_summary[i,t,6]=self._T[t_,y_core,x_core].max()
                    obs_summary[i,t,7]=self._T_diff[t_,y_core,x_core].max()
                    if self._SST is not None:
                        obs_summary[i,t,8]=self._SST[t_,y,x]
                    obs_summary[i,t,9]=self._land_mask[y,x]

        obs_summary=obs_summary[:,np.isfinite(np.nanmean(obs_summary,axis=(0,-1))),:]
        self._obs_track_info=da.DimArray(obs_summary,axes=[self._tc_sel.storm,range(obs_summary.shape[1]),['cat','obs_wind','obs_pres','VO','MSLP','Wind10','T','T_diff','SST','land']],dims=['storm','time','variable'])

        da.Dataset({'obs_track_info':self._obs_track_info}).write_nc(out_file)

        # print summary
        sys.stdout.write('Category:\t0\t\t1\t\t2\t\t3\t\t4\t\t5'); sys.stdout.flush()
        for vari,name in zip(range(2,7),['VO','MSLP','Wind10','T850','T500','SST']):
            sys.stdout.write('\n'+name+'\t\t'); sys.stdout.flush()
            for cat in range(6):
                pos=np.where(obs_summary==cat)
                sys.stdout.write(str(np.nanmean(obs_summary[pos[0],pos[1],vari]))+'\t'); sys.stdout.flush()

        return self._obs_track_info

    # plotting
    def init_map(self,ax,transform):
        """
        All plot functions except :meth:`plot_surrounding` plot on the map initialized in this function

        Parameters
        ----------
            ax: cartopy.mpl.geoaxes.GeoAxesSubplot
                an axis with cartopy projection. This ax object can be created with a command as `ax = plt.axes(projection=ccrs.PlateCarree())`
            transform: cartopy.crs, default=ccrs.PlateCarree()
                this argument will be added to all plot functions to  interprete given lat-lon coordinates
        """
        self._ax=ax
        self._transform=transform

        # for storm in range(len(self._tc_sel.storm)):
        #     self._m.plot(self._tc_lon[storm,:],self._tc_lat[storm,:],color='gray')

    def plot_on_map(self,ax,x_in,y_in,cat=None,latlon=False,**kwargs):
        """
        function used to plot tracks or positions on a map (usually the map initialized in :meth:`init_map`)

        Parameters
        ----------
            ax: cartopy.mpl.geoaxes.GeoAxesSubplot
                the axes on which to plot. This is usually `self._ax` except when called from :meth:`plot_surrounding`
            x_in: int or array
                x-location(s) to be ploted
            y_in: int or array
                y-location(s) to be ploted
            cat: int or array, default=None
                category of the TC-scale defining the plot-color
            latlon: bool, default=False
                if True, `x_in` and `y_in` are lons and lats. If False, they are indices of the `self._lons` and `self._lats` arrays initialized in :meth:`__init__`
            **kwargs
                key word arguments used in plot functions

        Returns
        -------
             ax.plot():
                the output of the plot command
        """
        if latlon:
            x,y=x_in,y_in
        if latlon==False:
            if isinstance(x_in,np.ndarray) or isinstance(x_in,list) or isinstance(x_in,da.core.dimarraycls.DimArray):
                x=self._lons[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
                y=self._lats[[int(yy) for yy in y_in],[int(xx) for xx in x_in]]
            elif isinstance(x_in,np.float64) or isinstance(x_in,int) or isinstance(x_in,float):
                x=self._lons[int(y_in),int(x_in)]
                y=self._lats[int(y_in),int(x_in)]
        if cat is not None:
            tmp=[]
            for i in range(len(x)-1):
                if np.isfinite(x[i+1]):
                    tmp.append(ax.plot(x[i:i+2],y[i:i+2],color=self._cat_colors[cat[i]],transform=self._transform,**kwargs))
            return tmp
        else:
            return ax.plot(x,y,transform=self._transform,**kwargs)

    def plot_track_path(self,track):
        """
        Plots an individual track and saves it in `self._working_dir+'track_path/'`

        Parameters
        ----------
            track: dimarray
                track saved by :meth:`combine_tracks`
        """

        t=int(track.ix[np.nanargmin(track[:,'MSLP'].values),0])
        tmp,txt=[],[]

        if self._obs_tc:
            storms=np.where(abs(self._tc_time-self._yr_frac[t])<0.002)[0]
            for storm in set(storms):
                tmp+=self.plot_on_map(self._ax,self._tc_lon[storm,:],self._tc_lat[storm,:],cat=self.tc_cat(self._tc_intens[storm,:],method='wind'),latlon=True)
                last_pos=np.where(np.isfinite(self._tc_lon[storm,:]))[0][-1]
                txt.append(self._ax.text(self._tc_lon[storm,last_pos],self._tc_lat[storm,last_pos],''.join(self._tc_sel['name'].ix[storm,:]),transform=self._transform))

        tmp.append(self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],c='orange'))
        tmp+=self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],cat=self.tc_cat(track[:,'MSLP'].values),marker='.',linestyle='')
        self._ax.set_title(str(self._dates[t]))

        #tmp+=self.plot_on_map(self._m,self._detected[:,'x'],self._detected[:,'y'],marker='.',linestyle='',color='m')

        plt.savefig(self._working_dir+'track_path/'+str(self._identifier)+'_'+self._add_name+'_'+str(t)+'_'+str(self._id)+'.png')

        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
        for element in txt:
            element.remove()

    def plot_all_tracks(self,out_name=None,start_point=True,facecolor='none'):
        """
        Plot all tracks found in the dataset

        Parameters
        ----------
            out_name: str, default=None
                filename under which the plot is saved. If None, it is saved in `self.working_dir`
            start_point: bool, default=True
                if True, start point of the track is highlighted
            facecolor: str, default="none"
                background color of the saved plot

        """

        tmp,txt=[],[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_tracks_'+self._add_name+'.png'

        self._ax.set_title('season '+self._identifier)#

        summary={0:[],1:[],2:[],3:[],4:[],5:[]}
        for id_,track in self._tcs.items():
            track=track[np.isfinite(track[:,'t']),:]
            if start_point:
                tmp.append(self.plot_on_map(self._ax,track.ix[0,2],track.ix[0,1],linestyle='',marker='o',c='r'))
                txt.append(self._ax.text(self._lons[int(track.ix[0,1]),int(track.ix[0,2])]-1,self._lats[int(track.ix[0,1]),int(track.ix[0,2])],str(track.ix[0,0]),color='red',va='bottom',ha='right',fontsize=7,transform=self._transform))
            tmp.append(self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],linestyle='-',linewidth=0.5,c='r'))
            tmp+=self.plot_on_map(self._ax,track[:,'x'],track[:,'y'],cat=self.tc_cat(track[:,'MSLP'].values),marker='.',linestyle='')
            summary[max(self.tc_cat(track[:,'MSLP'].values))].append(id_)

        if self._obs_tc:
            for storm in range(len(self._tc_sel.storm)):
                tmp+=self.plot_on_map(self._ax,self._tc_lon[storm,:],self._tc_lat[storm,:],cat=self.tc_cat(self._tc_intens[storm,:],method='wind'),latlon=True)


        summary.pop(0)
        for cat,y in zip(summary.keys(),[0.99,0.95,0.91,0.87,0.83]):
            txt.append(self._ax.text(0.005,y,self._cat_names[cat]+': '+''.join(['X']*len(summary[cat])),transform=self._ax.transAxes,color=self._cat_colors[cat],va='top',ha='left',fontsize=12))
        plt.tight_layout()
        plt.savefig(out_name,facecolor=facecolor)

        # clean map
        for element in txt:
            element.remove()
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def plot_detect_summary(self,thr_wind=17.5,out_name=None):
        """
        Plots all storm snapshots detected in the dataset. Highlights warm core snapshots with red triangles and strong storm snapshots with yellwo triangles

        Parameters
        ----------
            thr_wind: float [m/s], default=17.5
                threshold used to differentiate between strom and weak storms
            out_name: str, default=None
                filename under which the plot is saved. If None, it is saved in `self.working_dir`
        """

        tmp=[]
        if out_name is None:
            out_name=self._working_dir+'season_'+str(self._identifier)+'_found_positions_'+self._add_name+'.png'

        self._ax.set_title('season '+self._identifier)#

        detect=self._detected.copy()
        tmp.append(self.plot_on_map(self._ax,detect[:,'x'],detect[:,'y'],linestyle='',marker='o',c='g'))
        warm_core=detect[detect[:,'warm_core']==1]
        tmp.append(self.plot_on_map(self._ax,warm_core[:,'x'],warm_core[:,'y'],linestyle='',marker='v',c='r'))
        strong_wind=detect[detect[:,'Wind10']>=thr_wind]
        tmp.append(self.plot_on_map(self._ax,strong_wind[:,'x'],strong_wind[:,'y'],linestyle='',marker='^',c='y'))

        plt.tight_layout()
        plt.savefig(out_name)
        # clean map
        for element in tmp:
            l = element.pop(0); wl = weakref.ref(l); l.remove(); del l

    def plot_surrounding(self,axes,time_steps=None):
        """
        Function plotting the four fields MSLP, Ta, WS and Vort for selected time steps. The function also highlights fulfilled warm core and strong storm conditions for detected snapshots. If observed data is loaded, information about observed storms is also added.

        Parameters
        ----------
            axes: array of 4 cartopy.mpl.geoaxes.GeoAxesSubplot
                an array with 4 axes with cartopy projection. This ax object can be created with a command as `fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(11.5,5),subplot_kw={'projection': ccrs.PlateCarree()}); axes=axes.flatten()`
            time_steps: array of ints [indices of time axis], default=None
                time steps for which plots are created. If None, plots are created for all time steps

        """

        if time_steps is None:
            time_steps=self._time_i

        #plt.tight_layout()
        for t in time_steps:
            tmp,txt=[],[]
            ax=axes[0]; ax.set_title('mean sea level pressure')
            im=ax.pcolormesh(self._lons,self._lats,self._MSLP[t,:,:],vmin=980,vmax=1020,transform=self._transform)
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[1]; ax.set_title('temperature')
            im=ax.pcolormesh(self._lons,self._lats,self._T[t,:,:],transform=self._transform)
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            ax=axes[2]; ax.set_title('wind speed')
            im=ax.pcolormesh(self._lons,self._lats,self._Wind10[t,:,:],vmin=0,vmax=30,transform=self._transform)
            im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')

            if self._VO is not None:
                ax=axes[3]; ax.set_title('rel. Vorticity')
                im=ax.pcolormesh(self._lons,self._lats,self._VO[t,:,:],vmin=-9.5*10**(-5),vmax=0.0002,transform=self._transform)
                im.set_cmap('bone'); ax.autoscale(False); ax.axis('off')
            else:
                ax=axes[3]; ax.set_title('wind speed [m/s] and mslp [mbar]')

            for point in self._detected[self._detected[:,'t']==t].values.tolist():
                if point[3]==1:
                    tmp.append(self.plot_on_map(axes[0],int(point[2]),int(point[1]),c='b',marker='.'))
                    stats='wind: '+str(round(point[6],01))+'\nmslp: '+str(round(point[5],01))
                    cat=self.tc_cat(point[5])
                    if cat>0:
                        txt.append(axes[3].text(self._lons[int(point[1]),int(point[2])],self._lats[int(point[1]),int(point[2])],stats,color=self._cat_colors[cat],va='bottom',ha='right',fontsize=8,transform=self._transform))
                if point[4]==1:
                    tmp.append(self.plot_on_map(axes[1],int(point[2]),int(point[1]),c='g',marker='*'))


            if self._obs_tc:
                obs_tc=np.where(abs(self._tc_time-self._yr_frac[t])<0.0004)
                if len(obs_tc[0])>0:
                    for oo in range(len(obs_tc[0])):
                        if np.isfinite(self._tc_sel['source_wind'].ix[obs_tc[0][oo],obs_tc[1][oo],0]):
                            tmp.append(axes[3].plot(self._tc_lon[obs_tc[0][oo],obs_tc[1][oo]],self._tc_lat[obs_tc[0][oo],obs_tc[1][oo]],color=self._cat_colors[self.tc_cat(self._tc_intens[obs_tc[0][oo],obs_tc[1][oo]],method='wind')],marker='.'))


            plt.suptitle(str(self._dates[t]))
            plt.savefig(self._working_dir+'track_surrounding/'+self._add_name+'_'+str(t)+'.png', bbox_inches = 'tight')

            # clean map
            for ax in axes:
                for imm in ax.images:
                    ax.images.remove(imm)
            for element in tmp:
                l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
            for element in txt:
                element.remove()

    # analyze fields
    def get_box(self,y,x,window):
        """
        Finds square around position

        Parameters
        ----------
            y: index
                center of square as an index of `self._lats` and `self._lons`
            x: index
                center of square as an index of `self._lats` and `self._lons`
            window:
                half of the length of the square

        Returns
        -------
            corners: tuple (y_min,y_max,x_min,x_max)
                corner coordinates of square
        """
        y_min=int(max(0,y-window))
        y_max=int(min(self._lats.shape[0],y+window+1))
        x_min=int(max(0,x-window))
        x_max=int(min(self._lats.shape[1],x+window+1))
        return (y_min,y_max,x_min,x_max)

    def area_around(self,y,x,radius):
        """
        Finds area around position as ponints within a circle defining the area

        Parameters
        ----------
            y: index
                center of circle as an index of `self._lats` and `self._lons`
            x: index
                center of circle as an index of `self._lats` and `self._lons`
            window:
                radius of the circle

        Returns
        -------
            indices:
                ya and x indices of points within circle defining the area
        """

        box=self.get_box(y,x,radius)
        y_,x_=[],[]
        for i in range(box[0],box[1]):
            for j in range(box[2],box[3]):
                if ((y-i)**2+(x-j)**2)**0.5<=radius:
                    y_.append(i)
                    x_.append(j)
        return y_,x_

    def circle_around(self,y,x,radius):
        """
        Finds circle around position

        Parameters
        ----------
            y: index
                center of circle as an index of `self._lats` and `self._lons`
            x: index
                center of circle as an index of `self._lats` and `self._lons`
            window:
                radius of the circle

        Returns
        -------
            indices:
                y and x indices of points on the circle
        """
        box=self.get_box(y,x,radius)
        y_,x_=[],[]
        for i in range(box[0],box[1]):
            for j in range(box[2],box[3]):
                if radius-1<=((y-i)**2+(x-j)**2)**0.5<=radius:
                    y_.append(i)
                    x_.append(j)
        return y_,x_

    def find_closed_contours(self,field,y,x,step=2,search_radius=30,n_contours=None,method='min'):
        '''
        This function uses opencv to detect closed contours around a given point

        Parameters
        ----------
            field: np.array
                field in which contours are searched
            y:  int [index]
                latitude index within field of point around which contours are searched
            x:  int [index]
                longitude index within field of point around which contours are searched
            step: float,default=2
                step is added to the search threshold after each identified contour
            search_radius: int [lat/lon steps],
                maximal allowed radius for contours. If a contour touches this limit, the contour isn't considered as closed. Given in lat/lon steps
            n_contours: int, default=None
                If this is not None, the algorithm stops searching after `n_contours` are found
            method: str ["min" or "max"],default="min"
                defines whether closed contours are searched around a minimum or a maximum

        Returns
        -------
            cont: np.array
                array containing a lat-lon mask for each found contour. Cells in contained in the contour are 1, other are nans.
            ncont: int
                number of found contours
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

    # combine detected positions
    def combine_tracks(self, search_radius=6, total_steps=12, warm_steps=8, strong_steps=0, consecutive_warm_strong_steps=6, thr_wind=17.5, smooth_path_angle=np.pi*2, velocity_jump=3., lat_formation_cutoff=40, overwrite=False, plot=True):

        """
        Combines detected storm snapshots from :meth:`detect_contours` or :meth:`detect_knutson2007` into storm tracks.

        In a first step, starting with the storm snapshot with lowest MSLP, storm snapshots in the following and preceding time steps are added to the track if they are located within the `search_radius`. If multiple snapshots are found within the `search_radius` the one with lowest MSLP is selected. If a landmask was provided, all storm snapshots located over land are deleted from the beginning of the track. Finally if the track fullfills a set of conditions, it is saved.

        Parameters
        ----------
            search_radius: float [deg],default=6
                maximal allowed distance between two storm snapshots of consecutive time steps
            total_steps: int, default=12
                minimal number of steps required for the track to be saved
            warm_steps: int, default=8
                minimal number of warm core steps required for the track to be saved
            strong_steps: int, default=0
                minimal number of time steps with wind speeds above `thr_wind` required for the track to be saved
            consecutive_warm_strong_steps: int, default=6
                minimal number of consecutive warm core steps with wind speeds above `thr_wind` required for the track to be saved
            thr_wind:   float [m/s],default=17.5
                threshold above which a storm snapshot is considered to be a strong storm (relevant for `strong_steps`)
            lat_formation_cutoff, float [deg], default=40
                maximal allowed latitude for storm formation
            overwrite: bool, default=False
                if False, existing outputfiles of the function are loaded, if True they are overwritten
            plot:  bool, default=True
                create a plot of each saved track

        Returns
        -------
            self._tcs: dimarray
                a dimarray containg time, position and characteristics of each track

        """

        out_file=self._working_dir+'track_info_'+self._add_name+'.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._tcs=da.read_nc(out_file)
            return self._tcs

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
        self._tcs={}

        postions=self._detected.copy().values
        used_pos=[]
        # # sort for decreasing 10m Wind speed
        # for p in postions[postions[:,-1].argsort()[::-1],:].tolist():
        # sort for minimal mslp
        for p in postions[postions[:,-2].argsort()[::-1],:].tolist():
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
                    candidates=[]
                    for p_1 in postions[postions[:,0]==p[0]-1,:].tolist():
                        new_step=np.array((p_1[1]-p[1],p_1[2]-p[2],0))
                        if v_len(new_step-prev_step)<search_radius: # and self._land_mask[int(p_1[1]),int(p_1[2])]
                            # # exclude strong direction changes if TC isn't extremely slow
                            # if v_len(prev_step)/v_len(new_step-prev_step)<1/velocity_jump and angle_between(new_step,prev_step)>smooth_path_angle:
                            #     # strange movement
                            #     pass
                            # else:
                            candidates.append(p_1)
                            end=False
                    if len(candidates)>0:
                        # select snapshot with lowest mslp
                        track=[candidates[np.array(candidates)[:,-2].argmin()]]+track
                    else:
                        break

                #go forewards
                while True:
                    p=track[-1]
                    if len(track)==1:
                        prev_step=np.array((0,0,0))
                    else:
                        prev_step=np.array((p[1]-track[-2][1],p[2]-track[-2][2],0))
                    candidates=[]
                    for p_1 in postions[postions[:,0]==p[0]+1,:].tolist():
                        new_step=np.array((p_1[1]-p[1],p_1[2]-p[2],0))
                        if v_len(new_step-prev_step)<search_radius:
                            # # exclude strong direction changes if TC isn't extremely slow
                            # if v_len(prev_step)/v_len(new_step-prev_step)<1/velocity_jump and angle_between(new_step,prev_step)>smooth_path_angle:
                            #     # strange movement
                            #     pass
                            # else:
                            candidates.append(p_1)
                            end=False
                    if len(candidates)>0:
                        # select snapshot with lowest mslp
                        track=track+[candidates[np.array(candidates)[:,-2].argmin()]]
                    else:
                        break

                # delete land steps in the beginning of track
                while len(track)>0:
                    if self._land_mask[int(track[0][1]),int(track[0][2])]:
                        track=track[1:]
                    else:
                        break

                # check whether the track should be saved
                if len(track)>0:
                    if sum([pp in used_pos for pp in track])/float(len(track))<0.3:
                        used_pos+=track

                        track=da.DimArray(track,axes=[np.array(track)[:,0],['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['time','z'])
                        save_track=True
                        start_pos=track.values[0,1:3]
                        if track.shape[0]<total_steps:
                            save_track=False
                        if track[track[:,'warm_core']==1].shape[0]<warm_steps:
                            save_track=False
                        elif track[track[:,'warm_core']==1].shape[0]!=0:
                            start_pos=track[track[:,'warm_core']==1].values[0,1:3]
                        if track[track[:,'Wind10']>=thr_wind].shape[0]<strong_steps:
                            save_track=False
                        elif track[track[:,'Wind10']>=thr_wind].shape[0]!=0:
                            start_pos=track[track[:,'Wind10']>=thr_wind].values[0,1:3]
                        if consecutive_warm_strong_steps>0:
                            warm_strong=track[(track[:,'Wind10']>=thr_wind) & (track[:,'warm_core']==1)]
                            consecutive=np.diff(warm_strong[:,'t'],1)==1
                            consec_info=consecutive_sequence(consecutive)
                            if max(consec_info[:,1])<consecutive_warm_strong_steps:
                                save_track=False
                            else:
                                first_of_consec=consec_info[np.argmax(consec_info[:,1]),0]
                                start_pos=warm_strong.values[first_of_consec,1:3]

                        if self._lats[int(start_pos[0]),int(start_pos[1])]>=lat_formation_cutoff:
                            save_track=False

                        if save_track:
                            self._tcs[self._identifier+'_'+str(self._id)]=track
                            if plot:    self.plot_track_path(track)
                            self._id+=1

        self._tcs=da.Dataset(self._tcs)
        self._tcs.write_nc(out_file,mode='w')
        return self._tcs

    # detect positions
    def detect_contours(self,overwrite=False,p_radius=27,p_inc_step=1,warm_core_size=3,T_drop_step=1,dis_cores=1,dis_mslp_min=3):
        '''
        This detection algorithm is a adapted from Murakami 2015.
        The method is based on the detection mean sea level pressure (MSLP) minima. If closed contours in MSLP are found around a minimum (using :meth:`find_closed_contours`), the minimum is taken as storm center. Within a radius of `dis_cores` the maximum in temperature anomaly (Ta) is suspected to be the warm core of the storm. If a closed  contour in Ta is found around the maximum in Ta, the storm is treated as a warm core storm snapshot.

        Parameters
        ----------
            dis_mslp_min:   float [deg],default=3
                minimal allowed distance between detected mslp minima
            p_radius:       float [deg],defualt=27
                maximum radius for mslp closed contours
            p_inc_step:     float [hPa],default=1
                step for mslp contours
            warm_core_size: float [deg],default=3
                maximal radius of warm core contours
            T_drop_step:    float [K],default=1
                step for Ta contours
            dis_cores:      float [deg],default=1
                maximal distances between vort max and mslp min (mslp min and warm core)
            overwrite:      bool,default=False
                if False, existing outputfiles of the function are loaded, if True they are overwritten

        Returns
        -------
            self._detected: dimarray
                time, position and characteristics of detected snapshots
        '''

        self._add_name='contours'
        out_file=self._working_dir+'detected_positions_'+self._add_name+'.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected

        # convert distances from degrees into grid-cells
        dis_mslp_min=self.degree_to_step(dis_mslp_min)
        p_radius=self.degree_to_step(p_radius)
        dis_cores=self.degree_to_step(dis_cores)
        warm_core_size=self.degree_to_step(warm_core_size)**2*2*np.pi

        detect=np.array([[np.nan]*7])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()
            coords=peak_local_max(-self._MSLP_smoothed[t,:,:], min_distance=int(dis_mslp_min))
            for y_p,x_p in zip(coords[:,0],coords[:,1]):
                tc_area,ncont=self.find_closed_contours(self._MSLP_smoothed[t,:,:],y_p,x_p,step=p_inc_step,search_radius=p_radius,method='min')
                if ncont>0:
                    tmp=[t,y_p,x_p,1,0,0,0]
                    # have to check boundary issues here
                    box=self.get_box(y_p,x_p,dis_cores)
                    y_,x_=np.where(self._T[t,box[0]:box[1],box[2]:box[3]]==self._T[t,box[0]:box[1],box[2]:box[3]].max())
                    y_t,x_t=box[0]+y_[0],box[2]+x_[0]
                    warm_core_area,ncont=self.find_closed_contours(self._T[t,:,:],y_t,x_t,step=T_drop_step,search_radius=p_radius,n_contours=1,method='max')
                    yy,xx=np.where(warm_core_area==1)
                    if len(np.where(warm_core_area==1)[0])<warm_core_size and ncont==1:
                        tmp[4]=1

                    tmp[5]=self._MSLP[t,y_p,x_p]
                    tmp[6]=np.nanmax(tc_area*self._Wind10[t,:,:])
                    if np.isnan(tmp[6]):
                        print(ncont)
                    detect=np.concatenate((detect,np.array([tmp])))

        self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('done')
        return self._detected

    def detect_knutson2007(self,thr_vort=3.5*10**(-5),dis_vort_max=4,dis_cores=2,thr_MSLP_inc=4,dis_MSLP_inc=5,thr_T_drop=0.8,dis_T_drop=5,tc_size=7,overwrite=False):
        '''
        Detection based on thresholds

        Parameters
        ----------
            thr_vort:       float [1/s],default=3.5*10**(-5)
                threshold for relative vorticity maxima
            dis_vort_max:   float [deg],default=4
                minimal distance between vorticity maxima
            dis_cores:      float [deg],default=2
                maximal distances between vort max and mslp min (mslp min and warm core)
            thr_MSLP_inc:   float [hPa],default=4
                increase in MSLP from center over dis_MSLP_inc
            dis_MSLP_inc:   float [deg],default=5
                distance over which MSLP should increase by thr_MSLP_inc
            thr_T_drop:     float [K],default=0.8
                temperature drop from warm core center over dis_T_decrease
            dis_T_drop:     float [deg],default=5
                distance over which T should decrease by thr_T_decrease
            tc_size:        float [deg],default=7
                radius within which maximal Wind10 is searched for
            overwrite:      bool,default=False
                if False, existing outputfiles of the function are loaded, if True they are overwritten

        Returns
        -------
            self._detected: dimarray
                time, position and characteristics of detected snapshots
        '''

        self._add_name='knutson2007'
        out_file=self._working_dir+'detected_positions_thresholds_'+self._add_name+'.nc'
        if overwrite and os.path.isfile(out_file):
            os.system('rm '+out_file)
        elif overwrite==False and os.path.isfile(out_file):
            self._detected=da.read_nc(out_file)['detected']
            return self._detected

        # convert distances from degrees into grid-cells
        dis_vort_max=int(self.degree_to_step(dis_vort_max))
        dis_cores=int(self.degree_to_step(dis_cores))
        dis_MSLP_inc=int(self.degree_to_step(dis_MSLP_inc))
        dis_T_drop=int(self.degree_to_step(dis_T_drop))
        tc_size=int(self.degree_to_step(tc_size))

        detect=np.array([[np.nan]*7])
        print('detecting\n10------50-------100')
        for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
            sys.stdout.write(progress); sys.stdout.flush()
            # i vort max
            vo_=self._VO[t,:,:].copy()
            vo_[vo_<thr_vort]=0
            coords = peak_local_max(vo_, min_distance=int(dis_vort_max))
            if coords.shape[0]>0:
                for y_v,x_v in zip(coords[:,0],coords[:,1]):
                    y_circ,x_circ=self.area_around(y_v,x_v,dis_cores)
                    p_window=self._MSLP[t,y_circ,x_circ].flatten()
                    i=np.where(p_window==p_window.min())[0][0]; y_p,x_p=y_circ[i],x_circ[i]
                    for dis in range(1,dis_MSLP_inc+1):
                        y_circ,x_circ=self.circle_around(y_p,x_p,dis)
                        # ii relative pressure min
                        # could also be self._MSLP[t,y_circ,x_circ].mean()
                        if self._MSLP[t,y_circ,x_circ].min()-self._MSLP[t,y_p,x_p]>thr_MSLP_inc:
                            tmp=[t,y_p,x_p,1,0,0,0]
                            y_circ,x_circ=self.area_around(y_p,x_p,dis_cores)
                            # iv warm core
                            t_window=self._T[t,y_circ,x_circ].flatten()
                            i=np.where(t_window==t_window.max())[0][0]; y_t,x_t=y_circ[i],x_circ[i]
                            for dis in range(1,dis_T_drop+1):
                                y_circ,x_circ=self.circle_around(y_t,x_t,dis)
                                if self._T[t,y_t,x_t]-self._T[t,y_circ,x_circ].max()>thr_T_drop:
                                    tmp[4]=1
                            # iii wind speed
                            tmp[5]=self._MSLP[t,y_p,x_p]
                            y_circ,x_circ=self.area_around(y_p,x_p,tc_size)
                            tmp[6]=self._Wind10[t,y_circ,x_circ].max()
                            detect=np.concatenate((detect,np.array([tmp])))
                            break

        self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
        da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
        print('done')
        return self._detected

    # pld stuff
    # def detect_thresholds_simple(self,thr_vort=3.5*10**(-5),dis_vort_max=4,dis_cores=2,thr_MSLP_inc=2,dis_MSLP_inc=5,thr_T_drop=1,dis_T_drop=3,tc_size=5,overwrite=False):
    #     self._add_name='thresholds'
    #     out_file=self._working_dir+'detected_positions_thresholds_'+self._add_name+'.nc'
    #     if overwrite and os.path.isfile(out_file):
    #         os.system('rm '+out_file)
    #     elif overwrite==False and os.path.isfile(out_file):
    #         self._detected=da.read_nc(out_file)['detected']
    #         return self._detected
    #
    #     '''
    #     Required Thresholds:
    #     thr_vort:       float [1/s]:    threshold for relative vorticity maxima
    #     dis_vort_max:   float [deg]:    minimal distance between vorticity maxima
    #     dis_cores:      float [deg]:    maximal distances between vort max and mslp min (mslp min and warm core)
    #     dis_MSLP_inc:   float [deg]:    distance over which MSLP should increase by thr_MSLP_inc
    #     thr_T_drop:     float [K]:      temperature drop from warm core center over dis_T_decrease
    #     dis_T_drop:     float [deg]:    distance over which T should decrease by thr_T_decrease
    #     tc_size:        float [deg]:    radius within which maximal Wind10 is searched for
    #     '''
    #
    #     # convert distances from degrees into grid-cells
    #     for distance in [dis_vort_max,dis_cores,dis_MSLP_inc,dis_T_drop,tc_size]:
    #         distance=self.degree_to_step(distance)
    #
    #     detect=np.array([[np.nan]*7])
    #     print('detecting\n10------50-------100')
    #     for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
    #         sys.stdout.write(progress); sys.stdout.flush()
    #         # i vort max
    #         vo_=self._VO[t,:,:].copy()
    #         vo_[vo_<thr_vort]=0
    #         coords = peak_local_max(vo_, min_distance=int(dis_vort_max))
    #         if coords.shape[0]>0:
    #             for y_v,x_v in zip(coords[:,0],coords[:,1]):
    #                 y_circ,x_circ=self.area_around(y_v,x_v,dis_cores)
    #                 p_window=self._MSLP[t,y_circ,x_circ].flatten()
    #                 i=np.where(p_window==p_window.min())[0][0]; y_p,x_p=y_circ[i],x_circ[i]
    #                 y_area,x_area=self.area_around(y_p,x_p,dis_MSLP_inc)
    #                 # ii relative pressure min
    #                 if self._MSLP[t,y_p,x_p]==self._MSLP[t,y_area,x_area].min():
    #                     box_1=self.get_box(y_p,x_p,cores_distance)
    #                     box_2=self.get_box(y_p,x_p,tc_size)
    #                     tmp=[t,y_p,x_p,1,0,0,0]
    #                     # iv warm core
    #                     if self._T is None:
    #                         tmp[4]=1
    #                     elif self._T[t,1,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()-self._T[t,1,box_2[0]:box_2[1],box_2[2]:box_2[3]].mean()>thr_T_drop:
    #                         tmp[4]=1
    #                     # iii wind speed
    #                     tmp[5]=self._MSLP[t,y_p,x_p]
    #                     tmp[6]=self._Wind10[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()
    #                     detect=np.concatenate((detect,np.array([tmp])))
    #
    #     self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
    #     da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
    #     print('done')
    #     return self._detected
    #
    # def detect_hybrid(self,overwrite=False,thr_vort=5*10**(-5),dis_vort_max=4,p_radius=27,p_inc_step=1,warm_core_size=3,T_drop_step=1,dis_cores=1,dis_mslp_min=3):
    #     self._add_name='hybrid'
    #     out_file=self._working_dir+'detected_positions_'+self._add_name+'.nc'
    #     if overwrite and os.path.isfile(out_file):
    #         os.system('rm '+out_file)
    #     elif overwrite==False and os.path.isfile(out_file):
    #         self._detected=da.read_nc(out_file)['detected']
    #         return self._detected
    #
    #     '''
    #     Required Thresholds:
    #     dis_vort_max:   float [deg]:    minimal distance between vort maxima
    #     p_radius:       float [deg]:    maximum radius for mslp closed contours
    #     p_inc_step:     float [hPa]:    step for mslp contours
    #     warm_core_size: float [deg]:    maximal radius of warm core contours
    #     T_drop_step:    float [K]:      step for Ta contours
    #     dis_cores:      float [deg]:    maximal distances between vort max and mslp min (mslp min and warm core)
    #     '''
    #
    #     # convert distances from degrees into grid-cells
    #     p_radius=self.degree_to_step(p_radius)
    #     dis_cores=self.degree_to_step(dis_cores)
    #     warm_core_size=self.degree_to_step(warm_core_size)**2*2*np.pi
    #     dis_vort_max=int(self.degree_to_step(dis_vort_max))
    #
    #
    #     detect=np.array([[np.nan]*7])
    #     print('detecting\n10------50-------100')
    #     for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
    #         sys.stdout.write(progress); sys.stdout.flush()
    #
    #         vo_=self._VO[t,:,:].copy()
    #         vo_[vo_<thr_vort]=0
    #         coords = peak_local_max(vo_, min_distance=int(dis_vort_max))
    #         if coords.shape[0]>0:
    #             for y_v,x_v in zip(coords[:,0],coords[:,1]):
    #                 y_circ,x_circ=self.area_around(y_v,x_v,dis_cores)
    #                 # get small area around vort max
    #                 p_window=self._MSLP_smoothed[t,y_circ,x_circ].flatten()
    #                 # MSLP min next to vort max is center
    #                 i=np.where(p_window==p_window.min())[0][0]; y_p,x_p=y_circ[i],x_circ[i]
    #                 # find closed contours around MSLP min
    #                 tc_area,ncont=self.find_closed_contours(self._MSLP_smoothed[t,:,:],y_p,x_p,step=p_inc_step,search_radius=p_radius,method='min')
    #                 if ncont>0:
    #                     tmp=[t,y_p,x_p,1,0,0,0]
    #                     # find warm core center
    #                     box=self.get_box(y_p,x_p,dis_cores)
    #                     y_,x_=np.where(self._T[t,box[0]:box[1],box[2]:box[3]]==self._T[t,box[0]:box[1],box[2]:box[3]].max())
    #                     y_t,x_t=box[0]+y_[0],box[2]+x_[0]
    #                     # check if warm core exists
    #                     warm_core_area,ncont=self.find_closed_contours(self._T[t,:,:],y_t,x_t,step=T_drop_step,search_radius=p_radius,n_contours=1,method='max')
    #                     yy,xx=np.where(warm_core_area==1)
    #                     if len(np.where(warm_core_area==1)[0])<warm_core_size and ncont==1:
    #                         tmp[4]=1
    #
    #                     tmp[5]=self._MSLP[t,y_p,x_p]
    #                     # tc area is defined by closed contours
    #                     tmp[6]=np.nanmax(tc_area*self._Wind10[t,:,:])
    #                     detect=np.concatenate((detect,np.array([tmp])))
    #
    #     self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
    #     da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
    #     print('done')
    #     return self._detected

    # def set_thresholds(self,thr_wind,thr_mslp,p_radius,neighborhood_size,warm_core_size,cores_distance,search_radius,min_time_steps):
    #     self._thr_wind=thr_wind
    #     self._thr_mslp=thr_mslp
    #     self._min_time_steps=min_time_steps
    #
    #     self._p_radius=self.degree_to_step(p_radius)
    #     self._neighborhood_size=self.degree_to_step(neighborhood_size)
    #     self._warm_core_size=self.degree_to_step(warm_core_size)**2
    #     self._cores_distance=self.degree_to_step(cores_distance)
    #     self._search_radius=self.degree_to_step(search_radius)
