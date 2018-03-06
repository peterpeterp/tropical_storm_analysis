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

from cdo import *   # python version

sns.set_palette(sns.color_palette("plasma"))

try:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    data_path='data/CAM25/'
    cdo = Cdo()
except:
    os.chdir('/p/projects/tumble/carls/shared_folder/TC_detection/')
    data_path='/p/projects/tumble/carls/shared_folder/CPDN/data/batch_659/region/'
sys.path.append('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/tc_detection')
sys.path.append('/p/projects/tumble/carls/shared_folder/TC_detection/tc_detection')
from TC_support import *


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--verbosity",'-v', help="increase output verbosity",action="store_true")
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
parser.add_argument("--plot",'-p', help="plot output",action="store_true")
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False

if args.plot:
    plot=True
else:
    plot=False

if overwrite:
    # check for duplicates
    xxx,storms=[],[]
    found_tracks={}
    longest_track=0
    for identifier in identifiers:
        tmp=da.read_nc('detection/CAM25/'+str(identifier)+'_CAM25/track_info.nc')
        if len(tmp.values())>0:
            tmp_example=tmp
            for id_,track in tmp.items():
                track=np.array(track[np.isfinite(track[:,'t']),:])
                unique=True
                x_=[int(xx) for xx in track[:,2]]
                if x_ in xxx:
                    used=storms[xxx.index(x_)]
                    if [int(yy) for yy in track[:,1]]==[int(yy) for yy in found_tracks[used][:,1]]:
                        if [int(yy) for yy in track[:,0]]==[int(yy) for yy in found_tracks[used][:,0]]:
                            unique=False

                if unique:
                    xxx.append(x_)
                    storms.append(id_)
                    found_tracks[id_]=track
                    if track.shape[0]>longest_track:
                        longest_track=track.shape[0]


    all_tracks=da.DimArray(np.zeros([len(found_tracks.keys()),longest_track,13])*np.nan,axes=[found_tracks.keys(),range(longest_track),tmp_example.z],dims=['ID','time','z'])
    for id_,track in found_tracks.items():
        all_tracks[id_,0:track.shape[0]-1,:]=track
    da.Dataset({'all_tracks':all_tracks}).write_nc('detection/CAM25/CAM25_all_tracks.nc',mode='w')

else:
    all_tracks=da.read_nc('detection/CAM25/CAM25_all_tracks.nc')['all_tracks']


if plot:
    summary={'cat':[],'duration':[]}
    for id_ in all_tracks.ID:
        summary['cat'].append(np.nanmax(all_tracks[id_,:,'cat']))
        summary['duration'].append(len(np.where((all_tracks[id_,:,'cat']>0) & (np.isfinite(all_tracks[id_,:,'cat'])))[0]))

    runs=set([id_.split('_')[0] for id_ in all_tracks.ID])
    n_runs=len(runs)

    # intensity relations
    plt.close('all')
    fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(6,2.5))
    tmp=[]
    for id_ in all_tracks.ID:
        track=all_tracks[id_]
        track=track[np.isfinite(track[:,'t']),:]
        tmp+=[axes[0].plot(range(track.shape[0]),track[:,'MSLP'],c='k')]
        tmp+=[axes[1].plot(range(track.shape[0]),track[:,'Wind10'],c='k')]


    axes[0].set_ylim((91000,101000))
    axes[1].set_ylim((5,55))
    #
    axes[0].set_xlim((0,14))
    axes[1].set_xlim((0,14))

    plt.tight_layout(rect=(0,0,1,0.95))

    for storm in tmp:
        storm[0].set_alpha(0.1)
        storm[0].set_linewidth(1)

    plt.savefig('detection/CAM25/CAM25_intensities.png')


    # show summary plots
    plt.close('all')
    fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(7,4))
    ax=axes[0]
    for i in range(1,6):
        ax.bar([i-0.5],len(np.where(np.array(summary['cat'])==i)[0])/float(len(all_tracks.ID)),width=1)
    ax.set_xlim(0.5,5.5)
    ax.set_ylabel('probability density')
    ax.set_yscale('log')
    ax.set_xlabel('TC Category')

    ax=axes[1]
    tmp=np.histogram(summary['duration'])
    ax.hist(summary['duration'],10,normed=True)
    ax.set_ylabel('probability density')
    ax.set_xlabel('TC duration')

    plt.suptitle(str(len(all_tracks.ID))+' Tropical cyclones detected in CAM25')
    plt.tight_layout(rect=(0,0,1,0.95))
    plt.savefig('detection/CAM25/CAM25_summary.png')

    # show tracks on map
    nc = Dataset(data_path+'/item16222_daily_mean/item16222_daily_mean_o6uk_2017-06_2017-10.nc')
    lats = nc.variables['global_latitude0'][:]
    lons = nc.variables['global_longitude0'][:]
    rlats = nc.variables['latitude0'][:]
    rlons = nc.variables['longitude0'][:]
    rlons, rlats = np.meshgrid(rlons, rlats)
    data = nc.variables['item16222_daily_mean'][:,0,:,:].squeeze()
    data = np.ma.masked_values(data,-999.)
    rotpole = nc.variables['rotated_pole0']

    lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
    o_lon_p = rotpole.grid_north_pole_longitude
    o_lat_p = rotpole.grid_north_pole_latitude

    plt.close('all')
    tmp,txt=[],[]
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    m= Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
               llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
               llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
    x,y = m(lons,lats)
    m.drawmapboundary(fill_color='1.')
    m.drawmapboundary(fill_color='darkblue')
    m.fillcontinents(color='darkgreen',lake_color='darkblue')
    m.drawcoastlines(linewidth=0.3)
    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
    ax.invert_yaxis()
    ax.set_title('CAM25 '+str(n_runs)+' runs')

    for id_ in all_tracks.ID:
        track=all_tracks[id_]
        track=track[np.isfinite(track[:,'t']),:]
        x_=[int(xx) for xx in track[:,'x']]
        y_=[int(yy) for yy in track[:,'y']]
        x,y = m(lons[y_,x_],lats[y_,x_])
        tmp+=tc_plot(m,x,y,track[:,'MSLP'].values,color_type=tc_pressure_color,alpha=0.4)


    plt.tight_layout()
    # for storm in tmp:
    #     storm[0].set_alpha(0.1)
    #     #tmp.set_linewidth(0.001)

    plt.savefig('detection/CAM25/CAM25_summary_map.png')





    ############################
    # check for real duplicates
    ############################
        # # check for duplicates
        # xxx=[]
        # storms=[]
        # useful_runs=[]
        # not_unique={}
        # for identifier in identifiers:
        #     tmp=da.read_nc('detection/'+str(identifier)+'_CAM25/track_info.nc')
        #     if len(tmp.values())>0:
        #         for track in tmp.values():
        #             track=np.array(track[np.isfinite(track[:,'t']),:])
        #             x_=[int(xx) for xx in track[:,2]]
        #             if x_ in xxx:
        #                 used=storms[xxx.index(x_)]
        #                 cdo_diff=cdo.diff(input=data_path+'/item16222_daily_mean/item16222_daily_mean_'+used+'_2017-06_2017-10.nc'+' '+data_path+'/item16222_daily_mean/item16222_daily_mean_'+identifier+'_2017-06_2017-10.nc')
        #                 if len(cdo_diff)<150:
        #                     if len(cdo_diff)==72:
        #                         print('             ',identifier,used,len(cdo_diff))
        #                         break
        #                     elif len(cdo_diff)==0:
        #                         print('------------>',identifier,used,len(cdo_diff))
        #                         if used in not_unique.keys():
        #                             not_unique[used].append(identifier)
        #                         if used not in not_unique.keys():
        #                             not_unique[used]=[identifier]
        #                         break
        #
        #
        #                 if len(cdo_diff)>150:
        #                     xxx.append(x_)
        #                     storms.append(identifier)
        #                     useful_runs.append(identifier)
        #             if x_ not in xxx:
        #                 xxx.append(x_)
        #                 storms.append(identifier)
        #                 useful_runs.append(identifier)
        #
        # print(not_unique)
        # not_unique_summary=open('detection/CAM25_not_unique.txt','w')
        # for used,identic in not_unique.items():
        #     not_unique_summary.write(used+' '+' '.join(identic)+'\n')
        # not_unique_summary.close()
