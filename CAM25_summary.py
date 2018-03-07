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
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False


identifiers=[ff.split('_')[-3] for ff in glob.glob(data_path+'/item3225_daily_mean/item3225_daily*')]


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
