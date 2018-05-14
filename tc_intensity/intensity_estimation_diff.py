
from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
from scipy import stats
import numpy as np
from scipy.optimize import leastsq
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

obs_tracks=da.read_nc('detection/JRA55/JRA55_obs_track_info.nc')['obs_tracks']
obs_tracks_diff=np.diff(obs_tracks,axis=1)
obs_tracks_diff=da.DimArray(obs_tracks_diff,axes=[obs_tracks.unnamed,obs_tracks.time[1:],obs_tracks.variable],dims=['unnamed','time','variable'])

obs_stats_diff=obs_tracks_diff.values.reshape((obs_tracks_diff.shape[0]*obs_tracks_diff.shape[1],obs_tracks_diff.shape[2]))
obs_stats_diff=obs_stats_diff[np.isfinite(obs_stats_diff[:,0]),:]
obs_stats_diff=obs_stats_diff[np.isfinite(obs_stats_diff[:,-1]),:]

obs_stats_diff=da.DimArray(obs_stats_diff,axes=[range(obs_stats_diff.shape[0]),obs_tracks.variable],dims=['id','variable'])
obs_stats_diff=obs_stats_diff[np.isfinite(obs_stats_diff[:,'obs_pres']),:]

obs_stats_norm=da.DimArray(obs_stats_diff.copy(),axes=[range(obs_stats_diff.shape[0]),obs_tracks.variable],dims=['id','variable'])
obs_tracks_norm=da.DimArray(obs_tracks_diff.copy(),axes=[obs_tracks.unnamed,obs_tracks.time[1:],obs_tracks.variable],dims=['unnamed','time','variable'])
for var in obs_stats_diff.variable:
    tmp=obs_stats_diff[:,var].copy()
    obs_stats_norm[:,var]=(tmp-tmp.mean())/tmp.std()
    obs_tracks_norm[:,:,var]=(obs_tracks_diff[:,:,var]-tmp.mean())/tmp.std()

plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(8,8))
for ax,vari in zip(axes.flatten(),['MSLP','VO','Wind10','T']):
    x=obs_stats_norm[:,'obs_pres'].values
    y=obs_stats_norm[:,vari].values
    xmin,xmax=np.percentile(x,[0,100])
    ymin,ymax=np.percentile(y,[0,100])
    ax.plot(x, y,color='gray',marker='.',linestyle='',alpha=0.5, markersize=1)
    ax.set_xlim([-1, 1])
    # ax.set_ylim([ymin, ymax])

    ax.set_ylabel(vari)
    ax.set_xlabel('observed pressure')

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_diff_norm.png')


coefficient = obs_stats_norm[:,['MSLP','VO','Wind10']].values
dependent = obs_stats_norm[:,'obs_pres'].values

def linear_model(p,x):
    a,b,c,d = p
    mslp = x[:,0]
    vo = x[:,1]
    wind = x[:,2]
    return (mslp*a+vo*b+wind*c+d)


def linear_residuals(p, y, x):
    err = y - linear_model(p,x)
    return err

p_linear = leastsq(linear_residuals, np.array([1,1,1,1]), args=(dependent, coefficient))[0]


def new_model(p,x):
    a1,a2,b1,b2,c1,c2,d = p
    mslp = x[:,0]
    vo = x[:,1]
    wind = x[:,2]
    return (a1*mslp**2+a2*mslp+b1*vo**2+b2*vo+c1*wind**2+c2*wind+d)

def new_residuals(p, y, x):
    err = y - new_model(p,x)
    return err

p_new = leastsq(new_residuals, np.array([1,1,1,1,1,1,1]), args=(dependent, coefficient))[0]




obs_pres_mean=obs_stats[:,'obs_pres'].mean()
obs_pres_std=obs_stats[:,'obs_pres'].std()

cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}

strong_tcs=np.where(np.nanmax(obs_tracks[:,:,'cat'],axis=1)>3)[0]
medium_tcs=np.where((np.nanmax(obs_tracks[:,:,'cat'],axis=1)<4) & (np.nanmax(obs_tracks[:,:,'cat'],axis=1)>1))[0]
weak_tcs=np.where(np.nanmax(obs_tracks[:,:,'cat'],axis=1)<2)[0]
plt.close('all')
fig,axes=plt.subplots(nrows=18,ncols=10,figsize=(30,20))
for axx,tcs,title in zip([axes.flatten()[0:60],axes.flatten()[60:120],axes.flatten()[120:180]],[strong_tcs,medium_tcs,weak_tcs],['Cat 4 and 5','Cat 2 and 3','Cat 0 and 1']):
    for ax,storm in zip(axx,[tcs[int(rr*len(tcs))] for rr in np.random.sample(60)]):
        track=obs_tracks.ix[storm,:,:]
        track_norm=obs_tracks_norm.ix[storm,:,:]
        ax.plot(track[:,'obs_pres'],'k-')
        ax.plot(track[:,'MSLP'],'k--')
        ax.plot(linear_model(p_linear,track_norm[:,['MSLP','VO','Wind10']].values)*obs_pres_std+obs_pres_mean,'k:')
        #ax.axhline(y=0,color='gray')
        ax.set_facecolor(cat_colors[np.nanmax(track[:,'cat'])])

plt.tight_layout(h_pad=0,w_pad=0)
plt.savefig('detection/JRA55/JRA55_obs_linear')







###
