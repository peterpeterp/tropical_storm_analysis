
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
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

obs_tracks=da.read_nc('detection/JRA55/JRA55_obs_track_info.nc')['obs_tracks']


categories={}
obs_stats=obs_tracks.values.reshape((obs_tracks.shape[0]*obs_tracks.shape[1],obs_tracks.shape[2]))
obs_stats=obs_stats[np.isfinite(obs_stats[:,0]),:]
obs_stats=obs_stats[np.isfinite(obs_stats[:,-1]),:]

obs_stats=da.DimArray(obs_stats,axes=[range(obs_stats.shape[0]),obs_tracks.variable],dims=['id','variable'])
obs_stats[obs_stats[:,'obs_pres']<500,'obs_pres']=np.nan

pdfs={}
for ax_j,var_obs in zip(range(2),['obs_wind','obs_pres']):
    pdfs[var_obs]={}
    for ax_i,vari in zip(range(5),['MSLP','VO','Wind10','T','T_diff']):
        y=obs_stats[:,vari].values
        x=obs_stats[:,var_obs].values
        finite=np.where(np.isfinite(x) & np.isfinite(y))
        y=y[finite]
        x=x[finite]
        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        xmin,xmax=np.percentile(x,[0,100])
        ymin,ymax=np.percentile(y,[0,100])
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel(positions).T, X.shape)
        pdfs[var_obs][vari]={'X':X,'Y':Y,'Z':Z,'y':y,'x':x}

plt.close('all')
fig,axes=plt.subplots(nrows=1,ncols=2,figsize=(8,4))
for ax_j,var_obs,vari,xlabel,ylabel in zip(range(2),['obs_wind','obs_pres'],['Wind10','MSLP'],['observed 10min sustained wind speed [m/s]','observed minimum central pressure [hPa]'],['6hourly mean 10m wind speed [m/s]','6hourly mean sea level pressure [hPa]']):
    ax=axes[ax_j]
    pdf=pdfs[var_obs][vari]
    xmin,xmax=np.percentile(pdf['x'],[0,100])
    ymin,ymax=np.percentile(pdf['y'],[0,100])
    ax.plot(pdf['x'], pdf['y'],color='gray',marker='.',linestyle='',alpha=0.8, markersize=0.2)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.pcolormesh(pdf['X'],pdf['Y'],pdf['Z'], cmap=plt.cm.gist_earth_r)

    if (vari=='MSLP' and var_obs=='obs_pres') or (vari=='Wind10' and var_obs=='obs_wind'):
        ax.plot([0,99999],[0,99999],'k')

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_intensity.png',dpi=200)
