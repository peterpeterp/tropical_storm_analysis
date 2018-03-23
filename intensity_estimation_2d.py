
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
    for ax_i,vari in zip(range(4),['MSLP','VO','Wind10','T']):
        y=obs_stats[:,vari].values
        x=obs_stats[:,var_obs].values

        y=y[np.isfinite(x)]
        x=x[np.isfinite(x)]

        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)

        xmin,xmax=np.percentile(x,[0,100])
        ymin,ymax=np.percentile(y,[0,100])
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])

        Z = np.reshape(kernel(positions).T, X.shape)

        pdfs[var_obs][vari]={'X':X,'Y':Y,'Z':Z,'y':y,'x':x}

plt.close('all')
fig,axes=plt.subplots(nrows=4,ncols=2,figsize=(8,10))
for ax_j,var_obs in zip(range(2),['obs_wind','obs_pres']):
    for ax_i,vari in zip(range(4),['MSLP','VO','Wind10','T']):
        ax=axes[ax_i,ax_j]
        pdf=pdfs[var_obs][vari]
        xmin,xmax=np.percentile(pdf['x'],[0,100])
        ymin,ymax=np.percentile(pdf['y'],[0,100])
        ax.plot(pdf['x'], pdf['y'],color='gray',marker='.',linestyle='',alpha=0.5, markersize=1)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.pcolormesh(pdf['X'],pdf['Y'],pdf['Z'], cmap=plt.cm.gist_earth_r)

        if (vari=='MSLP' and var_obs=='obs_pres') or (vari=='Wind10' and var_obs=='obs_wind'):
            ax.plot([0,99999],[0,99999],'k')

        ax.set_ylabel(vari)
        ax.set_xlabel(var_obs)

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_hists_2D.png')





cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
cat_thr={'obs_wind':{0:0,1:64,2:82,3:95,4:112,5:136,6:900},
        'obs_pres':{0:2000,1:1020,2:980,3:965,4:945,5:920,6:0}}


cat4_tracks=np.where(np.nanmax(obs_tracks[:,:,'cat'],axis=1)==4)[0]
track=obs_tracks.ix[cat4_tracks[0],:,:]
track=track[np.isfinite(track[:,'VO']),:]

# detail plot
#for var_obs,var_primitiv in zip(['obs_pres','obs_wind'],['MSLP','Wind10']):
for var_obs,var_primitiv in zip(['obs_pres'],['MSLP']):
    plt.close('all')
    fig,axes=plt.subplots(nrows=4,ncols=3,figsize=(10,10))
    for ax in axes[1:,2]:
        ax.axis('off')

    for ax in [axes[0,1],axes[1,1],axes[2,1],axes[3,1],axes[0,2]]:
        for cat in range(0,6):
            ax.axvspan(cat_thr[var_obs][cat], cat_thr[var_obs][cat+1], alpha=0.7, color=cat_colors[cat])

    for ax_i,vari in zip(range(4),['MSLP','VO','Wind10','T']):
        # left column
        ax=axes[ax_i,0]
        pdf=pdfs[var_obs][vari]
        ymin,ymax=np.percentile(pdf['y'],[0,100])
        xmin,xmax=np.percentile(pdf['x'],[0,100])
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.pcolormesh(pdf['X'],pdf['Y'],pdf['Z'], cmap=plt.cm.gist_earth_r)
        ax.set_ylabel(vari)
        ax.set_xlabel(var_obs)

    for t in track.time:
        x=pdfs[var_obs]['VO']['X'][:,0]
        xmin,xmax=np.percentile(x,[0,100])
        result=np.zeros([x.shape[0],4])*np.nan
        for ax_i,vari in zip(range(4),['MSLP','VO','Wind10','T']):
            pdf=pdfs[var_obs][vari]
            # left column
            axes[ax_i,0].axhline(y=track[t,vari],color='k')

            # right column
            ax=axes[ax_i,1]
            yy=np.argmin(abs(pdf['Y'][0,:]-track[t,vari]))
            ax.plot(x,pdf['Z'][:,yy],color='k')
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([pdf['Z'][:,yy].min(), pdf['Z'][:,yy].max()])

            ax.axvline(x[np.argmax(pdf['Z'][:,yy])],color='k')
            ax.axvline(track[t,var_obs],color='k')
            ax.axvline(track[t,var_primitiv],color='k')
            ax.set_xlabel('wind estimate')
            ax.set_ylabel(vari)

            result[:,ax_i]=pdf['Z'][:,yy]

        # result
        ax=axes[0,2]
        #yyy=np.sum(result,axis=1)/np.sum(result)
        yyy=np.nanmean(result,axis=1)/np.sum(np.nanmean(result,axis=1))
        ax.plot(x,yyy,color='k')
        ax.axvline(x[np.argmax(yyy)],color='k')
        ax.axvline(track[t,var_obs],color='k')
        ax.axvline(track[t,var_primitiv],color='k')
        ax.set_ylim([yyy.min(),yyy.max()])
        ax.set_xlim([xmin,xmax])

        # probabilities
        ax=axes[1,2];   ax.clear(); ax.axis('off')
        for cat in range(1,6):
            ax.bar([cat],[np.sum(yyy[(cat_thr[var_obs][cat+1]<x) & (x<cat_thr[var_obs][cat])])/np.sum(yyy)],color='m')
        ax.plot(track[t,'cat']+1,0.8,'*g')
        ax.plot(np.where([track[t,var_primitiv]<cat_thr[var_obs][ccc] for ccc in range(1,6)])[0][-1]+1,0.7,'^b')
        ax.set_xlim((0,6))
        ax.set_ylim((0,1))

        plt.tight_layout()
        plt.savefig('detection/JRA55/intens_'+var_obs+'/JRA55_obs_hists_'+str(t)+'.png')
        for ax in axes.flatten():
            ax.lines=[]
