
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





plt.close('all')
tc_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
fig,axes=plt.subplots(nrows=2,ncols=2)
axes=axes.flatten()
pdfs={}
for ax,vari in zip(axes,['MSLP','VO','Wind10','T']):
    pdfs[vari]={}
    tmp=obs_stats[obs_stats[:,'cat']>0,vari]
    tmp=tmp[np.isfinite(tmp)]
    kde = stats.gaussian_kde(tmp)
    start,stop=np.percentile(tmp,[0,100])
    x=np.linspace(start,stop,1000)
    y=kde.evaluate(x)/np.sum(kde.evaluate(x))
    pdfs[vari]['all']={'pd':y,'x':x,'dx':np.diff(x,1)[0]}
    #ax.plot(x,y,color='green')
    for cat in range(1,6):
        tmp=obs_stats[obs_stats[:,'cat']==cat,vari]
        tmp=tmp[np.isfinite(tmp)]
        kde = stats.gaussian_kde(tmp,bw_method='silverman')#'silverman'
        start,stop=np.percentile(tmp,[0,100])
        #x=np.linspace(start,stop,10)
        y=kde.evaluate(x)/np.sum(kde.evaluate(x))
        pdfs[vari][cat]={'pd':y,'x':x,'dx':np.diff(x,1)[0]}
        # hist=np.histogram(tmp,bins=10,density=True)
        # ax.plot(hist[1][:-1],hist[0],color=tc_colors[cat])
        ax.plot(x,y,color=tc_colors[cat])
    ax.set_xlabel(vari)
    ax.set_facecolor('xkcd:black')

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_hists_1D.png')

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
            print('***')
            ax.plot([0,99999],[0,99999],'k')

        ax.set_ylabel(vari)
        ax.set_xlabel(var_obs)

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_hists_2D.png')



cat4_tracks=np.where(np.nanmax(obs_tracks[:,:,'cat'],axis=1)==4)[0]
track=obs_tracks.ix[cat4_tracks[0],:,:]
track=track[np.isfinite(track[:,'VO']),:]

fcat=np.array([0.0]*6)

for t in track.time:
    vo_=np.argmin(abs(pdfs['VO'][cat]['pd']-track[t,'VO']))
    mslp_=np.argmin(abs(pdfs['MSLP'][cat]['pd']-track[t,'MSLP']))
    axes[0].axvline(x=track[t,'MSLP'],color='cyan')
    axes[1].axvline(x=track[t,'VO'],color='cyan')
    for cat in range(1,6):
        p_vo_cat = pdfs['VO'][cat]['pd'][vo_] * pdfs['MSLP'][cat]['pd'][mslp_]
        normalize = sum([pdfs['VO'][ccc]['pd'][np.argmin(abs(pdfs['VO'][ccc]['pd']-track[t,'VO']))] * pdfs['MSLP'][ccc]['pd'][np.argmin(abs(pdfs['MSLP'][ccc]['pd']-track[t,'MSLP']))] for ccc in range(1,6)])
        #print(cat,p_vo_cat,normalize,p_vo_cat / normalize,fcat[cat])
        fcat[cat] = p_vo_cat / normalize
    print(fcat,np.argmax(fcat))
    plt.savefig('detection/JRA55/intensity/JRA55_obs_hists_'+str(t)+'.png')
    for ax in axes:
        ax.lines=ax.lines[:5]
#
# vo=track[:,'MSLP']
# dvo=20
# fcat=np.array([0.2]*6)
#
# for t in track.time:
#     #print('*******************')
#     for cat in range(1,6):
#         p_vo_cat = np.sum((obs_stats[obs_stats[:,'cat']==cat,'MSLP']>vo[t]-dvo) & (obs_stats[obs_stats[:,'cat']==cat,'MSLP']<vo[t]+dvo))/float(obs_stats[obs_stats[:,'cat']==cat,'MSLP'].shape[0])
#         normalize = sum([np.sum((obs_stats[obs_stats[:,'cat']==ccc,'MSLP']>vo[t]-dvo) & (obs_stats[obs_stats[:,'cat']==ccc,'MSLP']<vo[t]+dvo))/float(obs_stats[obs_stats[:,'cat']==ccc,'MSLP'].shape[0]) for ccc in range(1,6)])
#         #print(cat,p_vo_cat,normalize,p_vo_cat / normalize,fcat[cat])
#         fcat[cat] = p_vo_cat / normalize
#     print(fcat)
#     print(np.argmax(fcat),track[t,'cat'])
#
