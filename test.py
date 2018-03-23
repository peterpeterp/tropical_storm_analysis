
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
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

obs_tracks=da.read_nc('detection/JRA55/JRA55_obs_track_info.nc')['obs_tracks']


categories={}
obs_stats=obs_tracks.values.reshape((obs_tracks.shape[0]*obs_tracks.shape[1],obs_tracks.shape[2]))
obs_stats[np.isfinite(obs_stats[:,0]),:]

obs_stats=da.DimArray(obs_stats,axes=[range(obs_stats.shape[0]),obs_tracks.variable],dims=['id','variable'])

from scipy import stats



plt.close('all')
tc_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}
fig,axes=plt.subplots(nrows=2,ncols=2)
axes=axes.flatten()
pdfs={}
for ax,vari in zip(axes,['MSLP','VO','Wind10','T500']):
    pdfs[vari]={}
    tmp=obs_stats[obs_stats[:,'cat']>0,vari]
    tmp=tmp[np.isfinite(tmp)]
    kde = stats.gaussian_kde(tmp)
    start,stop=np.percentile(tmp,[10,90])
    x=np.linspace(start,stop,1000)
    y=kde.evaluate(x)/np.sum(kde.evaluate(x))
    pdfs[vari]['all']={'pd':y,'x':x,'dx':np.diff(x,1)[0]}
    #ax.plot(x,y,color='green')
    for cat in range(1,6):
        tmp=obs_stats[obs_stats[:,'cat']==cat,vari]
        tmp=tmp[np.isfinite(tmp)]
        kde = stats.gaussian_kde(tmp,bw_method=1)#'silverman'
        start,stop=np.percentile(tmp,[0,100])
        x=np.linspace(start,stop,10)
        y=kde.evaluate(x)#/np.sum(kde.evaluate(x))
        pdfs[vari][cat]={'pd':y,'x':x,'dx':np.diff(x,1)[0]}
        hist=np.histogram(tmp,bins=10,density=True)
        ax.plot(hist[1][:-1],hist[0],color=tc_colors[cat])
        ax.plot(x,y,color=tc_colors[cat])
    ax.set_xlabel(vari)
    ax.set_facecolor('xkcd:black')

plt.tight_layout()
plt.savefig('detection/JRA55/JRA55_obs_hists__.png')



cat4_tracks=np.where(np.nanmax(obs_tracks[:,:,'cat'],axis=1)==4)[0]
track=obs_tracks.ix[cat4_tracks[0],:,:]
track=track[np.isfinite(track[:,'VO']),:]

fcat=np.array([0.2]*6)

plot_tmp=[]
for t in track.time:
    vo_=np.argmin(abs(pdfs['VO'][cat]['pd']-track[t,'VO']))
    mslp_=np.argmin(abs(pdfs['MSLP'][cat]['pd']-track[t,'MSLP']))
    plot_tmp.append(axes[0].axvline(x=track[t,'MSLP'],color='cyan'))
    plot_tmp.append(axes[1].axvline(x=track[t,'VO'],color='cyan'))
    for cat in range(1,6):
        p_vo_cat = pdfs['VO'][cat]['pd'][vo_] * pdfs['MSLP'][cat]['pd'][mslp_]
        normalize = sum([pdfs['VO'][ccc]['pd'][np.argmin(abs(pdfs['VO'][ccc]['pd']-track[t,'VO']))] * pdfs['MSLP'][ccc]['pd'][np.argmin(abs(pdfs['MSLP'][ccc]['pd']-track[t,'MSLP']))] for ccc in range(1,6)])
        #print(cat,p_vo_cat,normalize,p_vo_cat / normalize,fcat[cat])
        fcat[cat] = p_vo_cat / normalize
    print(fcat)
    print(np.argmax(fcat))
    plt.savefig('detection/JRA55/JRA55_obs_hists_'+str(t)+'.png')
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
