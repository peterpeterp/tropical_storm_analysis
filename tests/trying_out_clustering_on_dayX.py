# import os,sys,glob,time,collections,gc,calendar
# from subprocess import Popen
#
#
# mslp=da.read_nc('data/ERA_6hourly/atl_MSL_'+str(year)+'.nc')['MSL']
# vort=da.read_nc('data/ERA_6hourly/atl_atmos_'+str(year)+'.nc')['VO']
#
# plt.close('all')
# fig,axes = plt.subplots(nrows=4,ncols=1,figsize=(5,10))
# axes.flatten()
# for i in range(4):
#     ax=axes[i]
#     data=vort.ix[93*4+i,0,:,:]
#
#     im=ax.imshow(data,vmin=-9.5*10**(-5),vmax=0.0002,interpolation='none')
#     im.set_cmap('bone')
#     ax.autoscale(False)
#     ax.axis('off')
# plt.tight_layout()
# plt.savefig('plots/test_6hour.png')
#

from sklearn.cluster import KMeans

def max_clusters__(data,threshold,neighborhood_size):
    y,x=np.where(data>threshold)
    points=[]
    y__,x__=[],[]
    groups=0
    for y_,x_ in zip(y,x):
        box=get_box(y_,x_,neighborhood_size)
        if len(np.where(data[box[0]:box[1],box[2]:box[3]]>threshold)[0])>5:
            points.append([y_,x_])
            y__.append(y_)
            x__.append(x_)
            if data[y_,x_]==np.max(data[box[0]:box[1],box[2]:box[3]]):
                groups+=1

    return KMeans(n_clusters=groups, random_state=0).fit(np.array(points)).labels_,np.array(y__),np.array(x__)


def get_box(y,x,window):
    y_min=max(0,y-window)
    y_max=min(len(vort.lat),y+window+1)
    x_min=max(0,x-window)
    x_max=min(len(vort.lon),x+window+1)
    return (y_min,y_max,x_min,x_max)

def get_corners(y,x,window):
    y_min,y_max,x_min,x_max=get_box(y,x,window)
    return ((y_min,x_min),(y_min,x_max),(y_max,x_max),(y_max,x_min))


def get_surrounding(y,x,window):
    y_min,y_max,x_min,x_max=get_box(y,x,window)
    points=[]
    for y in range(y_min,y_max):
        for x in range(x_min,x_max):
            points.append((y,x))
    return points

def get_neihgbors(y,x):
    return [y+nei[0] for nei in ((-1,0),(+1,0),(0,-1),(0,+1))],[x+nei[1] for nei in ((-1,0),(+1,0),(0,-1),(0,+1))]

import numpy as np
def path_length(points):
    return np.sum(np.sqrt(np.sum(np.diff(np.array(points), axis=0)**2, axis=1)))
path_length(points=[(-4,-4,-4),(4,4,1),(1,1,1)])

def max_clusters(data,threshold,neighborhood_size):
    y,x=np.where(data>threshold)
    points=[(y_,x_) for y_,x_ in zip(y,x)]
    group=0
    clusters={}
    for y_,x_ in zip(y,x):
        box_1=get_box(y_,x_,5)
        box_2=get_box(y_,x_,3)
        if data[y_,x_]==np.max(data[box_1[0]:box_1[1],box_1[2]:box_1[3]]):
            print '****************'
            candidates=points[:]
            found=[(y_,x_)]
            candidates.remove((y_,x_))
            for p in found:
                print p,found
                for nei in get_surrounding(p[0],p[1],1):
                    if nei in candidates:
                        found.append(nei)
                        candidates.remove(nei)

            clusters[group]={'y':[pp[0] for pp in found],'x':[pp[1] for pp in found]}
            group+=1

    return clusters


from scipy import stats

# mslp=da.read_nc('data/ERA_dayX/atl_MSL_'+str(year)+'.nc')['MSL']
# vort=da.read_nc('data/ERA_dayX/atl_atmos_'+str(year)+'.nc')['VO']
wind=da.read_nc('data/ERA_dayX/atl_10mWind_'+str(year)+'.nc')['ws']

plt.close('all')
fig,axes = plt.subplots(nrows=3,ncols=1,figsize=(5,7.5))

ax=axes[0]
data=vort.ix[95,0,:,:].values
im=ax.imshow(data,vmin=-9.5*10**(-5),vmax=0.0002,interpolation='none')
clustered=max_clusters(data,1*10**(-4),5)
y,x=np.where(data>1*10**(-4))
ax.plot(x,y,'.g',markersize=3)
for group in clustered.values():
    if len(group['x'])>4:
        ax.plot(group['x'],group['y'],'.',markersize=3)
        lr=stats.linregress(group['x'],group['y'])
        xxx=np.array(list(set(group['x'])))
        ax.plot(xxx,xxx*lr.slope+lr.intercept)
        for step in range(4):
            xx=min(xxx)+len(xxx)/4.*step
            ax.plot([xx],[xx*lr.slope+lr.intercept],'*r',markersize=5)
im.set_cmap('bone');ax.autoscale(False);ax.axis('off')



ax=axes[1]
data=mslp.ix[93,:,:].values
im=ax.imshow(data,vmin=100360,vmax=103000,interpolation='none')
# for i in range(len(set(kmeans))):
#     lr=stats.linregress(x_[kmeans==i], y_[kmeans==i])
#     xxx=np.array(list(set(x_[kmeans==i])))
#     ax.plot(xxx,xxx*lr.slope+lr.intercept)
#     for step in range(4):
#         xx=xxx[0]+len(xxx)/4.*step
#         ax.plot([xx],[xx*lr.slope+lr.intercept],'*r')
im.set_cmap('bone')
ax.autoscale(False)
ax.axis('off')


ax=axes[2]
data=wind.ix[93,:,:].values
im=ax.imshow(data,vmin=0,vmax=13,interpolation='none')
# for i in range(len(set(kmeans))):
#     lr=stats.linregress(x_[kmeans==i], y_[kmeans==i])
#     xxx=np.array(list(set(x_[kmeans==i])))
#     ax.plot(xxx,xxx*lr.slope+lr.intercept)
#     for step in range(4):
#         xx=xxx[0]+len(xxx)/4.*step
#         ax.plot([xx],[xx*lr.slope+lr.intercept],'*b')
im.set_cmap('bone')
ax.autoscale(False)
ax.axis('off')

plt.savefig('plots/test.png')
plt.savefig('plots/test.pdf')
