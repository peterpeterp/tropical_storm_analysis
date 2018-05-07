from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy


os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
u_wind='easterly_waves/u_2008_sample.nc'
v_wind='easterly_waves/v_2008_sample.nc'
time_steps=12
level=1
x1,x2,y1,y2=-100,0,40,0
out_file='test'

# Read u and v wind components from file.
u = da.read_nc(u_wind)['var33'].ix[time_steps,level,:,:].squeeze()
v = da.read_nc(v_wind)['var34'].ix[time_steps,level,:,:].squeeze()

u.lon[u.lon>180]-=360
v.lon[v.lon>180]-=360
lon_rolling=len(u.lon)-np.where(u.lon<0)[0][0]
u.lon=np.roll(u.lon,lon_rolling,axis=-1)
v.lon=np.roll(v.lon,lon_rolling,axis=-1)
u.values=np.roll(u,lon_rolling,axis=-1)
v.values=np.roll(v,lon_rolling,axis=-1)

u=u[y1:y2,x1:x2]
v=v[y1:y2,x1:x2]
lons,lats=np.meshgrid(u.lon,u.lat)


#VO=TC_support.rel_vort(u.values[:,:],v.values[:,:],u.lat,u.lon)

u[np.isnan(u)]=0
v[np.isnan(v)]=0

xx,yy = np.meshgrid(u.lon,u.lat)
dx,dy = np.meshgrid(u.lon.copy()*0+np.mean(np.diff(u.lon,1)),u.lat.copy()*0+np.mean(np.diff(u.lat,1)))
dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
dy*=6371000*2*np.pi/360.

W = (u**2+v**2)**0.5

du_dx = (u-np.roll(u,1,axis=1))/dx
du_dy = (u-np.roll(u,1,axis=0))/dy

dv_dx = (v-np.roll(v,1,axis=1))/dx
dv_dy = (v-np.roll(v,1,axis=0))/dy

dW_dx = (W-np.roll(W,1,axis=1))/dx
dW_dy = (W-np.roll(W,1,axis=0))/dy

vo=dv_dx-du_dy

shear_vort=v/W*dW_dx-u/W*dW_dy
curv_vort=vo-shear_vort

dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=1))/dx
dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=0))/dy
curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

curv_vort_coarse=ndimage.filters.gaussian_filter(curv_vort,3)
curv_vort_advect_coarse=ndimage.filters.gaussian_filter(curv_vort_advect,3)

kernel_size=5
kernel = np.ones([kernel_size,kernel_size])/float(kernel_size*kernel_size)

curv_vort_smoo=ndimage.convolve(curv_vort_coarse,kernel)
curv_vort_smoo=ndimage.convolve(curv_vort_smoo,kernel)

dcurv_vort_dx = (curv_vort_smoo-np.roll(curv_vort_smoo,1,axis=1))/dx
dcurv_vort_dy = (curv_vort_smoo-np.roll(curv_vort_smoo,1,axis=0))/dy
curv_vort_advect_smoo=-(u*dcurv_vort_dx+v*dcurv_vort_dy)

from itertools import combinations

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
        dist=np.sqrt((u.lat[pair[0][0]]-u.lat[pair[1][0]])**2+(u.lon[pair[0][1]]-u.lon[pair[1][1]])**2)
        if dist > max_square_distance:
            max_square_distance = dist
            max_pair = pair
    return(max_square_distance,max_pair)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(10,5),subplot_kw={'projection': plate_carree})
for ax in axes.flatten():
    #ax.axis('off')
    ax.set_global()
    ax.coastlines(edgecolor='magenta')
    ax.set_extent([np.min(u.lon),np.max(u.lon),np.min(u.lat),np.max(u.lat)],crs=plate_carree)


axes[0,0].contour(lons,lats,vo,np.arange(0,0.0001,10**(-5)))
axes[1,0].contourf(lons,lats,curv_vort_smoo,np.arange(-3,3,0.5)*10**(-5),cmap=plt.cm.RdYlBu_r)
axes[1,0].contour(u.lon,u.lat,curv_vort_advect_smoo,[-1,0,1],color='green')
axes[2,0].contourf(lons,lats,curv_vort_smoo,np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd)

coords=peak_local_max(curv_vort_smoo, min_distance=10)
for y,x in coords:
    if 5<u.lat[y]<35:
        if curv_vort_smoo[y,x]>1*10**(-5):
            if u.ix[y,x]<2.5:
                print(x,y,curv_vort_advect_smoo.ix[y,x],curv_vort_advect_smoo.ix[y,x-1],curv_vort_advect_smoo.ix[y,x+1])
                if np.sign(curv_vort_advect_smoo.ix[y,x])!=np.sign(curv_vort_advect_smoo.ix[y,x+1]):
                    axes[2,0].plot(u.lon[x],u.lat[y],'*b')

axes[0,1].contourf(lons,lats,curv_vort_smoo,np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd)

asign = np.sign(curv_vort_advect_smoo.values)
signchange_y = ((np.roll(asign, 1,axis=0) - asign) != 0).astype(int)
signchange_x = ((np.roll(asign, 1,axis=1) - asign) != 0).astype(int)
y,x=np.where((signchange_y+signchange_x>0))
axes[0,1].plot(u.lon[x],u.lat[y],'.g',linestyle='')

y,x=np.where((signchange_y+signchange_x>0) & (u.values<2.5) & (curv_vort_smoo>5*10**(-6)))

points=np.vstack((y,x)).T.tolist()
groups=group_points(points)

while True:
    for group in groups:
        done=True
        if group_extend(group)[0]>10:
            done=False
            groups.remove(group)
            thresh=5*10**(-6)
            while True:
                thresh*=1.5
                for pp in group:
                    if curv_vort_smoo[pp[0],pp[1]]<thresh:
                        group.remove(pp)
                sub_groups=group_points(group)
                stop=True
                for sub in sub_groups:
                    if group_extend(sub)[0]>10:
                        group=sub
                        stop=False
                    else:
                        groups.append(sub)
                if stop:
                    break
    if done:
        break


for group in groups:
    y,x=[p[0] for p in group],[p[1] for p in group]
    axes[1,1].plot(u.lon[x],u.lat[y],'.',linestyle='')

axes[2,1].contourf(lons,lats,curv_vort_smoo,np.arange(0.5,3,0.5)*10**(-5),cmap=plt.cm.YlOrRd)

for group in groups:
    if len(group)>1:
        dist,pair=group_extend(group)
        wave_y=[pair[0][0],int(round(np.array(group)[:,0].mean())),pair[1][0]]
        wave_x=[pair[0][1],int(round(np.array(group)[:,1].mean())),pair[1][1]]
        axes[2,1].plot(u.lon[wave_x],u.lat[wave_y],'.',linestyle='-')


plt.tight_layout()
plt.savefig(out_file)






# points=np.vstack((y,x)).tolist()
# used_pos=[]
# i=0
# while len(used_pos)<len(points):
#     p=points[i]
#     i+=1
#     if p not in used_pos:
#         used_pos.append(p)
#         group=[p]
#         while True:
#             found=False
#             yy,xx=group[-1][0],group[-1][1]
#             candidates=[(yy-1,xx),(yy+1,xx),(yy,xx-1),(yy,xx+1)]
#             for pp in candidates:
#                 if pp in points and pp not in used_pos:
#                     found=True
#                     group.append(pp)
#                     used_pos.append(pp)







#
