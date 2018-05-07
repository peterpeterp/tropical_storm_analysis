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

if False:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--u_wind','-u',help='u wind',required=True)
    parser.add_argument('--v_wind','-v',help='v wind',required=True)

    parser.add_argument('--out_file','-o',help='',default='test.png')

    parser.add_argument('--x1',help='',type=int,required=False)
    parser.add_argument('--x2',help='',type=int,required=False)
    parser.add_argument('--y1',help='',type=int,required=False)
    parser.add_argument('--y2',help='',type=int,required=False)
    parser.add_argument('--level',help='',type=int,required=False)
    parser.add_argument('--time_steps','-ts',help='time steps',nargs='+',type=int,required=False)
    args = parser.parse_args()
    print(args)
    x1,x2,y1,y2=args.x1,args.x2,args.y1,args.y2
    u_wind,v_wind=args.u_wind,args.v_wind
    time_steps,level=args.time_steps,args.level
    out_file=args.out_file

else:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    u_wind='easterly_waves/u_2008_sample.nc'
    v_wind='easterly_waves/v_2008_sample.nc'
    time_steps=12
    level=1
    x1,x2,y1,y2=-100,50,30,0
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

kernel = np.ones([9,9])
vo_smoo=ndimage.convolve(vo,kernel)
vo_smoo=ndimage.convolve(vo_smoo,kernel)

shear_vort_smoo=ndimage.convolve(shear_vort,kernel)
shear_vort_smoo=ndimage.convolve(shear_vort_smoo,kernel)

curv_vort_smoo=ndimage.convolve(curv_vort,kernel)
curv_vort_smoo=ndimage.convolve(curv_vort_smoo,kernel)


dcurv_vort_dx = (curv_vort-np.roll(curv_vort,1,axis=1))/dx
dcurv_vort_dy = (curv_vort-np.roll(curv_vort,1,axis=0))/dy

curv_vort_advect=-(u*dcurv_vort_dx+v*dcurv_vort_dy)
curv_vort_advect_smoo=ndimage.convolve(curv_vort_advect,kernel)
curv_vort_advect_smoo=ndimage.convolve(curv_vort_advect_smoo,kernel)

asp=abs(x1-x2)*3/(abs(y1-y2)*3)

plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(5*asp,5),subplot_kw={'projection': plate_carree})
for ax in axes.flatten():
    #ax.axis('off')
    ax.set_global()
    ax.coastlines(edgecolor='magenta')
    ax.set_extent([np.min(u.lon),np.max(u.lon),np.min(u.lat),np.max(u.lat)],crs=plate_carree)

# axes[0,0].pcolormesh(lons,lats,u,vmin=-30,vmax=30); axes[0,0].annotate('U', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
# axes[1,0].pcolormesh(lons,lats,v,vmin=-30,vmax=30); axes[1,0].annotate('V', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
# axes[2,0].pcolormesh(lons,lats,W,vmin=-30,vmax=30); axes[2,0].annotate('W', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')

axes[0,0].pcolormesh(lons,lats,vo,vmin=-3*10**(-5),vmax=3*10**(-5)); axes[0,1].annotate('rel. vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[1,0].pcolormesh(lons,lats,shear_vort,vmin=-3*10**(-5),vmax=3*10**(-5)); axes[1,1].annotate('shear vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[2,0].pcolormesh(lons,lats,curv_vort,vmin=-3*10**(-5),vmax=3*10**(-5)); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')

axes[0,1].pcolormesh(lons,lats,vo_smoo); axes[0,1].annotate('rel. vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
axes[1,1].pcolormesh(lons,lats,shear_vort_smoo); axes[1,1].annotate('shear vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
to_plot=curv_vort_smoo.copy()
to_plot[to_plot<0]=np.nan
axes[2,1].pcolormesh(lons,lats,to_plot); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')

to_plot=vo.copy()
to_plot[to_plot<0.5*10**(-5)]=np.nan
axes[0,2].pcolormesh(lons,lats,to_plot,vmin=-3*10**(-5),vmax=3*10**(-5),cmap=plt.cm.RdYlBu_r); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
to_plot=curv_vort.copy()
to_plot[to_plot<0.5*10**(-5)]=np.nan
axes[1,2].pcolormesh(lons,lats,to_plot,vmin=-3*10**(-5),vmax=3*10**(-5),cmap=plt.cm.RdYlBu_r); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
CS=axes[2,1].contour(u.lon,u.lat,curv_vort_advect_smoo,[-1,0,1],color='green'); axes[2,1].annotate('curvature vorticity', xy=(0.01, 0.01), xycoords='axes fraction',fontsize=10,color='w')
#axes[2,2].clabel(CS, fontsize=9, inline=1)

plt.tight_layout()
plt.savefig(out_file)
