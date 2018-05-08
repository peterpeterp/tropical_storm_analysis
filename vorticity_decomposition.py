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

if True:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--u_wind','-u',help='u wind',required=False)
    parser.add_argument('--v_wind','-v',help='v wind',required=False)
    parser.add_argument('--u_v_wind','-u_v',help='u v wind',required=False)
    parser.add_argument('--v_var',help='v variable name',required=True)
    parser.add_argument('--u_var',help='u variable name',required=True)

    parser.add_argument('--out_file','-o',help='',default='test.nc')

    parser.add_argument('--lev_var',help='',default='lev')
    parser.add_argument('--time_var',help='',default='time')
    parser.add_argument('--lat_var',help='',default='lat')
    parser.add_argument('--lon_var',help='',default='lon')
    args = parser.parse_args()
    print(args)
else:
    os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
    u_wind='easterly_waves/u_2008_sample.nc'
    v_wind='easterly_waves/v_2008_sample.nc'
    time_steps=12
    level=1
    x1,x2,y1,y2=-100,50,30,0
    out_file='test'

# Read u and v wind components from file.
if args.u_v_wind is not None:
    u = da.read_nc(args.u_v_wind)[args.u_var].ix[:,:,:,:]
    v = da.read_nc(args.u_v_wind)[args.v_var].ix[:,:,:,:]
    nc=da.read_nc(args.u_v_wind)

else:
    u = da.read_nc(args.u_wind)[args.u_var].ix[:,:,:,:]
    v = da.read_nc(args.v_wind)[args.v_var].ix[:,:,:,:]
    nc=da.read_nc(args.u_wind)

lon,lat=nc[args.lon_var].values,nc[args.lat_var].values

xx,yy = np.meshgrid(lon,lat)
dx,dy = np.meshgrid(lon.copy()*0+np.mean(np.diff(lon,1)),lat.copy()*0+np.mean(np.diff(lat,1)))
dx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
dy*=6371000*2*np.pi/360.

W = (u**2+v**2)**0.5

du_dx = (u-np.roll(u,1,axis=-1))/dx
du_dy = (u-np.roll(u,1,axis=-2))/dy

dv_dx = (v-np.roll(v,1,axis=-1))/dx
dv_dy = (v-np.roll(v,1,axis=-2))/dy

dW_dx = (W-np.roll(W,1,axis=-1))/dx
dW_dy = (W-np.roll(W,1,axis=-2))/dy

vo=dv_dx-du_dy

shear_vort=v/W*dW_dx-u/W*dW_dy
curv_vort=vo-shear_vort

curv_vort.ix[:,:,-1,:]=np.nan
curv_vort.ix[:,:,0,:]=np.nan
curv_vort.ix[:,:,:,-1]=np.nan
curv_vort.ix[:,:,:,0]=np.nan


curv_vort=da.DimArray(curv_vort,axes=[nc[args.time_var].values,nc[args.lev_var].values,lat,lon],dims=[args.time_var,args.lev_var,args.lat_var,args.lon_var])
time_=nc[args.time_var]
level=nc[args.lev_var]
lon=nc[args.lon_var]
lat=nc[args.lat_var]
# time_.units=nc['time'].units
# time_.calendar=nc['time'].calendar
ds=da.Dataset({args.time_var:time_,args.lon_var:lon,args.lat_var:lat,args.lev_var:level,'curv_vort':curv_vort})
ds.write_nc(args.out_file)
