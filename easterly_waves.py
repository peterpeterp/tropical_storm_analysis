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

from scipy import fftpack

from windspharm.standard import VectorWind
from skimage.feature import peak_local_max

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones')


from scipy.signal import butter, lfilter


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5,axis=0):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data,axis=axis)
    return y


# Read u and v wind components from file.
u = da.read_nc('easterly_waves/atl_U_2017_bp.nc')['var33'].ix[:,0,:,:].squeeze()
#u=u.ix[:,:,100:160,range(100)+range(500,640)]

v = da.read_nc('easterly_waves/atl_V_2017_bp.nc')['var34'].ix[:,0,:,:].squeeze()
#v=v.ix[:,:,100:160,range(100)+range(500,640)]

# fs=1/(6.*60*60)
# lowcut=1/(144.*60*60)
# highcut=1/(48.*60*60)
#
# u.values=butter_bandpass_filter(u.values,lowcut,highcut,fs)
# v.values=butter_bandpass_filter(v.values,lowcut,highcut,fs)

lon,lat=u.lon,u.lat
lon[lon>180]-=360
lon_rolling=len(lon)-np.where(lon<0)[0][0]
lon=np.roll(lon,lon_rolling,axis=-1)
u.values=np.roll(u,lon_rolling,axis=-1)
v.values=np.roll(v,lon_rolling,axis=-1)

lons,lats=np.meshgrid(lon,lat)

VO=TC_support.rel_vort(u.values[:,:,:],v.values[:,:,:],lat,lon)

u[np.isnan(u)]=0
v[np.isnan(v)]=0


plt.close('all')
plate_carree = ccrs.PlateCarree()
fig,axes = plt.subplots(nrows=6,ncols=1,figsize=(8,6),subplot_kw={'projection': plate_carree})
for ax in axes.flatten():
    #ax.axis('off')
    ax.set_global()
    ax.coastlines(edgecolor='magenta')
    ax.add_feature(cartopy.feature.LAND, facecolor='darkgreen')
    ax.add_feature(cartopy.feature.OCEAN,facecolor='darkblue')
    ax.set_extent([np.min(lon),np.max(lon),np.min(lat),np.max(lat)],crs=plate_carree)

for ax,t in zip(axes.flatten(),range(100,106)):
    w = VectorWind(u.values[t,:,:], v.values[t,:,:])
    vo=w.vorticity(truncation=None)
    vo=VO[t,:,:]
    vo[vo<1*10**(-5)]=0
    ax.pcolormesh(lons,lats,vo)
    coords=peak_local_max(vo, min_distance=6)
    ax.plot(lon[coords[:,1]],lat[coords[:,0]],'*m',linestyle='')

plt.tight_layout()
plt.savefig('test.png')



# def fft_filter(data,keep_fraction=0.1):
#     im_fft = fftpack.fft2(data)
#
#     # Call ff a copy of the original transform. Numpy arrays have a copy
#     # method for this purpose.
#     im_fft2 = im_fft.copy()
#
#     # Set r and c to be the number of rows and columns of the array.
#     r, c = im_fft2.shape
#
#     # Set to zero all rows with indices between r*keep_fraction and
#     # r*(1-keep_fraction):
#     im_fft2[int(r*keep_fraction):int(r*(1-keep_fraction))] = 0
#
#     # Similarly with the columns:
#     im_fft2[:, int(c*keep_fraction):int(c*(1-keep_fraction))] = 0
#     return fftpack.ifft2(im_fft2).real
#
# def fft_filter_truncate(data,k_low=3,k_high=63):
#     im_fft = fftpack.fft2(data)
#     im_fft2 = im_fft.copy()
#
#     r, c = im_fft2.shape
#
#     freq=(im_fft2.real**2+im_fft2.imag**2)**0.5/0.56*360
#     k=2*np.pi/freq
#
#     im_fft2[k<k_low] = 0
#     im_fft2[k>k_high] = 0
#
#     return fftpack.ifft2(im_fft2).real
