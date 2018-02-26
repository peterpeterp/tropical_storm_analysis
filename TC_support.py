import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

#------------ Met stuff
def rel_vort(U,V,y,x):
    """adapted from wrftools: github.com/keltonhalbert/wrftools"""
    xx,yy = np.meshgrid(x,y)
    xx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
    yy*=6371000*2*np.pi/360.
    dy=np.gradient(yy)
    dx=np.gradient(xx)

    vort=U.copy()*np.nan
    for i in range(U.shape[0]):
        du = np.gradient( U[i,:,:] )
        dv = np.gradient( V[i,:,:] )
        vort[i,:,:]= dv[-1]/dx[-1] - du[-2]/dy[-2]
    return vort

def wind_speed(U,V):
    return np.sqrt(U**2+V**2)

# ----------- Time conversions

def toYearFraction(date):
    date=datetime(year=date.year,month=date.month,day=date.day,hour=date.hour)
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def yearFractiontoDate(yr_fr):
    year=int(yr_fr)
    fr=yr_fr-year
    startOfEpoch = datetime(year=1970, month=1, day=1)
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    seconds=fr*yearDuration+(startOfThisYear - startOfEpoch).total_seconds()
    dt=time.localtime(seconds)
    return datetime(year=dt.tm_year, month=dt.tm_mon, day=dt.tm_mday, hour=dt.tm_hour)

def running_mean_func(xx,N):
	if N==1:
		return xx
	if N!=1:
	    x=np.ma.masked_invalid(xx.copy())
	    ru_mean=x.copy()*np.nan
	    for t in range(int(N/2),len(x)-int(N/2)):
	        ru_mean[t]=np.nanmean(x[t-int(N/2):t+int(N/2)])
	    return ru_mean

#------------ TC stuff

def normalize180(lon):
    """Normalize lon to range [180, 180)"""
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

def tc_wind_color(z):
    if z<=64: color= 'lightblue'
    if z>64: color= '#ffffcc'
    if z>82: color= '#ffe775'
    if z>95: color= '#ffc148'
    if z>112: color= '#ff8f20'
    if z>136: color= '#ff6060'
    if np.isnan(z): color= 'lightblue'
    return color

def tc_pressure_color(z):
    if z>=1020000: color= 'lightblue'
    if z<1020000: color= '#ffffcc'
    if z<98000: color= '#ffe775'
    if z<96500: color= '#ffc148'
    if z<94500: color= '#ff8f20'
    if z<92000: color= '#ff6060'
    if np.isnan(z): color= 'lightblue'
    return color

def tc_plot(m, x, y, z, color_type=tc_wind_color,**kwargs):
    tmp=[]
    for i in range(len(x)-1):
        if np.isfinite(x[i+1]):
            tmp.append(m.plot(x[i:i+2],y[i:i+2],color=color_type(z[i]),**kwargs))
    return tmp

def latlon_to_meshgrid(lat,lon):
    x,y=lon.copy(),lat.copy()
    #x-=np.diff(x,1)[0]/2.
    #y-=np.diff(y,1)[0]/2.
    x=np.append(x,[x[-1]+np.diff(x,1)[0]])
    y=np.append(y,[y[-1]+np.diff(y,1)[0]])
    return np.meshgrid(x,y)

def local_extrem_filters(data,method,threshold=13,neighborhood_size=30):
    data_max = filters.maximum_filter(data, neighborhood_size)
    data_min = filters.minimum_filter(data, neighborhood_size)
    if method=='max':
        extreme = (data == data_max)
        diff = ((data_max - data_min) > threshold)
    if method=='min':
        extreme = (data == data_min)
        diff = ((data_max - data_min) < threshold)
    extreme[diff == 0] = 0

    labeled, num_objects = ndimage.label(extreme)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(int(x_center))
        y_center = (dy.start + dy.stop - 1)/2
        y.append(int(y_center))

    return y,x


def local_extrem(data,method='max',threshold=0,neighborhood_size=3):
    ex_x,ex_y=[],[]
    if method=='max':
        y,x=np.where(data>threshold)
        for y_,x_ in zip(y,x):
            if data[y_,x_]==np.max(data[max([0,y_-neighborhood_size]):min([y_+neighborhood_size,len(lat)-1]),max([0,x_-neighborhood_size]):min([x_+neighborhood_size,len(lon)-1])]):
                ex_y.append(y_)
                ex_x.append(x_)

    if method=='min':
        y,x=np.where(data<threshold)
        for y_,x_ in zip(y,x):
            if data[y_,x_]==np.min(data[max([0,y_-neighborhood_size]):min([y_+neighborhood_size,len(lat)-1]),max([0,x_-neighborhood_size]):min([x_+neighborhood_size,len(lon)-1])]):
                ex_y.append(y_)
                ex_x.append(x_)

    return ex_y,ex_x

# -------------- other

# unique / set for a list of lists
def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

def unique(l):
    return list(uniq(sorted(l, reverse=True)))
