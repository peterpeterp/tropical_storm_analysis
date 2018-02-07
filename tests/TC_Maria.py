import os,sys,glob,time,collections,gc
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import pandas as pd
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import webcolors
import matplotlib
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')

def projected_polygon(map,points):
    verts=[]
    for lon,lat in points:
        verts.append(m(lon,lat))

    codes2 = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    path = Path(verts, codes2)
    return path


plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
m = Basemap(projection='ortho',lon_0=-50,lat_0=20,resolution='l',ax=ax)
m.bluemarble(scale=0.2)

# Maria on 800 AM AST Tue Sep 19 2017 http://www.nhc.noaa.gov/archive/2017/al15/al152017.public_a.013.shtml?
m.tissot(-62.8, 16.2, 205/40000.*360, 50, color='white',alpha=0.4)
m.tissot(-62.8, 16.2, 95/40000.*360, 50, color='white',alpha=0.4)

for lon in np.arange(-80,-40,2):
    for lat in np.arange(0,30,2):
        path = projected_polygon(m,[(lon,lat),(lon+2,lat),(lon+2,lat+2),(lon,lat+2),(lon,lat)])
        patch = patches.PathPatch(path, facecolor='r',lw=1.5,edgecolor='r',alpha=0.2)
        ax.add_patch(patch)

plt.tight_layout()
plt.savefig('plots/GCM_Maria_.png')

plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
m = Basemap(lat_0=16.2,lon_0=-62.8,width=1400000,height=1400000,projection='lcc',lat_1=10.,lat_2=20.,resolution ='l',area_thresh=1000.,ax=ax)
m.bluemarble(scale=0.2)

# Maria on 800 AM AST Tue Sep 19 2017 http://www.nhc.noaa.gov/archive/2017/al15/al152017.public_a.013.shtml?
m.tissot(-62.8, 16.2, 205/40000.*360, 50, color='white',alpha=0.4)
m.tissot(-62.8, 16.2, 95/40000.*360, 50, color='white',alpha=0.4)

for lon in np.arange(-80,-40,2):
    for lat in np.arange(0,30,2):
        path = projected_polygon(m,[(lon,lat),(lon+2,lat),(lon+2,lat+2),(lon,lat+2),(lon,lat)])
        patch = patches.PathPatch(path, facecolor='r',lw=1.5,edgecolor='r',alpha=0.2)
        ax.add_patch(patch)

plt.tight_layout()
plt.savefig('plots/GCM_Maria.png')

plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
m = Basemap(lat_0=16.2,lon_0=-62.8,width=1400000,height=1400000,projection='lcc',lat_1=10.,lat_2=20.,resolution ='l',area_thresh=1000.,ax=ax)
m.bluemarble()

# Maria on 800 AM AST Tue Sep 19 2017 http://www.nhc.noaa.gov/archive/2017/al15/al152017.public_a.013.shtml?
m.tissot(-62.8, 16.2, 205/40000.*360, 50, color='white',alpha=0.4)
m.tissot(-62.8, 16.2, 95/40000.*360, 50, color='white',alpha=0.4)

for lon in np.arange(-80,-50,0.5):
    for lat in np.arange(5,25,0.5):
        path = projected_polygon(m,[(lon,lat),(lon+0.5,lat),(lon+0.5,lat+0.5),(lon,lat+0.5),(lon,lat)])
        patch = patches.PathPatch(path, facecolor='r',lw=1,edgecolor='r',alpha=0.2)
        ax.add_patch(patch)

plt.tight_layout()
plt.savefig('plots/RCM_Maria.png')

# # CHIPS coordinates
# def potential_radius(r,V,lat):
#     f=7.2921/10**5*2*np.sin(np.radians(lat))
#     R=(2*r*V/f+r**2)**0.5
#     return R
#
# def potential_radius_back(R,V,lat):
#     f=7.2921/10**5*2*np.sin(np.radians(lat))
#     r=-V/f+(v**2/f**2+R**2)**0.5
#     return r
#
# r=1000
# V=70

print potential_radius(r,V,16.2)/1000.

plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
m = Basemap(lat_0=16.2,lon_0=-62.8,width=1400000,height=1400000,projection='lcc',lat_1=10.,lat_2=20.,resolution ='l',area_thresh=1000.,ax=ax)
m.bluemarble()

# Maria on 800 AM AST Tue Sep 19 2017 http://www.nhc.noaa.gov/archive/2017/al15/al152017.public_a.013.shtml?
#m.tissot(-62.8, 16.2, 1000/40000.*360, 50, color='blue',alpha=0.4)
for i in np.arange(51,1,-1):
    rad=i**2/float(50**2)*1000-1
    print rad
    m.tissot(-62.8, 16.2, rad/40000.*360, 50, facecolor='white',alpha=0.1,edgecolor='red',lw=2)
    #m.tissot(-62.8, 16.2, (i-1)**2*10/40000.*360, 50, color='blue',alpha=1)


plt.tight_layout()
plt.savefig('plots/CHIPS_Maria.png')



plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
m = Basemap(lat_0=16.,lon_0=-62.8,width=500000,height=500000,projection='lcc',lat_1=12.,lat_2=18.,resolution ='l',area_thresh=1000.,ax=ax)
m.bluemarble()

# Maria on 800 AM AST Tue Sep 19 2017 http://www.nhc.noaa.gov/archive/2017/al15/al152017.public_a.013.shtml?
m.tissot(-62.8, 16.2, 205/40000.*360, 50, color='white',alpha=0.4)
m.tissot(-62.8, 16.2, 95/40000.*360, 50, color='white',alpha=0.4)
m.tissot(-62.8, 16.2, 10/40000.*360, 50, color='blue',alpha=1)

for lon in np.arange(-80,-50,0.5):
    for lat in np.arange(5,25,0.5):
        path = projected_polygon(m,[(lon,lat),(lon+0.5,lat),(lon+0.5,lat+0.5),(lon,lat+0.5),(lon,lat)])
        patch = patches.PathPatch(path, facecolor='r',lw=1,edgecolor='r',alpha=0.2)
        ax.add_patch(patch)

plt.tight_layout()
plt.savefig('plots/RCM_Maria_zoom.png')







plt.close()
# lon_0, lat_0 are the center point of the projection.
# resolution = 'l' means use low resolution coastlines.
m = Basemap(projection='ortho',lon_0=-50,lat_0=40,resolution='l')
# m.drawcoastlines()
# m.drawmapboundary(fill_color='darkcyan')
# m.fillcontinents(color='white',lake_color='darkcyan')
# # draw parallels and meridians.
# m.drawparallels(np.arange(-90.,120.,30.))
# m.drawmeridians(np.arange(0.,420.,60.))

m.bluemarble(scale=0.2)   # full scale will be overkill
#m.drawcoastlines(color='white', linewidth=0.3)  # add coastlines


for id_ in TC.storm:
    if TC['basin'][id_,0]==0.0 and TC['season'][id_]>2014:
        x,y=TC['lon_for_mapping'][id_,:].values,TC['lat_for_mapping'][id_,:].values
        x1,y1=m(x[np.isfinite(x)],y[np.isfinite(y)])
        msw=TC['source_wind'][id_,:,0].values[np.isfinite(x)]
        color=np.array(['#ffffcc']*msw.shape[0],dtype='|S7')
        color[msw>82]='#ffe775'
        color[msw>95]='#ffc148'
        color[msw>112]='#ff8f20'
        color[msw>136   ]='#ff6060'
        color=[tuple(np.array(webcolors.hex_to_rgb(cc))/255.0) for cc in color]
        m.scatter(x1,y1,marker='o',color=color)

plt.title("Tropical Cyclones")
plt.savefig('plots/test.png')
