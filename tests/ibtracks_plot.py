import os,sys,glob,time,collections,gc,calendar,weakref
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import seaborn as sn
import webcolors

from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry.polygon import Polygon

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

#TC=da.read_nc('data/Allstorms.ibtracs_all.v03r10.nc')
# tc_sel=TC.ix[np.where(TC['basin'][:,0]==0)[0]]
# tc_lat=tc_sel['lat_for_mapping']
# tc_lon=tc_sel['lon_for_mapping']
# tc_lon[tc_lon<0]+=360
# tc_wind=tc_sel['source_wind']
# tc_wind[np.isnan(tc_wind)]=-999
#
plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(10,5))
axes=axes.flatten()
maps=[]
for ax in axes:
    m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
    m.drawmapboundary(fill_color='1.')
    m.drawmapboundary(fill_color='darkblue')
    m.fillcontinents(color='black',lake_color='darkblue')
    m.drawcoastlines(linewidth=0.3)
    m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
    m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
    maps.append(m)


tmp=[{},{},{},{}]
polys=[]
for storm in tc_sel.storm:
    tmp[0][storm]=maps[0].plot(tc_lon[storm,:],tc_lat[storm,:],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=0.05)


caribbean=Polygon([(300,19),(282,28),(275,22),(280,13),(300,9),(300,19)])
for storm in tc_sel.storm:
    inter=0
    for i in range(tc_lat.shape[1]):
        if caribbean.contains(Point(tc_lon[storm,i],tc_lat[storm,i])):
            inter+=1
            tmp[1][storm]=maps[1].plot(tc_lon[storm,:],tc_lat[storm,:],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=0.05)
            break
    if inter==0:
        tmp[2][storm]=maps[2].plot(tc_lon[storm,:],tc_lat[storm,:],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=0.05)


xs,ys=caribbean.exterior.xy
polys.append(axes[1].fill(xs, ys, alpha=0.5, fc='r', ec='k'))

axes[0].set_title('All TCs: '+str(len(tmp[0])))
axes[1].set_title('TCs passing through the Caribbean: '+str(len(tmp[1])))
axes[2].set_title('TCs not passing through the Caribbean: '+str(len(tmp[2])))

for i in range(3):
    for storm in tmp[i].keys():
        tmp[i][storm][0].set_alpha(0.2)
        tmp[i][storm][0].set_linewidth(1)

plt.savefig('plots/ibtracks_atlantic_dark.png',dpi=300)

for i in range(3):
    for storm in tmp[i].keys():
        tmp[i][storm][0].set_alpha(1)
        tmp[i][storm][0].set_linewidth(0.001)

plt.savefig('plots/ibtracks_atlantic_dark.pdf')


for element in polys:
    l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
    polys.remove(element)


angle=83
alpha=np.radians(angle)
beta=np.radians(90-angle)
domain=Polygon([(320,5),
                (320+np.cos(alpha)*17,5+np.sin(alpha)*17),
                (320+np.cos(alpha)*17,5+np.sin(alpha)*17),
                (320+np.cos(alpha)*17-np.cos(beta)*50,5+np.sin(alpha)*17+np.sin(beta)*50),
                (320-np.cos(beta)*50,5+np.sin(beta)*50),
                (320,5)])
xs,ys=domain.exterior.xy
polys.append(axes[3].fill(xs, ys, alpha=0.5, fc='c', ec='w'))

for i in range(1,4):
    xs,ys=caribbean.exterior.xy
    polys.append(axes[i].fill(xs, ys, alpha=0.5, fc='r', ec='w'))

plt.savefig('plots/ibtracks_atlantic_dark_domain.pdf')


for id_ in tmp[3]:
    element=tmp[3][id_]
    if len(element)==1:
        l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
tmp[3]={}

for storm in tc_sel.storm:
    inter=0
    for i in range(tc_lat.shape[1]):
        if caribbean.contains(Point(tc_lon[storm,i],tc_lat[storm,i])):
            inter+=1
            break
    if inter!=0:
        for j in range(tc_lat.shape[1])[0:i]:
            if domain.contains(Point(tc_lon[storm,j],tc_lat[storm,j]))==False:
                if np.isfinite(tc_lon[storm,j]+tc_lat[storm,j]):
                    print tc_lon[storm,j],tc_lat[storm,j]
                    tmp[3][storm]=maps[3].plot(tc_lon[storm,0:i],tc_lat[storm,0:i],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=1,linewidth=0.001)
                    break

axes[3].set_title('TCs passing through the Caribbean\n starting outside the domain: '+str(len(tmp[3])))

plt.savefig('plots/ibtracks_atlantic_domain.pdf')


# for storm in tc_sel.storm:
#     tmp[storm][0].set_color(get_tc_color(tc_wind[storm,:,0].max()))
