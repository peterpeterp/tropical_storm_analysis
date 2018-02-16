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
fig,axes=plt.subplots(nrows=2,ncols=2)
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
for storm in tc_sel.storm:
    tmp[0][storm]=maps[0].plot(tc_lon[storm,:],tc_lat[storm,:],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=1,linewidth=0.001)


caribbean=Polygon([(300,19),(282,28),(275,22),(280,13),(300,9),(300,19)])
for storm in tc_sel.storm:
    inter=0
    for i in range(tc_lat.shape[1]):
        if caribbean.contains(Point(tc_lon[storm,i],tc_lat[storm,i])):
            inter+=1
            tmp[1][storm]=maps[1].plot(tc_lon[storm,:],tc_lat[storm,:],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=1,linewidth=0.001)
            break
    if inter==0:
        pass

polys=[]
xs,ys=caribbean.exterior.xy
for i in range(1,4):
    polys.append(axes[i].fill(xs, ys, alpha=0.5, fc='r', ec='k'))

axes[0].set_title('All TCs: '+str(len(tmp[0])),fontsize=8)
axes[1].set_title('TCs passing through the Caribbean: '+str(len(tmp[1])),fontsize=8)

angle,height,width=86,25,83
y,x=340,5
alpha=np.radians(angle)
beta=np.radians(90-angle)
domain=Polygon([(y,x),
                (y+np.cos(alpha)*height,x+np.sin(alpha)*height),
                (y+np.cos(alpha)*height,x+np.sin(alpha)*height),
                (y+np.cos(alpha)*height-np.cos(beta)*width,x+np.sin(alpha)*height+np.sin(beta)*width),
                (y-np.cos(beta)*width,x+np.sin(beta)*width),
                (y,x)])
xs,ys=domain.exterior.xy
domain_2=axes[2].fill(xs, ys, alpha=0.5, fc='c', ec='w')

for storm in tc_sel.storm:
    inter=0
    for i in range(tc_lat.shape[1]):
        if domain.contains(Point(tc_lon[storm,i],tc_lat[storm,i])):
            inter+=1
            break
    if inter!=0:
        for j in range(tc_lat.shape[1])[0:i]:
            if domain.contains(Point(tc_lon[storm,j],tc_lat[storm,j]))==False:
                if np.isfinite(tc_lon[storm,j]+tc_lat[storm,j]):
                    print tc_lon[storm,j],tc_lat[storm,j]
                    tmp[2][storm]=maps[2].plot(tc_lon[storm,0:i],tc_lat[storm,0:i],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=1,linewidth=0.001)
                    break

axes[2].set_title('TCs starting outside the domain: '+str(len(tmp[2])),fontsize=8)

angle,height,width=86,25,60
x,y=317,6.6
alpha=np.radians(angle)
beta=np.radians(90-angle)
domain=Polygon([(x,y),
                (x+np.cos(alpha)*height,y+np.sin(alpha)*height),
                (x+np.cos(alpha)*height,y+np.sin(alpha)*height),
                (x+np.cos(alpha)*height-np.cos(beta)*width,y+np.sin(alpha)*height+np.sin(beta)*width),
                (x-np.cos(beta)*width,y+np.sin(beta)*width),
                (x,y)])
xs,ys=domain.exterior.xy
domain_3=axes[3].fill(xs, ys, alpha=0.5, fc='c', ec='w')

for storm in tc_sel.storm:
    inter=0
    for i in range(tc_lat.shape[1]):
        if domain.contains(Point(tc_lon[storm,i],tc_lat[storm,i])):
            inter+=1
            break
    if inter!=0:
        for j in range(tc_lat.shape[1])[0:i]:
            if domain.contains(Point(tc_lon[storm,j],tc_lat[storm,j]))==False:
                if np.isfinite(tc_lon[storm,j]+tc_lat[storm,j]):
                    print tc_lon[storm,j],tc_lat[storm,j]
                    tmp[3][storm]=maps[3].plot(tc_lon[storm,0:i],tc_lat[storm,0:i],color=get_tc_color(tc_wind[storm,:,0].max()),alpha=1,linewidth=0.001)
                    break

axes[3].set_title('TCs starting outside the domain: '+str(len(tmp[3])),fontsize=8)



plt.savefig('plots/TC_atlantic.pdf')


for i in range(4):
    for storm in tmp[i].keys():
        tmp[i][storm][0].set_alpha(0.2)
        tmp[i][storm][0].set_linewidth(1)


plt.savefig('plots/TC_atlantic.png',dpi=300)
