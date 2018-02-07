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



def method_plot(option=5):
    plt.close()
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    m = Basemap(projection='ortho',lon_0=-50,lat_0=20,resolution='l',ax=ax)
    m.bluemarble(scale=0.2)   # full scale will be overkill


    for id_ in TC.storm:
        if TC['basin'][id_,0]==0.0 and TC['season'][id_]>2014:
            x,y=TC['lon_for_mapping'][id_,:].values,TC['lat_for_mapping'][id_,:].values
            x1,y1=m(x[np.isfinite(x)],y[np.isfinite(y)])

            m.scatter(x1[0],y1[0],marker='*',color='magenta',s=25)
            if option==1:
                m.plot(x1,y1,color='cyan',lw=0.7)

            if option>1:
                msw=TC['source_wind'][id_,:,0].values[np.isfinite(x)]
                color=np.array(['#ffffcc']*msw.shape[0],dtype='|S7')
                color[msw>82]='#ffe775'
                color[msw>95]='#ffc148'
                color[msw>112]='#ff8f20'
                color[msw>136   ]='#ff6060'
                color=[tuple(np.array(webcolors.hex_to_rgb(cc))/255.0) for cc in color]
                m.scatter(x1,y1,marker='o',color=color,s=15)

    plt.tight_layout
    plt.savefig('plots/method_'+str(option)+'.png')

for i in range(3):
    method_plot(i)
