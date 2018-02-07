import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import weak


plt.close()
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,5))
m = Basemap(ax=ax,llcrnrlon=min(lon),urcrnrlon=max(lon),llcrnrlat=min(lat),urcrnrlat=max(lat),resolution="l",projection='cyl')
m.drawmapboundary(fill_color='1.')
m.drawmapboundary(fill_color='#99ffff')
m.fillcontinents(color='#cc9966',lake_color='#99ffff')
m.drawcoastlines(linewidth=0.3)
m.drawparallels(np.arange(-60,100,30),labels=[0,0,0,0],color='grey',linewidth=0.5)
m.drawmeridians([-120,0,120],labels=[0,0,0,0],color='grey',linewidth=0.5)
lines=m.plot([300,320,300,290,280],[10,20,30,20,30],'-k')
title=ax.set_title('bla')
plt.savefig('plots/detected/track_1.png')

l = lines.pop(0); wl = weakref.ref(l); l.remove(); del l

lines=m.plot([300,280],[10,30],'-r')
title=ax.set_title('asdasd')
plt.savefig('plots/detected/track_2.png')
