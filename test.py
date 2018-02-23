
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

nc = Dataset('data/CAM25/item16222_daily_mean/item16222_daily_mean_o6uk_2017-06_2017-10.nc')
lats = nc.variables['global_latitude0'][:]
lons = nc.variables['global_longitude0'][:]
rlats = nc.variables['latitude0'][:]
rlons = nc.variables['longitude0'][:]
rlons, rlats = np.meshgrid(rlons, rlats)
data = nc.variables['item16222_daily_mean'][0,0,:,:].squeeze()
data = np.ma.masked_values(data,-999.)
rotpole = nc.variables['rotated_pole0']


def normalize180(lon):
    """Normalize lon to range [180, 180)"""
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
o_lon_p = rotpole.grid_north_pole_longitude
o_lat_p = rotpole.grid_north_pole_latitude
print( rotpole )
print( 'lon_0,o_lon_p,o_lat_p=',lon_0,o_lon_p,o_lat_p)

plt.close('all')
tmp,txt=[],[]
m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
           llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
           llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
x,y = m(lons,lats)
m.drawcoastlines()
#m.contourf(x,y,data,20)
im=m.imshow(data,interpolation='none')
im.set_cmap('bone'); ax.autoscale(False)
m.drawmeridians(np.arange(-180,180,20))
m.drawparallels(np.arange(20,80,20))
m.colorbar()
ax=plt.gca()
ax.invert_yaxis()
txt.append(ax.text(80,-3,'bla'))
xx,yy=m([-100,-90],[20,30])
tmp.append(m.plot(xx,yy,'r'))
tmp.append(m.plot(xx,yy,linestyle='',marker='*'))
plt.title('rotated pole data in native map using real sphere corner lat/lons' )
plt.savefig('plots/test.png')


for imm in ax.images:
    ax.images.remove(imm)
for element in tmp:
    l = element.pop(0); wl = weakref.ref(l); l.remove(); del l
for element in txt:
    element.remove()

plt.savefig('plots/test2.png')

plt.figure()
m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
           llcrnry = rlats[0,0], urcrnry = rlats[-1,-1],\
           llcrnrx = rlons[0,0], urcrnrx = rlons[-1,-1],resolution='c')
x,y = m(lons,lats)
m.drawcoastlines()
m.contourf(x,y,data,20)
m.drawmeridians(np.arange(-180,180,20))
m.drawparallels(np.arange(20,80,20))
m.colorbar()

x,y=m([-60,-30],[10,10])
m.plot(x,y,'-*r')
plt.title('rotated pole data in native map using rotated sphere corner lat/lons' )

plt.savefig('plots/test.png')
