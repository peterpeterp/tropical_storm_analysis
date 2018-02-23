
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
data = nc.variables['item16222_daily_mean'][:,0,:,:].squeeze()
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
fig,axes=plt.subplots(nrows=2,ncols=3,figsize=(10,5))
for ax,filt in zip(axes.flatten(),range(6)):
    m= Basemap(ax=ax,projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
               llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
               llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
    x,y = m(lons,lats)
    m.drawcoastlines()
    #m.contourf(x,y,data,20)
    if filt>0:
        im=m.imshow(ndimage.gaussian_filter(data,sigma=(0,filt,filt))[0,:,:],interpolation='none')
        #im=m.imshow(ndimage.gaussian_filter(data[0,:,:],sigma=(filt,filt)),interpolation='none')
    else:
        im=m.imshow(data[0,:,:],interpolation='none')

    im.set_cmap('bone'); ax.autoscale(False)
    m.drawmeridians(np.arange(-180,180,20))
    m.drawparallels(np.arange(20,80,20))
    ax.invert_yaxis()
    ax.set_title('filter: '+str(filt))
plt.tight_layout()
plt.savefig('plots/filters.png')
