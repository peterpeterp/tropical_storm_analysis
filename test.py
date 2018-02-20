


def wrf_vort(U,V,x,y ):
    """adapted from wrftools: github.com/keltonhalbert/wrftools"""
    xx,yy = np.meshgrid(x,y)
    xx*=np.cos(np.radians(yy))*6371000*2*np.pi/360.
    yy*=6371000*2*np.pi/360.
    dy=np.gradient(yy)
    dx=np.gradient(xx)

    vort=U.copy*np.nan
    for i in range(U.shape[0]):
        du = np.gradient( U )
        dv = np.gradient( V )
        vort[i,:,:]= dv[-1]/dx[-1] - du[-2]/dy[-2]

    return vort

# u=da.read_nc('data/ERA5/atl_2016_UV.nc')['U']
# v=da.read_nc('data/ERA5/atl_2016_UV.nc')['V']
#
# VO=da.read_nc('data/ERA5/atl_2016_VO.nc')['vo']

vo=wrf_vort(u.ix[0,0,:,:],v.ix[0,0,:,:],u.lon,u.lat)


plt.close('all')
fig,axes = plt.subplots(nrows=1,ncols=2)
axes[0].imshow(vo)
axes[1].imshow(VO.ix[0,0,:,:])
plt.savefig('plots/test.png')
