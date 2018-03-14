

from scipy import interpolate

y_i=range(60,80)
x_i=range(80,100)
data=T[612,50:70,65:85]

X=lons[50:70,65:85]
Y=lats[50:70,65:85]

f = interpolate.interp2d(X[0,:], Y[:,0], data, kind='cubic')
zoom_factor=100
zoomed_lat_spacing=abs(np.diff(lats[:,0],1)[0]/zoom_factor)
zoomed_lon_spacing=abs(np.diff(lons[0,:],1)[0]/zoom_factor)

xnew = np.arange(X.min(),X.max(), zoomed_lon_spacing)
ynew = np.arange(Y.min(),Y.max(), zoomed_lat_spacing)
znew = f(xnew, ynew)


plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=2)
axes[0,0].pcolormesh(X,Y,data)
axes[0,1].pcolormesh(xnew,ynew,znew)

y_max,x_max=np.where(data==np.max(data))[0][0],np.where(data==np.max(data))[1][0]
y_zoomed_max=np.argmin(abs(ynew-Y[y_max,0]))
axes[1,0].plot(X[y_max,:], data[y_max, :], 'ro-', xnew, znew[y_zoomed_max, :], 'b-')
x_zoomed_max=np.argmin(abs(xnew-X[0,x_max]))
axes[1,1].plot(Y[:,x_max], data[:,x_max], 'ro-', ynew, znew[:,x_zoomed_max], 'b-')

print(data.min(),znew.min())
print(data.max(),znew.max())

plt.savefig('plots/test_bicubic.png')

#
# plt.close('all')
# fig,axes=plt.subplots(nrows=2,ncols=2)
# zoomed_lat_spacing=abs(np.diff(lats.mean(axis=1),1)[0]/float(20))
# zoomed_lon_spacing=abs(np.diff(lons.mean(axis=0),1)[0]/float(20))
# y_interpol=lats.mean(axis=1)
# x_interpol=lons.mean(axis=0)
# y_zoomed=np.arange(lats.mean(axis=1).min(),lats.mean(axis=1).max(), zoomed_lat_spacing)
# x_zoomed=np.arange(lons.mean(axis=0).min(),lons.mean(axis=0).max(), zoomed_lon_spacing)
#
#
# f = interpolate.interp2d(x_interpol,y_interpol,MSLP[613,:,:], kind='cubic')
# test=f(x_zoomed,y_zoomed)
#
#
# axes[0,0].pcolormesh(x_interpol,y_interpol,MSLP[613,:,:])
# axes[0,1].pcolormesh(x_zoomed,y_zoomed,test)
# print(MSLP[613,:,:].min(),test.min())
# plt.savefig('plots/detect_test_zoom.png')
#
#
#
# #
