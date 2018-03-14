from scipy import interpolate
x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**2)
f = interpolate.interp2d(x, y, z, kind='cubic')



plt.close('all')
data=MSLP[3,20:40,20:40]

X=lons[20:40,20:40]
Y=lats[20:40,20:40]

y_zoom=np.arange(0,100,1)
x_zoom=np.arange(0,100,1)
yy,xx=np.meshgrid(y_zoom,x_zoom)


new_grid=scipy.interpolate.griddata((X.flatten(),Y.flatten()),data.flatten() , (xx,yy),method='cubic')
plt.close('all')
fig,axes=plt.subplots(nrows=1,ncols=2)

axes[0].imshow(data)


plt.savefig('plots/test_bicubic.png')






from scipy import interpolate

y_i=range(60,80)
x_i=range(80,100)
data=MSLP[613,65:85,65:85]

X=lons[65:85,65:85]
Y=lats[65:85,65:85]

f = interpolate.interp2d(X[0,:], Y[:,0], data, kind='cubic')

xnew = np.arange(X.min(),X.max(), 1e-2)
ynew = np.arange(Y.min(),Y.max(), 1e-2)
znew = f(xnew, ynew)


plt.close('all')
fig,axes=plt.subplots(nrows=1,ncols=2)
axes[0].pcolormesh(X,Y,data)
axes[1].pcolormesh(xnew,ynew,znew)
#axes[1].plot(X[0,:], data[0, :], 'ro-', xnew, znew[0, :], 'b-')

print(data.min(),znew.min())

plt.savefig('plots/test_bicubic.png')






#
