import os,sys,glob,time,collections,gc,calendar
from subprocess import Popen


plt.close()
fig,ax = plt.subplots(nrows=1,ncols=1)

ax.set_title('10m wind speed')
data=wind.values[i,:,:]
y_w,x_w = local_extrem(data,'max',threshold=12,neighborhood_size=15)

im=ax.imshow(data,vmin=0,vmax=13)
im.set_cmap('bone')
ax.autoscale(False)
ax.axis('off')


y,x=np.where(data>12)
ax.scatter(x,y)

ex_x=[]
ex_y=[]
for y_,x_ in zip(y,x):
    if wind.ix[i,y_,x_]==np.max(wind.ix[i,max([0,y_-3]):min([y_+3,len(lat)-1]),max([0,x_-3]):min([x_+3,len(lon)-1])]):
        ex_y.append(y_)
        ex_x.append(x_)

ax.scatter(x_w,y_w,color='r',s=60,facecolors='none')
ax.scatter(ex_x,ex_y,marker='*',color='m')


plt.savefig('plots/test.png', bbox_inches = 'tight')
