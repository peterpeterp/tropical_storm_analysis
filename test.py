

import cv2
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
mslp=da.read_nc('data/ERA5/atl_2016_surface.nc')['MSL'].values[483,:,:]
ta=da.read_nc('data/ERA5/atl_2016_T.nc')['T'].values[483,1,:,:]

from scipy import ndimage as ndi
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from skimage import data, img_as_float

im = img_as_float(data.coins())

# image_max is the dilation of im with a 20*20 structuring element
# It is used within peak_local_max function
image_max = ndi.maximum_filter(im, size=20, mode='constant')

# Comparison between image_max and im to find the coordinates of local maxima
coordinates = peak_local_max(im, min_distance=20)

# display results
fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharex=True, sharey=True,
                         subplot_kw={'adjustable': 'box-forced'})
ax = axes.ravel()
ax[0].imshow(im, cmap=plt.cm.gray)
ax[0].axis('off')
ax[0].set_title('Original')

ax[1].imshow(image_max, cmap=plt.cm.gray)
ax[1].axis('off')
ax[1].set_title('Maximum filter')

ax[2].imshow(im, cmap=plt.cm.gray)
ax[2].autoscale(False)
ax[2].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
ax[2].axis('off')
ax[2].set_title('Peak local max')

fig.tight_layout()

plt.savefig('plots/loc_min_test.png')


plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=4,figsize=(10,5))

y,x=np.where(mslp==mslp.min())
y_p,x_p=y[0],x[0]

start = time.time()
axes[0,0].imshow(mslp[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
cont_p=find_closed_contours(mslp,y_p,x_p,radius=20)
axes[1,0].imshow(cont_p[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
y,x=np.where(ta[y_p-1:y_p+1,x_p-1:x_p+1]==ta[y_p-1:y_p+1,x_p-1:x_p+1].max())
y_t,x_t=y_p-1+y[0],x_p-1+x[0]
axes[0,1].imshow(ta[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
cont_t=find_closed_contours(ta,y_t,x_t,step=-1,radius=10,n_contours=1)
axes[1,1].imshow(cont_t[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
print(time.time() - start)

start = time.time()
axes[0,2].imshow(mslp[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
cont_p=find_closed_contours_opencv(mslp,y_p,x_p,radius=20,step=-200)
axes[1,2].imshow(cont_p[y_p-50:y_p+50,x_p-50:x_p+50])
y,x=np.where(ta[y_p-1:y_p+1,x_p-1:x_p+1]==ta[y_p-1:y_p+1,x_p-1:x_p+1].max())
y_t,x_t=y_p-1+y[0],x_p-1+x[0]
axes[0,3].imshow(ta[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
cont_t=find_closed_contours(ta,y_t,x_t,step=-1,radius=10,n_contours=1)
axes[1,3].imshow(cont_t[y_p-50:y_p+50,x_p-50:x_p+50],interpolation='none')
print(time.time() - start)


plt.savefig('plots/contour_test.png')
