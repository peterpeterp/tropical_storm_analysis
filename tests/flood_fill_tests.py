from __future__ import print_function
import os,sys,glob,time,collections,gc,calendar,weakref,resource,datetime
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.ndimage as ndimage
import cv2
from skimage.feature import peak_local_max

def find_closed_contours_opencv(field,y,x,tabus=[],step=2,search_radius=30,n_contours=None,method='min'):
    '''
    finds closed contours around center
    center is supposed to be a minimum
    '''
    any_in = lambda a, b: any(i in b for i in a)

    y_min=int(max(0,y-search_radius))
    y_max=int(min(field.shape[0],y+search_radius+1))
    x_min=int(max(0,x-search_radius))
    x_max=int(min(field.shape[1],x+search_radius+1))

    if method=='min':
        im=-field.copy()
        threshold=-field[y,x]
    else:
        im=field.copy()
        threshold=field[y,x]

    ny,nx=im.shape
    running=True
    ncont=0
    cont=np.zeros((ny, nx))
    while n_contours is None or ncont<n_contours:
        threshold-=step
        th, im_floodfill = cv2.threshold(im, threshold, 255, cv2.THRESH_BINARY_INV);
        mask = np.zeros((ny+2, nx+2), np.uint8)
        cv2.floodFill(im_floodfill, mask, (x,y), 1);
        y_,x_=np.where(mask[1:-1,1:-1]==1)
        if y_min in y_ or x_min in x_ or y_max-1 in y_ or x_max-1 in x_:
            cont[cont==0]=np.nan
            return cont,ncont
        # if any_in(tabus,[[yy,xx] for yy,xx in zip(y_,x_)]):
        #     print('tabu')
        #     cont[cont==0]=np.nan
        #     return cont,ncont
        cont=np.array(mask[1:-1,1:-1],np.float)
        ncont+=1

    cont[cont==0]=np.nan
    return cont,ncont

identifier='2016'
os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
data_path='data/ERA5/'
# MSLP=da.read_nc(data_path+'atl_'+identifier+'_surface.nc')['MSL'].values/100.
# T=da.read_nc(data_path+'atl_'+identifier+'_T.nc')['T'].values
for t in range(700,710):
    mslp=MSLP[t,:,:]
    ta=T[t,1,:,:]

    mslp_smoo=ndimage.gaussian_filter(mslp,sigma=(3,3))

    plt.close('all')
    fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(10,5))
    axes[0,0].imshow(mslp_smoo,interpolation='none')
    axes[1,0].imshow(ta,interpolation='none')
    axes[0,1].imshow(mslp,interpolation='none')
    axes[1,1].imshow(ta,interpolation='none')

    coords=peak_local_max(-mslp_smoo, min_distance=7)
    for y_p,x_p in zip(coords[:,0],coords[:,1]):
        axes[0,0].plot(x_p,y_p,'.r',linestyle='')
        cont_p,ncont=find_closed_contours_opencv(mslp_smoo,y_p,x_p,search_radius=30,tabus=coords.tolist(),step=1,method='min')
        yy,xx=np.where(cont_p==1)
        axes[1,1].plot(x_t,y_t,'*c')
        if len(xx)>0:
            axes[0,1].plot(xx,yy,'.g',linestyle='')
            y,x=np.where(ta[y_p-2:y_p+2,x_p-2:x_p+2]==ta[y_p-2:y_p+2,x_p-2:x_p+2].max())
            y_t,x_t=y_p-1+y[0],x_p-1+x[0]
            cont_t,ncont=find_closed_contours_opencv(ta,y_t,x_t,step=2,search_radius=96,n_contours=1,method='max')
            yy,xx=np.where(cont_t==1)
            if len(xx)>0:
                axes[1,1].plot(xx,yy,'.r',linestyle='')


    for ax in axes.flatten():
        ax.axis('off')

    plt.tight_layout()
    plt.savefig('plots/contour_test'+str(t)+'.png')
