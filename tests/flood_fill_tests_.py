


def find_closed_contours(field,y,x,step=2,search_radius=30,n_contours=None):
    '''
    finds closed contours around center
    center is supposed to be a minimum
    '''

    y_min=int(max(0,y-search_radius))
    y_max=int(min(field.shape[0],y+search_radius+1))
    x_min=int(max(0,x-search_radius))
    x_max=int(min(field.shape[1],x+search_radius+1))

    loc=field[y_min:y_max,x_min:x_max].copy()

    if step<0:
        threshold=loc.max()
    else:
        threshold=loc.min()

    ny,nx=field.shape
    cont=np.zeros([1,ny,nx])*np.nan
    cont[0,y,x]=1
    running=True
    ncont=0
    while running and (n_contours is None or ncont<n_contours):
        threshold+=step
        cont_tmp=cont[-1,:,:].copy()
        p_=[(y,x) for y_,x_ in zip(np.where(cont_tmp==1)[0],np.where(cont_tmp==1)[1])]
        p_searched=p_[:]
        for pp in p_:
            y_,x_=pp
            if y_==y_min or y_==y_max or x_==x_min or x_==x_max:
                return cont[-1,:,:].copy(),ncont
            elif (step<0 and field[y_,x_]>threshold) or (step>0 and field[y_,x_]<threshold):
                cont_tmp[y_,x_]=1
                p_searched+=[(y_,x_)]
                if (y_-1,x_) not in p_searched+p_: p_.append((y_-1,x_))
                if (y_+1,x_) not in p_searched+p_: p_.append((y_+1,x_))
                if (y_,x_-1) not in p_searched+p_: p_.append((y_,x_-1))
                if (y_,x_+1) not in p_searched+p_: p_.append((y_,x_+1))
        cont=np.concatenate((cont,cont_tmp.reshape(1,ny,nx)))
        ncont+=1
    return cont[-1,:,:].copy(),ncont

def find_closed_contours_opencv(field,y,x,step=2,search_radius=30,n_contours=None,method='min'):
    '''
    finds closed contours around center
    center is supposed to be a minimum
    '''
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
            break
        cont=np.array(mask[1:-1,1:-1],np.float)
        ncont+=1

    cont[cont==0]=np.nan
    return cont,ncont
#
# y_p,x_p=coords[5,0],coords[5,1]
# plt.close('all')
# fig,axes=plt.subplots(nrows=9,ncols=3,figsize=(10,20)); axes=axes.flatten()
# nn=0
# for radius in [10,20,30,40,50,60,70,80,90]:
#     print('____________',y_p,x_p)
#     cont_p=find_closed_contours(mslp,y_p,x_p,search_radius=radius,step=2)
#     yy,xx=np.where(cont_p==1)
#     if len(xx)>1:
#         axes[nn].plot(xx,yy,'.g',linestyle='')
#     nn+=1
#     cont_p=find_closed_contours_opencv(mslp,y_p,x_p,search_radius=radius,step=-2)
#     nn+=1
#     yy,xx=np.where(cont_p==1)
#     if len(xx)>1:
#         axes[nn].plot(xx,yy,'.g',linestyle='')
#     nn+=1
#
# for ax in axes.flatten():
#     ax.axis('off')
# plt.tight_layout(h_pad=0, w_pad=0)
# plt.savefig('plots/contour_test.png')
# #


# import cv2
# os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
# mslp=da.read_nc('data/ERA5/atl_2016_surface.nc')['MSL'].values[483,:,:]
# ta=da.read_nc('data/ERA5/atl_2016_T.nc')['T'].values[483,1,:,:]

mslp=MSLP[708,:,:]

mslp_smoo=ndimage.gaussian_filter(MSLP[704,:,:],sigma=(3,3))

plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(10,5))
axes[0,0].imshow(mslp,interpolation='none')
axes[1,0].imshow(ta,interpolation='none')
axes[0,1].imshow(mslp,interpolation='none')
axes[1,1].imshow(ta,interpolation='none')

coords=peak_local_max(-mslp_smoo, min_distance=10)
for y_p,x_p in zip(coords[:,0],coords[:,1]):

    start = time.time()
    cont_p,ncont=find_closed_contours_opencv(mslp,y_p,x_p,search_radius=20,step=2)
    yy,xx=np.where(cont_p==1)
    #if len(xx)>1:
        #axes[0,0].plot(xx,yy,'.g',linestyle='')
    y,x=np.where(ta[y_p-1:y_p+1,x_p-1:x_p+1]==ta[y_p-1:y_p+1,x_p-1:x_p+1].max())
    y_t,x_t=y_p-1+y[0],x_p-1+x[0]
    cont_t,ncont=find_closed_contours(ta,y_t,x_t,step=-1,search_radius=30,n_contours=1)
    #print(time.time() - start)

    start = time.time()
    cont_p,ncont=find_closed_contours_opencv(mslp_smoo,y_p,x_p,search_radius=96,step=1,method='min')
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

    #print(time.time() - start)

for ax in axes.flatten():
    ax.axis('off')

plt.tight_layout()
plt.savefig('plots/contour_test.png')
