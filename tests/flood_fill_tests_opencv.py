

# import cv2
# os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
# mslp=da.read_nc('data/ERA5/atl_2016_surface.nc')['MSL'].values[483,:,:]
# ta=da.read_nc('data/ERA5/atl_2016_T.nc')['T'].values[483,1,:,:]

def find_closed_contours(field,y,x,step=2,radius=30,n_contours=None):
    '''
    finds closed contours around center
    center is supposed to be a minimum
    '''

    y_min=int(max(0,y-radius))
    y_max=int(min(field.shape[0],y+radius+1))
    x_min=int(max(0,x-radius))
    x_max=int(min(field.shape[1],x+radius+1))

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
                return cont[-1,:,:].copy()
            elif (step<0 and field[y_,x_]>threshold) or (step>0 and field[y_,x_]<threshold):
                cont_tmp[y_,x_]=1
                p_searched+=[(y_,x_)]
                if (y_-1,x_) not in p_searched+p_: p_.append((y_-1,x_))
                if (y_+1,x_) not in p_searched+p_: p_.append((y_+1,x_))
                if (y_,x_-1) not in p_searched+p_: p_.append((y_,x_-1))
                if (y_,x_+1) not in p_searched+p_: p_.append((y_,x_+1))
        cont=np.concatenate((cont,cont_tmp.reshape(1,ny,nx)))
        ncont+=1
    return cont[-1,:,:].copy()

def find_closed_contours_opencv(field,y,x,step=2,search_radius=30,n_contours=None):
    '''
    finds closed contours around center
    center is supposed to be a minimum
    '''
    y_min=int(max(0,y-search_radius))
    y_max=int(min(field.shape[0],y+search_radius+1))
    x_min=int(max(0,x-search_radius))
    x_max=int(min(field.shape[1],x+search_radius+1))

    if step<0:
        loc=-field[y_min:y_max,x_min:x_max].copy()
    else:
        loc=field[y_min:y_max,x_min:x_max].copy()

    threshold=field[y,x]

    ny,nx=loc.shape
    running=True
    ncont=0
    cont=np.zeros((ny, nx))
    while n_contours is None or ncont<n_contours:
        threshold+=step
        print('*** '+str(threshold))
        th, im_floodfill = cv2.threshold(loc, threshold, 255, cv2.THRESH_BINARY_INV);
        mask = np.zeros((ny+2, nx+2), np.uint8)
        cv2.floodFill(im_floodfill, mask, (y-y_min,x-x_min), 1);
        y_,x_=np.where(mask[1:-1,1:-1]==1)
        print(y_,x_)
        if 0 in list(y_)+list(x_) or ny-1 in y_ or nx-1 in x_:
            cont_out=np.zeros(field.shape)*np.nan
            #cont=np.array(mask[1:-1,1:-1],np.float)
            cont[cont==0]=np.nan
            cont_out[y_min:y_max,x_min:x_max]=cont
            return cont_out
            break
        cont=np.array(mask[1:-1,1:-1],np.float)
        ncont+=1

    cont_out=np.zeros(field.shape)*np.nan
    cont[cont==0]=np.nan
    cont_out[y_min:y_max,x_min:x_max]=cont
    return cont_out


mslp=MSLP[483,:,:]

plt.close('all')
fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(10,5))
axes[0,0].imshow(mslp,interpolation='none')
axes[1,0].imshow(ta,interpolation='none')
axes[0,1].imshow(mslp,interpolation='none')
axes[1,1].imshow(ta,interpolation='none')

coords=peak_local_max(-mslp, min_distance=10)
start = time.time()
for y_p,x_p in zip(coords[:,0],coords[:,1]):
    cont_p=find_closed_contours(mslp,y_p,x_p,radius=20)
    yy,xx=np.where(cont_p==1)
    if len(xx)>1:
        axes[0,0].plot(xx,yy,'.g',linestyle='')
    y,x=np.where(ta[y_p-1:y_p+1,x_p-1:x_p+1]==ta[y_p-1:y_p+1,x_p-1:x_p+1].max())
    y_t,x_t=y_p-1+y[0],x_p-1+x[0]
    cont_t=find_closed_contours(ta,y_t,x_t,step=-1,radius=96,n_contours=1)
print(time.time() - start)

coords=peak_local_max(-mslp, min_distance=10)
start = time.time()
for y_p,x_p in zip(coords[:,0],coords[:,1]):
    start = time.time()
    print(mslp[y_p,x_p])
    cont_p=find_closed_contours_opencv(mslp,y_p,x_p,search_radius=20,step=-2)
    yy,xx=np.where(cont_p==1)
    if len(xx)>0:
        axes[0,1].plot(xx,yy,'.g',linestyle='')

print(time.time() - start)

for ax in axes.flatten():
    ax.axis('off')

plt.tight_layout()
plt.savefig('plots/contour_test.png')
