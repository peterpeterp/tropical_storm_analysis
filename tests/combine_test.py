
self=found_tcs


self._add_name='thresholds'
out_file=self._working_dir+'detected_positions_thresholds_'+self._add_name+'.nc'
if overwrite and os.path.isfile(out_file):
    os.system('rm '+out_file)
elif overwrite==False and os.path.isfile(out_file):
    self._detected=da.read_nc(out_file)['detected']
    return self._detected

'''
Required Thresholds:
thr_vort:       float [1/s]:    threshold for relative vorticity maxima
dis_vort_max:   float [deg]:    minimal distance between vorticity maxima
dis_cores:      float [deg]:    maximal distances between vort max and mslp min (mslp min and warm core)
thr_MSLP_inc:   float [hPa]:    increase in MSLP from center over dis_MSLP_inc
dis_MSLP_inc:   float [deg]:    distance over which MSLP should increase by thr_MSLP_inc
thr_T_drop:     float [K]:      temperature drop from warm core center over dis_T_decrease
dis_T_drop:     float [deg]:    distance over which T should decrease by thr_T_decrease
tc_size:        float [deg]:    radius within which maximal Wind10 is searched for
'''

# convert distances from degrees into grid-cells
for distance in [dis_vort_max,dis_cores,dis_MSLP_inc,dis_T_drop,tc_size]:
    distance=self.degree_to_step(distance)

detect=np.array([[np.nan]*7])
print('detecting\n10------50-------100')
for t,progress in zip(self._time_i,np.array([['-']+['']*(len(self._time_i)/20+1)]*20).flatten()[0:len(self._time_i)]):
    sys.stdout.write(progress); sys.stdout.flush()
    # i vort max
    vo_=self._VO[t,:,:]
    vo_[vo_<thr_vort]=0
    coords = peak_local_max(vo_, min_distance=int(dis_vort_max))
    if coords.shape[0]>0:
        for y_v,x_v in zip(coords[:,0],coords[:,1]):
            print(y_v,x_v)
            y_,x_=self.circle_around(y_v,x_v,dis_cores)
            print(y_,x_)
            tmp=self._MSLP[t,y_,x_].flatten()
            print(np.where(tmp==tmp.min()))
            print(tmp)
            i=np.where(tmp==tmp.min())[0][0]; y_p,x_p=y_[i],x_[i]
            print(tmp)
            print(i)
            print(y_p,x_p)
            asdas
            y_circ,x_cirx=self.circle_around(y_p,x_p,dis_MSLP_inc)
            # ii relative pressure min
            if self._MSLP[t,y_circ,x_circ].min()-self._MSLP[t,y_p,x_p]>thr_MSLP_inc:
                box_1=self.get_box(y_p,x_p,cores_distance)
                box_2=self.get_box(y_p,x_p,tc_size)
                tmp=[t,y_p,x_p,1,0,0,0]
                # iv warm core
                if self._T is None:
                    tmp[4]=1
                elif self._T[t,1,box_1[0]:box_1[1],box_1[2]:box_1[3]].max()-self._T[t,1,box_2[0]:box_2[1],box_2[2]:box_2[3]].mean()>thr_T_decrease:
                    tmp[4]=1
                # iii wind speed
                tmp[5]=self._MSLP[t,y_p,x_p]
                tmp[6]=self._Wind10[t,box_2[0]:box_2[1],box_2[2]:box_2[3]].max()
                detect=np.concatenate((detect,np.array([tmp])))

self._detected=da.DimArray(np.array(detect[1:,:]),axes=[range(detect.shape[0]-1),['t','y','x','pressure_low','warm_core','MSLP','Wind10']],dims=['ID','z'])
da.Dataset({'detected':self._detected}).write_nc(out_file,mode='w')
print('done')
return self._detected
