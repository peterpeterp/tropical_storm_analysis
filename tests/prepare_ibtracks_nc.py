import os,sys,glob,time,collections,gc
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import matplotlib.pylab as plt
import dimarray as da
import pandas as pd

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')


dat=pd.read_csv('data/Allstorms.ibtracs_all.v03r10.csv',header=1)

dat=dat.drop(dat.index[[0]])

names=np.asarray(sorted(set(dat['Serial_Num'])))

#lenghts=[len(dat.loc[dat['Serial_Num'] == id_]) for id_ in names]
#print max(lengths)

TC={'Sub_basin':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time'],dtype='|S5'),
    'ISO_time':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time'],dtype='|S20'),
    'Nature':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time'],dtype='|S5'),
    'Latitude':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Longitude':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Wind(WMO)':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Pres(WMO)':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Latitude_for_mapping':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Longitude_for_mapping':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Name':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time'],dtype='|S10'),
    'Basin':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time'],dtype='|S5'),
    'Season':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
    'Num':da.DimArray(axes=[names,np.arange(0,200,1)],dims=['ID','time']),
}

for var in ['Sub_basin','ISO_time','Nature','Name','Basin']:
    TC[var].values[:,:]=''


for id_,i in zip(names,range(len(names))):
    rows=dat.loc[dat['Serial_Num'] == id_]
    print i,list(rows['Name'])[0]

    for var in TC.keys():
        tmp=list(rows[var])
        TC[var][id_,0:len(tmp)-1]=tmp


ds=da.Dataset(TC)
ds.write_nc('data/TC_ibtracs.nc', mode='w')
