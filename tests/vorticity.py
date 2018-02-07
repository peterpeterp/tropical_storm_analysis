import os,sys,glob,time,collections,gc,calendar
from datetime import datetime as datetime
from datetime import timedelta
from netCDF4 import Dataset,netcdftime,num2date
import dimarray as da
import numpy as np
import wrf

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

nc=Dataset('detection/ERA_interim/atl_interim_ocean_1979.nc')
data=da.read_nc('detection/ERA_interim/atl_interim_atmos_1979.nc')



wrf.vorticity.get_avo(nc, timeidx=0, method=u'cat', squeeze=True, cache=None, meta=True, _key=None)
