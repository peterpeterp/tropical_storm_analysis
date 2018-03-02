
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import dimarray as da

hybrid=[(0,1),\
        (0.00000000000000E+00,9.97000000000000E-01),\
        (0.00000000000000E+00,9.94000000000000E-01),\
        (0.00000000000000E+00,9.89000000000000E-01),\
        (0.00000000000000E+00,9.82000000000000E-01),\
        (0.00000000000000E+00,9.72000000000000E-01),\
        (0.00000000000000E+00,9.60000000000000E-01),\
        (0.00000000000000E+00,9.46000000000000E-01),\
        (1.33051011276943E+02,9.26669489887231E-01),\
        (3.64904148871589E+02,9.04350958511284E-01),\
        (6.34602716447362E+02,8.79653972835526E-01),\
        (9.59797167291774E+02,8.51402028327082E-01),\
        (1.34768004165515E+03,8.19523199583449E-01),\
        (1.79090739595110E+03,7.85090926040489E-01),\
        (2.29484168994850E+03,7.48051583100515E-01),\
        (2.84748477771176E+03,7.09525152222882E-01),\
        (3.46887148811864E+03,6.68311285118814E-01),\
        (4.16295646296916E+03,6.24370435370308E-01),\
        (4.89188083250491E+03,5.80081191674951E-01),\
        (5.67182423980408E+03,5.34281757601959E-01),\
        (6.47671299638532E+03,4.88232870036147E-01),\
        (7.29746989472049E+03,4.42025301052795E-01),\
        (8.12215979124915E+03,3.95778402087509E-01),\
        (8.91408220106234E+03,3.50859177989377E-01),\
        (9.65618191050164E+03,3.07438180894984E-01),\
        (1.03294361777746E+04,2.65705638222254E-01),\
        (1.09126384442387E+04,2.25873615557613E-01),\
        (1.13696478308432E+04,1.89303521691568E-01),\
        (1.16953715974700E+04,1.55046284025300E-01),\
        (1.18612530873948E+04,1.24387469126052E-01),\
        (1.18554343163493E+04,9.64456568365075E-02),\
        (1.16633553655803E+04,7.23664463441966E-02),\
        (1.12854040644942E+04,5.21459593550578E-02),\
        (1.07299494055679E+04,3.57005059443214E-02),\
        (1.00146150535107E+04,2.28538494648935E-02),\
        (9.16724703583310E+03,1.33275296416689E-02),\
        (8.22624490770442E+03,6.73755092295582E-03),\
        (7.20156898029828E+03,2.48431019701722E-03),\
        (6.08867300853392E+03,1.13269914660783E-04),\
        (4.95000000000000E+03,0.00000000000000E+00),\
        (4.00000000000000E+03,0.00000000000000E+00),\
        (3.23000000000000E+03,0.00000000000000E+00),\
        (2.61000000000000E+03,0.00000000000000E+00),\
        (2.10500000000000E+03,0.00000000000000E+00),\
        (1.70000000000000E+03,0.00000000000000E+00),\
        (1.37000000000000E+03,0.00000000000000E+00),\
        (1.10500000000000E+03,0.00000000000000E+00),\
        (8.93000000000000E+02,0.00000000000000E+00),\
        (7.20000000000000E+02,0.00000000000000E+00),\
        (5.81000000000000E+02,0.00000000000000E+00),\
        (4.69000000000000E+02,0.00000000000000E+00),\
        (3.77000000000000E+02,0.00000000000000E+00),\
        (3.01000000000000E+02,0.00000000000000E+00),\
        (2.37000000000000E+02,0.00000000000000E+00),\
        (1.82000000000000E+02,0.00000000000000E+00),\
        (1.36000000000000E+02,0.00000000000000E+00),\
        (9.70000000000000E+01,0.00000000000000E+00),\
        (6.50000000000000E+01,0.00000000000000E+00),\
        (3.90000000000000E+01,0.00000000000000E+00),\
        (2.00000000000000E+01,0.00000000000000E+00),\
        (0.00000000000000E+00,0.00000000000000E+00)]

def pressure(i,ps):
    a1=hybrid[i][0]
    b1=hybrid[i][1]
    a2=hybrid[i+1][0]
    b2=hybrid[i+1][1]

    pk1=a1+b1*ps
    pk2=a2+b2*ps

    pk=np.exp(1/(pk1-pk2)*(pk1*np.log(pk1)-pk2*np.log(pk2))-1)
    return pk

os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')

nc = da.read_nc('data/JRA55/test_with_hybrid_levels/atl_U.nc')
u=nc['var33']

nc = da.read_nc('data/JRA55/test_with_hybrid_levels/atl_MSLP.nc')
ps=nc['var1']


u_out=np.zeros([u.shape[0],2,u.shape[2],u.shape[3]])*np.nan

for t_i,t in enumerate(u.time):
    print(t_i,t)


    p_=np.zeros([60,ps.shape[1],ps.shape[2]])
    for i in range(60):
        p_[i,:,:]=pressure(i,ps[t,:,:])

    for i,y in enumerate(ps.lat):
        for j,x in enumerate(ps.lon):
            tmp=p_[:,i,j]
            for p_i,pres in zip(range(2),[85000,100]):
                pos=np.where(tmp>pres)[0]
                if len(pos)!=0:
                    p1,p2=pos[-1],pos[-1]+1
                    u_out[t_i,p_i,i,j]=(u.ix[t_i,p1,i,j]+u.ix[t_i,p2,i,j])*0.5


uuu=da.DimArray(u_out,axes=[u.time,[85000,100],u.lat,u.lon],dims=['time','plev','lat','lon'])
uuu.axes['plev'].units='Pa'
uuu.axes['plev'].standard_name='air_pressure'
uuu.axes['plev'].longname='pressure'
uuu.axes['plev'].positive='down'
uuu.axes['plev'].axis='Z'

uuu.axes['lat'].standard_name = "latitude"
uuu.axes['lat'].long_name = "latitude"
uuu.axes['lat'].units = "degrees_north"
uuu.axes['lat'].axis = "Y"

uuu.axes['lon'].standard_name = "longitude"
uuu.axes['lon'].long_name = "longitude"
uuu.axes['lon'].units = "degrees_east"
uuu.axes['lon'].axis = "X"

da.Dataset({'u':uuu}).write_nc('data/JRA55/test_with_hybrid_levels/atl_U__.nc')


# os.chdir('/Users/peterpfleiderer/Documents/Projects/tropical_cyclones/')
#
# nc = da.read_nc('data/JRA55/test_with_hybrid_levels/U.nc')
# u=nc['var33']
# u.table=int(3)
#
# nc = da.read_nc('data/JRA55/test_with_hybrid_levels/MSLP.nc')
# ps=nc['var1']
# ps.table=int(3)
#
# topo=da.read_nc('data/JRA55/test_with_hybrid_levels/topo.nc')['var8']
#
# hyai=da.DimArray([hybrid[i][0] for i in range(61)],axes=[range(61)],dims=['nhyi'])
# hyam=da.DimArray([(hybrid[i][0]+hybrid[i+1][0])*0.5 for i in range(60)],axes=[range(60)],dims=['nhym'])
# hybi=da.DimArray([hybrid[i][1] for i in range(61)],axes=[range(61)],dims=['nhyi'])
# hybm=da.DimArray([(hybrid[i][1]+hybrid[i+1][1])*0.5 for i in range(60)],axes=[range(60)],dims=['nhym'])
#
# ds=da.Dataset({'var1':ps[u.time,:,:],'var8':topo,'var11':u,'hyai':hyai,'hyam':hyam,'hybi':hybi,'hybm':hybm})
# ds.write_nc('data/JRA55/test_with_hybrid_levels/U__.nc')
