from mpl_toolkits.basemap import Basemap, cm
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime
from netCDF4 import Dataset
import glob, os
from scipy import stats
import numpy.ma as ma
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Polygon
from math import pi, cos, radians
import numpy.matlib
from pylab import *
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import interpolate
from datetime import datetime


'''

Code to prepare Kazu's 3-D CO production fields for inversion

'''

#(surface)                                                                                               
Ap =np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
              3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
              7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
              1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
              2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
              2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
              7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 1.979160e+01, 9.292942e+00,
              4.076571e+00, 1.650790e+00, 6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])
#(top of atmosphere)                                                                                     
#(surface)                                                                                               
Bp = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
               8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
               7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
               5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
               2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
               6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])
#(top of atmosphere) 

Ap1 = np.repeat(Ap[:, np.newaxis], 91, axis=1)
Ap2 = np.repeat(Ap1[:,:, np.newaxis], 144, axis=2)

Bp1 = np.repeat(Bp[:, np.newaxis], 91, axis=1)
Bp2 = np.repeat(Bp1[:,:, np.newaxis], 144, axis=2)

days_in_months = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

for yy in range(1):
    year = yy + 2019


    ncfile='/u/bbyrne1/prod_co_6hr_'+str(year).zfill(4)+'_p_20221114.nc'
    f=Dataset(ncfile,mode='r')
    lon1o=f.variables['lon'][:]
    lat1=f.variables['lat'][:]
    PS1=f.variables['lev'][:] # pressure
    time_day1=f.variables['time'][:]/(60.*24.) # minutes since 2003-12-11 00:00
    co1o=f.variables['co'][:] # ppb/day #(time,lev,lat,lon)
    f.close()
    
    co1o[np.where(co1o.mask)]=-1e10
    lon1 = np.append(lon1o[72:144]-360,lon1o[0:72])
    co1 = np.append(co1o[:,:,:,72:144].data,co1o[:,:,:,0:72].data,axis=3)
    
    day_of_year = 0
    for month in range(12):
        for day in range(days_in_months[month]):
            ncfile='/nobackup/mlee7/GEOS/GEOS_2x2.5/MERRA2/'+str(year).zfill(4)+'/'+str(month+1).zfill(2)+'/MERRA2.'+str(year).zfill(4)+str(month+1).zfill(2)+str(day+1).zfill(2)+'.I3.2x25.nc4'
            print ncfile
            f=Dataset(ncfile,mode='r')
            time_day2=f.variables['time'][:]/(60.*24.)
            lon2=f.variables['lon'][:]
            lat2=f.variables['lat'][:]
            lev2=f.variables['lev'][:]
            PS2=f.variables['PS'][:]/100. # hPa
            f.close()
            
            CO_regrid = np.zeros((24,47,91,144))
            for loop in range(4):
                i=loop*2
                PStemp = PS2[i,:,:]
                Pedge = Ap2 + ( Bp2 * np.repeat(PStemp[np.newaxis,:,:],48,axis=0) )
                Pmid = (Pedge[0:47,:,:]+Pedge[1:48,:,:])/2.
                
                co1_now = co1[day_of_year*4+loop,:,:,:]
                
                for x in range(91):
                    for y in range(144):
                        Pmid_temp = Pmid[:,x,y]
                        co1_temp = co1_now[:,x,y]
                        
                        co_interp = np.zeros(47)
                        
                        IND = np.where(co1_temp>0)
                        if np.size(IND[0]>0):
                            f = interpolate.interp1d(PS1[IND], co1_temp[IND], fill_value='extrapolate')
                            
                            It = np.where(Pmid_temp>45)
                            co_interp[It] = f(Pmid_temp[It])
                            co_interp[It[0][-1]:47] = co_interp[It[0][-1]]
                            
                            CO_regrid[loop*6+0,:,x,y] = co_interp
                            CO_regrid[loop*6+1,:,x,y] = co_interp
                            CO_regrid[loop*6+2,:,x,y] = co_interp
                            CO_regrid[loop*6+3,:,x,y] = co_interp
                            CO_regrid[loop*6+4,:,x,y] = co_interp
                            CO_regrid[loop*6+5,:,x,y] = co_interp
                            
                            # Convert [ppb CO / day] to [kg C / kg dry air /s]
                            
                            # (1day/60*60*24s)*(1e-9(molCO/molair)/ppb)*(12.01g/mol C)*(1mol/28.97g dry air) 

                            #  (1/60*60*24s) * (1e-9molCO/molair) * (12.01g/mol C) * (1mol/28.97g dry air) 
                            
            CO_regrid_out = CO_regrid * ( (12.01 * 1e-9) / (60. * 60. * 24. * 28.97) )
                            
            nc_out='/nobackupp17/bbyrne1/CO_production_2x25/'+str(year).zfill(4)+'/'+str(month+1).zfill(2)+'/'+str(day+1).zfill(2)+'.nc'
            dataset = Dataset(nc_out,'w')
            times = dataset.createDimension('time',24)
            lats = dataset.createDimension('lat',91)
            lons = dataset.createDimension('lon',144)
            levs = dataset.createDimension('lev',47)
            FCO_prods = dataset.createVariable('CO_prodDens', np.float64, ('time','lev','lat','lon'))
            FCO_prods[:,:,:,:] = CO_regrid_out
            FCO_prods.unit = 'kgC/kgAir/s'
            dataset.close()
            
            day_of_year = day_of_year+1
            
