from mpl_toolkits.basemap import Basemap, cm
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime
from netCDF4 import Dataset
import glob, os
import numpy as np
import h5py # if this creates an error please make sure you have the h5py library
from datetime import datetime, timedelta
import numpy.ma as ma
import math



for yy in range(14):
    year = yy + 2009
    for mm in range(12):
        month = mm + 1
        b = glob.glob("/nobackup/bbyrne1/MOPITT_v9_xCO_Log_20221216/"+str(year).zfill(4)+"/"+str(month).zfill(2)+"/*.nc")
        a = np.sort(b)
        for file_in in a:
            f=Dataset(file_in,mode='r')
            lon_superobN=f.variables['longitude'][:]#.data
            lat_superobN=f.variables['latitude'][:]#.data
            hour_superobN=f.variables['time'][:]#.data
            Pprof=f.variables['pressure'][:].data
            xCOprior=f.variables['xCO'][:]#.data
            PWF=f.variables['xCO-pressureWeight'][:]#.data
            error_superobsN=f.variables['xCO-uncertainty'][:]#.data
            AK=f.variables['xCO-averagingKernel'][:].data
            prior_CO=f.variables['CO-apriori'][:]#.data
            xCO=f.variables['xCO'][:]#.data
            xCOapriori=f.variables['xCO-apriori'][:]#.data
            f.close()            


            III = np.where( np.logical_and(np.isfinite(lon_superobN) , np.logical_and(np.logical_and(np.isfinite(lat_superobN),np.isfinite(xCOprior)),np.logical_and(np.isfinite(error_superobsN),np.isfinite(xCO)))) )

            lon_superobNTo = lon_superobN[III] 
            lat_superobNTo = lat_superobN[III]
            hour_superobNTo = hour_superobN[III]
            PprofTo = Pprof[III[0],:]
            xCOpriorTo = xCOprior[III]
            PWFTo = PWF[III[0],:]
            error_superobsNTo = error_superobsN[III]
            AKTo = AK[III[0],:]
            prior_COTo = prior_CO[III[0],:]
            xCOTo = xCO[III]
            xCOaprioriTo = xCOapriori[III]


            III = np.where( np.logical_and( np.logical_and(lon_superobNTo>-179,lon_superobNTo<179) , np.logical_and(lat_superobNTo>-89,lat_superobNTo<89) ) )

            print np.shape(III)

            lon_superobN_out = lon_superobNTo[III] 
            lat_superobN_out = lat_superobNTo[III]
            hour_superobN_out = hour_superobNTo[III]
            Pprof_out = PprofTo[III[0],:]
            xCOprior_out = xCOpriorTo[III]
            PWF_out = PWFTo[III[0],:]
            error_superobsN_out = error_superobsNTo[III]
            AK_out = AKTo[III[0],:]
            prior_CO_out = prior_COTo[III[0],:]
            xCO_out = xCOTo[III]
            xCOapriori_out = xCOaprioriTo[III]

            IIIoo = np.where(AK_out == -999)
            AK_out[IIIoo] = 0

            print np.max(lon_superobN_out)


            #### FIX ERROR
            # error = ln(1+delta/x)
            error_fixed_t1 = np.log( 1 + np.exp(error_superobsN_out)/np.exp(xCO_out) )
            error_fixed_t2 = np.exp(error_superobsN_out)/np.exp(xCO_out)

            cost_t1 = ( xCO_out - xCOapriori_out ) / error_fixed_t1
            cost_t2 = ( xCO_out - xCOapriori_out ) / error_fixed_t2


            IIIt = np.where(np.isfinite(cost_t2) == 0)
            if np.size(IIIt[0]) > 0:
                stophere


            file_out='/nobackup/bbyrne1/MOPITT_v9_xCO_LogF/'+file_in[-13:-3]+'.nc'
            print file_out
            dataset = Dataset(file_out,'w')
            nSamples = dataset.createDimension('nSamples',np.size(lon_superobN_out))
            maxLevels = dataset.createDimension('maxLevels',10)
            #
            longitudes = dataset.createVariable('longitude', np.float64, ('nSamples',))
            longitudes[:]=lon_superobN_out
            #
            latitudes = dataset.createVariable('latitude', np.float64, ('nSamples',))
            latitudes[:]=lat_superobN_out
            #
            modes = dataset.createVariable('mode', np.float64, ('nSamples',))
            modes[:]=hour_superobN_out*0
            #
            times = dataset.createVariable('time', np.float64, ('nSamples',))
            times[:]=hour_superobN_out
            #
            pressures = dataset.createVariable('pressure', np.float64, ('nSamples','maxLevels'))
            pressures[:,:]=Pprof_out
            #
            xCOs = dataset.createVariable('xCO', np.float64, ('nSamples',))
            xCOs[:]=xCO_out
            #
            xCOaprioris = dataset.createVariable('xCO-apriori', np.float64, ('nSamples',))
            xCOaprioris[:]=xCOapriori_out
            #
            xCOpressureWeights = dataset.createVariable('xCO-pressureWeight', np.float64, ('nSamples','maxLevels'))
            xCOpressureWeights[:,:]=PWF_out
            #
            xCOuncertaintys = dataset.createVariable('xCO-uncertainty', np.float64, ('nSamples',))
            xCOuncertaintys[:] = np.abs(error_fixed_t2) #error_superobsN_out
            #
            xCOaveragingKernels = dataset.createVariable('xCO-averagingKernel', np.float64, ('nSamples','maxLevels'))
            xCOaveragingKernels[:,:]=AK_out
            #
            COaprioris = dataset.createVariable('CO-apriori', np.float64, ('nSamples','maxLevels'))
            COaprioris[:,:]=prior_CO_out
            #
            dataset.AKSpace = "LogVMR"
            dataset.close()
            
            
            
            




