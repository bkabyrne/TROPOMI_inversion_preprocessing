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

'''                                         
This program writes posterior fire emissions for MOPITT and TROPOMI inversions (7-day opt)
                                            
'''

def calc_fluxes(SF,nc_CO2_fire,nc_CO_fire,prior_model,year, lat, lon):
    #
    # ================================
    #
    # This function reads in posterior scale factors and daily prior fluxes
    # during Apr-Sep then calculates timeseries of prior and posterior fluxes
    #
    # inputs:
    #  - SF: array of posterior scale factors
    #  - nc_CO2_fire: path to prior CO fire flux directory
    #
    # outputs:
    #  - CO2_Prior_flux: timeseries of prior CO2 fire fluxes (time,lat,lon) 
    #  - CO2_Posterior_flux: timeseries of posterior CO2 fire fluxes (time,lat,lon) 
    #
    # ================================

    # ----------
    if np.logical_and(year % 4 == 0, year % 100 != 0):
        days_in_year = 366
        days_in_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        Jan1 = 0#31+29+31
    else:
        days_in_year = 365
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        Jan1 = 0#31+28+31
    # ----------

    # =============================================    
    day_of_year_inv = np.arange(days_in_year)+1
    year_inv = np.zeros(days_in_year)+(year)
    #                                                  
    CURRENT_GROUP = np.zeros(days_in_year)
    B_Y = year-1
    for i in range(np.size(day_of_year_inv)):
        if (day_of_year_inv[i] < 360):
            CURRENT_GROUP[i] = (np.floor((day_of_year_inv[i]-1.)/7.0)-38) + (year_inv[i]-B_Y)*52 - 1. # python uses zero indexing
        else:
            CURRENT_GROUP[i] = 52 + (year_inv[i]-B_Y)*52 - 40 + 1 - 1. # python uses zero indexing           
    # =============================================
    
    
    #
    # Track day of year
    days_in_month_cum = np.zeros(13)
    for i in range(13):
        days_in_month_cum[i] = np.sum(days_in_month[0:i])
    #
    # Apply scale factors to fluxes
    CO_Prior_flux = np.zeros((days_in_year,np.size(lat),np.size(lon)))
    CO_Posterior_flux = np.zeros((days_in_year,np.size(lat),np.size(lon)))
    CO2_Prior_flux = np.zeros((days_in_year,np.size(lat),np.size(lon)))
    CO2_Posterior_flux = np.zeros((days_in_year,np.size(lat),np.size(lon)))
    for nn in range(days_in_year):
        #
        SF_index = int(CURRENT_GROUP[nn+Jan1])
        #
        month = np.argmax( nn < days_in_month_cum)
        day = int(nn-days_in_month_cum[month-1])
        #
        print(str(month).zfill(2)+'/'+str(day+1).zfill(2))
        file_in = nc_CO_fire+str(month).zfill(2)+'/'+str(day+1).zfill(2)+'.nc'
        f=Dataset(file_in,mode='r')
        if prior_model == 'GFED':
            CO_Prior_flux[nn+Jan1,:,:] = np.mean(f.variables['CO_Flux'][:],0)  * (60.*60.*24.)/1000. # kgC/km2/s -> gC/m2/d
        else:
            CO_Prior_flux[nn+Jan1,:,:] = f.variables['CO_Flux'][:]  * (60.*60.*24.)/1000. # kgC/km2/s -> gC/m2/d
        CO_Posterior_flux[nn+Jan1,:,:] = CO_Prior_flux[nn+Jan1,:,:] * SF[SF_index,:,:]
        #
        file_in = nc_CO2_fire+str(month).zfill(2)+'/'+str(day+1).zfill(2)+'.nc'
        f=Dataset(file_in,mode='r')
        if prior_model == 'GFED':
            CO2_Prior_flux[nn+Jan1,:,:] = np.mean(f.variables['CO2_Flux'][:],0)  * (60.*60.*24.)/1000. # kgC/km2/s -> gC/m2/d
        else:
            CO2_Prior_flux[nn+Jan1,:,:] = f.variables['CO2_Flux'][:]  * (60.*60.*24.)/1000. # kgC/km2/s -> gC/m2/d 
        CO2_Posterior_flux[nn+Jan1,:,:] = CO2_Prior_flux[nn+Jan1,:,:] * SF[SF_index,:,:]
    #
    return CO2_Prior_flux, CO2_Posterior_flux, CO2_Prior_flux, CO2_Posterior_flux


def calculate_2x25_grid_area():
    #
    # =============================
    # Returns grid area (lat,lon) in m2
    # =============================
    #
    grid_area_2x25 = np.array([2.70084e+08,  2.16024e+09,  4.31787e+09,  6.47023e+09,  8.61471e+09,
                               1.07487e+10,  1.28696e+10,  1.49748e+10,  1.70617e+10,  1.91279e+10,
                               2.11708e+10,  2.31879e+10,  2.51767e+10,  2.71348e+10,  2.90599e+10,
                               3.09496e+10,  3.28016e+10,  3.46136e+10,  3.63835e+10,  3.81090e+10,
                               3.97881e+10,  4.14187e+10,  4.29988e+10,  4.45266e+10,  4.60001e+10,
                               4.74175e+10,  4.87772e+10,  5.00775e+10,  5.13168e+10,  5.24935e+10,
                               5.36063e+10,  5.46538e+10,  5.56346e+10,  5.65477e+10,  5.73920e+10,
                               5.81662e+10,  5.88696e+10,  5.95014e+10,  6.00606e+10,  6.05466e+10,
                               6.09588e+10,  6.12968e+10,  6.15601e+10,  6.17484e+10,  6.18615e+10,
                               6.18992e+10,  6.18615e+10,  6.17484e+10,  6.15601e+10,  6.12968e+10,
                               6.09588e+10,  6.05466e+10,  6.00606e+10,  5.95014e+10,  5.88696e+10,
                               5.81662e+10,  5.73920e+10,  5.65477e+10,  5.56346e+10,  5.46538e+10,
                               5.36063e+10,  5.24935e+10,  5.13168e+10,  5.00775e+10,  4.87772e+10,
                               4.74175e+10,  4.60001e+10,  4.45266e+10,  4.29988e+10,  4.14187e+10,
                               3.97881e+10,  3.81090e+10,  3.63835e+10,  3.46136e+10,  3.28016e+10,
                               3.09496e+10,  2.90599e+10,  2.71348e+10,  2.51767e+10,  2.31879e+10,
                               2.11708e+10,  1.91279e+10,  1.70617e+10,  1.49748e+10,  1.28696e+10,
                               1.07487e+10,  8.61471e+09,  6.47023e+09,  4.31787e+09,  2.16024e+09,
                               2.70084e+08])
    #
    grid_area_2x25_arr = np.zeros((91,144))
    for ii in range(144):
        grid_area_2x25_arr[:,ii] = grid_area_2x25
    #
    return grid_area_2x25_arr


def write_dataset(nc_out, CO2_Flux_prior, CO2_Flux_post, CO_Flux_prior, CO_Flux_post, lat, lon, year):
    #
    # =============================
    # Write prior/posterior fluxes to netcdf
    # =============================
    #
    if np.logical_and(year % 4 == 0, year % 100 != 0):
        days_in_year = 366
    else:
        days_in_year = 365
    # Read grid to write out area (m2)
    grid_area_2x25 = calculate_2x25_grid_area()
    #
    # Write out data
    dataset = Dataset(nc_out,'w')
    print(nc_out)
    times = dataset.createDimension('time',days_in_year)
    lats = dataset.createDimension('lat',91)
    lons = dataset.createDimension('lon',144)
    gridareas = dataset.createVariable('grid_area', np.float64, ('lat','lon'))
    gridareas[:,:] = grid_area_2x25
    gridareas.units = 'm2'
    latss = dataset.createVariable('latitude', np.float64, ('lat',))
    latss[:] = lat
    lonss = dataset.createVariable('longitude', np.float64, ('lon',))
    lonss[:] = lon
    CO2_priors = dataset.createVariable('CO2_prior', np.float64, ('time','lat','lon'))
    CO2_priors[:,:,:] = CO2_Flux_prior
    CO2_priors.units = 'gC/m2/day'
    CO2_posts = dataset.createVariable('CO2_post', np.float64, ('time','lat','lon'))
    CO2_posts[:,:,:] = CO2_Flux_post
    CO2_posts.units = 'gC/m2/day'
    CO_priors = dataset.createVariable('CO_prior', np.float64, ('time','lat','lon'))
    CO_priors[:,:,:] = CO_Flux_prior
    CO_priors.units = 'gC/m2/day'
    CO_posts = dataset.createVariable('CO_post', np.float64, ('time','lat','lon'))
    CO_posts[:,:,:] = CO_Flux_post
    CO_posts.units = 'gC/m2/day'
    dataset.close()


def run_optimized_fluxes(yyyy, sfnum, obsType, dir_out):
    """
    Run optimized flux calculations based on observational data type and write the results to a netCDF file.

    Parameters:
    -----------
    yyyy : int
        The year for which the fluxes are being calculated.
    sfnum : int
        Iteration number of optimization
    obsType : str
        Type of observational data, either 'MOPITT' or 'TROPOMI'.
    dir_out : str
        Directory where the output netCDF file will be saved.

    Raises:
    -------
    ValueError
        If obsType is not one of 'MOPITT' or 'TROPOMI'.
    FileNotFoundError
        If the specified netCDF file for scaling factors or flux directories does not exist.
    """
    
    # Ensure the observation type is valid
    if obsType == 'MOPITT':
        nc_file = f'/nobackup/bbyrne1/GHGF-CMS-7day-COinv-{str(yyyy).zfill(4)}/Run_COinvMOPITT2/GDT-EMS/EMS-sf-{str(sfnum).zfill(2)}.nc'
    elif obsType == 'TROPOMI':
        nc_file = f'/nobackup/bbyrne1/GHGF-CMS-7day-COinv-{str(yyyy).zfill(4)}/Run_COinv/GDT-EMS/EMS-sf-{str(sfnum).zfill(2)}.nc'
    else:
        raise ValueError(f"Invalid obsType '{obsType}'. Must be either 'MOPITT' or 'TROPOMI'.")

    # Check if the netCDF file exists
    if not os.path.isfile(nc_file):
        raise FileNotFoundError(f"NetCDF file not found: {nc_file}")

    try:
        # Open the netCDF file
        with Dataset(nc_file, mode='r') as f:
            lon = f.variables['lon'][:]
            lat = f.variables['lat'][:]
            SF = f.variables['EMS-01'][:]
    except Exception as e:
        raise RuntimeError(f"Error reading NetCDF file: {nc_file}. Details: {e}")
    
    # Define directories of prior fluxes (GFED is the only prior model)
    ncdir_CO_fire = f'/nobackup/bbyrne1/Flux_2x25_CO/BiomassBurn/GFED41s/{str(yyyy).zfill(4)}/'
    ncdir_CO2_fire = f'/nobackup/bbyrne1/GFED41s_2x25/{str(yyyy).zfill(4)}/'

    # Check if prior flux directories exist
    if not os.path.isdir(ncdir_CO_fire):
        raise FileNotFoundError(f"CO fire flux directory not found: {ncdir_CO_fire}")
    if not os.path.isdir(ncdir_CO2_fire):
        raise FileNotFoundError(f"CO2 fire flux directory not found: {ncdir_CO2_fire}")

    try:
        # Read and calculate prior and posterior fluxes using GFED as the prior model
        CO2_Flux_prior, CO2_Flux_post, CO_Flux_prior, CO_Flux_post = calc_fluxes(SF, ncdir_CO2_fire, ncdir_CO_fire, 'GFED', yyyy, lat, lon)
    except Exception as e:
        raise RuntimeError(f"Error calculating fluxes. Details: {e}")
    
    # Prepare the output netCDF file path
    ncfile_out = os.path.join(dir_out, f'Inversion_{obsType}_{str(yyyy).zfill(4)}.nc')
    print(f"Writing output to {ncfile_out}")

    try:
        # Write the calculated fluxes to a new netCDF file
        write_dataset(ncfile_out, CO2_Flux_prior, CO2_Flux_post, CO_Flux_prior, CO_Flux_post, lat, lon, yyyy)
    except Exception as e:
        raise RuntimeError(f"Error writing dataset to {ncfile_out}. Details: {e}")
    

        
if __name__ == "__main__":

    dir_out = '/nobackupp19/bbyrne1/Global_CO_inversion_results/'

    run_optimized_fluxes(2019,39,'TROPOMI',dir_out)
    run_optimized_fluxes(2020,39,'TROPOMI',dir_out)
    run_optimized_fluxes(2022,32,'TROPOMI',dir_out)
    run_optimized_fluxes(2023,29,'TROPOMI',dir_out)
    
    run_optimized_fluxes(2010,39,'MOPITT',dir_out)
    run_optimized_fluxes(2011,39,'MOPITT',dir_out)
    run_optimized_fluxes(2012,39,'MOPITT',dir_out)
    run_optimized_fluxes(2013,39,'MOPITT',dir_out)
    run_optimized_fluxes(2014,39,'MOPITT',dir_out)
    run_optimized_fluxes(2015,39,'MOPITT',dir_out)
    run_optimized_fluxes(2016,39,'MOPITT',dir_out)
    run_optimized_fluxes(2017,39,'MOPITT',dir_out)
    run_optimized_fluxes(2018,38,'MOPITT',dir_out)
    run_optimized_fluxes(2019,38,'MOPITT',dir_out)
    run_optimized_fluxes(2020,38,'MOPITT',dir_out)
    run_optimized_fluxes(2021,38,'MOPITT',dir_out)
