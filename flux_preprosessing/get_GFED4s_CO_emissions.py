from mpl_toolkits.basemap import Basemap, cm
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime
from netCDF4 import Dataset
import glob, os
import numpy as np
import h5py # if this creates an error please make sure you have the h5py library

months       = '01','02','03','04','05','06','07','08','09','10','11','12'
sources      = 'SAVA','BORF','TEMF','DEFO','PEAT','AGRI'

# in this example we will calculate annual CO emissions for the 14 GFED 
# basisregions over 1997-2014. Please adjust the code to calculate emissions
# for your own specie, region, and time period of interest. Please
# first download the GFED4.1s files and the GFED4_Emission_Factors.txt
# to your computer and adjust the directory where you placed them below

directory    = '/u/bbyrne1/fire_emissions_v4_R1_1293/data'


"""
Read in emission factors
"""
species = [] # names of the different gas and aerosol species
EFs     = np.zeros((41, 6)) # 41 species, 6 sources

k = 0
f = open(directory+'/GFED4_Emission_Factors.txt')
while 1:
    line = f.readline()
    if line == "":
        break
        
    if line[0] != '#':
        contents = line.split()
        species.append(contents[0])
        EFs[k,:] = contents[1:]
        k += 1
                
f.close()

# we are interested in CO for this example (4th row):
EF_CO = EFs[3,:]
##print EF_CO
#stophere
start_year = 2023
end_year   = 2024


lat025_temp=-1.*np.arange(720)*0.25+89.875
lon025=np.arange(1440)*0.25-179.875

lat025=lat025_temp[::-1]

"""
make table with summed DM emissions for each region, year, and source
"""
CO_table = np.zeros((15, end_year - start_year + 1)) # region, year

CO_emiss_vecArr = np.zeros((700,720,1440))
for year in range(start_year, end_year+1):

    day_per_month = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    CO_emissions2 = np.zeros((365,8,720, 1440))    
    diurnal2 = np.zeros((365,8,720, 1440))    
    if np.logical_or(np.logical_or(np.logical_or(np.logical_or(np.logical_or(np.logical_or(year == 2000, year == 2004), year == 2008), year == 2012), year == 2016), year == 2020), year == 2024):
        day_per_month = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
        CO_emissions2 = np.zeros((366,8,720, 1440))    
        diurnal2 = np.zeros((366,8,720, 1440))    

    print('/u/bbyrne1/fire_emissions_v4_R1_1293/data/GFED4.1s_2023_beta_20240607.hdf5')
    string = directory+'/GFED4.1s_'+str(year)+'_beta_20240607.hdf5'
    print(string)
    f = h5py.File(string, 'r')
    
    #print year
    
    #if year == start_year: # these are time invariable    
    #    basis_regions = f['/ancill/basis_regions'][:]
    #    grid_area     = f['/ancill/grid_cell_area'][:]
    n=-1

    CO_emissions1 = np.zeros((12,720, 1440))
    for month in range(12):
        # read in DM emissions
        string = '/emissions/'+months[month]+'/DM'
        DM_emissions = f[string][:]
        #
        daily_fraction=np.zeros((day_per_month[month],720,1440))
        for i in range(day_per_month[month]):
            #print str(i+1)
            string2 = '/emissions/'+months[month]+'/daily_fraction/day_'+str(i+1)
            print(string2)
            daily_fraction[i,:,:] = f[string2][:] * day_per_month[month]
        #
        diurnal = np.zeros((8,720,1440))
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_0-3h'
        diurnal[0,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_3-6h'
        diurnal[1,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_6-9h'
        diurnal[2,:,:]= f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_9-12h'
        diurnal[3,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_12-15h'
        diurnal[4,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_15-18h'
        diurnal[5,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_18-21h'
        diurnal[6,:,:] = f[string1][:]*8.
        string1 = '/emissions/'+months[month]+'/diurnal_cycle/UTC_21-24h'
        diurnal[7,:,:] = f[string1][:]*8.
        #
        CO_emissions_all = np.zeros((720,1440))
        for source in range(6):
            # read in the fractional contribution of each source
            string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
            contribution = f[string][:]
            # calculate CO emissions as the product of DM emissions (kg DM per 
            # m2 per month), the fraction the specific source contributes to 
            # this (unitless), and the emission factor (g CO per kg DM burned
            #
            #  kgDM m-2 month-1 * gCO kgDM-1 = gCO m-2 month-1
            #
            # gCO/m2/month
            # (12.01gC/28.01gCO) * (1kgC/1000gC) * (1000m/1km)^2 * ( 1 month / ( 60 * 60 * 24 * xxx day/month ))  
            # kgC km-2 s-1
            CO_emissions_all += DM_emissions * contribution * EF_CO[source] * ( (12.01*1000.)/(28.01*60.*60.*24.*day_per_month[month]) )  
        CO_emissions1[month,:,:] = CO_emissions_all
        #
        for i in range(day_per_month[month]):
            n=n+1
            #print 'month: '+str(month+1)
            #print 'day of year: '+str(n+1)
            for j in range(8):
                CO_emissions2[n,j,:,:] = CO_emissions1[month,:,:]*daily_fraction[i,:,:]*diurnal[j,:,:]
                diurnal2[n,j,:,:] = diurnal[j,:,:]
            #
            CO_emissions_out = np.squeeze(CO_emissions2[n,:,::-1,:])
            #
            file_out='/nobackupp17/bbyrne1/GFED41s_025_CO/'+str(year).zfill(4)+'/'+str(month+1).zfill(2)+'/'+str(i+1).zfill(2)+'.nc'
            dataset = Dataset(file_out,'w')
            lon = dataset.createDimension('lon',1440)
            lat = dataset.createDimension('lat',720)
            time = dataset.createDimension('time',8)
            lons = dataset.createVariable('lon', np.float64, ('lon',))
            lons[:]=lon025
            lats = dataset.createVariable('lat', np.float64, ('lat',))
            lats[:]=lat025
            CO_fluxes = dataset.createVariable('CO_Flux',np.float64,('time','lat','lon'))
            CO_fluxes[:,:,:]=CO_emissions_out
            dataset.close()
            #
            #diurnal_out = np.squeeze(diurnal2[n,:,::-1,:])
            #
            #file_out='/nobackupp17/bbyrne1/GFED41s_025_diurnal_scale/'+str(year).zfill(4)+'/'+str(month+1).zfill(2)+'/'+str(i+1).zfill(2)+'.nc'
            #dataset = Dataset(file_out,'w')
            #lon = dataset.createDimension('lon',1440)
            #lat = dataset.createDimension('lat',720)
            #time = dataset.createDimension('time',8)
            #lons = dataset.createVariable('lon', np.float64, ('lon',))
            #lons[:]=lon025
            #lats = dataset.createVariable('lat', np.float64, ('lat',))
            #lats[:]=lat025
            #CO_fluxes = dataset.createVariable('diurnal',np.float64,('time','lat','lon'))
            #CO_fluxes[:,:,:]=diurnal_out
            #dataset.close()

            CO_emiss_vecArr[n,:,:]=np.mean(CO_emissions_out,0)
    
    ## fill table with total values for the globe (row 15) or basisregion (1-14)
    #for region in range(15):
    #    if region == 14:
    #        mask = np.ones((720, 1440))
    #    else:
    #        mask = basis_regions == (region + 1)            
    # 
    #    CO_table2[region, year-start_year] = np.sum(grid_area * mask * CO_emissions)
        
    #print year


# kgCO/m2/month
# (12.01gC/28.01gCO) * (1000m/1km)^2 * ( 1 month / ( 60 * 60 * 24 * xxx day/month )) 


#lon025[1295:1365]
#lat025[195:232]

#CO_emiss_vecArr[:,195:232,1295:1365]
CO_vec=np.mean(np.mean(CO_emiss_vecArr[:,195:232,1295:1365],1),1)
plt.figure(2)
plt.plot(CO_vec)
plt.xlim([366-31-30,366+31+29])
plt.show()
stophere

#CO_emissions = CO_emissions2[:,:,::-1,:]

#lat025_temp=-1.*np.arange(720)*0.25+89.875
#lon025=np.arange(1440)*0.25-179.875

#lat025=lat025_temp[::-1]


#365,8,720, 1440



fig = plt.figure(1)
fig.set_size_inches(20,5)
m = Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180.,urcrnrlat=90.,projection='cyl',lon_0=0)
ax1 = fig.add_axes([0.05, 0.05, .9, 0.9])
xx,yy=m(lon025,lat025)
m.pcolormesh(xx,yy,CO_emissions[0,0,:,:],cmap='hot_r',shading='flat',vmin=0.0, vmax=0.0001)#,vmin=387, vmax=398)
#m.pcolormesh(xx,yy,np.sum(daily_fraction,0),cmap='hot_r',shading='flat',vmin=0.0, vmax=2)#,vmin=387, vmax=398)
m.drawcoastlines()
plt.colorbar()
plt.show()
stophere
        
# convert to Tg CO 
#CO_table2 = CO_table2 / 1E12
##print CO_table

# please compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
