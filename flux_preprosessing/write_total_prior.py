from datetime import datetime, timedelta
from netCDF4 import Dataset
import numpy as np
import os

'''

This program calculated and writes the total prior CO flux, which is the combined Fossil Fuel, biomass burning and biogenic CO emissions


'''


def get_flux_data(filepath, variable_name='CO_Flux'):
    """
    Reads the flux data from a netCDF file and ensures it has the shape (8, lat, lon).

    Parameters:
    filepath (str): The path to the netCDF file.
    variable_name (str): The name of the variable to read from the netCDF file.

    Returns:
    numpy.ndarray: The flux data with shape (8, lat, lon).
    """
    with Dataset(filepath, mode='r') as dataset:
        flux_data = dataset.variables[variable_name][:]
        if flux_data.ndim == 2:  # shape (lat, lon)
            flux_data = np.repeat(flux_data[np.newaxis, :, :], 8, axis=0)
    return flux_data

def calc_total_flux(model, date):
    """
    Calculates the total flux by summing biomass burning, fossil fuel, and biogenic fluxes.

    Parameters:
    model (str): The biomass burning model to use ('GFED', 'GFAS', or 'QFED').
    date (datetime): The date for which to calculate the flux.

    Returns:
    numpy.ndarray: The total flux data with shape (8, lat, lon).

    Raises:
    ValueError: If the date is not a datetime object or if the model is invalid.
    """
    if not isinstance(date, datetime):
        raise ValueError("The 'date' parameter must be a datetime object.")
    
    year, month, day = date.year, date.month, date.day

    model_paths = {
        'GFED': f'/nobackup/bbyrne1/Flux_2x25_CO/BiomassBurn/GFED41s/{year:04}/{month:02}/{day:02}.nc',
        'GFAS': f'/nobackup/bbyrne1/Flux_2x25_CO/BiomassBurn/QFED/{year:04}/{month:02}/{day:02}.nc',
        'QFED': f'/nobackup/bbyrne1/Flux_2x25_CO/BiomassBurn/QFED/{year:04}/{month:02}/{day:02}.nc'
    }

    if model not in model_paths:
        raise ValueError("Invalid model specified. Choose 'GFED', 'GFAS', or 'QFED'.")

    CO_flux_BB = get_flux_data(model_paths[model])
    CO_flux_FF = get_flux_data(f'/nobackup/bbyrne1/Flux_2x25_CO/FossilFuel/CEDSdaily/{year:04}/{month:02}/{day:02}.nc')
    CO_flux_Biogenic = get_flux_data(f'/nobackup/bbyrne1/Flux_2x25_CO/Biogenic_units/{year:04}/{month:02}/{day:02}.nc')
    
    return CO_flux_BB + CO_flux_FF + CO_flux_Biogenic

def write_total_flux(total_flux, date, model):
    """
    Writes the total flux data to a netCDF file.

    Parameters:
    total_flux (numpy.ndarray): The total flux data with shape (8, lat, lon).
    date (datetime): The date for which the flux data is written.
    model (str): The biomass burning model used ('GFED', 'GFAS', 'QFED').

    Raises:
    ValueError: If the date is not a datetime object or if the model is invalid.
    """
    if not isinstance(date, datetime):
        raise ValueError("The 'date' parameter must be a datetime object.")

    year, month, day = date.year, date.month, date.day

    model_subdirs = {
        'GFED': 'FF_BB_Bio',
        'GFAS': 'FF_GFAS_Bio',
        'QFED': 'FF_QFED_Bio'
    }

    if model not in model_subdirs:
        raise ValueError("Invalid model specified. Choose 'GFED', 'GFAS', or 'QFED'.")

    nc_out_dir = f'/nobackup/bbyrne1/Flux_2x25_CO/Combined/{model_subdirs[model]}/{year:04}/{month:02}'
    os.makedirs(nc_out_dir, exist_ok=True)
    nc_out = f'{nc_out_dir}/{day:02}.nc'

    print(f"Writing to {nc_out}")
    
    with Dataset(nc_out, 'w') as dataset:
        dataset.createDimension('time', 8)
        dataset.createDimension('lat', total_flux.shape[1])
        dataset.createDimension('lon', total_flux.shape[2])
        
        postBBs = dataset.createVariable('CO_Flux', np.float64, ('time', 'lat', 'lon'))
        postBBs[:, :, :] = total_flux
        postBBs.units = 'kgC/km2/s'

def process_year(model, year):
    """
    Processes the total flux for each day of a given year and writes the output to netCDF files.

    Parameters:
    model (str): The biomass burning model to use ('GFED', 'GFAS', or 'QFED').
    year (int): The year to process.
    """
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    current_date = start_date

    while current_date <= end_date:
        try:
            total_flux = calc_total_flux(model, current_date)
            write_total_flux(total_flux, current_date, model)
        except Exception as e:
            print(f"Error processing {current_date}: {e}")
        current_date += timedelta(days=1)

if __name__ == "__main__":
    for model in ['GFAS','GFED','QFED']:
        for year in range(2018,2024):
            process_year(model, year)
