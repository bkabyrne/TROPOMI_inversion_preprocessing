import os
import datetime
from subprocess import call
import calc_TROPOMI_mole_fractions_2x25

def directory_exists(path):
    """
    Check if a directory exists at the given path.
    
    :param path: The path to the directory.
    :return: True if the directory exists, False otherwise.
    """
    return os.path.isdir(path)

def create_directory(path):
    """
    Create a directory at the given path if it doesn't exist.
    
    :param path: The path to the directory.
    """
    if not directory_exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")

def directory_is_empty(path):
    """
    Check if a directory is empty.
    
    :param path: The path to the directory.
    :return: True if the directory is empty, False otherwise.
    """
    return not any(os.scandir(path))

def date_to_path(path_in,date):
    """
    Convert a date to the corresponding directory path.
    
    :param date: The date object to convert.
    :return: The path string in the format /path_in/YYYY/MM/DD.
    """
    return os.path.join(path_in, date.strftime("%Y/%m/%d"))

def run_command(command):
    """
    Run an external command using the shell.
    
    :param command: The command string to run.
    :return: The exit status of the command.
    """
    result = call(command, shell=True)
    if result != 0:
        print(f"Command failed: {command}")
    return result



# Define the base path for the directories and the start date
path_raw_TROPOMI='/nobackupp19/bbyrne1/Test_TROPOMI_download/PRODUCT/'
path_processed_TROPOMI='/nobackup/bbyrne1/TROPOMI_XCO_2x25/'
start_date = datetime.date(2023, 1, 1)
end_date = datetime.date.today() - datetime.timedelta(weeks=1)

# Load the CDO module
run_command('module load cdo')

# Loop through each day from start_date to end_date
date = start_date
while date <= end_date:

    # File path for TROPOMI super-obs on given date
    file_out = date_to_path(path_processed_TROPOMI,date)+'.nc'

    # If super-obs don't exist
    if not os.path.isfile(file_out):

        # Path to raw TROPOMI data
        dir_path = date_to_path(path_raw_TROPOMI,date)
    
        # Create the directory if it doesn't exist
        create_directory(dir_path)

        # If data hasn't been downloaded
        if directory_is_empty(dir_path):
            
            # Directory exists but is empty, run the necessary commands
            print(f"Directory {dir_path} is empty. Running commands.")
            
            # Run download_tropomi_data from download_TROPOMI.sh
            download_command = f"bash download_TROPOMI.sh {date.year} {date.month:02d} {date.day:02d}"
            run_command(download_command)

        calc_TROPOMI_mole_fractions_2x25.make_TROPOMI_obs(path_raw_TROPOMI, path_processed_TROPOMI, date, '_OFFL_')

    else:
        print(f" -- {file_out} exists")

    # Move to the next day
    date += datetime.timedelta(days=1)
