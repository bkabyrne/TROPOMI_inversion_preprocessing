import os
import datetime
import subprocess
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

def date_to_path(path_in, date):
    """
    Convert a date to the corresponding directory path.
    
    :param date: The date object to convert.
    :return: The path string in the format /path_in/YYYY/MM/DD.
    """
    return os.path.join(path_in, date.strftime("%Y/%m/%d"))

def run_command(command, use_shell=False):
    """
    Run an external command safely.
    
    :param command: List of command arguments or a string (if use_shell=True).
    :param use_shell: Whether to use the shell to run the command.
    :return: The exit status of the command.
    """
    try:
        if use_shell:
            result = subprocess.run(command, shell=True, check=True)
        else:
            result = subprocess.run(command, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")
        return e.returncode


# Define the base path for the directories and the start date
path_raw_TROPOMI = '/nobackupp19/bbyrne1/Test_TROPOMI_download/PRODUCT/'
path_processed_TROPOMI = '/nobackup/bbyrne1/TROPOMI_XCO_2x25/'
start_date = datetime.date(2024, 1, 1)
end_date = datetime.date.today() - datetime.timedelta(weeks=1)

# Load the CDO module safely (this requires shell=True to work)
run_command('module load cdo', use_shell=True)

# Loop through each day from start_date to end_date
date = start_date
while date <= end_date:

    # File path for TROPOMI super-obs on the given date
    file_out = date_to_path(path_processed_TROPOMI, date) + '.nc'

    # If super-obs don't exist
    if not os.path.isfile(file_out):

        # Path to raw TROPOMI data
        dir_path = date_to_path(path_raw_TROPOMI, date)
    
        # Create the directory if it doesn't exist
        create_directory(dir_path)

        # If data hasn't been downloaded
        if directory_is_empty(dir_path):
            
            # Directory exists but is empty, run the necessary commands
            print(f"Directory {dir_path} is empty. Running commands.")
            
            # Ensure year, month, and day are sanitized as they are integers.
            year = f"{date.year}"
            month = f"{date.month:02d}"
            day = f"{date.day:02d}"
            
            # Run download_tropomi_data from download_TROPOMI.sh safely using a list
            download_command = [
                'bash', 
                '/nobackupp19/bbyrne1/TROPOMI_inversion_preprocessing/download_TROPOMI.sh', 
                year, month, day
            ]
            run_command(download_command)

        calc_TROPOMI_mole_fractions_2x25.make_TROPOMI_obs(path_raw_TROPOMI, path_processed_TROPOMI, date, '_OFFL_')

    else:
        print(f" -- {file_out} exists")

    # Move to the next day
    date += datetime.timedelta(days=1)
