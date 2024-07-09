import os
from datetime import datetime, timedelta
import shutil

def create_folder_for_date(base_dir, date):
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    
    dir_path = os.path.join(base_dir, year, month, day)
    os.makedirs(dir_path, exist_ok=True)
    print(f"Folder created: {dir_path}")

def create_folders_for_range(base_dir, start_date, end_date):
    current_date = start_date
    while current_date <= end_date:
        create_folder_for_date(base_dir, current_date)
        current_date += timedelta(days=1)
    
    print("All folders created successfully.")

def move_files_to_subdirectories(base_dir):
    """
    Moves files in the specified base directory into subdirectories based on the date in their filenames.
    
    Parameters:
    base_dir (str): The path to the base directory containing the files.
    
    Raises:
    FileNotFoundError: If the base directory does not exist.
    OSError: If the target directory cannot be created.
    """
    if not os.path.exists(base_dir):
        raise FileNotFoundError(f"The directory '{base_dir}' does not exist.")
    
    # Iterate over all files in the base directory
    for filename in os.listdir(base_dir):
        if filename.endswith(".nc"):
            # Extract the date parts from the filename
            YYYY = filename[20:24]
            MM = filename[24:26]
            DD = filename[26:28]
            
            # Create the target directory path
            target_dir = os.path.join(base_dir, YYYY, MM, DD)
            
            # Try to create the directories if they don't exist
            try:
                os.makedirs(target_dir, exist_ok=False)
            except OSError as e:
                if not os.path.exists(target_dir):
                    raise OSError(f"Could not create directory '{target_dir}': {e}")
            
            # Move the file to the target directory
            src_path = os.path.join(base_dir, filename)
            dest_path = os.path.join(target_dir, filename)
            shutil.move(src_path, dest_path)
            
            print(f"Moved {filename} to {dest_path}")

    print("File reorganization complete.")



    
# Example usage
base_dir = '/nobackupp19/bbyrne1/Test_TROPOMI_download/PRODUCT'
move_files_to_subdirectories(base_dir)

## Define base directory and date range
#start_date = datetime(2018, 1, 1)
#end_date = datetime(2024, 12, 31)

## Create folders for the specified date range
#create_folders_for_range(base_dir, start_date, end_date)


