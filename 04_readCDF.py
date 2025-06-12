from netCDF4 import Dataset

# # Open the NetCDF file
# nc_file = 'netcdfs\ITP_2.nc'
# dataset = Dataset(nc_file, mode='r')  # 'r' for read-only

# # Explore the contents
# print(dataset.variables.keys())  # List all variable names
# print(dataset.dimensions.keys())  # List dimensions

# # Access a specific variable
# # lat = dataset.variables['lat'][:]  # Replace with your variable name
# # print(lat.shape)

# depth = dataset.variables['depth'][:]  # Replace with your variable name
# print(depth.shape)


# # Close when done
# dataset.close()


import xarray as xr
import numpy as np

# Open dataset
ds = xr.open_dataset('netcdfs\ITP_2.nc')

# List all variables and check for NaNs or missing values
for var in ds.data_vars:
    data = ds[var]
    missing_count = data.isnull().sum().item()
    print(f"Variable '{var}' has {missing_count} missing values.")


# Open the NetCDF file
nc_file = 'netcdfs\ITP_2.nc'
dataset = Dataset(nc_file, mode='r')  # 'r' for read-only

# Explore the contents
print(dataset.variables.keys())  # List all variable names
print(dataset.dimensions.keys())  # List dimensions

# Access a specific variable
# lat = dataset.variables['lat'][:]  # Replace with your variable name
# print(lat.shape)

depth = dataset.variables['depth'][:]  # Replace with your variable name
print(depth.shape)


# Close when done
dataset.close()
