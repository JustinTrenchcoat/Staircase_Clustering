"""
Author: Mikhail Schee
Created: 2022-03-24

This script will load data from various Arctic Ocean profile datasets and
reformat them into netcdf files, one for each instrument.
This script expects the following file structure of the raw data files:

science_data_file_path/
    ITPs/
        itp1/
            itp1cormat/
                cor0001.mat
                cor0002.mat
                ...
            itp1final/
                itp1grd0001.dat
                itp1grd0002.dat
                ...
        itp2/
            ...
        ...

"""

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
# For searching and listing directories
import os
# For reading the ITP `cormat` files
import mat73
from scipy import io
# For reading netcdf files
import netCDF4 as netcdf
# For matching regular expressions
import re
# Import the Thermodynamic Equation of Seawater 2010 (TEOS-10) from GSW
# For converting from depth to pressure
import gsw

# For test plots
import matplotlib.pyplot as plt
plt.style.use('dark_background')

# Define attribute strings
science_data_file_path = '/Users/Grey/Documents/Research/Science_Data/'
me_string = 'Mikhail Schee, University of Toronto'
doi = 'No DOI yet'

# Sizes of netcdf dimensions (because netcdfs can't do ragged arrays)
max_ml_count = 200              # Maximum number of mixed layers stored

################################################################################
# Loading in data

def list_data_files(file_path):
    """
    Creates a list of all data files within a containing folder

    file_path       The path to the folder containing the data files
    """
    # Search the provided file path
    if os.path.isdir(file_path):
        data_files = os.listdir(file_path)
        # Remove .DS_Store directory from list
        if '.DS_Store' in data_files: data_files.remove('.DS_Store')
        # Restrict to just unique instruments (there shouldn't be any duplicates hopefully)
        data_files = np.unique(data_files)
    else:
        print(file_path, " is not a directory")
        data_files = None
    #
    return data_files

################################################################################

def make_all_ITP_netcdfs(science_data_file_path, format='cormat'):
    """
    Finds ITP data files for all instruments available and formats them into netcdfs

    science_data_file_path      string of the filepath where the data is stored
    format                      which version of the data files to use
                                    either 'cormat' or 'final'
    """
    # Declare file path
    main_dir = science_data_file_path+'ITPs/'
    # Search the provided file path
    if os.path.isdir(main_dir):
        ITP_dirs = os.listdir(main_dir)
        # Remove .DS_Store directory from list
        if '.DS_Store' in ITP_dirs: ITP_dirs.remove('.DS_Store')
    else:
        print(main_dir, " is not a directory")
        exit(0)
    # Read in data for all ITPs available
    for itp in ITP_dirs:
        # Get just the number for the itp
        itp_number = ''.join(filter(str.isdigit, itp))
        read_instrmt('ITP', itp_number, main_dir+itp+'/'+itp+format, 'netcdfs/ITP_'+itp_number.zfill(3)+'.nc')
    #

################################################################################

def read_instrmt(source, instrmt_name, instrmt_dir, out_file):
    """
    Reads in all the data for the specified instrument and formats it into a
    single netcdf

    source              string of the name of the data source (ex: 'AIDJEX', 'ITP')
    instrmt_name        string of the name of this instrmt
    instrmt_dir         string of a file path to this instrmt's directory
    out_file            string of the file path in which to save the netcdf
    """
    print('Reading',source,instrmt_name)
    # Select the corresponding read function for the provided data source
    if source == 'ITP':
        read_data_file  = read_ITP_data_file
        attribution    = 'The Ice-Tethered Profiler data were collected and made available by the Ice-Tethered Profiler Program (Toole et al., 2011; Krishfield et al., 2008) based at the Woods Hole Oceanographic Institution'
        source_url      = 'https://www.whoi.edu/itp'
        og_vert         = 'press'
        og_temp         = 'iT'
        og_salt         = 'SP'
    else:
        print(source,'is not a valid source')
        exit(0)
    # Make blank arrays for data
    list_of_entries         = []
    list_of_pf_nos          = []
    list_of_black_list      = []
    list_of_datetimes_start = []
    list_of_datetimes_end   = []
    list_of_lons            = []
    list_of_lats            = []
    list_of_regs            = []
    list_of_up_casts        = []
    list_of_press_arrs      = []
    list_of_depth_arrs      = []
    list_of_iT_arrs         = []
    list_of_CT_arrs         = []
    list_of_SP_arrs         = []
    list_of_SA_arrs         = []
    # Keep track of the maximum number of vertical measurements per profile
    max_vert_count = 0
    # Loop through the data files in this instrmt's directory
    data_files = list_data_files(instrmt_dir)
    if isinstance(data_files, type(None)):
        print('Did not find any files for',instrmt_name)
        exit(0)
    else:
        # Keep track of the entry number with i
        i = 0
        for file in data_files:
            if os.path.isfile(instrmt_dir+'/'+file):
                # Read in the data file for this profile
                out_dict = read_data_file(instrmt_dir, file, instrmt_name)
                if not isinstance(out_dict, type(None)):
                    # Append that data
                    list_of_entries.append(i)
                    list_of_pf_nos.append(out_dict['prof_no'])
                    list_of_black_list.append(out_dict['black_list'])
                    list_of_datetimes_start.append(out_dict['dt_start'])
                    list_of_datetimes_end.append(out_dict['dt_end'])
                    list_of_lons.append(out_dict['lon'])
                    list_of_lats.append(out_dict['lat'])
                    list_of_regs.append(out_dict['region'])
                    list_of_up_casts.append(out_dict['up_cast'])
                    list_of_press_arrs.append(out_dict['press'])
                    list_of_depth_arrs.append(out_dict['depth'])
                    list_of_iT_arrs.append(out_dict['iT'])
                    list_of_CT_arrs.append(out_dict['CT'])
                    list_of_SP_arrs.append(out_dict['SP'])
                    list_of_SA_arrs.append(out_dict['SA'])
                    # Check for a new maximum vertical dimension length
                    max_vert_count = max(max_vert_count, len(out_dict['depth']))
                    # Increase entry number
                    i += 1
                #
            #
        #
    #
    # Make sure the vertical data are all the same length lists
    for i in range(len(list_of_depth_arrs)):
        list_of_press_arrs[i] = list(list_of_press_arrs[i]) +[None]*(max_vert_count - len(list_of_press_arrs[i]))
        list_of_depth_arrs[i] = list(list_of_depth_arrs[i]) +[None]*(max_vert_count - len(list_of_depth_arrs[i]))
        list_of_iT_arrs[i] = list(list_of_iT_arrs[i]) + [None]*(max_vert_count - len(list_of_iT_arrs[i]))
        list_of_SP_arrs[i] = list(list_of_SP_arrs[i]) + [None]*(max_vert_count - len(list_of_SP_arrs[i]))
        list_of_CT_arrs[i] = list(list_of_CT_arrs[i]) + [None]*(max_vert_count - len(list_of_CT_arrs[i]))
        list_of_SA_arrs[i] = list(list_of_SA_arrs[i]) + [None]*(max_vert_count - len(list_of_SA_arrs[i]))
    #
    # Make a blank array for each dimension
    Time_blank = [None]*len(list_of_datetimes_start)
    Vertical_blank = [[None]*max_vert_count]*len(list_of_datetimes_start)
    Layer_blank = [[None]*max_ml_count]*len(list_of_datetimes_start)

    # Define variables with data and attributes
    nc_vars = {
                'entry':(
                        ['Time'],
                        list_of_entries,
                        {
                            'units':'N/A',
                            'label':'Entry (index)',
                            'long_name':'Sequential entry number',
                            'dtype':'int64'
                        }
                ),
                'prof_no':(
                        ['Time'],
                        list_of_pf_nos,
                        {
                            'units':'N/A',
                            'label':'Profile number',
                            'long_name':'Profile number',
                            'dtype':'int64'
                        }
                ),
                'BL_yn':(
                        ['Time'],
                        list_of_black_list,
                        {
                            'units':'True/False',
                            'label':'On Black List',
                            'long_name':'Profile is on Black List, True(1) or False(0)'
                        }
                ),
                'dt_start':(
                        ['Time'],
                        list_of_datetimes_start,
                        {
                            'units':'YYYY-MM-DD HH:MM:SS',
                            'label':'Datetime of profile start',
                            'long_name':'Datetime when data for this entry started'
                        }
                ),
                'dt_end':(
                        ['Time'],
                        list_of_datetimes_end,
                        {
                            'units':'YYYY-MM-DD HH:MM:SS',
                            'label':'Datetime of profile end',
                            'long_name':'Datetime when data for this entry ended'
                        }
                ),
                'lon':(
                        ['Time'],
                        list_of_lons,
                        {
                            'units':'Degrees',
                            'label':'Longitude ($^\circ$E+)',
                            'long_name':'Longitude'
                        }
                ),
                'lat':(
                        ['Time'],
                        list_of_lats,
                        {
                            'units':'Degrees',
                            'label':'Latitude ($^\circ$N)',
                            'long_name':'Latitude'
                        }
                ),
                'region':(
                        ['Time'],
                        list_of_regs,
                        {
                            'units':'N/A',
                            'label':'Geographical region',
                            'long_name':'Geographical region in which this entry was measured'
                        }
                ),
                'up_cast':(
                        ['Time'],
                        list_of_up_casts,
                        {
                            'units':'N/A',
                            'label':'Upcast or Downcast',
                            'long_name':'Profile taken with instrument going up(True) or down(False)'
                        }
                ),
                'R_rho':(
                        ['Time'],
                        Time_blank,
                        {
                            'units':'N/A',
                            'label':'$R_\\rho$ Density Ratio',
                            'long_name':'Density ratio (defined how?)'
                        }
                ),
                'press':(
                        ['Time','Vertical'],
                        list_of_press_arrs,
                        {
                            'units':'dbar',
                            'label':'Pressure (dbar)',
                            'long_name':'Pressure'
                        }
                ),
                'depth':(
                        ['Time','Vertical'],
                        list_of_depth_arrs,
                        {
                            'units':'m',
                            'label':'Depth (m)',
                            'long_name':'Depth'
                        }
                ),
                'iT':(
                        ['Time','Vertical'],
                        list_of_iT_arrs,
                        {
                            'units':'degrees Celcius',
                            'label':'in-situ Temperature ($^\circ$C)',
                            'long_name':'in-situ Temperature'
                        }
                ),
                'CT':(
                        ['Time','Vertical'],
                        list_of_CT_arrs,
                        {
                            'units':'degrees Celcius',
                            'label':'Conservative Temperature ($^\circ$C)',
                            'long_name':'Conservative Temperature'
                        }
                ),
                'SP':(
                        ['Time','Vertical'],
                        list_of_SP_arrs,
                        {
                            'units':'g/kg',
                            'label':'Practical Salinity (g/kg)',
                            'long_name':'Practical Salinity'
                        }
                ),
                'SA':(
                        ['Time','Vertical'],
                        list_of_SA_arrs,
                        {
                            'units':'g/kg',
                            'label':'Absolute Salinity (g/kg)',
                            'long_name':'Absolute Salinity'
                        }
                ),
                'sigma':(
                        ['Time','Vertical'],
                        gsw.sigma1(list_of_SA_arrs, list_of_CT_arrs),
                        {
                            'units':'kg/m^3',
                            'label':'Density anomaly (kg/m$^3$)',
                            'long_name':'Density anomaly referenced to 1000 dbar'
                        }
                ),
                'alpha':(
                        ['Time','Vertical'],
                        gsw.alpha(list_of_SA_arrs, list_of_CT_arrs, list_of_press_arrs),
                        {
                            'units':'1/(degrees Celcius)',
                            'label':'$\\alpha$ (1/$^\circ$C)',
                            'long_name':'Thermal expansion coefficient alpha'
                        }
                ),
                'beta':(
                        ['Time','Vertical'],
                        gsw.beta(list_of_SA_arrs, list_of_CT_arrs, list_of_press_arrs),
                        {
                            'units':'1/(g/kg)',
                            'label':'$\\beta$ (kg/g)',
                            'long_name':'Saline contraction coefficient beta'
                        }
                ),
                'ss_mask':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'N/A',
                            'label':'Sub-sample mask',
                            'long_name':'A mask to apply the sub-sample scheme'
                        }
                ),
                'ma_iT':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'degrees Celcius',
                            'label':'Moving average in-situ temperature ($^\circ$C)',
                            'long_name':'Moving average in-situ temperature'
                        }
                ),
                'ma_CT':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'degrees Celcius',
                            'label':'Moving average conservative temperature ($^\circ$C)',
                            'long_name':'Moving average conservative temperature'
                        }
                ),
                'ma_SP':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'g/kg',
                            'label':'Moving average practical salinity (g/kg)',
                            'long_name':'Moving average practical salinity'
                        }
                ),
                'ma_SA':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'g/kg',
                            'label':'Moving average absolute salinity (g/kg)',
                            'long_name':'Moving average absolute salinity'
                        }
                ),
                'ma_sigma':(
                        ['Time','Vertical'],
                        Vertical_blank,
                        {
                            'units':'kg/m^3',
                            'label':'Moving average density anomaly (kg/m$^3$)',
                            'long_name':'Moving average density anomaly'
                        }
                )
    }
    #
    # Define coordinates
    nc_coords = {
               'Time':(
                        ['Time'],
                        list_of_datetimes_start,
                        {
                            'units':'YYYY-MM-DD HH:MM:SS',
                            'long_name':'Datetime of the start of the measurement',
                            'dtype':'datetime64[ns]'
                        }
                ),
               'Vertical':(
                        ['Vertical'],
                        np.arange(max_vert_count),
                        {
                            'units':'N/A',
                            'long_name':'An index of vertical position (depth or pressure)',
                            'dtype':'int64'
                        }
                ),
               'Layer':(
                        ['Layer'],
                        np.arange(max_ml_count),
                        {
                            'units':'N/A',
                            'long_name':'An index for detected mixed/gradient layers',
                            'dtype':'int64'
                        }
                )
    }
    #
    # Define global attributes
    nc_attrs = {        # Note: can't store datetime objects in netcdfs
                'Title':'Arctic Ocean Profile Data',
                'Creation date':str(datetime.now()),
                'Expedition':source,
                'Instrument':instrmt_name,
                'Data Attribution':attribution,
                'Source':source_url,
                'Data Organizer':me_string,
                'Reference':doi,
                'Last modified':str(datetime.now()),
                'Last modification':'File was created',
                'Original vertical measure':og_vert,
                'Original temperature measure':og_temp,
                'Original salinity measure':og_salt,
                'Sub-sample scheme':'None',
                'Moving average window':'None'
    }

    # Convert into a dataset
    ds = xr.Dataset(data_vars=nc_vars, coords=nc_coords, attrs=nc_attrs)
    # Write out to netcdf
    print('Writing data to',out_file)
    ds.to_netcdf(out_file, 'w')

################################################################################

black_list = {'BigBear': [531, 535, 537, 539, 541, 543, 545, 547, 549],
              'BlueFox': [94, 308, 310],
              'Caribou': [],
              'Snowbird': [443],
              'Seacat': ['SH15200.UP', 'SH03200', 'SH30500', 'SH34100']}

################################################################################

def read_ITP_data_file(file_path, file_name, instrmt):
    """
    Reads certain data from an ITP profile file
    Returns a dictionary with that specific information formatted in a manageable way

    file_path           string of a file path to the containing directory
    file_name           string of the file name of a specific file
    instrmt             string of the name of this instrmt (just an unpadded number here)
    """
    # Make sure it isn't a 'sami' file instead of a 'grd' file
    if 'sami' in file_name:
        print('Skipping',file_name)
        return
    # Get just the proper file name, after the slash
    filename2 = file_name.split('/')[-1]
    # Assuming filename format itpXgrdYYYY.dat where YYYY is always 4 digits
    #   works for `final` format above or `cormat` format corYYYY.mat
    try:
        prof_no = int(filename2[-8:-4])
    except:
        print('Skipping',instrmt,'file',filename2)
        return None
    #
    # Check to see whether loading `final` or `cormat` file format
    if 'cormat' in file_path:
        return read_ITP_cormat(file_path, file_name, instrmt, prof_no)
    elif 'final' in file_path:
        return read_ITP_final(file_path, file_name, instrmt, prof_no)
    #

def read_ITP_cormat(file_path, file_name, instrmt, prof_no):
    """
    Loads the data from an ITP profile file in the `cormat` format
    Returns a dictionary with specific information

    file_path           string of a file path to the containing directory
    file_name           string of the file name of a specific file
    instrmt             The number of the ITP that took the profile measurement
    prof_no             The number identifying this specific profile
    """
    # Check to make sure this one isn't on the black list
    if instrmt in black_list.keys():
        if prof_no in black_list[instrmt]:
            on_black_list = True
        else:
            on_black_list = False
    else:
        on_black_list = False
    #
    # Load cormat file into dictionary with mat73
    #   (specific to version of MATLAB used to make cormat files)
    try:
        dat = mat73.loadmat(file_path+'/'+file_name)
    except:
        dat = io.loadmat(file_path+'/'+file_name)
    # Extract certain data from the object, specific to how the files are formatted
    #   The latitude and longitude values where the profile was taken
    lon = float(dat['longitude'])
    lat = float(dat['latitude'])
    #   The date this profile was taken, psdate: profile start or pedate: profile end
    date_MMDDYY_start = dat['psdate']
    date_MMDDYY_end   = dat['pedate']
    #   The time this profile was taken, pstart: profile start or pstop: profile end
    time_HHMMSS_start = dat['pstart']
    time_HHMMSS_end   = dat['pstop']
    # Sometimes, the date comes out as an array. In that case, take the first index
    if not isinstance(date_MMDDYY_start, str):
        date_MMDDYY_start = date_MMDDYY_start[0]
    if not isinstance(date_MMDDYY_end, str):
        date_MMDDYY_end   = date_MMDDYY_end[0]
    if not isinstance(time_HHMMSS_start, str):
        time_HHMMSS_start = time_HHMMSS_start[0]
    if not isinstance(time_HHMMSS_end, str):
        time_HHMMSS_end   = time_HHMMSS_end[0]
    # Attempt to format the datetimes
    try:
        dt_start = str(datetime.strptime(date_MMDDYY_start+' '+time_HHMMSS_start, r'%m/%d/%y %H:%M:%S'))
    except:
        dt_start = None
    try:
        dt_end = str(datetime.strptime(date_MMDDYY_end+' '+time_HHMMSS_end, r'%m/%d/%y %H:%M:%S'))
    except:
        dt_end = None
    # Determine the region
    reg = find_geo_region(lon, lat)
    #
    # If it finds the correct column headers, put data into arrays
    if 'te_adj' in dat and 'sa_adj' in dat and 'pr_filt' in dat:
        press0 = dat['pr_filt'].flatten()
        iT0  = dat['te_adj'].flatten()
        SP0  = dat['sa_adj'].flatten()
        # Convert to SA and CT to follow vdB21a
        pf_vert_len = len(press0)
        SA1 = gsw.SA_from_SP(SP0, press0, [lon]*pf_vert_len, [lat]*pf_vert_len)
        CT1 = gsw.CT_from_t(SA1, iT0, press0)
        # Convert to depth from pressure,
        #   take the negative so that depth increases downward
        depth = -gsw.z_from_p(press0, [lat]*pf_vert_len)
        # Down-casts have an issue with the profiler wake, so note whether the
        #   profile was taken going up or down
        if press0[0] < press0[-1]:
            up_cast = False
        else:
            up_cast = True
        # Create output dictionary for this profile
        out_dict = {'prof_no': prof_no,
                    'black_list': on_black_list,
                    'dt_start': dt_start,
                    'dt_end': dt_end,
                    'lon': lon,
                    'lat': lat,
                    'region': reg,
                    'up_cast': up_cast,
                    'press': press0,
                    'depth': depth,
                    'iT': iT0,
                    'CT': CT1,
                    'SP': SP0,
                    'SA': SA1
                    }
        #
        # Return all the relevant values
        return out_dict
    else:
        return None

def read_ITP_final(file_path, file_name, instrmt, prof_no):
    """
    Loads the data from an ITP profile file in the `final` format
    Returns a dictionary with specific information

    file_path           string of a file path to the containing directory
    file_name           string of the file name of a specific file
    instrmt             The number of the ITP that took the profile measurement
    prof_no             The number identifying this specific profile
    """
    # Check to make sure this one isn't on the black list
    if instrmt in black_list.keys():
        if prof_no in black_list[instrmt]:
            on_black_list = True
        else:
            on_black_list = False
    else:
        on_black_list = False
    #
    # Read data file in as a pandas data object
    #   The arguments used are specific to how the ITP files are formatted
    #   Reading in the file one line at a time because the number of items
    #       on each line is inconsistent between files
    dat0 = pd.read_table(file_path+'/'+file_name, header=None, nrows=2, engine='python', delim_whitespace=True).iloc[1]
    # Extract certain data from the object, specific to how the files are formatted
    #   The date this profile was taken
    try:
        date = str(pd.to_datetime(float(dat0[1])-1, unit='D', origin=str(dat0[0])))
    except:
        date = None
    #   The latitude and longitude values where the profile was taken
    lon = float(dat0[2])
    lat = float(dat0[3])
    # Determine the region
    reg = find_geo_region(lon, lat)
    #
    # Read in data from the file
    dat = pd.read_table(file_path+'/'+file_name,header=2,skipfooter=1,engine='python',delim_whitespace=True)
    # If it finds the correct column headers, put data into arrays
    if 'temperature(C)' and 'salinity' and '%pressure(dbar)' in dat.columns:
        press0 = dat['%pressure(dbar)'][:].values
        iT0  = dat['temperature(C)'][:].values
        SP0  = dat['salinity'][:].values
        # Convert to SA and CT to follow vdB21a
        pf_vert_len = len(press0)
        SA1 = gsw.SA_from_SP(SP0, press0, [lon]*pf_vert_len, [lat]*pf_vert_len)
        CT1 = gsw.CT_from_t(SA1, iT0, press0)
        # `final` formatted profiles are sorted, so no way to tell which direction
        #   They were taken in. So, just mark all as up-casts
        up_cast = True
        # Create output dictionary for this profile
        out_dict = {'prof_no': prof_no,
                    'black_list': on_black_list,
                    'dt_start': date,
                    'dt_end': None,
                    'lon': lon,
                    'lat': lat,
                    'region': reg,
                    'up_cast': up_cast,
                    'press': press0,
                    'depth': [None]*len(press0),
                    'iT': iT0,
                    'CT': CT1,
                    'SP': SP0,
                    'SA': SA1
                    }
        #
        # Return all the relevant values
        return out_dict
    else:
        return None

################################################################################

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def find_geo_region(lon, lat):
    """
    Returns a string of the geographical region of the Arctic Ocean in which
    the coordinates are located. Region boundaries based on:
    Peralta-Ferriz and Woodgate (2015), doi:10.1016/j.pocean.2014.12.005

    lon         longitude value
    lat         latitude value
    """
    if isinstance(lon, type(None)) or isinstance(lat, type(None)):
        return 'error'
    # Chukchi Sea
    elif (lon > -180) & (lon < -155) & (lat > 68) & (lat < 76):
        return 'CS'
    # Southern Beaufort Sea
    elif (lon > -155) & (lon < -120) & (lat > 68) & (lat < 72):
        return 'SBS'
    # Canada Basin
    elif (lon > -155) & (lon < -130) & (lat > 72) & (lat < 84):
        return 'CB'
    # Makarov Basin part 1
    elif (lon > 50) & (lon < 180) & (lat > 83.5) & (lat < 90):
        return 'MB'
    # Makarov Basin part 2
    elif (lon > 141) & (lon < 180) & (lat > 78) & (lat < 90):
        return 'MB'
    # Eurasian Basin part 1
    elif (lon > 30) & (lon < 140) & (lat > 82) & (lat < 90):
        return 'EB'
    # Eurasian Basin part 2
    elif (lon > 110) & (lon < 140) & (lat > 78) & (lat < 82):
        return 'EB'
    # Barents Sea part 1
    elif (lon > 15) & (lon < 60) & (lat > 75) & (lat < 80):
        return 'BS'
    # Barents Sea part 2
    elif (lon > 15) & (lon < 55) & (lat > 67) & (lat < 75):
        return 'BS'
    else:
        return False

################################################################################

## Read instrument makes a netcdf for just the given instrument
read_instrmt('ITP', '1', science_data_file_path+'ITPs/itp1/itp1cormat', 'netcdfs/ITP_1.nc')
# read_instrmt('ITP', '2', science_data_file_path+'ITPs/itp2/itp2cormat', 'netcdfs/ITP_2.nc')
# read_instrmt('ITP', '3', science_data_file_path+'ITPs/itp3/itp3cormat', 'netcdfs/ITP_3.nc')
# read_instrmt('ITP', '13', science_data_file_path+'ITPs/itp13/itp13cormat', 'netcdfs/ITP_13.nc')

## These will make all the netcdfs for a certain source (takes a long time)
# make_all_ITP_netcdfs(science_data_file_path)

exit(0)

## Tests for how the netcdfs turned out

ds2 = xr.load_dataset('netcdfs/ITP_2.nc')
# ds2['Time'] = pd.DatetimeIndex(ds2['Time'].values)
# # print ds.info() to only see the structure
# print(ds2.info())

# This gets all the data for a profile that has the given `prof_no`
# ds_test = ds2.where(ds2.prof_no==185, drop=True).squeeze()
# print(ds_test)
# ds_test.iT.plot()
# plt.show()

# This gets all the TS data from all profiles between certain depths
# ds_test = ds2[['press','iT','SP']]
# ds_test = ds_test.where(ds_test.press>200, drop=True).squeeze()
# ds_test = ds_test.where(ds_test.press<300, drop=True).squeeze()
# ds_test.plot.scatter(x='SP', y='iT', marker='.', s=0.5)
# plt.show()
