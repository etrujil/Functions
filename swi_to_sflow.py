'''
2018-08-10 Ernesto Trujillo

This code generates necessary inputs to StreamFlow (Hydrology model / WSL SLF)
from iSnobal outputs, namely, surface water input (SWI snowmelt/rain output)
into grid ascii format of
'''
import netCDF4 as nc
import numpy as np
import pandas as pd
import os

from subprocess import Popen, check_output, PIPE
from datetime import datetime
from datetime import timedelta

# ------------------------------------------------------------------------------
# Basic Functions


def parse_extent(fname, cellsize_return=0):
    """
    Author: Micah Johnson, mod. by Ernesto Trujillo
    Uses ogr to parse the information of some GIS file and returns a list of the
    response of the things important to this script.
    Args:
        fname: Full path point to file containing GIS information
        cellsize_return: Optional (deafault = 0). 1 will add cellsize as the
                         last element of the return list. Option only for '.asc'
    Returns:
        extent: containing images extent in list type
        [x_ll, y_ll, x_ur, y_ur, cellsize (optional)]
    """
    file_type = fname.split('.')[-1]
    if file_type == 'shp':
        basin_shp_info = check_output(['ogrinfo', '-al', fname],
                                      universal_newlines=True)
        parse_list = basin_shp_info.split('\n')

        # Parse extents from basin info
        for l in parse_list:
            if 'extent' in l.lower():
                k, v = l.split(':')
                parseable = ''.join(c for c in v if c not in ' ()\n')
                parseable = parseable.replace('-', ',')
                extent = [i for i in parseable.split(',')]
                break

    elif file_type == 'tif':
        basin_shp_info = check_output(['gdalinfo', fname],
                                      universal_newlines=True)
        parse_list = basin_shp_info.split('\n')
        extent = []
        for l in parse_list:
            if 'lower left' in l.lower() or 'upper right' in l.lower():
                for w in l.split(' '):
                    try:
                        if len(extent) <= 4:
                            parseable = \
                                ''.join(c for c in w if c not in ' ,()\n')
                            extent.append(float(parseable))
                    except:
                        pass

    elif file_type == 'asc':
        ascii_file = open(fname, 'r')
        ascii_headlines = []
        ascii_headlines = [ascii_file.readline().strip('\n')
                           for i_line in range(6)]
        ascii_file.close()

        parse_list = ascii_headlines
        extent = []
        n_rows = 0
        n_cols = 0
        x_ll = 0
        y_ll = 0
        cellsize = 0

        for l in parse_list:
            if 'xllcorner' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                x_ll = float(w[-1])
            elif 'yllcorner' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                y_ll = float(w[-1])
            elif 'ncols' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                n_cols = float(w[-1])
            elif 'nrows' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                n_rows = float(w[-1])
            elif 'cellsize' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                cellsize = float(w[-1])

            if cellsize_return == 1:
                extent = [x_ll, y_ll,
                          x_ll + (n_cols + 1) * cellsize,
                          y_ll + (n_rows + 1) * cellsize,
                          cellsize]

            else:

                extent = [x_ll, y_ll,
                          x_ll + n_cols * cellsize,
                          y_ll + n_rows * cellsize]

    else:
        raise IOError("File type .{0} not recoginizeable for parsing extent"
                      "".format(file_type))

    return extent

# ------------------------------------------------------------------------------
# File preparation functions


def swi_to_ascii(swinc, swi_x_field, swi_y_field, swi_time_field, swi_field,
                 date_zero, date_ini, date_end, convert_factor=1, utc_in=0,
                 utc_out=0):
    '''
    Converts the SWI output from iSnobal netcdf files to ascii grids
    One grid file per time step with name following the Alpine3D convention
    "YYYY-MM-DDTHH.MM.SS_ROT.asc" as in "2014-08-01T01.00.00_ROT.asc"
    ROT denotes Runoff Total in Alpine3D
    Args:
        swinc: string with name of netcdf file generated from iSnobal/AWSM
               in local folder "./". All output will be written in "./"
        swi_x_field: name of field in "swinc" containing the x dimension
        swi_y_field: name of field in "swinc" containing the y dimension
        swi_time_field: name of field in "swinc" containing the time stamp
        swi_field: name of the field containing SWI (Surface Water Input)
                   SWI from iSNOBAL/AWSM in [kg / m^2 / timestep]
                   (or [mm / m^2 / timestep])
        date_zero: string with the zero date (see swinc
                   fields descriptions for the datestamp field)
                   e.g., in AWSM swi files, description states:
                   float time(time) ;
                                time:units = "hours since 2014-10-01 00:00:00" ;
                                time:time_zone = "utc" ;
                                time:calendar = "standard" ;
                    zero date should be '2014-10-01 00:00:00'
        date_ini: start date for extraction from "swinc"
        date_end: final date of extraction from "swinc"
            all dates in format "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DDTHH:MM:SS"
                                 As in iSnobal/AWSM   or  Alpine3D/MeteoIO
        utc_in: Optional (default = 0): UTC of the input ncfile timestamp
                (e.g.,0 for UTC+0)
        utc_out: Optional (default = 0): desired UTC for the output files
                (e.g.,-7 for UTC-7)
        convert_factor: convertion factor from input units to mm / m^2
                        Typically, SWI input from iSNOBAL/AWSM is in
                        [mm m^-2 timestep^-1] or [kg m^-2 timestep^-1]
                        This convertion factor only applies to depth
                        e.g., default = 1 (from mm to mm). Use 0.001 from m to mm
        nans will be replaced with zero values
        units needed for Alpine3D/StremFlow should be [kg/m^2/hr] (or mm / hr).
    '''

    # Reading in data from files and extracting said data
    swi_dataset = os.path.abspath('./' + swinc)
    ncfile = nc.Dataset(swi_dataset, 'r')

    # Extract fields of interest
    x_vector = ncfile.variables[swi_x_field][:]
    y_vector = ncfile.variables[swi_y_field][:]
    time_field = ncfile.variables[swi_time_field][:]  # time field comes in hrs

    # Determine time steps to export
    zero_timestamp = pd.to_datetime(date_zero, infer_datetime_format=True)
    ini_timestamp = pd.to_datetime(date_ini, infer_datetime_format=True)
    end_timestamp = pd.to_datetime(date_end, infer_datetime_format=True)
    delta_hours = time_field[1] - time_field[0]
    delta_timestamp = timedelta(hours=int(delta_hours))
    delta_timestamp_utc = timedelta(hours=int(utc_out - utc_in))

    # Indices for data extraction
    ini_index = np.floor((ini_timestamp - zero_timestamp) / delta_timestamp)
    if ini_index < 0:
        print('Initial date prior to initial date on file, '
              + 'date reset to initial date on file')
        ini_index = 0
    elif ini_index > len(time_field):
        print('Initial date after final date on file, '
              + 'check initial date and re-run')
        return None

    end_index = np.ceil((end_timestamp - zero_timestamp) / delta_timestamp)
    if end_index < ini_index:
        print('Final date prior to initial date, '
              + 'date reset to match initial date')
        end_index = ini_index
    elif end_index >= len(time_field):
        print('Final date after final date of file, '
              + 'date reset to match final date on file')
        end_index = len(time_field) - 1

    n_x = len(x_vector)
    n_y = len(y_vector)
    # Be careful if coordinate system is lat-lon and southern hemisphere, etc.
    # Input fields should be in projected coords (in meters)
    x_ll = x_vector.min()
    y_ll = y_vector.min()
    cell_size = (x_vector[1]-x_vector[0])

    # Create header for file
    header = "ncols %s\n" % n_x
    header += "nrows %s\n" % n_y
    header += "xllcorner %s\n" % x_ll
    header += "yllcorner %s\n" % y_ll
    header += "cellsize %s\n" % cell_size
    header += "NODATA_value -999"

    # Loop through dates and create files
    for i_files in range(int(ini_index), int(end_index+1)):
        swi_matrix = ncfile.variables[swi_field][i_files, :, :]
        # in case timestep is more than one hr
        swi_matrix *= (convert_factor/delta_hours)
        swi_matrix[np.isnan(swi_matrix)] = 0.0
        i_timestamp = (zero_timestamp
                       + i_files*delta_timestamp
                       + delta_timestamp_utc)

        # Alpine3D name convention for Runoff Total (ROT)
        file_name_swi = i_timestamp.strftime('%Y-%m-%dT%H.%M.%S_ROT.asc')
        print(file_name_swi)
        np.savetxt(file_name_swi, swi_matrix,
                   header=header, fmt="%1.3e", comments='')

    ncfile.close()

# ------------------------------------------------------------------------------


def swi_to_cactmentswi(watershed_file, swinc,
                       swi_x_field,
                       swi_y_field, swi_time_field, swi_field, date_zero,
                       date_ini, date_end, convert_factor=0.001, utc_in=0,
                       utc_out=0, output_option=1,
                       filename_out='watershed_swi_flow.csv'):
    '''
    Converts the SWI output from iSnobal netcdf files to time series per
    catchment in a txt file. The file is organized to have the area of each
    catchment and a time series of SWI per catchment (see description below)
    Args:
        watershed_file: full path of ascii grid file defining the watersheds over
                        which to calculate total SWI,
                        grid in file needs to coincide with grid
                        in SWI files. Units of gridcell in meters
                        (projected coord. sys.. Must use same projection
                         as 'swinc')
        swinc: full path of netcdf file generated from iSnobal/AWSM.
               All output will be written in "./"
               swi is assumed to be in kg m^-2 per timestep (or mm water per
               timestep)

               output is given in m^3 s^-1

        swi_x_field: name of field in "swinc" containing the x dimension
        swi_y_field: name of field in "swinc" containing the y dimension
        swi_time_field: name of field in "swinc" containing the time stamp
        swi_field: name of the field containing SWI (Surface Water Input)
                   SWI from iSNOBAL/AWSM in [kg / m^2 / timestep]
                   (or [mm / m^2 / timestep])
        date_zero: string with the zero date (see swinc
                   fields descriptions for the datestamp field)
                   e.g., in AWSM swi files, description states:
                   float time(time) ;
                                time:units = "hours since 2014-10-01 00:00:00" ;
                                time:time_zone = "utc" ;
                                time:calendar = "standard" ;
                    zero date should be '2014-10-01 00:00:00'
        date_ini: start date for extraction from "swinc"
        date_end: final date of extraction from "swinc"
            all dates in format "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DDTHH:MM:SS"
                                 As in iSnobal/AWSM   or  Alpine3D/MeteoIO
        convert_factor: Optional (default = 0.001) convertion factor from input
                        units to desired output units
                        The default of 0.001 is from mm-SWI to m-SWI
                        nans will be replaced with zero values
        utc_in: Optional (default = 0): UTC of the input ncfile timestamp
                (e.g.,0 for UTC+0)
        utc_out: Optional (default = 0): desired UTC for the output files
                (e.g.,-7 for UTC-7)
        output_option: Optional (default = 1): 0 for only a return pd.DataFrame
                                               1 if csv file and pd.DataFrame
                                                    are to be saved and returned
                                               2 if only csv file is to be
                                                    generated
        filename_out: Optional if no csv output is to be created
                      (default = "watershed_swi_flow.csv")
                      (see output_option below) . Otherwise, provide output csv
                      file name - include '.csv' extension at the end
    Returns:
        pd.DataFrame with watershed areas and time series of total swi volume
        per unit time (m^3 s-1) if output_option = 0 or 1.
    '''

    # Reading in data from files and extracting said data
    # watershed_path = os.path.abspath('./' + watershed_file)
    # Uncomment line above and comment line below to use local path
    watershed_path = watershed_file
    watershed_grid = np.loadtxt(watershed_path, skiprows=6)
    # swi_dataset = os.path.abspath('./' + swinc)
    swi_dataset = swinc
    # Uncomment line above and comment line below to use local path
    ncfile = nc.Dataset(swi_dataset, 'r')

    # Extract fields of interest
    x_vector = ncfile.variables[swi_x_field][:]
    y_vector = ncfile.variables[swi_y_field][:]
    time_field = ncfile.variables[swi_time_field][:]

    # Determine time steps to export
    zero_timestamp = pd.to_datetime(date_zero, infer_datetime_format=True)
    ini_timestamp = pd.to_datetime(date_ini, infer_datetime_format=True)
    end_timestamp = pd.to_datetime(date_end, infer_datetime_format=True)
    delta_hours = time_field[1] - time_field[0]
    delta_timestamp = timedelta(hours=int(delta_hours))
    delta_timestamp_utc = timedelta(hours=int(utc_out - utc_in))

    # Indices for data extraction
    ini_index = np.floor((ini_timestamp - zero_timestamp) / delta_timestamp)
    if ini_index < 0:
        print('Initial date prior to initial date on file, '
              + 'date reset to initial date on file')
        ini_index = 0
    elif ini_index > len(time_field):
        print('Initial date after final date on file, '
              + 'check initial date and re-run')
        return None

    end_index = np.ceil((end_timestamp - zero_timestamp) / delta_timestamp)
    if end_index < ini_index:
        print('Final date prior to initial date, '
              + 'date reset to match initial date')
        end_index = ini_index
    elif end_index >= len(time_field):
        print('Final date after final date of file, '
              + 'date reset to match final date on file')
        end_index = len(time_field) - 1

    n_x = len(x_vector)
    n_y = len(y_vector)
    # Be careful if coordinate system is lat-lon and southern hemisphere, etc.
    # Should be in meters (projected coordinate system)
    x_ll = x_vector.min()
    y_ll = y_vector.min()
    cell_size = (x_vector[1]-x_vector[0])  # in meters

    time_zone = utc_out - utc_in

    # determine watershed indices
    watershed_unique = np.unique(watershed_grid, return_counts=True)
    watershed_area = pd.DataFrame(watershed_unique[1][watershed_unique[0] > 0],
                                  index=watershed_unique[0][watershed_unique[0] > 0],
                                  columns=['area'])
    watershed_numbers = watershed_unique[0][watershed_unique[0] > 0]
    watershed_numbers = watershed_numbers.astype(int)
    watershed_area.area = watershed_area.area * cell_size * cell_size

    swi_timeseries = pd.DataFrame(
        watershed_area.area.values, index=watershed_numbers, columns=['area'])

    # Loop through dates and create files
    for i_files in range(int(ini_index), int(end_index+1)):

        swi_matrix = ncfile.variables[swi_field][i_files, :, :]
        swi_matrix *= convert_factor  # should give m SWI
        # (or m m^-2 per time step)
        swi_matrix[np.isnan(swi_matrix)] = 0.0
        i_timestamp = (zero_timestamp
                       + i_files*delta_timestamp
                       + delta_timestamp_utc)

        swi_volume = pd.DataFrame(np.array(
            [np.nan] * len(watershed_numbers)), index=watershed_numbers,
            columns=[i_timestamp])

        for i_catch in watershed_numbers:

            swi_volume.at[i_catch, i_timestamp] = \
                (swi_matrix[watershed_grid == i_catch].sum()
                 * cell_size * cell_size
                 / (delta_hours * 3600.))
            # Output in m^3 s^-1

        swi_timeseries[i_timestamp] = swi_volume.loc[:, i_timestamp].values

    swi_timeseries = swi_timeseries.T

    # Output
    if output_option == 0:
        return swi_timeseries
    elif output_option == 1:
        swi_timeseries.to_csv(filename_out)
        return swi_timeseries
    elif output_option == 2:
        swi_timeseries.to_csv(filename_out)
        return None
    else:
        print("Invalid output selection, no output or files were generated")
        return None

    ncfile.close()

# ------------------------------------------------------------------------------
# UNFINISHED
# the idea here is to be able to provide swi and watershed defined in grids with
# different cell sizes and extents


def swi_to_cactmentswi_diff_grid(watershed_file, swinc,
                                 swi_x_field,
                                 swi_y_field, swi_time_field, swi_field, date_zero,
                                 date_ini, date_end, convert_factor=0.001, utc_in=0,
                                 utc_out=0, output_option=1,
                                 filename_out='watershed_swi_flow.csv'):
    '''
    Converts the SWI output from iSnobal netcdf files to time series per
    catchment in a txt file. The file is organized to have the area of each
    catchment and a time series of SWI per catchment (see description below)
    Args:
        watershed_file: full path of ascii grid file defining the watersheds over
                        which to calculate total SWI,
                        grid in file does not need to coincide with grid
                        in SWI files. Units of gridcell in meters
                        (projected coord. sys.. Must use same projection
                         as 'swinc')
        swinc: full path of netcdf file generated from iSnobal/AWSM.
               All output will be written in "./"
               swi is assumed to be in kg m^-2 per timestep (or mm water per
               timestep)

               output is given in m^3 s^-1

        swi_x_field: name of field in "swinc" containing the x dimension
        swi_y_field: name of field in "swinc" containing the y dimension
        swi_time_field: name of field in "swinc" containing the time stamp
        swi_field: name of the field containing SWI (Surface Water Input)
                   SWI from iSNOBAL/AWSM in [kg / m^2 / timestep]
                   (or [mm / m^2 / timestep])
        date_zero: string with the zero date (see swinc
                   fields descriptions for the datestamp field)
                   e.g., in AWSM swi files, description states:
                   float time(time) ;
                                time:units = "hours since 2014-10-01 00:00:00" ;
                                time:time_zone = "utc" ;
                                time:calendar = "standard" ;
                    zero date should be '2014-10-01 00:00:00'
        date_ini: start date for extraction from "swinc"
        date_end: final date of extraction from "swinc"
            all dates in format "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DDTHH:MM:SS"
                                 As in iSnobal/AWSM   or  Alpine3D/MeteoIO
        convert_factor: Optional (default = 0.001) convertion factor from input
                        units to desired output units
                        The default of 0.001 is from mm-SWI to m-SWI
                        nans in input will be replaced with zero values
        utc_in: Optional (default = 0): UTC of the input ncfile timestamp
                (e.g.,0 for UTC+0)
        utc_out: Optional (default = 0): desired UTC for the output files
                (e.g.,-7 for UTC-7)
        output_option: Optional (default = 1): 0 for only a return pd.DataFrame
                                               1 if csv file and pd.DataFrame
                                                    are to be saved and returned
                                               2 if only csv file is to be
                                                    generated
        filename_out: Optional if no csv output is to be created
                      (default = "watershed_swi_flow.csv")
                      (see output_option below) . Otherwise, provide output csv
                      file name - include '.csv' extension at the end
    Returns:
        pd.DataFrame with watershed areas and time series of total swi volume
        per unit time (m^3 s-1) if output_option = 0 or 1.
    '''

    # Reading in data from files and extracting said data
    # watershed_path = os.path.abspath('./' + watershed_file)
    # Uncomment line above and comment line below to use local path
    watershed_path = watershed_file
    watershed_grid = np.loadtxt(watershed_path, skiprows=6)
    # swi_dataset = os.path.abspath('./' + swinc)
    swi_dataset = swinc
    # Uncomment line above and comment line below to use local path
    ncfile = nc.Dataset(swi_dataset, 'r')

    # Extract fields of interest
    x_vector = ncfile.variables[swi_x_field][:]
    y_vector = ncfile.variables[swi_y_field][:]
    time_field = ncfile.variables[swi_time_field][:]

    # Determine time steps to export
    zero_timestamp = pd.to_datetime(date_zero, infer_datetime_format=True)
    ini_timestamp = pd.to_datetime(date_ini, infer_datetime_format=True)
    end_timestamp = pd.to_datetime(date_end, infer_datetime_format=True)
    delta_hours = time_field[1] - time_field[0]
    delta_timestamp = timedelta(hours=int(delta_hours))
    delta_timestamp_utc = timedelta(hours=int(utc_out - utc_in))

    # Indices for data extraction
    ini_index = np.floor((ini_timestamp - zero_timestamp) / delta_timestamp)
    if ini_index < 0:
        print('Initial date prior to initial date on file, '
              + 'date reset to initial date on file')
        ini_index = 0
    elif ini_index > len(time_field):
        print('Initial date after final date on file, '
              + 'check initial date and re-run')
        return None

    end_index = np.ceil((end_timestamp - zero_timestamp) / delta_timestamp)
    if end_index < ini_index:
        print('Final date prior to initial date, '
              + 'date reset to match initial date')
        end_index = ini_index
    elif end_index >= len(time_field):
        print('Final date after final date of file, '
              + 'date reset to match final date on file')
        end_index = len(time_field) - 1

    # Determine extents of input SWI (swinc) and Watershed definition file
    # (watershed_file)
    swi_n_x = len(x_vector)
    swi_n_y = len(y_vector)
    # Be careful if coordinate system is lat-lon and southern hemisphere, etc.
    # Should be in meters (projected coordinate system)
    swi_xll = x_vector.min()
    swi_yll = y_vector.min()
    swi_dx = abs(x_vector[1]-x_vector[0])  # in meters
    swi_dy = abs(x_vector[1]-x_vector[0])  # in meters
    watershed_extent = parse_extent(watershed_file, cellsize_return=1)

    time_zone = utc_out - utc_in

    # determine watershed indices
    watershed_unique = np.unique(watershed_grid, return_counts=True)
    watershed_area = pd.DataFrame(watershed_unique[1][watershed_unique[0] > 0],
                                  index=watershed_unique[0][watershed_unique[0] > 0], columns=['area'])
    watershed_numbers = watershed_unique[0][watershed_unique[0] > 0]
    watershed_numbers = watershed_numbers.astype(int)
    watershed_area.area = watershed_area.area * swi_dx * swi_dy

    swi_timeseries = pd.DataFrame(
        watershed_area.area.values, index=watershed_numbers, columns=['area'])

    # Loop through dates and create files
    for i_files in range(int(ini_index), int(end_index+1)):

        swi_matrix = ncfile.variables[swi_field][i_files, :, :]
        swi_matrix *= convert_factor  # should give m SWI
        # (or m m^-2 per time step)
        swi_matrix[np.isnan(swi_matrix)] = 0.0
        i_timestamp = (zero_timestamp
                       + i_files*delta_timestamp
                       + delta_timestamp_utc)

        swi_volume = pd.DataFrame(np.array(
            [np.nan] * len(watershed_numbers)), index=watershed_numbers, columns=[i_timestamp])

        for i_catch in watershed_numbers:

            swi_volume.at[i_catch, i_timestamp] = (swi_matrix[watershed_grid == i_catch].sum()
                                                   * swi_dx * swi_dy
                                                   / (delta_hours * 3600.))
            # Output in m^3 s^-1

        swi_timeseries[i_timestamp] = swi_volume.loc[:, i_timestamp].values

    swi_timeseries = swi_timeseries.T
    swi_timeseries.to_csv(filename_out)

    ncfile.close()


# ------------------------------------------------------------------------------

# UNFINISHED!
#
#
# def cactmentswi_to_smet(watershed_swi_file_names, date_ini, date_end, freq_out='60min', fill_nans=True):
#     '''
#     Generates SMET files (Alpine3D format) containing the time series of SWI
#     for each subwatershed (one SMET file per subwatershed) to be used as input
#     for StreamFlow. The function compiles the data in the SWI files listed in
#     'watershed_swi_file_names' to fill in the period between 'date_ini' and
#     'date_end'
#
#     Args:
#         watershed_swi_file_names: name of text file containing the names of the
#                                   individual files to combine in the SMET
#                                   files. The file must contain one filename per
#                                   line. E.g.,
#                                   'SWI_2007_tollgate_wfile_150.csv
#                                   SWI_2008_tollgate_wfile_150.csv
#                                   SWI_2009_tollgate_wfile_150.csv'
#
#                                   input is in m^3 s^-1
#
#                                   Provide full path if files are not in current
#                                   directory
#
#         date_ini: start date of extraction
#         date_end: final date of extraction
#             all dates in format "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DDTHH:MM:SS"
#                                  As in iSnobal/AWSM   or  Alpine3D/MeteoIO
#
#         [freq_out]: time step of output (and input), e.g., '15min', '60min'
#                     , '1440min'
#             (default = '60min', adjust to match frequency of the input SWI)
#
#         [fill_nans]: option to replace nans with zeros, default: True
#     '''
#
#     # Read files and store data
#     n_files = 0
#     swi_n_fields = []  # list
#     swi_dict = {}  # dictionary
#
#     with open(watershed_swi_file_names) as f_swi:
#         file_names = f_swi.read().split('\n')
#
#     for i_files in range(len(file_names)):
#
#         if os.path.isfile(file_names[i_files]) == True:
#
#             n_files += 1
#             swi_data = pd.read_csv(
#                 file_names[i_files], index_col=0)
#             # First row includes watershed areas
#             if n_files == 1:
#                 swi_dict['area'] = swi_data.loc['area']
#
#             swi_data = swi_data[1:]
#             swi_data.index = pd.to_datetime(swi_data.index)
#             # the option parse_dates does not read the csv file dates to
#             # datetime, likely because the first row has an index named 'area'
#             swi_dict[n_files] = swi_data
#             # used to keep track of fields in swi_dict
#             swi_n_fields.append(n_files)
#
#     ini_timestamp = pd.to_datetime(date_ini, infer_datetime_format=True)
#     end_timestamp = pd.to_datetime(date_end, infer_datetime_format=True)
#
#     index_series = pd.period_range(start=ini_timestamp, end=end_timestamp,
#                                    freq=freq_out)
#
#     index_series = index_series.to_timestamp()
#
#     swi_df = pd.DataFrame(np.array([np.nan] * len(index_series)),
#                           index=index_series, columns=['SWI'])
#
#     for i_catch in swi_dict['area'].index:
#         file_name = 'catch' + str(int(i_catch)) + '.smet'
#         # print(file_name)
#
#         swi_out = swi_df.copy()  # contains nans
#
#         for i_file in swi_n_fields:
#
#             index_zero = max(swi_out.index[0], swi_dict[i_file].index[0])
#             index_end = min(swi_out.index[-1], swi_dict[i_file].index[-1])
#
#             print(index_zero)
#             print(index_end)
#
#             print(swi_out.loc[index_zero:index_end])
#             print(swi_dict[i_file][i_catch].loc[index_zero:index_end])
#             swi_out['SWI'].loc[index_zero:index_end] = swi_dict[i_file][i_catch].loc[index_zero:index_end]
#             # print(i_file)
#             # swi_out.to_csv(file_name)
