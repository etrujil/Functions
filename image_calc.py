#!/usr/bin/env python3

import argparse
import netCDF4 as nc
import numpy as np
import xarray as xr
from subprocess import check_output
from skimage import io
from os.path import isfile


# Function to get extent


def get_extent(fname, x_field='x', y_field='y'):
    """Uses gdal or python to parse the information of a raster or netcdf file and returns a list of the
    response of the things important to this script.

    Args:
        fname (str): Full path to file
        x_field (str, optional): variable name for x coords in netcdf file. Defaults to 'x'.
        y_field (str, optional): variable name for y coords in netcdf file. Defaults to 'y'.

    Raises:
        IOError: File type not recoginizeable for parsing extent

    Returns:
        extent (list): contains file extent in list type
        [x_ll, y_ll, x_ur, y_ur, cellsize, nodata_value(not for '.nc' files)]
    """

    file_type = fname.split('.')[-1].lower()

    if file_type == 'tif':
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

        for l in parse_list:
            if 'pixel size' in l.lower():
                
                # two-step parsing
                l_lst = l.split()
                l_parsed = []
                for w in l_lst:
                    l_parsed += w.strip('()').split(',')
                    
                try:
                    extent.append(abs(float(l_parsed[-2])))
                except:
                    extent.append(np.nan)

        for l in parse_list:
            if 'nodata value' in l.lower():
                try:
                    extent.append(float(l.split('=')[-1].strip(' ')))
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
        nodata_value = -9999

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
                cellsize = abs(float(w[-1]))
            elif 'nodata_value' in l.lower():
                w = l.split(' ')
                w = [w[i_w] for i_w in range(len(w)) if w[i_w] != '']
                nodata_value = float(w[-1])

            extent = [x_ll, y_ll,
                      x_ll + (n_cols) * cellsize,
                      y_ll + (n_rows) * cellsize,
                      cellsize, nodata_value]

    elif file_type == 'nc':

        ncfile = nc.Dataset(fname, 'r')

        # Extract fields of interest
        try:
            x_vector = ncfile.variables[x_field][:]
        except KeyError:
            print('KeyError: x_field key not found')
            return None

        try:
            y_vector = ncfile.variables[y_field][:]
        except KeyError:
            print('KeyError: y_field key not found')
            return None

        # Determine extents of input netCDF file

        n_cols = len(x_vector)
        n_rows = len(y_vector)
        # Be careful if coordinate system is lat-lon and southern hemisphere, etc.
        # Should be in meters (projected coordinate system)
        dx = abs(x_vector[1]-x_vector[0])  # in meters
        dy = abs(y_vector[1]-y_vector[0])  # in meters
        x_ll = x_vector.min() - dx/2  # the nc_file uses center of cell coords
        y_ll = y_vector.min() - dy/2  # Change if not cell center coords

        extent = [x_ll, y_ll,
                  x_ll + (n_cols) * dx,
                  y_ll + (n_rows) * dy,
                  dx, np.nan]

        ncfile.close()

    else:
        raise IOError("File type .{0} not recoginizeable for parsing extent"
                      "".format(file_type))

    return extent


def load_img(fpath, nc_field=None, timestep=None):
    """Function to read an image and return as a numpy array

    Args:
        fpath (str): Path to file (.asc, .tif, or .nc)
        nc_field (str, optional): Name of variable in netcdf file. Defaults to None.
        timestep (str, optional): Date to extract in netcdf file (e.g., '2020-10-01', or '2020-10-01 23:00'). Defaults to None.

    Raises:
        IOError: File not found
        KeyError: Time coordinate does not exist in netcdf file
        IOError: File type not recognized

    Returns:
        numpy.array: Array read in input file
    """

    # Check files exist
    if not isfile(fpath):
        raise IOError("{} not found, please check and rerun".format(fpath))

    # load fpath
    file_type = fpath.split('.')[-1].lower()
    if file_type == 'asc':
        np_out = np.loadtxt(fpath, skiprows=6)
        extent = get_extent(fpath)
        if len(extent) == 6:
            np_out[np_out == extent[5]] = np.nan

    elif file_type == 'tif':
        np_out = io.imread(fpath)
        extent = get_extent(fpath)
        if len(extent) == 6:
            np_out[np_out == extent[5]] = np.nan

    elif file_type == 'nc':
        nc_ds = xr.open_dataset(fpath)
        extent = get_extent(fpath, x_field='x', y_field='y')
        try:
            np_out = nc_ds[nc_field].sel(time=timestep).values[:]
            if len(np_out.shape) == 3:
                # e.g., when user provide day date: '2013-04-01',
                # not to hour precision: '2013-04-01 23:00'
                np_out = nc_ds[nc_field].sel(time=timestep).values[-1, :, :]
        except (ValueError, KeyError):
            raise KeyError(
                'Check that the file has a time coordinate and try again...')

    else:
        raise IOError("File type .{} not recoginizeable for parsing extent"
                      "".format(file_type))

    return np_out


def image_calc(fpath_a, fpath_b, fpath_mask_a=None, fpath_mask_b=None,
                   nc_field_a=None, time_a=None, convert_factor_a=1.0,
                   nc_field_b=None, time_b=None, convert_factor_b=1.0,
                   calc='*', fname_out='img_out.asc'):
    """Performs operation on input rasters

    Args:
        fpath_a (str): path to file A  (tif, ascii-grid, netcdf)
        fpath_b (str): path to file B (tif, ascii-grid, netcdf)
        fpath_mask_a (str, optional): path to mask file A  (tif, ascii-grid). Defaults to None.
        fpath_mask_b (str, optional): path to mask file B  (tif, ascii-grid). Defaults to None.
        nc_field_a (str, optional): Name of netcdf variable in file A to calculate over, required for netcdf. Defaults to None.
        time_a (str, optional): date-time to calculate over in file A, required for netcdf (e.g., '2020-10-01', or '2020-10-01 23:00'). Defaults to None.
        convert_factor_a (float, optional): multiplier for field A for unit conversion. Defaults to 1.0.
        nc_field_b (str, optional): Name of netcdf variable in file B to calculate over, required for netcdf. Defaults to None.
        time_b (str, optional): date-time to calculate over in file B, required for netcdf (e.g., '2020-10-01', or '2020-10-01 23:00'). Defaults to None.
        convert_factor_b (float, optional): multiplier for file B for unit conversion. Defaults to 1.0.
        calc (str, optional): python operator, use (*, +, -, /, %, **, //). Defaults to '*'.
                              Note: Zero division returns a no value
        fname_out (str, optional): Name of output ascid-grid desired - include '.asc'. Defaults to 'img_out.asc'.

    Raises:
        IOError: File A not found or file type not recognized
        KeyError: Key or time dimension not found in file A
        IOError: Mask A file not found or file type not recognized
        ValueError: Extents of file A and mask A are not the same
        IOError: File A not found or file type not recognized
        KeyError: Key or time dimension not found in file B
        IOError: Mask B file not found or file type not recognized
        ValueError: Extents of file B and mask B are not the same
        ValueError: Unrecognized operator
        ValueError: Extents of file A and file B are not the same
        IOError: Saving output file failed

    Returns:
        saves ascii-grid with A calc B

    """

    # load files

    # fpath_a
    try:
        np_in_a = load_img(fpath_a, nc_field=nc_field_a, timestep=time_a)
        extent_a = get_extent(fpath_a)
    except IOError:
        raise IOError(
            'File not found or file type not recognized in {}, please check and try again...'.format(fpath_a))
    except KeyError:
        raise KeyError(
            'Key or time dimension not found in {}, please check and try again...'.format(fpath_a))

    # mask_a
    if fpath_mask_a is not None:
        try:
            mask_a = load_img(fpath_mask_a)
            extent_mask_a = get_extent(fpath_mask_a)
        except IOError:
            raise IOError('File not found or file type not recognized in {}, please check and try again...'.format(fpath_mask_a))
        
        if all([ex1==ex2 for ex1, ex2 in zip(extent_a[:4], extent_mask_a[:4])]):
            np_in_a[mask_a < 1.0] = np.nan
        else:
            raise ValueError('Extents of {} and {} are not the same'.format(fpath_a, fpath_mask_a))

    # fpath_b
    try:
        np_in_b = load_img(fpath_b, nc_field=nc_field_b, timestep=time_b)
        extent_b = get_extent(fpath_b)
    except IOError:
        raise IOError('File not found or file type not recognized in {}, please check and try again...'.format(fpath_b))
    except KeyError:
        raise KeyError('Key or time dimension not found in {}, please check and try again...'.format(fpath_b))
    
    # mask_b
    if fpath_mask_b is not None:
        try:
            mask_b = load_img(fpath_mask_b)
            extent_mask_b = get_extent(fpath_mask_b)
        except IOError:
            raise IOError('File not found or file type not recognized in {}, please check and try again...'.format(fpath_mask_b))
        
        if all([ex1==ex2 for ex1, ex2 in zip(extent_b[:4], extent_mask_b[:4])]):
            np_in_b[mask_b < 1.0] = np.nan
        else:
            raise ValueError('Extents of {} and {} are not the same'.format(fpath_b, fpath_mask_b))
    
    # operation
    if all([ex1==ex2 for ex1, ex2 in zip(extent_a[:5], extent_b[:5])]):
        np_in_a = np_in_a * convert_factor_a
        np_in_b = np_in_b * convert_factor_b
        if calc == '*':
            np_out = np_in_a * np_in_b
        elif calc == '+':
            np_out = np_in_a + np_in_b
        elif calc == '-':
            np_out = np_in_a - np_in_b
        elif calc == '/':
            np_in_b[np_in_b == 0] = np.nan # added to avoid zero division errors
            np_out = np_in_a / np_in_b
        elif calc == '%':
            np_out = np_in_a % np_in_b
        elif calc == '**':
            np_out = np_in_a ** np_in_b
        elif calc == '//':
            np_out = np_in_a // np_in_b
        else:
            raise ValueError('Unrecognized operator {}, use (*, +, -, /, %, **, //)'.format(calc))
        
    else:
        raise ValueError('Extents of {} and {} are not the same'.format(fpath_b, fpath_b))

    nan_flag = -9999
    np_out[np.isnan(np_out)] = nan_flag

    # file output
    # Create header for file
    header = "ncols {}\n".format(np_out.shape[1])
    header += "nrows {}\n".format(np_out.shape[0])
    header += "xllcorner {}\n".format(extent_a[0])
    header += "yllcorner {}\n".format(extent_a[1])
    header += "cellsize {}\n".format(extent_a[4])
    header += "NODATA_value {}".format(int(nan_flag))

    try:
        np.savetxt(fname_out, np_out,
                header=header, fmt="%.3f", comments='')
        print('file {} saved in current directory'.format(fname_out))
    except:
        raise IOError("Saving file {} failed...")
    
def main():
    p = argparse.ArgumentParser(description='Operates over two fields in specified files - '
                                ' optimized for iSnobal snow densities and depths')

    p.add_argument("-A", "--pathA", dest="fpath_a",
                   required=True, type=str,
                   help="Path to file A (tif, ascii-grid, netcdf)")
    p.add_argument("-B", "--pathB", dest="fpath_b",
                   required=True, type=str,
                   help="Path to file B (tif, ascii-grid, netcdf)")
    p.add_argument("-mA", "--maskA", dest="fpath_mask_a",
                   required=False, type=str, default=None,
                   help="Path to mask file A (tif, ascii-grid)")
    p.add_argument("-mB", "--maskB", dest="fpath_mask_b",
                   required=False, type=str, default=None,
                   help="Path to mask file B (tif, ascii-grid, netcdf)")
    p.add_argument("-fA", "--fieldA", dest="nc_field_a",
                   required=False, type=str, default=None,
                   help="Name of netcdf variable in file A to calculate over (optional, required for netcdf)")
    p.add_argument("-tA", "--timeA", dest="time_a",
                   required=False, type=str, default=None,
                   help="string with date-time to calculate over in file A (e.g., '2012-10-01 23:00)")
    p.add_argument("-cfA", "--convfactA", dest="convert_factor_a",
                   required=False, type=float, default=1.0,
                   help="multiplier for file A for unit conversion. Default: 1.0")
    p.add_argument("-fB", "--fieldB", dest="nc_field_b",
                   required=False, type=str, default=None,
                   help="Name of netcdf variable in file B to calculate over (optional, required for netcdf)")
    p.add_argument("-tB", "--timeB", dest="time_b",
                   required=False, type=str, default=None,
                   help="string with date-time to calculate over in file B (e.g., '2012-10-01 23:00)")
    p.add_argument("-cfB", "--convfactB", dest="convert_factor_b",
                   required=False, type=float, default=1.0,
                   help="multiplier for file B for unit conversion. Default: 1.0")
    p.add_argument("-c", "--calc", dest="calc",
                   required=False, type=str, default="*",
                   help="operator (e.g., '*' - default)")
    p.add_argument("-fno", "--fnout", dest="fname_out",
                   required=False, type=str, default="img_out.asc",
                   help="Name of output ascid-grid desired - include '.asc' (optional, defaul: 'img_out.asc')")

    args = p.parse_args()

    image_calc(args.fpath_a, args.fpath_b , args.fpath_mask_a, args.fpath_mask_b, 
                   args.nc_field_a, args.time_a, args.convert_factor_a,
                   args.nc_field_b, args.time_b, args.convert_factor_b,
                   args.calc, args.fname_out)

if __name__ == '__main__':
    main()
