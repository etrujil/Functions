#!/usr/bin/env python3

import argparse
import netCDF4 as nc
import numpy as np
from skimage import io
from os.path import isfile

def replace_nc_field(fpath_nc, field_name, field_array, mask=None):
    """
    Parameters: fpath_nc: str
                    full path to the netCDF file with the field to
                    replace
                field_name: str
                    name of field in netCDF file to be replaced
                field_array: array_like or str
                    array to replace the field with (e.g., ndarray), or
                    full path to the file containing the new array
                    (e.g., ascii-grid, or tif)
                mask: array_like or str, Optional
                    numpy array with mask or
                    full path to the file containing a mask array (1 or True for the values inside the mask), 
                    accepts .asc or .tif files
                    With a mask, only the values inside the mask are replaced

    Returns:    bool
                True if field_name replaced and modified netCDF file in fapath_nc

    Note: this method is written for 2-dimensional arrays, and it might not work properly for 3- or n- dimensional arrays
    """
    # Check files exist
    if isfile(fpath_nc):
        nc_in = nc.Dataset(fpath_nc, 'r+')
        if not(field_name in nc_in.variables.keys()):
            print("field_name not found in netCDF file")
            nc_in.close()
            return False

    else:
        print('File ' + fpath_nc +' not found')
        return False

    if type(field_array) is str:
        if isfile(field_array):
            parse_path = field_array.split('/')
            if '.asc' in parse_path[-1].lower():
                np_field = np.loadtxt(field_array, skiprows=6)
            elif '.tif' in parse_path[-1].lower():
                np_field = io.imread(field_array)
            else:
                print('Error: unrecognized format of field_array, use .asc or .tif')
                nc_in.close()
                return False
        else:
            print('File ' + field_array +' not found')
            nc_in.close()
            return False
    elif isinstance(field_array, np.ndarray):
        np_field = field_array
    else:
        print('Unrecognized field_array data type')

    if mask is not None:
        if type(mask) is str:
            if isfile(mask):
                parse_path = mask.split('/')
                if '.asc' in parse_path[-1].lower():
                    np_mask = np.loadtxt(mask, skiprows=6)
                elif '.tif' in parse_path[-1].lower():
                    np_mask = io.imread(mask)
                else:
                    print('Error: unrecognized format of mask, use .asc or .tif')
                    nc_in.close()
                    return False
            else:
                print('File ' + mask +' not found')
                nc_in.close()
                return False
        elif isinstance(mask, np.ndarray):
            np_mask = mask

# Replace variables
    if mask is None:
        if nc_in.variables[field_name][:].shape == np_field.shape:
            nc_in.variables[field_name][:] = np_field
            nc_in.close()
            return True
        else:
            print('The input field_array does not match the dimensions of ' + field_name)
            nc_in.close()
            return False
    else:
        if (nc_in.variables[field_name][:].shape == np_field.shape) and \
                (nc_in.variables[field_name][:].shape == np_mask.shape):
            
            np_aux = nc_in.variables[field_name][:]
            np_aux[np_mask > 0] = np_field[np_mask > 0]
            nc_in.variables[field_name][:] = np_aux
            nc_in.close()
            return True
        else:
            print('The input field_array does not match the dimensions of ' + field_name)
            nc_in.close()
            return False

    nc_in.close()

def main():
    p = argparse.ArgumentParser(description='Replaces a variable in a netcdf file - '
                                ' optimized for the topo.nc files in iSnobal')
    p.add_argument("-nc", "--ncpath", dest="fpath_nc",
                   required=True, type=str,
                   help="Path to netcdf file to modify")
    p.add_argument("-fn", "--fieldname", dest="field_name",
                   required=True, type=str,
                   help="Name of field/variable to replace")
    p.add_argument("-fs", "--fieldarray", dest="field_array",
                   required=True, type=str,
                   help="Path to file (tif or ascii-grid) with array that will replace the field")
    p.add_argument("-m", "--mask", dest="mask",
                   required=False, type=str,
                   help="Path to file (tif or ascii-grid) with mask array, only values inside mask are replaced")

    args = p.parse_args()

    replace_nc_field(args.fpath_nc, args.field_name,
                       args.field_array, args.mask)

if __name__ == '__main__':
    main()
