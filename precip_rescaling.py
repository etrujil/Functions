#!/usr/bin/env python3

import argparse
import netCDF4 as nc
import numpy as np
from subprocess import run
from skimage import io


def snowfall_rescaling_timestep(precip_np, snow_fraction_np, snow_map_np):
    '''
    Parameters:
        Precip_np: 2D numpy array with the precipitation field for a given timestep
        snow_fraction_np: 2D numpy array with the snow fraction (0.0 - 1.0) field for a given timestep
        snow_map_np: 2D numpy array with the snow field to use for the rescaling
            of solid precip

        All input grids must have the same projection, extent, and resolution

    Return:
        precip_np_out: 2D numpy array with rescaled total precipitation (rain + snow) mass
                       for the time interval
        snow_fraction_np_out: 2D numpy array with the snow fraction (0.0 - 1.0)

    Note: This function assumes that there are not no-data flags, just values with a range between 0-inf
    '''

    precip_np_out = precip_np.copy()
    snow_fraction_np_out = snow_fraction_np.copy()
    return_flag = False  # True if rescaled

    if np.sum(snow_fraction_np) > 0:

        snowfall = precip_np * snow_fraction_np
        rainfall = precip_np - snowfall
        # The rainfall matrix does not change

        # Select the area that has a snowfall value for the time step
        snow_map_select = np.zeros(np.shape(snowfall), dtype=float)
        snow_map_select[(snowfall > 0) & (snow_map_np > 0)
                        ] = snow_map_np[(snowfall > 0) & (snow_map_np > 0)]
        # With the lines above we make sure that we only select the area with
        # positive values on both.
        # Also, this way we mantain the no-data flags in the original NetCDF4 file
        # remember both grids must have the same extent, dim., res.

        # won't enter if no area intersect between snowfall and snow_map_np
        if np.sum(snow_map_select > 0) > 0:

            mean_snowfall = np.mean(snowfall[snow_map_select > 0])
            # we use snow_map_select because we want the area with positive values in both snowfall and snow_map_np
            mean_snow_map = np.mean(
                snow_map_select[snow_map_select > 0])  # non-zero mean

            snowfall_rescaled = snowfall.copy()
            snowfall_rescaled[snow_map_select > 0] = snow_map_select[snow_map_select >
                                                                     0] * mean_snowfall / mean_snow_map
            # The line above only modifies areas with positive values in snow_map_select
            # which means it only modifies the areas where we know the snow distribution and we have snowfall

            # output
            precip_np_out = rainfall + snowfall_rescaled  # only the snowfall was modified
            snow_fraction_np_out[snowfall_rescaled > 0] = snowfall_rescaled[snowfall_rescaled >
                                                                            0] / precip_np_out[snowfall_rescaled > 0]
            return_flag = True

    return precip_np_out, snow_fraction_np_out, return_flag

# Function to rescale a full NetCDF4 precip and snow_percent (which should be fraction, nor percent) files
# using the 'snowfall_rescaling_timestep' function above


def snowfall_rescaling(precip_path, percent_snow_path, snow_map_path, precip_field='precip', snow_percent_field='percent_snow'):
    '''
    Parameters:
        precip_path: full path of NetCDF4 file with the precipitation information
        percent_snow_path: full path of NetCDF4 file with the snow fraction information
        precip_field: name of the precipitation field in the nc file (str)
        snow_percent_field: name of the field containing the
            percentage of precip in the form of snow (str)

        All input grids must have the same projection, extent, and resolution

    Return:
        precip_np_out: 3D numpy array with total precipitation (rain + snow) mass per time interval
        snow_fraction_np_out: 3D numpy array with fraction (0-1)

    Note: Assuming that there are not no-data flags, just values with a range between 0-inf
    '''

    run('cp ' + precip_path + ' precip_rescaled.nc', shell=True)
    run('cp ' + percent_snow_path + ' percent_snow_rescaled.nc', shell=True)

    snow_map_np = io.imread(snow_map_path)

    precip_nc = nc.Dataset('precip_rescaled.nc', 'r+')
    percent_nc = nc.Dataset('percent_snow_rescaled.nc', 'r+')

    precip_rescaled = np.zeros(
        precip_nc.variables[precip_field].shape, dtype=float)
    percent_rescaled = np.zeros(
        percent_nc.variables[snow_percent_field].shape, dtype=float)
    rescaling_flag = np.empty(
        precip_nc.variables[precip_field].shape[0], dtype=bool)

    for i_time in range(len(rescaling_flag)):

        precip_rescaled[i_time, :, :], percent_rescaled[i_time, :, :], rescaling_flag[i_time] = \
            snowfall_rescaling_timestep(precip_nc.variables[precip_field][i_time, :, :],
                                        percent_nc.variables[snow_percent_field][i_time, :, :],
                                        snow_map_np)

    precip_nc.variables[precip_field][:] = precip_rescaled
    percent_nc.variables[snow_percent_field][:] = percent_rescaled

    precip_nc.close()
    percent_nc.close()

    # return None
    # precip_rescaled, percent_rescaled, rescaling_flag


def main():
    p = argparse.ArgumentParser(description='Generates rescaled precipitation'
                                ' files for SMRF/AWSM')

    p.add_argument("-p", "--precip", dest="precip_path",
                   required=True, type=str,
                   help="Path to precip.nc")
    p.add_argument("-ps", "--percentsnow", dest="percent_snow_path",
                   required=True, type=str,
                   help="Path to percent_snow.nc")
    p.add_argument("-sm", "--snowmap", dest="snow_map_path",
                   required=True, type=str,
                   help="Path to snow map/snow distribution .tif file")
    p.add_argument("-pf", "--precipfield", dest="precip_field",
                   default="precip", type=str,
                   help="name of precip field in precip.nc"
                   ", default=precip")
    p.add_argument("-spf", "--percentfield", dest="snow_percent_field",
                   default="percent_snow", type=str,
                   help="name of the snow percent field in snow_percent.nc"
                   ", default=percent_snow")
    args = p.parse_args()

    snowfall_rescaling(args.precip_path, args.percent_snow_path,
                       args.snow_map_path, args.precip_field,
                       args.snow_percent_field)


if __name__ == '__main__':
    main()
