import os
import nctoolkit as nc
import xarray as xr
import numpy as np
import subprocess

session = dict()


def get_extent(ff):
    # add docstring
    """ "
    Get the extent of a netcdf file

    Parameters
    ----------
    ff : str
        The path to the netcdf file

    Returns
    -------
    extent : list
        A list of the form [lon_min, lon_max, lat_min, lat_max]

    """

    ds = nc.open_data(ff)
    ds.subset(variables=ds.variables[0])
    ds.top()
    ds.subset(time=0)
    ds_xr = ds.to_xarray()
    lon_name = [x for x in ds_xr.coords if "lon" in x][0]
    lat_name = [x for x in ds_xr.coords if "lat" in x][0]
    lons = ds_xr[lon_name].values
    lons = lons.flatten()
    # as a unique, sorted list
    lons = list(set(lons))
    lats = ds_xr[lat_name].values
    lats = lats.flatten()
    # as a unique, sorted list
    lats = list(set(lats))
    lon_res = np.abs(lons[2] - lons[1])
    lat_res = np.abs(lats[2] - lats[1])
    ds = nc.open_data(ff)
    ds.subset(variables=ds.variables[0])
    ds.top()
    ds.tmax()
    ds.to_latlon(lon=[-180, 180], lat=[-90, 90], res=[lon_res, lat_res])
    df = ds.to_dataframe().dropna().reset_index()
    lon_min = df.lon.min()
    lon_max = df.lon.max()
    lat_min = df.lat.min()
    lat_max = df.lat.max()
    lons = [lon_min, lon_max]
    lats = [lat_min, lat_max]
    extent = [
        lons[0] - lon_res,
        lons[1] + lon_res,
        lats[0] - lat_res,
        lats[1] + lat_res,
    ]
    #
    return extent


def extension_of_directory(starting_directory, exclude=[]):
    levels = session["levels_down"]

    new_directory = ""
    for i in range(levels):
        new_directory = new_directory + "/**"
    return new_directory + "/"



def is_latlon(ff):
    ds = nc.open_data(ff, checks = False)

    cdo_result = subprocess.run(
            f"cdo griddes {ff}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return "lonlat" in cdo_result.stdout.decode("utf-8")

def get_resolution(ff):
    ds = nc.open_data(ff, checks = False)
    var = ds.variables[0]
    ds = xr.open_dataset(ff)
    lon_name = [x for x in list(ds.coords) if 'lon' in x][0]
    lat_name = [x for x in list(ds.coords) if 'lat' in x][0]
    var_dims =  ds[var].dims
    extent = get_extent(ff)
    if lon_name in var_dims:
        if lat_name in var_dims:
            n_lon = len(ds[lon_name])
            n_lat = len(ds[lat_name])
            lon_res = (extent[1] - extent[0]) / n_lon
            lat_res = (extent[3] - extent[2]) / n_lat
            return [lon_res, lat_res]
    else:
        # get the final two var_dims
        var_dims = var_dims[-2:]
        # lat should be the second
        n_lat = len(ds[var_dims[1]])
        n_lon = len(ds[var_dims[0]])
        lon_res = (extent[1] - extent[0]) / n_lon
        lat_res = (extent[3] - extent[2]) / n_lat
        return [lon_res, lat_res]
