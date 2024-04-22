

import os
import glob


session = dict()

import nctoolkit as nc
import numpy as np
def get_extent(ff):
    # add docstring
    """"
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
    ds.top()
    ds.subset(time = 0)
    ds.subset(variables = ds.variables[0])
    ds_xr = ds.to_xarray()
    lon_name = [x for x in ds_xr.coords if 'lon' in x][0]
    lat_name = [x for x in ds_xr.coords if 'lat' in x][0]
    lons = ds_xr[lon_name].values
    # as a unique, sorted list
    lons = list(set(lons))
    lats = ds_xr[lat_name].values
    # as a unique, sorted list
    lats = list(set(lats))
    lon_res = np.abs(lons[2] - lons[1])
    lat_res = np.abs(lats[2] - lats[1])
    ds = nc.open_data(ff)
    ds.top()
    ds.tmax()
    ds.to_latlon(lon = [-180, 180], lat = [-90, 90], res = [lon_res, lat_res])
    df = ds.to_dataframe().dropna().reset_index()
    lon_min = df.lon.min()
    lon_max = df.lon.max()
    lat_min = df.lat.min()
    lat_max = df.lat.max()
    lons = [lon_min, lon_max]
    lats = [lat_min, lat_max]
    extent = [lons[0]-lon_res, lons[1]+lon_res, lats[0] -lat_res, lats[1] + lat_res]
    # 
    return extent


def find_config(level = 0):
    # first look in the working directory
    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(ff):
            return ff

    if level == -2:
        for ff in ["../../.ecovalrc", "../../ecovalrc"]:
            if os.path.exists(ff):
                return ff

    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(ff):
            return ff

    # now look in the home directory....
    from os.path import expanduser

    home = expanduser("~")
    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(home + "/" + ff):
            return home + "/" + ff

    return None

def get_datadir(level = 0):

    data_dir = "/data/proteus1/scratch/rwi/evaldata/data/"
    config_file = find_config(level = level)

    if config_file is not None:
        # valid_keys = ["thread_safe", "lazy", "cores", "precision", "temp_dir"]

        file1 = open(config_file, "r")
        Lines = file1.readlines()

        count = 0
        # Strips the newline character
        for line in Lines:
            text = line.replace(" ", "").strip()
            if text.count(":") != 1:
                if len(text) > 0:
                    raise ValueError(f"Line in {config_file} is invalid: {line}")

        for line in Lines:
            text = line.replace(" ", "").strip()
            if len(text) > 0:
                terms = text.split(":")
                key = terms[0]
                value = None
                # print(key)
            data_path = terms[1].replace(" ", "")

            if os.path.exists(data_path):
                data_dir = data_path
            else:
                raise ValueError(f"{data_path} does not exist")
    return data_dir
    

def extension_of_directory(starting_directory, exclude = []):
    levels = session["levels"]

    new_directory = "" 
    for i in range(levels):
        new_directory = new_directory + "/**"
    return new_directory + "/"

    n_remove = len(starting_directory)

    # n_max = 0
    # for root, directories, files in os.walk(starting_directory):
    #     r_n = (len(root[n_remove:].split("/")))
    #     paths = root[n_remove+1:].split("/")
    #     paths = [x for x in paths if len(x) > 0]
    #     r_n = len(paths)
    #     if r_n > n_max:
    #         n_max = r_n 
    #     checks =  glob.glob(root +   "/*.nc")
    #     # remove anything in exclude from checks
    #     for exc in exclude:
    #         checks = [x for x in checks if exc not in x]
    #     checks = [x for x in checks if "part" not in os.path.basename(x)]
    #     checks = [x for x in checks if "restart" not in os.path.basename(x)]
    #     if len(checks) > 0:
    #         n_max = r_n
    #         break
    ## repeat "**" n_max
    return "/" + "/".join(["**" for i in range(n_max)]) + "/"