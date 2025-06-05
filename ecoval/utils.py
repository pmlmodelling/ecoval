import os
import nctoolkit as nc
import warnings
import xarray as xr
import numpy as np
import subprocess
import pandas as pd
from tqdm import tqdm
from ecoval.session import session_info



def bin_value(x, bin_res):
    return np.floor((x + bin_res / 2) / bin_res + 0.5) * bin_res - bin_res / 2



def extension_of_directory(starting_directory, exclude=[]):
    levels = session_info["levels_down"]

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
    ds.as_missing(0)
    ds_xr = ds.to_xarray()
    lon_name = [x for x in ds_xr.coords if "lon" in x][0]
    lat_name = [x for x in ds_xr.coords if "lat" in x][0]

    df = ds_xr.to_dataframe().reset_index()
    df = ds.to_dataframe().dropna().reset_index()
    df = df.rename(columns={lon_name: "lon", lat_name: "lat"})
    df = df.dropna()
    lon_min = float(df.lon.min())
    lon_max = float(df.lon.max())
    lat_min = float(df.lat.min())
    lat_max = float(df.lat.max())
    return [lon_min, lon_max, lat_min, lat_max]


    lons = ds_xr[lon_name].values
    lons = lons.flatten()
    # as a unique, sorted list
    lons = list(set(lons))
    lats = ds_xr[lat_name].values
    lats = lats.flatten()
    # as a unique, sorted list
    lats = list(set(lats))
    lats.sort()
    lons.sort()
    lon_res = np.abs(lons[2] - lons[1])
    lat_res = np.abs(lats[2] - lats[1])
    print(lon_res)
    print(lat_res)
    ds = nc.open_data(ff)
    ds.subset(variables=ds.variables[0])
    ds.top()
    ds.tmax()
    ds.to_latlon(lon=[-180, 180], lat=[-90, 90], res=[lon_res, lat_res])
    # rename
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
            lon_max = float(ds[lon_name].max())
            lon_min = float(ds[lon_name].min())
            lat_max = float(ds[lat_name].max())
            lat_min = float(ds[lat_name].min())
            lon_res = (lon_max - lon_min) / (n_lon - 1)
            lat_res = (lat_max - lat_min) / (n_lat - 1)
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




def fvcom_regrid(ff = None, new_grid = None, vv = None, lons = None, lats = None, res = None):
    # multiple = False
    # if len(vv) > 1:
    #     ds_all = nc.open_data()
    #     multiple = True
    # vars = vv
    # del vv
    # for vv in vars:
    if True:
        # print(vv)
        drop_variables = ["siglay", "siglev"]
        ds_xr = xr.open_dataset( ff, drop_variables=drop_variables, decode_times=False)
        ds1 = nc.from_xarray(ds_xr[vv])
        lon = ds1.to_xarray().lon.values
        lat = ds1.to_xarray().lat.values

        lon_min = float(lon.min())
        lon_max = float(lon.max())
        lat_min = float(lat.min())
        lat_max = float(lat.max())
        # handle longitudes over 180 appropriately
        if lon_min > 180:
            lon_min = lon_min - 360
        if lon_max > 180:
            lon_max = lon_max - 360
        extent = [lon_min, lon_max, lat_min, lat_max]
        session_info["extent"] = extent

        ds1.run()
        ds1.nco_command( "ncks -d siglay,0,0")
        ds_xr = ds1.to_xarray()
        try:
            ds_xr = ds_xr.squeeze("siglay")
        except:
            pass
        ds1 = nc.from_xarray(ds_xr)
        ds1.subset(variable = vv)
        #ds1.tmean(align = "left")
        grid = pd.DataFrame({"lon": lon, "lat": lat})
        ds1.run()
        out_grid = nc.generate_grid.generate_grid(grid)
        nc.session.append_safe(out_grid)

        ds2 = ds1.copy()
        ds2.run()
        ds2.cdo_command(f"setgrid,{out_grid}")
        if lons is not None:
            ds2.to_latlon(lon = lons, lat = lats, res = res, method = "nn")
        else:
            ds2.regrid(new_grid, method = "nn")


        ds2.as_missing(0)
        ds2.run()
        df_mask = grid.assign(value=1)
        bin_res = list(get_resolution(ds2[0]))
        df_mask["lon"] = bin_value(df_mask["lon"], 0.25)
        df_mask["lat"] = bin_value(df_mask["lat"], 0.25)
        df_mask = df_mask.groupby(["lon", "lat"]).sum().reset_index()
        df_mask = df_mask.set_index(["lat", "lon"])
        ds_mask = nc.from_xarray(df_mask.to_xarray())
        os.system(f"cdo griddes {ds_mask[0]} > /tmp/mygrid")
        # open the text file text.txt and replace the string "generic" with "lonlat"
        with open("/tmp/mygrid", "r") as f:
            lines = f.readlines()

        # write line by line to /tmp/newgrid
        with open("/tmp/newgrid", "w") as f:
            for ll in lines:
                f.write(ll.replace("generic", "lonlat"))

        ds_mask.cdo_command(f"setgrid,/tmp/newgrid")
        # ds_mask.to_nc("/tmp/mask.nc")
        ds_mask.regrid(ds2, method = "con")
        # ds_mask.to_nc("/tmp/mask1.nc")
        ds_mask > 0
        # ds_mask.to_nc("/tmp/mask2.nc")
        os.remove("/tmp/mygrid")
        os.remove("/tmp/newgrid")
        ds_mask.as_missing(0)
        ds_mask.set_fill(-9999)
        ds2 * ds_mask

        # if multiple is False:
        return ds2
    #     else:
    #         ds_all.append(ds2)
    # ds_all.merge("variables")
    # ds_all.run()
    # return ds_all





