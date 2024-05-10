import xarray as xr
import nctoolkit as nc


def mask_all(ds):


    try:

        data_path = pkg_resources.resource_filename("ecoval", "data/amm7_val_subdomains.nc")
        ds_regions = nc.open_data(data_path, checks = False)
        ds_xr = ds_regions.to_xarray()
        lon_size = len(ds_xr.lon)
        lat_size = len(ds_xr.lat)

        ds_xr = ds.to_xarray()
        lon_size_ds = len(ds_xr.lon)
        lat_size_ds = len(ds_xr.lat)

        if lon_size_ds == lon_size:
            if lat_size_ds == lat_size:
                ds_regions.subset(variables = ["Shelf", "Ocean"])
                ds_regions.sum_all()
                ds_regions.as_missing(0)
                ds * ds_regions
                ds.run()
    except Exception as e:
        x = "not maskable"



def mask_shelf(ds):
    try:

        data_path = pkg_resources.resource_filename("ecoval", "data/amm7_val_subdomains.nc")
        ds_regions = nc.open_data(data_path, checks = False)
        ds_xr = ds_regions.to_xarray()
        lon_size = len(ds_xr.lon)
        lat_size = len(ds_xr.lat)

        ds_xr = ds.to_xarray()
        lon_size_ds = len(ds_xr.lon)
        lat_size_ds = len(ds_xr.lat)

        if lon_size_ds == lon_size:
            if lat_size_ds == lat_size:
                ds_regions.subset(variables = ["Shelf"])
                ds_regions.sum_all()
                ds_regions.as_missing(0)
                ds * ds_regions
                ds.run()
    except Exception as e:
        x = "not maskable"


