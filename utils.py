import ncplot
import xarray as xr
def grid_resolution(ff):
    """
    Function to get the grid resolution. This probably only works on nemo files right now.
    Likely enough for the tools' purpose
    """
    ds = xr.open_dataset(ff)
    dims = ncplot.utils.get_dims(ff)
    lon_name = dims.longitude.values[0]
    lat_name = dims.latitude.values[0]
    lon_range = ds[lon_name].max() - ds[lon_name].min()
    lat_range = ds[lat_name].max() - ds[lat_name].min()
    lon_size = len(ds[lon_name].x)
    lat_size = len(ds[lat_name].x)
    res = [float((lon_range / lon_size).values), float((lat_range / lat_size).values)]
    return res
