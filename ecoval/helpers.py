import glob
import os
import warnings
import xarray as xr

def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
