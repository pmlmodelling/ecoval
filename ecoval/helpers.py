import glob
import os
import copy
import warnings
import xarray as xr
from ecoval.session import session_info


def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


def list_observational_data(domain = "nws", obs_dir = "default"):

    # first do the gridded_surface
    if obs_dir == "default":
        obs_dir = session_info["obs_dir"]
        print(obs_dir)
    print(session_info)
    valid_surface = []
    valid_surface += [ os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/global/*") ]
    if domain == "nws":
        valid_surface += [ os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/nws/*") ]
    gridded_surface = copy.deepcopy(valid_surface)
    # now do the point_surface
    if domain == "nws":
        bottom = [ os.path.basename(x) for x in glob.glob(obs_dir + "point/nws/all/*") ]
    else:
        valid_surface = []
    point_surface = copy.deepcopy(valid_surface)
    if domain == "nws":
        valid_surface += [ os.path.basename(x) for x in glob.glob(obs_dir + "/point/nws/*") ]
    surface = {"gridded": gridded_surface, "point": point_surface} 
    # now do benthic
    if domain == "nws":
        benthic = [ os.path.basename(x) for x in glob.glob(obs_dir + "point/nws/benthic/*") ]
    else:
        benthic = [] 
    # now do bottom
    if domain == "nws":
        bottom = [ os.path.basename(x) for x in glob.glob(obs_dir + "point/nws/bottom/*") ]
    else:
        bottom = [] 

    print("*****************")
    print("Valid surface variables:")
    print(surface)
    print("*****************")
    print("Valid benthic variables:")
    print(benthic)
    print("*****************")
    print("Valid bottom variables:")
    print(bottom)



