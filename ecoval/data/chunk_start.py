import nctoolkit as nc
nc.options(parallel=True)
import os
import pkg_resources
import pandas as pd

def fix_basename(x):
    #annualmean_nitrate_nsbc.nc
    # 3 part names like above need the final part removed
    x_fixed = x.split('_')
    if len(x_fixed) == 3:
        x_fixed = x_fixed[:-1] 
        return '_'.join(x_fixed) + "." + x.split('.')[1] 
    else:
        return x


def fix_variable_name(x):
    valid_vars = [
        "temperature", "salinity", "oxygen", "phosphate",
        "silicate", "nitrate", "ammonium", "alkalinity",
        "ph", "chlorophyll", "co2flux", "pco2",
        "doc", "poc", "carbon", "benbio",
        "benthic_carbon_flux", "mesozoo", "oxycons" ]
    if x.lower() == "temperature":
        return "temperature"
    x_lower = x.lower()
    if x_lower == "doc":
        x = "DOC concentration"
        return x
    if x_lower == "poc":
        x = "POC concentration"
        return x
    if x_lower == "oxygen":
        x = "oxygen concentration"
        return x
    if x_lower == "phosphate":
        x = "phosphate concentration"
        return x
    if x_lower == "silicate":
        x = "silicate concentration"
        return x
    if x_lower == "nitrate":
        x = "nitrate concentration"
        return x
    if x_lower == "ammonium":
        x = "ammonium concentration"
        return x
    if x_lower == "alkalinity":
        x = "alkalinity"
        return x
    if x_lower == "ph":
        x = "pH"
        return x
    if x_lower == "chlorophyll":
        x = "chlorophyll concentration"
        return x
    if x_lower == "co2flux":
        x = "air-sea carbon dioxide flux"
        return x
    if x_lower == "pco2":
        x = "pCO2"
        return x
    if x_lower == "carbon":
        x = "carbon concentration in sediments"
        return x
    if x_lower == "benbio":
        x = "macrobenthos biomass concentration"
        return x
    if x_lower == "benthic_carbon_flux":
        x = "carbon flux in sediments"
        return x
    if x_lower == "mesozoo":
        x = "mesozooplankton concentration"
        return x
    if x_lower == "oxycons":
        x = "benthic oxygen consumption"
        return x
    if x_lower == "temperature":
        x = "temperature"
        return x
    if x_lower == "salinity":
        x = "salinity"
        return x

    return x




def df_display(df):
    # only 2 decimal places
    df = df.round(2)
    # coerce numeric columns to str
    df = df.astype(str)
    # capitalize unit column name, if it exists
    if "unit" in df.columns:
        df = df.rename(columns={"unit": "Unit"})
    # capitalize variable column name, if it exists
    if "variable" in df.columns:
        df = df.rename(columns={"variable": "Variable"})
        # fix variable names
        df["Variable"] = df["Variable"].apply(fix_variable_name)
        # capitalize variable
        df["Variable"] = df["Variable"].str.capitalize()
        # ensure "Poc " is "POC "
        df["Variable"] = df["Variable"].str.replace("Poc ", "POC ")
        # ensure "Doc" is "DOC"
        df["Variable"] = df["Variable"].str.replace("Doc", "DOC")
    # convert nan to N/A
    df = df.replace("nan", "N/A")
    return df.style.hide(axis="index")

import warnings
warnings.filterwarnings('ignore')
from IPython.display import Markdown as md
%load_ext rpy2.ipython
import jellyfish
compact = False
import pickle
import ecoval
import copy
import calendar
import nctoolkit as nc
import hvplot.xarray
import geopandas as gpd
import xarray as xr
import holoviews as hv
from holteandtalley import HolteAndTalley
import numpy as np
#from mask import mask_all, mask_shelf
import pandas as pd
build = "book_build"
test_status = the_test_status
if build == "pdf":
    pd.set_option('styler.render.repr', 'latex')
chapter = "book_chapter"
if build == "pdf":
    chapter = f"{chapter}."
else:
    chapter = ""
import glob
def tidy_summary_paths(paths):
    """
    Function to tidy up the paths for the summary stats
    Right now this just ignores occci files if chlo and nsbc files are present
    """
    paths_new = paths
    paths_nsbc = [x for x in paths if "nsbc" in x]
    for ff in paths_nsbc:
        for ff_bad in glob.glob(ff.replace("nsbc", "**")):
            if "nsbc" not in os.path.basename(ff_bad):
                paths_new.remove(ff_bad)
    #if len([x for x in paths if "chlo" in x and "nsbc" in x]) > 0:
    #    paths_new = [x for x in paths if ("chlo" in x and "occci" in x) == False]
    return paths_new

import cmocean as cm
from tqdm import tqdm
import hvplot.pandas
from plotnine import *
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
i_figure = 1
i_table = 1
stamp = nc.session_info["stamp"]
out = ".trackers/" + stamp
fast_plot = fast_plot_value
if not os.path.exists(".trackers"):
    os.makedirs(".trackers")
# save out as empty file
with open(out, 'w') as f:
    f.write("")
nws = False
