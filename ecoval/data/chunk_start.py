import nctoolkit as nc
nc.options(parallel=True)
import os
import pkg_resources
import pandas as pd
import re

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
    return x




def df_display(df):
    # only 2 decimal places
    for col in df.columns:
        if "Number" in col:
            try:
                df[col] = df[col].astype(int)
                # put in commas
                df[col] = df[col].apply(lambda x: "{:,}".format(x))
            except:
                pass
    # if rsmd is a column title change to RSMD
    if "rsmd" in df.columns:
        df = df.rename(columns={"rsmd": "RSMD"})
    df = df.round(2)
    # coerce numeric columns to str
    # if number is in the column title, make sure the variable is int

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
    # if R2 is in the column name, make sure the 2 is superscript
    if "R2" in df.columns:
        df = df.rename(columns={"R2": "R²"}) 
    # convert nan to N/A
    df = df.replace("nan", "N/A")
    return df.style.hide(axis="index")

import warnings
warnings.filterwarnings('ignore')
from IPython.display import Markdown as md_markdown
def md(x):
    valid_vars = [
        "temperature", "salinity", "oxygen", "phosphate",
        "silicate", "nitrate", "ammonium", "alkalinity",
        "ph", "chlorophyll", "co2flux", "pco2",
        "doc", "poc", "carbon", "benbio",
        "benthic_carbon_flux", "mesozoo", "oxycons" ]
    x = x.replace(" 5 m ", " 5m ")
    x = x.replace(" 5 m.", " 5m.")
    if x.lower() == "temperature":
        return "temperature"
    x = x.replace(" doc ", " DOC concentration ")
    x = x.replace(" doc.", " DOC concentration.")
    x = x.replace(" poc ", " POC concentration ")
    x = x.replace(" poc.", " POC concentration.")
    x = x.replace(" oxygen ", " oxygen concentration ")
    x = x.replace(" oxygen.", " oxygen concentration.")
    x = x.replace(" phosphate ", " phosphate concentration ")
    x = x.replace(" phosphate.", " phosphate concentration.")
    x = x.replace(" silicate ", " silicate concentration ")
    x = x.replace(" silicate.", " silicate concentration.")
    x = x.replace(" nitrate ", " nitrate concentration ")
    x = x.replace(" nitrate.", " nitrate concentration.")
    x = x.replace(" ammonium ", " ammonium concentration ")
    x = x.replace(" ammonium.", " ammonium concentration.")
    x = x.replace(" ph ", " pH ")
    x = x.replace(" ph.", " pH.") 
    x = x.replace(" chlorophyll ", " chlorophyll concentration ")
    x = x.replace(" chlorophyll.", " chlorophyll concentration.")
    x = x.replace(" co2flux ", " air-sea carbon dioxide flux ")
    x = x.replace(" co2flux.", " air-sea carbon dioxide flux.")
    x = x.replace(" carbon ", " carbon concentration in sediments ")
    x = x.replace(" carbon.", " carbon concentration in sediments.")
    x = x.replace(" benbio ", " macrobenthos biomass concentration ")
    x = x.replace(" benbio.", " macrobenthos biomass concentration.")
    x = x.replace(" benthic_carbon_flux ", " carbon flux in sediments ")
    x = x.replace(" benthic_carbon_flux.", " carbon flux in sediments.")
    x = x.replace(" mesozoo ", " mesozooplankton concentration ")
    x = x.replace(" mesozoo.", " mesozooplankton concentration.")
    x = x.replace(" oxycons ", " benthic oxygen consumption ")
    x = x.replace(" oxycons.", " benthic oxygen consumption.")
    x = x.replace(" color", " colour")

    x = x.replace("pco2", "pCO<sub>2</sub>")
    # make CO2 subscript
    x = x.replace("CO2", "CO<sub>2</sub>")
    x = x.replace(" ph ", " pH ")
    x = x.replace(" R2 ", " R<sup>2</sup> ")
    x = x.replace(" R2.", " R<sup>2</sup>.")
    # get rid of double spaces
    x = x.replace("  ", " ")
    # use regex to ensure any numbers have commas
    if "**Figure" in x:
        # ensure the sentence ends with .
        if x[-1] != ".":
            x = x + "."
    if "**Table" in x:
        # ensure the sentence ends with .
        if x[-1] != ".":
            x = x + "."

    x = x.replace(" .", ".")
    x = x.replace(" ,", ",")
    x = x.replace(" :", ":")
    x = x.replace(" ;", ";")
    x = x.replace(" %", "%")
    # /m^3
    x = x.replace("/m3", "m<sup>-3</sup>")
    x = x.replace("/m^3", "m<sup>-3</sup>")

    if "year" not in x.lower(): 
        if "period" not in x.lower(): 
            x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)
    return md_markdown(x)
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
