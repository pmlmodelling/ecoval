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

def fix_unit(x, build = "book_build"):
    if build == "html":
        x = x.replace("/m^3", "m<sup>-3</sup>") 
        x = x.replace("/m3", "m<sup>-3</sup>")
        x = x.replace("m-3", "m<sup>-3</sup>")
        x = x.replace("/m^2", "m<sup>-2</sup>")
        x = x.replace("/m2", "m<sup>-2</sup>")
        x = x.replace("m2", "m<sup>2</sup>")
        x = x.replace("m3", "m<sup>3</sup>")
        #O_2
        x = x.replace("O_2", "O<sub>2</sub>")
        # fix CO2
        x = x.replace("CO2", "CO<sub>2</sub>")
        # fix /yr
        x = x.replace("/yr", "year<sup>-1</sup>")
        # degC
        x = x.replace("degC", "°C")
    if build == "pdf":
        x = x.replace("/m^3", "m$^{-3}$")
        x = x.replace("/m3", "m$^{-3}$")
        x = x.replace("m-3", "m$^{-3}$")
        x = x.replace("/m^2", "m$^{-2}$")
        x = x.replace("/m2", "m$^{-2}$")
        x = x.replace("m2", "m$^2$")
        x = x.replace("m3", "m$^3$")
        #O_2
        x = x.replace("O_2", "O$_2$")
        # fix CO2
        x = x.replace("CO2", "CO$_2$")
        # fix /yr
        x = x.replace("/yr", "year$^{-1}$")
        # fix /day
        x = x.replace("/day", "day$^{-1}$")
        # degC
        x = x.replace("degC", "$^\circ$C")


    return x

def fix_variable_name(x):
    if "concentration" in x:
        return x
    if "nitrate" in x:
        return "nitrate concentration"
    if "phosphate" in x:
        return "phosphate concentration"
    if "silicate" in x:
        return "silicate concentration"
    if "ammonium" in x:
        return "ammonium concentration"
    if "temperature" in x:
        return "temperature"
    if "salinity" in x:
        return "salinity"
    if "oxygen" in x and "enthic" not in x:
        return "oxygen concentration"
    if x =="ph":
        return "pH"
    if "chlorophyll" in x:
        return "chlorophyll concentration"

    return x




def df_display(df, build = "book_build"):
    # only 2 decimal places
    for col in df.columns:
        if "Number" in col:
            try:
                df[col] = df[col].astype(int)
                # put in commas
                df[col] = df[col].apply(lambda x: "{:,}".format(x))
            except:
                pass
    # if rmsd is a column title change to RMSD
    if "rmsd" in df.columns:
        df = df.rename(columns={"rmsd": "RMSD"})
    # round to 2 decimal places
    df = df.round(2)

    # coerce numeric columns to str
    # if number is in the column title, make sure the variable is int

    df = df.applymap(lambda x: f"{x:.2f}" if isinstance(x, float) else x)
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
    if "Variable" in df.columns:
        if build == "html":
            # pCO2
            df["Variable"] = df["Variable"].str.replace("pCO2", "pCO<sub>2</sub>")
        if build == "pdf":
            # pCO2
            df["Variable"] = df["Variable"].str.replace("pCO2", "pCO$_2$") 

    # if R2 is in the column name, make sure the 2 is superscript
    if "R2" in df.columns:
        df = df.rename(columns={"R2": "R²"}) 
    # convert nan to N/A
    df = df.replace("nan", "N/A")
    if "Region" in df.columns:
        if "Full Domain" in df["Region"].values:
            # ensure the Full Domain region is the first row
            # get the index of the row
            i_domain = df[df["Region"] == "Full Domain"].index[0]
            df1 = df.iloc[i_domain:i_domain+1]
            df2 = df.drop(i_domain)
            df = pd.concat([df1, df2]) 
            df = df.reset_index(drop = True)

    if "Unit" in df.columns:
        #format this appropriately. Markdown, superscripts etc.
        df["Unit"] = df["Unit"].apply(fix_unit)

    return df.style.hide(axis="index")

import warnings
warnings.filterwarnings('ignore')
from IPython.display import Markdown as md_markdown

def md_basic(x, build = "book_build"):
    x = x.replace(" 5 m ", " 5m ")
    x = x.replace(" 5 m.", " 5m.")

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
    # ensure there are spaces after commas, using regex
    x = re.sub(r",(\w)", r", \1", x)
    if build == "html":
        x = x.replace("CO", "CO<sub>2</sub>")
    if build == "pdf":
        x = x.replace("CO", "CO$_2$")
    # x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)

    # if "year" not in x.lower(): 
    #     if "period" not in x.lower(): 
    #         x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)
    return md_markdown(x)

def md(x, number = False, build = "book_build"):
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
    if "DOC conc" not in x:
        x = x.replace(" doc ", " DOC concentration ")
    if "DOC conc" not in x:
        x = x.replace(" doc.", " DOC concentration.")
    if "POC conc" not in x:
        x = x.replace(" poc ", " POC concentration ")
    if "POC conc" not in x:
        x = x.replace(" poc.", " POC concentration.")
    if "oxygen conc" not in x and "enthic" not in x:
        x = x.replace(" oxygen ", " oxygen concentration ")
    if "oxygen conc" not in x and "enthic" not in x:
        x = x.replace(" oxygen.", " oxygen concentration.")
    if "phosphate conc" not in x:
        x = x.replace(" phosphate ", " phosphate concentration ")
    if "phosphate conc" not in x:
        x = x.replace(" phosphate.", " phosphate concentration.")
    if "silicate conc" not in x:
        x = x.replace(" silicate ", " silicate concentration ")
    if "silicate conc" not in x:
        x = x.replace(" silicate.", " silicate concentration.")
    if "nitrate conc" not in x:
        x = x.replace(" nitrate ", " nitrate concentration ")
    if "nitrate conc" not in x:
        x = x.replace(" nitrate.", " nitrate concentration.")
    if "ammonium conc" not in x:
        x = x.replace(" ammonium ", " ammonium concentration ")
    if "ammonium conc" not in x:
        x = x.replace(" ammonium.", " ammonium concentration.")
    if "ph conc" not in x:
        x = x.replace(" ph ", " pH ")
    if "ph conc" not in x:
        x = x.replace(" ph.", " pH.")
    if "chlorophyll conc" not in x:
        x = x.replace(" chlorophyll ", " chlorophyll concentration ")
    if "chlorophyll conc" not in x:
        x = x.replace(" chlorophyll.", " chlorophyll concentration.")
    if "co2flux" not in x:
        x = x.replace(" co2flux ", " air-sea carbon dioxide flux ")
    if "co2flux" not in x:
        x = x.replace(" co2flux.", " air-sea carbon dioxide flux.")
    if "carbon" not in x:
        x = x.replace(" carbon ", " carbon concentration in sediments ")
    if "carbon" not in x:
        x = x.replace(" carbon.", " carbon concentration in sediments.")
    if "benbio" not in x:
        x = x.replace(" benbio ", " macrobenthos biomass concentration ")
    if "benbio" not in x:
        x = x.replace(" benbio.", " macrobenthos biomass concentration.")
    if "benthic_carbon_flux" not in x:
        x = x.replace(" benthic_carbon_flux ", " carbon flux in sediments ")
    if "benthic_carbon_flux" not in x:
        x = x.replace(" benthic_carbon_flux.", " carbon flux in sediments.")
    if "mesozoo" not in x:
        x = x.replace(" mesozoo ", " mesozooplankton concentration ")
    if "mesozoo" not in x:
        x = x.replace(" mesozoo.", " mesozooplankton concentration.")
    if "oxycons" not in x:
        x = x.replace(" oxycons ", " benthic oxygen consumption ")
    if "oxycons" not in x:
        x = x.replace(" oxycons.", " benthic oxygen consumption.")
    if "color" in x:
        x = x.replace(" color", " colour")
    if build == "html":
        if "pco2" in x:
            x = x.replace("pco2", "pCO<sub>2</sub>")
        # make CO2 subscript
        x = x.replace("CO2", "CO<sub>2</sub>")
        # fix O_2
        x = x.replace("O_2", "O<sub>2</sub>")
        x = x.replace(" ph ", " pH ")
        x = x.replace(" R2 ", " R<sup>2</sup> ")
        x = x.replace(" R2.", " R<sup>2</sup>.")
        # fix g/kg
        x = x.replace("g/kg", "g kg<sup>-1</sup>")

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
        # ensure there are spaces after commas, using regex
        x = re.sub(r",(\w)", r", \1", x)


        x = x.replace(" .", ".")
        x = x.replace(" ,", ",")
        x = x.replace(" :", ":")
        x = x.replace(" ;", ";")
        x = x.replace(" %", "%")
        # /m^3
        x = x.replace("/m3", "m<sup>-3</sup>")
        x = x.replace("/m^3", "m<sup>-3</sup>")
        # handle /m^2
        x = x.replace("/m2", "m<sup>-2</sup>")
        x = x.replace("/m^2", "m<sup>-2</sup>")
        # handl m-3
        x = x.replace(" m-3", " m<sup>-3</sup>")
        # fix /yr
        x = x.replace("/yr", "year<sup>-1</sup>")
        # fix /day
        x = x.replace("/day", "day<sup>-1</sup>")
        if "air-sea" in x:
            # no need for surface, it's redundant
            x = x.replace("surface", "")
            # ensure no double spaces
            x = x.replace("  ", " ")

        if number:
            if "year" not in x.lower(): 
                if "period" not in x.lower(): 
                    # do not use numbers between brackets ()
                    # x = re.sub(r"\((\d{1,3})(\d{3})\)", r"(\1,\2)", x)
                    x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)
                    # x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)
    if build == "pdf":
        # this is latex
        if "pco2" in x:
            x = x.replace("pco2", "pCO$_2$")
        # make CO2 subscript
        x = x.replace("CO2", "CO$_2$")
        # fix O_2
        x = x.replace("O_2", "O$_2$")
        x = x.replace(" ph ", " pH ")
        x = x.replace(" R2 ", " R$^2$ ")
        x = x.replace(" R2.", " R$^2$.")
        # fix g/kg
        x = x.replace("g/kg", "g kg$^{-1}$")
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
        # ensure there are spaces after commas, using regex
        x = re.sub(r",(\w)", r", \1", x)
        x = x.replace(" .", ".")
        x = x.replace(" ,", ",")
        x = x.replace(" :", ":")
        x = x.replace(" ;", ";")
        x = x.replace(" %", "%")
        # /m^3
        x = x.replace("/m3", "m$^{-3}$")
        x = x.replace("/m^3", "m$^{-3}$")
        # handle /m^2

        x = x.replace("/m2", "m$^{-2}$")
        x = x.replace("/m^2", "m$^{-2}$")
        # handl m-3
        x = x.replace(" m-3", " m$^{-3}$")
        # fix /yr
        x = x.replace("/yr", "year$^{-1}$")
        # fix /day
        x = x.replace("/day", "day$^{-1}$")
        if "air-sea" in x:
            # no need for surface, it's redundant
            x = x.replace("surface", "")
            # ensure no double spaces
            x = x.replace("  ", " ")
        if number:
            if "year" not in x.lower(): 
                if "period" not in x.lower(): 
                    # do not use numbers between brackets ()
                    # x = re.sub(r"\((\d{1,3})(\d{3})\)", r"(\1,\2)", x)
                    x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)
                    # x = re.sub(r"(\d{1,3})(\d{3})", r"\1,\2", x)




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
