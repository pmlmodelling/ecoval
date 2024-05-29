import nctoolkit as nc
import copy
import re
import glob
import multiprocessing
import os
import pandas as pd
import string
import random
import warnings
# A custom format for warnings.
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return "Warning: " + str(msg) + "\n"

warnings.formatwarning = custom_formatwarning
import numpy as np
import xarray as xr
from ecoval.session import session_info


import random
def extract_start(paths):
    random_files = random.sample(paths, 5)

    doable = True
    all_years = []
    all_indices = []
    while True:
        for x in random_files:
            print(x)
            ds = nc.open_data(x, checks = False)
            years = ds.years 
            if len(years) > 1:
                doable = False
                break
            if years[0] not in all_years:
                all_years += years
            # find indices where 2080 appear in the string x
            indices = [i for i in range(len(x)) if x.startswith(str(years[0]), i)]
            all_indices.append(indices)
            if len(all_years) > 4:
                break
        if len(all_years) > 4:
            break
    len_0 = len(all_indices[0])
    all_indices = list(set([item for sublist in all_indices for item in sublist]))
    if len(all_indices) == len_0:
        return all_indices[0]
    return None


from multiprocessing import Manager
session_warnings = Manager().list() 
from tqdm import tqdm
from ecoval.utils import session
from ecoval.utils import extension_of_directory
from ecoval.ices import generate_mapping
from ecoval.gridded import gridded_matchup
from ecoval.fixers import tidy_warnings

nc.options(parallel=True)
nc.options(progress=False)


def bin_value(x, bin_res):
    return np.floor((x + bin_res / 2) / bin_res + 0.5) * bin_res - bin_res / 2


ices_variables = [
    "temperature",
    "salinity",
    "oxygen",
    "chlorophyll",
    "phosphate",
    "silicate",
    "nitrate",
    "pH",
    "ammonium",
    "co2flux",
    "pco2",
    "doc",
    "poc",
    "alkalinity",
    "benbio"
]
data_dir = "/data/proteus1/scratch/rwi/evaldata/data/"


def find_config():
    # first look in the working directory
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


config_file = find_config()

if config_file is not None:
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
        data_path = terms[1].replace(" ", "")

        if os.path.exists(data_path):
            data_dir = data_path
        else:
            raise ValueError(f"{data_path} does not exist")


def matchup_wod(ff=None, variable=None, df_all=None, depths=None):
    data_dir = session_info["data_dir"] 
    ds = nc.open_data(ff, checks=False)
    years = list(set(ds.years))
    # read in the WOD data
    df_wod = []

    try:
        for year in years:
            ff_year = f"{data_dir}/wod/temp_csv/yearly/wod_temp_{year}.csv"
            df_wod.append(pd.read_csv(ff_year))
    except:
        return None

    df_wod = pd.concat(df_wod)

    df_locs = df_wod.drop(columns="temperature")
    df_wod = df_wod.rename(columns={"temperature": "observation"})

    depths_file = ds.levels

    with warnings.catch_warnings(record=True) as w:
        df_out = ds.match_points(
            df_locs, variables=variable, depths=depths_file, quiet=True
        )
    for ww in w:
        if str(ww.message) not in session_warnings:
            session_warnings.append(str(ww.message))

    df_out.rename(columns={variable: "model"}, inplace=True)
    df_out = df_out.merge(df_wod)

    df_all.append(df_out)




def get_out():
    # choose from all lowercase letter
    length = 8
    letters = string.ascii_lowercase
    result_str = "".join(random.choice(letters) for i in range(length))
    mapping = "mapping_" + result_str + ".csv"
    return mapping


def get_res(x, folder=None):
    if "_1d_" in x:
        return "d"
    if "_1m_" in x:
        return "m"

    final_extension = extension_of_directory(folder)

    path = glob.glob(folder + final_extension + x)[0]

    ds = nc.open_data(path)
    ds_times = ds.times
    months = [x.month for x in ds_times]
    days = [x.day for x in ds_times]
    years = [x.year for x in ds_times]
    df_times = pd.DataFrame({"month": months, "day": days, "year": years})

    n1 = len(
        df_times.loc[:, ["month", "year"]].drop_duplicates().reset_index(drop=True)
    )
    n2 = len(df_times)
    if n1 == n2:
        return "m"
    else:
        return "d"


def find_paths(folder, fvcom=False, exclude=[]):
    i = 1
    while True:

        levels = session["levels"]

        new_directory = folder + "/"
        for i in range(levels):
            dir_glob = glob.glob(new_directory + "/**")
            # randomize dir_glob
            import random
            random.shuffle(dir_glob)
            for x in dir_glob:
                # figure out if the the base directory is an integer
                try:
                    y = int(os.path.basename(x))
                    new_directory = x + "/"
                except:
                    blah = "blah"
        options = glob.glob(new_directory + "/**.nc")
        if not fvcom:
            options = [x for x in options if "part" not in os.path.basename(x)]
            options = [x for x in options if "restart" not in os.path.basename(x)]
        if fvcom:
            options = [x for x in options if "restart" not in os.path.basename(x)]

        if len([x for x in options if ".nc" in x]) > 0:
            break
        if fvcom:
            if len(options) == 1:
                break

    all_df = []
    print("********************************")
    print("Identifying variables in model output")
    print("********************************")

    # remove any files from options if parts of exclude are in them
    for exc in exclude:
        options = [x for x in options if f"{exc}" not in os.path.basename(x)]



    for ff in options:
        ds = nc.open_data(ff, checks=False)
        stop = True
        try:
            ds_dict = generate_mapping(ds, fvcom=fvcom)
            stop = False
        # output error and ff
        except:
            pass
        if stop:
            continue

        ds_vars = ds.variables
        # vosaline and votemper are special cases

        if "vosaline" in ds_vars:
            if ds_dict["salinity"] is None:
                ds_dict["salinity"] = "vosaline"

        if "votemper" in ds_vars:
            if ds_dict["temperature"] is None:
                ds_dict["temperature"] = "votemper"

        if len([x for x in ds_dict.values() if x is not None]) > 0:
            new_name = ""
            for x in os.path.basename(ff).split("_"):
                try:
                    y = int(x)
                    if len(new_name) > 0:
                        new_name = new_name + "_**"
                    else:
                        new_name = new_name + "**"
                except:
                    if len(new_name) > 0:
                        new_name = new_name + "_" + x
                    else:
                        new_name = x
            # replace integers in new_name with **

            new_dict = dict()
            for key in ds_dict:
                if ds_dict[key] is not None:
                    new_dict[ds_dict[key]] = [key]
            # new_name. Replace numbers between _ with **

            # replace integers with 4 or more digits with **
            new_name = re.sub(r"\d{4,}", "**", new_name)

            all_df.append(
                pd.DataFrame.from_dict(new_dict).melt().assign(pattern=new_name)
            )
    all_df = pd.concat(all_df).reset_index(drop=True)

    if fvcom is False:
        all_df["resolution"] = [get_res(x, folder) for x in all_df.pattern]
    else:
        all_df["resolution"] = "d"

    all_df = (
        all_df.sort_values("resolution").groupby("value").head(1).reset_index(drop=True)
    )
    all_df = all_df.rename(columns={"variable": "model_variable"})
    all_df = all_df.rename(columns={"value": "variable"})
    all_df = all_df.drop(columns="resolution")
    all_df = all_df.loc[:, ["variable", "model_variable", "pattern"]]

    print("********************************")
    print("Automatic variable identification complete")
    print("********************************")
    return all_df


def trends(
    folder=None,
    spinup=None,
    start=None,
    end=None,
    surface_level=None,
    surface = [
        "temperature",
        "salinity",
        "oxygen",
        "phosphate",
        "silicate",
        "nitrate",
        "ammonium",
        "alkalinity",
        "ph",
        "chlorophyll",
        "doc",
        "pco2",
        "co2flux",
        "poc",
    ],
    bottom= ["ph", "oxygen"],
    benthic = ["carbon", "benbio"],
    cores=None,
    e3t=None,
    mapping=None,
    exclude=[],
    levels=2,
    lon_lim = None,
    lat_lim = None, 
    fixed_format = True,
    **kwargs,
):
    """
    Match up model with observational data

    Parameters
    -------------
    folder: str
        Folder containing model output
    spinup: int
        Number of years to view as spinup. Default is None, which means the start year is used.
    surface_level: str
        Surface level of the model netCDF files. Either 'top' or 'bottom'. Default is None, so this must be supplied.
    surface: list
        List of surface variables. Internally, ecoval will decide which variables should matchup with gridded or point observations.
        Potential surface variables are temperature, salinity, oxygen, phosphate, silicate, nitrate, ammonium, alkalinity, ph, chlorophyll, doc, poc.
        If you want finer control, you can provide a dictionary with keys 'gridded' and 'point'.
        So surface = {'gridded': ['temperature'], 'point': ['ph', 'poc']} would matchup temperature with gridded data and pH and POC with point data.   
    bottom: list
        List of bottom variables to matchup with observational data.
        Potential bottom variables are temperature, salinity, oxygen, phosphate, silicate, nitrate, ammonium, alkalinity, ph, chlorophyll, doc, poc.    
    benthic: list
        List of benthic variables.
        Potential benthic variables are carbon. 
    cores: int
        Number of cores to use.
    e3t: str
        Path to e3t file, i.e. cell thickness. This only needs to be supplied if the variable is missing from the raw data.
    mapping: str
        Path to mapping file. This is a csv. A starting point can be generated by running `matchup` and saying you are not happy with the matchups.
    start: int
        Start year. First year of the simulations to matchup.
    end: int
        End year. Final year of the simulations to matchup.
    levels: int
        Number of levels down to look for netCDF files. Default is 2, ie. the files are of the format **/**/*.nc.
    exclude: list
        List of strings to exclude. This is useful if you have files in the directory that you do not want to include in the matchup.
    fvcom: bool
        Is the model fvcom? Default is False.
    strict: bool
        Strict temporal matching. Default is True. Set this to False if you want to match data up by month and day, but not year.
    mld: bool
        Matchup temperature data for mixed layer depth validation. Default is False.
    daily_match: bool
        Matchup data daily. Default is True. Set this to False if you want to match data up monthly.
    lon_lim: list
        List of two floats. The first is the minimum longitude, the second is the maximum longitude. Default is None.
    lat_lim: list
        List of two floats. The first is the minimum latitude, the second is the maximum latitude. Default is None.
    data_dir: str
        Path to data directory. Default is 'default'. If 'default', the data directory is taken from the session_info dictionary.
    
    kwargs: dict
        Additional arguments



    """

    # check that lon_lim and lat_lim and valid when either is not None

    if lon_lim is not None or lat_lim is not None:
        # check both are lists
        if not isinstance(lon_lim, list) or not isinstance(lat_lim, list):
            raise ValueError("lon_lim and lat_lim must be lists")

    if isinstance(benthic, str):
        benthic = [benthic]

    point_surface = []
    surf_dict = False

    surf_default = True
    if isinstance(surface, str):
        surface = [surface]
        surface = {"gridded": surface, "point": []}
    surf_all = False
    if isinstance(surface, list):
        surface = {"gridded": surface, "point": []}
        surf_all = True
    if isinstance(surface, dict):
        # throw error if gridded and point not in surface
        if "gridded" not in surface and "point" not in surface:
            raise ValueError("Please provide gridded or point variables")
        
        if "gridded" not in surface:
            surface["gridded"] = []
        else:
            if "point" not in surface:
                surface["point"] = []

        point_surface = surface["point"]
        surface = surface["gridded"]
        if isinstance(surface, str):
            surface = [surface]
        if isinstance(point_surface, str):
            point_surface = [point_surface]
    
    if len(surface) > 0:
        surf_default = False
    if len(point_surface) > 0:
        surf_default = False

    # fix benthic if something like "benthic biomass" is an element
    for pp in benthic:
        if "ben" in pp and "bio" in pp:
            benthic.remove(pp)
            benthic.append("benbio")

    session["levels"] = levels

    if isinstance(exclude, str):
        exclude = [exclude]

    # check if the folder exists
    if folder is None:
        raise ValueError("Please provide a folder")

    if not os.path.exists(folder):
        raise ValueError(f"{folder} does not exist")

    # loop through kwargs, if first three characters match arg and arg is None, set arg to value

    with open("matchup_report.md", "w") as f:
        pass

    def write_report(x):
        # append x to report
        with open("matchup_report.md", "a") as f:
            f.write(x + "\n")
            # add blank line
            f.write("\n")

    # add title
    write_report("## Summary of matchups")

    # convert matchup_report_md to pdf

    if surface_level is None:
        raise ValueError(
            "You need to specify if the surface is the top or the bottom level"
        )

    if surface_level not in ["top", "bottom"]:
        raise ValueError("surface_level must be top or bottom")

    sim_start = -1000
    sim_end = 10000
    for key in kwargs:
        if key[:3] == "fol":
            if folder is None:
                folder = kwargs[key]
        if key[:3] == "spi":
            if spinup is None:
                spinup = kwargs[key]
        if key[:3] == "map":
            if mapping is None:
                mapping = kwargs[key]
        if key[:3] == "cor":
            if cores is None:
                cores = kwargs[key]

    if start is not None:
        if spinup is not None:
            raise ValueError("You can only provide one of start or spinup")

    if end is not None:
        sim_end = end

    if start is not None:
        sim_start = start

    # check validity of variables chosen

    valid_vars = [
        "temperature",
        "salinity",
        "oxygen",
        "phosphate",
        "silicate",
        "nitrate",
        "ammonium",
        "alkalinity",
        "ph",
        "chlorophyll",
        "co2flux",
        "pco2",
        "doc",
        "poc",
        "carbon"
        "benbio"
    ]

    # type check for spinup

    if not isinstance(spinup, int):
        if start is None:
            raise ValueError("Please set spinup to int")

    if not isinstance(cores, int):
        raise ValueError("Please set cores to int")
    nc.options(cores=cores)

    all_df = None
    if isinstance(mapping, pd.DataFrame):
        all_df = mapping
    if isinstance(mapping, str):
        all_df = pd.read_csv(mapping)

    # create lists for working out which variables are needed for point matchups
    point_all = []
    point_bottom = []
    point_benthic = benthic
    
    if isinstance(bottom, str):
        bottom = [bottom]
    if bottom is None:
        bottom = []

    var_choice = surface + bottom + point_surface
    var_choice = list(set(var_choice))
    for vv in var_choice:
        if vv not in valid_vars and vv != "all":
            # suggest another variable based on similarity to valid_vars
            from difflib import get_close_matches

            close = get_close_matches(vv, valid_vars)
            if len(close) > 0:
                raise ValueError(
                    f"{vv} is not a valid variable. Did you mean {close[0]}?"
                )
            raise ValueError(
                f"{vv} is not a valid variable. Please choose from {valid_vars}"
            )
    if len(bottom) > 0:
        if bottom != "all":
            point_bottom = bottom

    if all_df is None:
        all_df = find_paths(folder, fvcom=False, exclude=exclude)

        # add in anything that is missing
        all_vars = [
            "temperature",
            "salinity",
            "oxygen",
            "chlorophyll",
            "phosphate",
            "silicate",
            "nitrate",
            "ph",
            "ammonium",
            "alkalinity",
            "co2flux",
            "pco2",
        ]
        missing_df = pd.DataFrame({"variable": all_vars}).assign(
            model_variable=None, pattern=None
        )

        all_df = (
            pd.concat([all_df, missing_df])
            .groupby("variable")
            .head(1)
            .reset_index(drop=True)
        )
        # add in poc
        df_poc = all_df.query("variable == 'chlorophyll'").reset_index(drop=True)
        if df_poc.model_variable[0] is not None:
            poc_mapping = df_poc.model_variable[0]
            # replace Chl with c
            poc_mapping = poc_mapping.replace("Chl", "c")
            poc_mapping = poc_mapping + "+Z5_c+Z6_c+R4_c+R6_c+R8_c"
            df_poc["model_variable"] = [poc_mapping]
            df_poc["variable"] = ["poc"]

            all_df = pd.concat([all_df, df_poc]).reset_index(drop=True)

        print(f"** Inferred mapping of model variable names from {folder}")
        print(all_df)
        print("Are you happy with these matchups? Y/N")
        x = input()

        if x.lower() not in ["y", "n"]:
            print("Provide Y or N")
            x = input()

        if x.lower() == "n":
            out = get_out()
            print(f"Inferred mapping saved as {out}")
            all_df.to_csv(out, index=False)
            return None

    out = "matched/mapping.csv"
    if not os.path.exists("matched"):
        os.mkdir("matched")
    df_out = all_df.dropna().reset_index(drop=True)
    final_extension = extension_of_directory(folder)
    df_out["pattern"] = [folder + final_extension + x for x in df_out.pattern]
    df_out.to_csv(out, index=False)


    vv = "temperature"

    model_variable = all_df.query(f"variable == '{vv}'").model_variable[0]
    model_pattern = all_df.query(f"variable == '{vv}'").pattern[0]
    print(model_variable)
    print(model_pattern)

    df = all_df.query("variable == @vv").reset_index(drop=True)
    mapping = dict()
    if len(df) > 0:
        mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

        selection = []
        try:
            selection += mapping[vv].split("+")
        except:
            selection = selection

        patterns = set(df.pattern)
        if len(patterns) > 1:
            raise ValueError(
                "Something strange going on in the string patterns. Unable to handle this. Bug fix time!"
            )
        pattern = list(patterns)[0]

        final_extension = extension_of_directory(folder)
        paths = glob.glob(folder + final_extension + pattern)

        for exc in exclude:
            paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

        new_paths = []
        ds = nc.open_data(paths, checks=False)
        # set up model_grid if it doesn't exist

        # This really should be a function....
        if not os.path.exists("matched/model_grid.csv"):
            ds_grid = nc.open_data(paths[0], checks=False)
            var = ds_grid.variables[0]
            ds_grid.subset(variables=var)
            if surface_level == "top":
                ds_grid.top()
            else:
                ds_grid.bottom()
            ds_grid.as_missing(0)
            # ds = nc.open_data(paths[0], checks=False)
            amm7 = False
            if max(ds_grid.contents.npoints) == 111375:
                amm7_out = "matched/amm7.txt"
                # create empty file
                with open(amm7_out, "w") as f:
                    f.write("")

                data_dir = "/data/proteus1/scratch/rwi/evaldata/data"

                ff_grid = f"{data_dir}/amm7_val_subdomains.nc"
                ds_grid.cdo_command(f"setgrid,{ff_grid}")
                # ds_grid.fix_amm7_grid()
                amm7 = True
            df_grid = ds_grid.to_dataframe().reset_index().dropna()
            columns = [x for x in df_grid.columns if "lon" in x or "lat" in x]
            df_grid = df_grid.loc[:, columns].drop_duplicates()
            if not os.path.exists("matched"):
                os.makedirs("matched")
            df_grid.to_csv("matched/model_grid.csv", index=False)

        yy_start = None
        all_years = []
        if fixed_format:
            yy_start = extract_start(paths)

        if yy_start is not None:
            for ff in tqdm(paths):
                ds_years = [int(ff[yy_start : yy_start + 4])]
                all_years += ds_years
        else:
            for ff in paths:
                try:
                    ds = nc.open_data(ff, checks=False)
                    ds_years = ds.years
                    all_years += ds_years
                except:
                    print(f"Unable to find relevant years in  {ff}")
        all_years = list(set(all_years)) 

        if spinup is not None:
            years = [x for x in all_years if x >= min(all_years) + spinup]

        sim_years = range(sim_start, sim_end + 1)
        if start is not None:
            years = [x for x in all_years if x in sim_years]

        yy_start = None
        if fixed_format:
            yy_start = extract_start(paths)

        if yy_start is not None:
            for ff in tqdm(paths):
                print(ff)
                ds_years = [int(ff[yy_start : yy_start + 4])]
                if len([x for x in ds_years if x in years]) > 0:
                    new_paths.append(ff)
        else:
            for ff in paths:
                print(ff)
                try:
                    ds = nc.open_data(ff, checks=False)
                    ds_years = ds.years
                    if len([x for x in ds_years if x in years]) > 0:
                        new_paths.append(ff)
                except:
                    print(f"Unable to find relevant years in  {ff}")

        paths = list(set(new_paths))
        paths.sort()

        ds = nc.open_data(new_paths, checks = False)
        ds.subset(variable = model_variable)
        if surface_level == "top":
            ds.top()
        else:
            ds.bottom()
        ds.as_missing(0)
        
        ds.tmean(["year", "month"])
        ds.merge("time")
        ds.tmean("year")
        #ds.rolling_mean(20)
        out = f"trends/{vv}.nc"
        # check if directory exists for out
        if os.path.exists(out):
            os.remove(out)
        if not os.path.exists("trends"):
            os.makedirs("trends")
        ds.to_nc(out, zip = True)


    print(paths)
