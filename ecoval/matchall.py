import copy
import nctoolkit as nc
import re
import glob
import multiprocessing
import pathlib
import os
import pandas as pd
import string
import random
import warnings
import pickle
import numpy as np
import xarray as xr
from ecoval.session import session_info
from multiprocessing import Manager
from tqdm import tqdm
from ecoval.utils import session
from ecoval.utils import extension_of_directory
from ecoval.parsers import generate_mapping
from ecoval.gridded import gridded_matchup


# A custom format for warnings.
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return "Warning: " + str(msg) + "\n"


warnings.formatwarning = custom_formatwarning

session_warnings = Manager().list()

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
    "benbio",
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


def mm_match(
    ff,
    ersem_variable,
    df,
    df_times,
    ds_depths,
    variable,
    df_all,
    top_layer=False,
    bottom_layer=False,
):
    """
    Parameters
    -------------
    ff: str
        Path to file
    ersem_variable: str
        Variable name in ERSEM
    df: pd.DataFrame
        Dataframe of observational data
    df_times: pd.DataFrame
        Dataframe of observational data with time information
    ds_depths: list
        Depths to match

    """

    if ds_depths is not None:
        nc.session.append_safe(ds_depths[0])
    try:
        with warnings.catch_warnings(record=True) as w:
            ds = nc.open_data(ff, checks=False)
            var_match = ersem_variable.split("+")
            ds.subset(variables=var_match)
            if top_layer:
                ds.top()
            if bottom_layer:
                ds.bottom()
            ds.as_missing(0)
            ds.run()
            if variable != "pft":
                if len(var_match) > 1:
                    ds.sum_all()
            valid_locs = ["lon", "lat", "year", "month", "day", "depth"]
            valid_locs = [x for x in valid_locs if x in df.columns]

            valid_times = (
                "year" in df.columns or "month" in df.columns or "day" in df.columns
            )

            if valid_times:
                df_locs = (
                    df_times.query("path == @ff")
                    .merge(df)
                    .loc[:, valid_locs]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )
            else:
                df_locs = df.loc[:, valid_locs]

            # idenify if the files have data from multiple days
            if "day" in df_locs.columns:
                if len(set(df_locs.day)) < 10:
                    df_locs = (
                        df_locs.drop(columns=["month"])
                        .drop_duplicates()
                        .reset_index(drop=True)
                    )

            if len(df_locs) > 0:
                if top_layer:
                    df_ff = ds.match_points(df_locs, quiet=True, top=top_layer)
                else:
                    df_ff = ds.match_points(
                        df_locs, depths=ds_depths, quiet=True, top=top_layer
                    )
                if df_ff is not None:
                    valid_vars = ["lon", "lat", "year", "month", "day", "depth"]
                    for vv in ds.variables:
                        valid_vars.append(vv)
                    valid_vars = [x for x in valid_vars if x in df_ff.columns]
                    df_ff = df_ff.loc[:, valid_vars]
                    df_all.append(df_ff)
        for ww in w:
            if str(ww.message) not in session_warnings:
                session_warnings.append(str(ww.message))

    except Exception as e:
        print(e)


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


def find_paths(folder, fvcom=False, exclude=[], n_check = None):
    while True:

        levels = session["levels"]

        new_directory = folder + "/"
        if levels > 0:
            for i in range(levels+1):
                dir_glob = glob.glob(new_directory + "/**")
                # randomize dir_glob

                random.shuffle(dir_glob)
                for x in dir_glob:
                    # figure out if the the base directory is an integer
                    try:
                        if levels != 0:
                            y = int(os.path.basename(x))
                        new_directory = x + "/"
                    except:
                        pass
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
    print("Parsing model information from netCDF files")

    # remove any files from options if parts of exclude are in them
    for exc in exclude:
        options = [x for x in options if f"{exc}" not in os.path.basename(x)]

    print("Searching through files in a random directory to identify variable mappings")
    # randomize options
    if n_check is not None:
        options = random.sample(options, n_check)
    for ff in tqdm(options):
        ds = nc.open_data(ff, checks=False)
        stop = True
        ds_dict = generate_mapping(ds, fvcom=fvcom)
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

    return all_df


def matchup(
    folder=None,
    start=None,
    end=None,
    surface_level=None,
    surface=[
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
    bottom=["ph", "oxygen"],
    benthic=["carbon", "benbio"],
    pft=False,
    cores=None,
    thickness=None,
    mapping=None,
    mld=False,
    exclude=[],
    levels=2,
    point_time_res=["year", "month", "day"],
    lon_lim=None,
    lat_lim=None,
    data_dir="default",
    n_check = None,
    **kwargs,
):
    """
    Match up model with observational data

    Parameters
    -------------

    folder : str
        Folder containing model output
    start : int
        Start year. First year of the simulations to matchup.
        This must be supplied
    end : int
        End year. Final year of the simulations to matchup.
        This must be supplied
    surface_level : str
        Surface level of the model netCDF files. Either 'top' or 'bottom'. Default is None, so this must be supplied.
    surface : list
        List of surface variables. Internally, ecoval will decide which variables should matchup with gridded or point observations.
        If you want finer control, you can provide a dictionary with keys 'gridded' and 'point', 
        So surface = {'gridded': ['temperature'], 'point': ['ph', 'poc']} would matchup temperature with gridded data and pH and POC with point data.
        Options for the NWS:
            gridded: ['ammonium', 'poc', 'doc', 'temperature', 'phosphate', 'salinity', 'nitrate', 'silicate', 'chlorophyll', 'oxygen']
            point: ['ammonium', 'poc', 'doc', 'alkalinity', 'temperature', 'pco2', 'phosphate', 'pft', 'nitrate', 'nitrogen', 'salinity', 'silicate', 'ph', 'chlorophyll', 'oxygen']
        Options for non-NWS models:
            gridded: ['chlorophyll', 'alkalinity', 'pco2', 'phosphate', 'co2flux', 'nitrate', 'salinity', 'silicate', 'ph', 'temperature', 'oxygen']
            point: []
    bottom : list
        List of bottom variables to matchup with observational data.
        Full list of options for NWS: ["temperature", "salinity", "oxygen", "phosphate", "silicate", "nitrate", "ammonium", "alkalinity", "ph", "chlorophyll", "doc", "poc"]
    benthic : list
        List of benthic variables.
        Full list of options for NWS: ["carbon", "benbio"]
    cores : int
        Number of cores to use.
    thickness : str
        Path to a thickness file, i.e. cellvertical  thickness. This only needs to be supplied if the variable is missing from the raw data.
        If the e3t variable is in the raw data, it will be used, and thickness does not need to be supplied.
    mapping : str
        Path to mapping file. This is a csv. A starting point can be generated by running `matchup` and saying you are not happy with the matchups.
    levels : int
        Number of levels down to look for netCDF files. Default is 2, ie. the files are of the format **/**/*.nc.
    point_time_res : list
        List of strings. Default is ['year', 'month', 'day']. This is the time resolution of the point data matchup.
    exclude : list
        List of strings to exclude. This is useful if you have files in the directory that you do not want to include in the matchup.
    mld : bool
        Matchup temperature data for mixed layer depth validation. Default is False.
    lon_lim : list
        List of two floats. The first is the minimum longitude, the second is the maximum longitude. Default is None.
    lat_lim : list
        List of two floats. The first is the minimum latitude, the second is the maximum latitude. Default is None.
    data_dir : str
        Path to validation data directory. Default is 'default'. If 'default', the data directory is taken from the session_info dictionary.
    n_check : int
        Number of files when identifying mapping. Default is None, which means all files are checked.
        The mapping is identified by looking at a random output file directory on the assumption on all directories have the same structure.
        In general each directory will only have a small number of files. Only set n_check if there are many files. 
    kwargs: dict
        Additional arguments

    """
    global_grid = None
    fvcom=False

    # make sure start and end are integers

    if start is None:
        raise ValueError("Please provide a start year")

    if end is None:
        raise ValueError("Please provide an end year")

    if isinstance(start, int) is False:
        raise ValueError("Start must be an integer")
    
    if isinstance(end, int) is False:
        raise ValueError("End must be an integer")

    # ensure time resolution is a list
    if isinstance(point_time_res, list) is False:
        raise ValueError("point_time_res must be a list")
    for x in point_time_res:
        if x not in ["year", "month", "day"]:
            raise ValueError("point_time_res must be year, month or day")
    if len(point_time_res) == 1:
        if point_time_res[0] == "day":
            raise ValueError(
                "This is not a sensible time resolution for point data matchups"
            )

    if isinstance(surface, dict):
        # make sure there are no more than 2 keys
        if len(surface.keys()) > 2:
            raise ValueError("surface dictionary can only have two keys")
        # loop through the keys

    # coerce bottom to list
    if isinstance(bottom, str):
        bottom = [bottom]

    if data_dir != "default":
        if not os.path.exists(data_dir):
            raise ValueError(f"{data_dir} does not exist")
        session_info["data_dir"] = data_dir
    else:
        data_dir = session_info["data_dir"]

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

    valid_points = list(set([x for x in glob.glob(data_dir + "/point/nws/all/*")]))
    # extract directory base name
    valid_points = [os.path.basename(x) for x in valid_points]
    for pp in point_surface:
        if pp not in valid_points:
            raise ValueError(f"{pp} is not a valid point dataset")

    valid_surface = [os.path.basename(x) for x in glob.glob(data_dir + "/gridded/*/*")]

    valid_benthic = [
        os.path.basename(x) for x in glob.glob(data_dir + "/point/nws/benthic/*")
    ]

    # fix benthic if something like "benthic biomass" is an element
    for pp in benthic:
        if "ben" in pp and "bio" in pp:
            benthic.remove(pp)
            benthic.append("benbio")

    for pp in benthic:
        if pp not in valid_benthic:
            raise ValueError(f"{pp} is not a valid benthic dataset")

    for pp in surface:
        if pp not in valid_surface:
            raise ValueError(f"{pp} is not a valid gridded dataset")

    valid_bottom = [
        os.path.basename(x) for x in glob.glob(data_dir + "/point/nws/bottom/*")
    ]
    for pp in bottom:
        if pp not in valid_bottom:
            raise ValueError(f"{pp} is not a valid bottom dataset")

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
        if key[:3] == "map":
            if mapping is None:
                mapping = kwargs[key]
        if key[:3] == "cor":
            if cores is None:
                cores = kwargs[key]

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
        "carbon" "benbio",
    ]

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
    if mld:
        point_all = ["temperature"]
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
        all_df = find_paths(folder, fvcom=fvcom, exclude=exclude, n_check = n_check)

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

            pattern = all_df.iloc[0, :].pattern

            df = all_df
            df = df.dropna()
            df = df.iloc[0:1, :]
            pattern = list(df.pattern)[0]
            pattern = pattern.replace("//", "/")

            final_extension = extension_of_directory(folder)

            if final_extension[0] == "/":
                final_extension = final_extension[1:]

            wild_card = final_extension + pattern
            wild_card = wild_card.replace("**", "*")
            for x in pathlib.Path(folder).glob(wild_card):
                path = x
                # convert to string
                path = str(path)
                break

            ds = nc.open_data(path, checks=False).to_xarray()
            lon_name = [x for x in ds.coords if "lon" in x]
            lat_name = [x for x in ds.coords if "lat" in x]
            lon = ds[lon_name[0]].values
            lat = ds[lat_name[0]].values
            lon_max = lon.max()
            lon_min = lon.min()
            lat_max = lat.max()
            lat_min = lat.min()

            global_grid = False
            if lon_max - lon_min > 350:
                global_grid = True
            if lat_max - lat_min > 170:
                global_grid = True
            if lon_max > 50:
                global_grid = True

            if global_grid:
                model_domain = "global"
            else:
                model_domain = "nws"
            print("********************************")

            # add the global checker here
            # sort all_df alphabetically by variable
            all_df = all_df.sort_values("variable").reset_index(drop=True)

            print(f"** Inferred mapping of model variable names from {folder}")
            print(all_df)
            print(
                "Are you happy with these matchups? Y/N \nNote: all possible variables are listed, not just those requested. Non-requested variables will be ignored if you answer yes."
            )
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

    # if fvcom:

    #     # matching up when fvcom
    #     print("Creating gridded data for NSBC matchups")

    #     vars = [
    #         "ammonium",
    #         "chlorophyll",
    #         "nitrate",
    #         "phosphate",
    #         "oxygen",
    #         "silicate",
    #         "temperature",
    #         "salinity",
    #     ]
    #     vars = [x for x in vars if x in var_choice]

    #     ds_total = nc.open_data()

    #     for vv in vars:
    #         pattern = all_df.query("variable == @vv").reset_index(drop=True).pattern[0]

    #         good_to_go = True

    #         if pattern is not None:
    #             final_extension = extension_of_directory(folder)
    #             ersem_paths = glob.glob(folder + final_extension + pattern)
    #             if len(ersem_paths) > 0:
    #                 good_to_go = True

    #         if good_to_go:
    #             final_extension = extension_of_directory(folder)
    #             ersem_paths = glob.glob(folder + final_extension + pattern)

    #             for exc in exclude:
    #                 ersem_paths = [
    #                     x for x in ersem_paths if f"{exc}" not in os.path.basename(x)
    #                 ]

    #             ds_all = nc.open_data()

    #             for ff in ersem_paths:
    #                 drop_variables = ["siglay", "siglev"]
    #                 ds_xr = xr.open_dataset(
    #                     ff, drop_variables=drop_variables, decode_times=False
    #                 )
    #                 model_variables = (
    #                     all_df.query("variable == @vv")
    #                     .reset_index(drop=True)
    #                     .model_variable
    #                 )
    #                 ds_xr = ds_xr[model_variables[0].split("+")]
    #                 ds1 = nc.from_xarray(ds_xr)
    #                 ds1.nco_command("ncks -d siglay,0,0")
    #                 if vv == "temp":
    #                     ds1.nco_command("ncks -O -C -v temp")
    #                 lon = ds1.to_xarray().lon.values
    #                 lat = ds1.to_xarray().lat.values
    #                 grid = pd.DataFrame({"lon": lon, "lat": lat})
    #                 out_grid = nc.generate_grid.generate_grid(grid)
    #                 ds1.subset(variables=model_variables[0].split("+"))
    #                 ds1.run()
    #                 out_grid = nc.generate_grid.generate_grid(grid)
    #                 nc.session.append_safe(out_grid)
    #                 os.path.exists(out_grid)
    #                 ds2 = ds1.copy()
    #                 ds2.run()
    #                 ds2.cdo_command(f"setgrid,{out_grid}")
    #                 ds2.as_missing(0)

    #                 if vv == "doc":
    #                     command = "-aexpr,doc=" + model_variables[0]
    #                     ds2.cdo_command(command)
    #                     drop_these = model_variables[0].split("+")
    #                     ds_contents = ds2.contents
    #                     ds_contents = ds_contents.query("variable in @drop_these")
    #                     doc_unit = ds_contents.unit[0]
    #                     ds2.set_units({"doc": doc_unit})
    #                     ds2.drop(variables=drop_these)

    #                 if vv == "chlorophyll":
    #                     command = "-aexpr,chlorophyll=" + model_variables[0]
    #                     ds2.cdo_command(command)
    #                     drop_these = model_variables[0].split("+")
    #                     ds_contents = ds2.contents
    #                     ds_contents = ds_contents.query("variable in @drop_these")
    #                     chl_unit = ds_contents.unit[0]
    #                     ds2.set_units({"chlorophyll": chl_unit})
    #                     ds2.drop(variables=drop_these)

    #                 ds_nsbc = nc.open_data(
    #                     f"{data_dir}/nsbc/level_3/climatological_monthly_mean/NSBC_Level3_phosphate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN_1960_2014.nc",
    #                     checks=False,
    #                 )
    #                 ds2.regrid(ds_nsbc, "nn")
    #                 # create a netcdf mask for the fvcom grid
    #                 df_mask = grid.assign(value=1)
    #                 bin_res = 0.25
    #                 df_mask["lon"] = bin_value(df_mask["lon"], bin_res)
    #                 df_mask["lat"] = bin_value(df_mask["lat"], bin_res)
    #                 df_mask = df_mask.groupby(["lon", "lat"]).sum().reset_index()
    #                 df_mask = df_mask.set_index(["lat", "lon"])
    #                 ds_mask = nc.from_xarray(df_mask.to_xarray())
    #                 os.system(f"cdo griddes {ds_mask[0]} > /tmp/mygrid")
    #                 # open the text file text.txt and replace the string "generic" with "lonlat"
    #                 with open("/tmp/mygrid", "r") as f:
    #                     lines = f.readlines()

    #                 # write line by line to /tmp/newgrid

    #                 with open("/tmp/newgrid", "w") as f:
    #                     for ll in lines:
    #                         f.write(ll.replace("generic", "lonlat"))
    #                 ds_mask.cdo_command(f"setgrid,/tmp/newgrid")
    #                 ds_mask.regrid(ds_nsbc, "bil")
    #                 ds_mask > 0

    #                 ds4 = ds2.copy()
    #                 # rename the variable to the correct name
    #                 ds4.rename({ds4.variables[0]: vv})
    #                 ds_all.append(ds4)
    #             ds_all.merge("time")

    #             out = "matched/gridded/nsbc/nsbc_" + vv + ".nc"
    #             if not os.path.exists(os.path.dirname(out)):
    #                 os.makedirs(os.path.dirname(out))

    #             ds_total.append(ds_all)

    #     ds_year = min(ds_total.year)

    #     ds_total.merge("variables", ["year", "month", "day"])
    #     ds_total.set_year(ds_year)

    #     out = "matched/gridded/nsbc/nsbc_model.nc"
    #     if not os.path.exists(os.path.dirname(out)):
    #         os.makedirs(os.path.dirname(out))

    #     ds_total.to_nc(out, zip=True, overwrite=True)

    #     return None

    if global_grid is None:
        final_extension = extension_of_directory(folder)
        path = glob.glob(folder + final_extension + all_df.pattern[0])[0]
        ds = nc.open_data(path, checks=False).to_xarray()
        lon_name = [x for x in ds.coords if "lon" in x]
        lat_name = [x for x in ds.coords if "lat" in x]
        lon = ds[lon_name[0]].values
        lat = ds[lat_name[0]].values
        lon_max = lon.max()
        lon_min = lon.min()
        lat_max = lat.max()
        lat_min = lat.min()

        global_grid = False
        if lon_max - lon_min > 350:
            global_grid = True
        if lat_max - lat_min > 170:
            global_grid = True
        if lon_max > 50:
            global_grid = True

        if global_grid:
            model_domain = "global"
        else:
            model_domain = "nws"

    if not surf_dict and surf_default:
        surf_all = False
        if surface == ["all"]:
            surface = copy.deepcopy(all_vars)
            surf_all = True
    

    if "ph" in surface and model_domain == "nws":
        surface.remove("ph")
        point_surface.append("ph")
    if "alkalinity" in surface and model_domain == "nws":
        surface.remove("alkalinity")
        point_surface.append("alkalinity")
    # do the same for alkalinity, poc and doc
    if "poc" in surface and model_domain == "nws":
        point_surface.append("poc")
    if "doc" in surface and model_domain == "nws":
        point_surface.append("doc")

    if pft:
        point_surface.append("pft")

    if type(surface) is str:
        surface = [surface]

    for vv in surface:
        if not os.path.exists(f"{data_dir}/gridded/{model_domain}/{vv}"):
            if not os.path.exists(f"{data_dir}/gridded/global/{vv}"):
                surface.remove(vv)

    if surf_all and surf_default:
        point_surface.append("ph")
        point_surface.append("poc")
        point_surface.append("doc")
        point_surface.append("alkalinity")

    point_surface = list(set(point_surface))

    # combine all variables into a list
    all_vars = surface + bottom + point_surface + benthic + point_all
    all_vars = list(set(all_vars))

    df_variables = all_df.query("variable in @all_vars").reset_index(drop=True)
    # remove rows where model_variable is None
    df_variables = df_variables.dropna().reset_index(drop=True)

    patterns = list(set(df_variables.pattern))


    times_dict = dict()

    print("*************************************")
    if thickness is None:
        for pattern in patterns:
            print("Identifying whether e3t exists in the files")
            final_extension = extension_of_directory(folder)
            ensemble = glob.glob(folder + final_extension + pattern)
            for exc in exclude:
                ensemble = [x for x in ensemble if f"{exc}" not in os.path.basename(x)]

            ds = nc.open_data(ensemble[0])
            if "e3t" in ds.variables:
                print(f"Extracting and saving thickness from {ensemble[0]} as matched/e3t.nc")
                ds.subset(variable = "e3t")
                ds.subset(time= 0)
                ds.as_missing(0)
                if os.path.exists("matched/e3t.nc"):
                    os.remove("matched/e3t.nc")
                ds.to_nc("matched/e3t.nc", zip=True, overwrite=True)
                thickness = "matched/e3t.nc"
                break

    print("*************************************")
    for pattern in patterns:
        print(f"Indexing file time information for {pattern} files")
        final_extension = extension_of_directory(folder)
        ensemble = glob.glob(folder + final_extension + pattern)
        for exc in exclude:
            ensemble = [x for x in ensemble if f"{exc}" not in os.path.basename(x)]

        ds = xr.open_dataset(ensemble[0])
        time_name = [x for x in list(ds.dims) if "time" in x][0]

        days = []
        for ff in tqdm(ensemble):
            ds = xr.open_dataset(ff)
            ff_month = [int(x.dt.month) for x in ds[time_name]]
            ff_year = [int(x.dt.year) for x in ds[time_name]]
            days = [int(x.dt.day) for x in ds[time_name]]
            df_ff = pd.DataFrame(
                {
                    "year": ff_year,
                    "month": ff_month,
                    "day": days,
                }
            )
            times_dict[ff] = df_ff

    # save this as a pickle
    with open("matched/times_dict.pkl", "wb") as f:
        pickle.dump(times_dict, f)

    print("********************************")
    if True:
        # figure out the lon/lat extent in the model
        lons = [lon_min, lon_max]
        lats = [lat_min, lat_max]
        # start of with the raw coords
        # This will not work with nemo, which outputs the grid incorrectly
        # so we will check if the step between the first lon/lat and the second lon/lat is
        # far bigger than the rest. If this is the case, the first should be ignored
        # get the lon/lat values
        lon_name = [x for x in ds.coords if "lon" in x]
        lat_name = [x for x in ds.coords if "lat" in x]
        lon_vals = ds[lon_name[0]].values
        lat_vals = ds[lat_name[0]].values
        # make them unique and ordered, and 1d
        lon_vals = np.unique(lon_vals)
        # make a list
        lon_vals = lon_vals.tolist()
        diff_1 = lon_vals[1] - lon_vals[0]
        diff_2 = lon_vals[2] - lon_vals[1]
        diff_3 = lon_vals[3] - lon_vals[2]
        if diff_1 / diff_2 > 10:
            if diff_1 / diff_3 > 10:
                lons[0] = lon_vals[1]
        # do it for lats
        lat_vals = np.unique(lat_vals)
        lat_vals = lat_vals.tolist()
        diff_1 = lat_vals[1] - lat_vals[0]
        diff_2 = lat_vals[2] - lat_vals[1]
        diff_3 = lat_vals[2] - lat_vals[1]
        if diff_1 / diff_2 > 10:
            if diff_1 / diff_3 > 10:
                lats[0] = lat_vals[1]

    all_df = all_df.dropna().reset_index(drop=True)
    df_mapping = all_df
    good_model_vars = [x for x in all_df.model_variable if x is not None]

    point_surface = list(set(point_surface))

    df_mapping = all_df
    if model_domain == "nws":

        if len(point_all) > 0 or len(point_bottom) > 0:
            print("Matching up with observational point data")
            print("********************************")

        # if model_variable is None remove from all_df

        for depths in ["bottom", "all", "surface", "benthic"]:
            the_vars = list(df_out.dropna().variable)
            var_choice = [x for x in var_choice if x in the_vars]
            if depths == "all":
                point_vars = point_all
            else:
                if depths == "bottom":
                    point_vars = point_bottom
                    if isinstance(point_bottom, str):
                        point_vars = [point_bottom]
                    if point_bottom is None:
                        point_bottom = []
                    # do the same for ices_all
                    if isinstance(point_all, str):
                        point_all = [point_all]
                    if point_all is None:
                        point_all = []
                if depths == "surface":
                    point_vars = point_surface
                    if isinstance(point_surface, str):
                        point_surface = [point_surface]
                    if point_surface is None:
                        point_surface = []
            if depths == "benthic":
                if isinstance(point_benthic, str):
                    point_benthic = [point_benthic]
                if point_benthic is None:
                    point_benthic = []
                point_vars = point_benthic

            for vv in point_vars:
                all_df = df_mapping
                all_df = all_df.query("model_variable in @good_model_vars").reset_index(
                    drop=True
                )

                all_df = all_df.dropna()
                if vv != "pft":
                    all_df = all_df.query("variable == @vv").reset_index(drop=True)
                else:
                    all_df = all_df.query("variable == 'chlorophyll'").reset_index(
                        drop=True
                    )
                patterns = list(set(all_df.pattern))

                for pattern in patterns:
                    final_extension = extension_of_directory(folder)
                    ensemble = glob.glob(folder + final_extension + pattern)
                    for exc in exclude:
                        ensemble = [
                            x for x in ensemble if f"{exc}" not in os.path.basename(x)
                        ]

                    df_times = []
                    days = []
                    for ff in ensemble:
                        df_ff = times_dict[ff]
                        df_times.append(
                            pd.DataFrame(
                                {
                                    "month": df_ff.month,
                                    "year": df_ff.year,
                                    "day": df_ff.day,
                                }
                            ).assign(path=ff)
                        )
                    df_times = pd.concat(df_times)

                    # figure out if it is monthly or daily data

                    df_times = df_times.query(
                        "year >= @sim_start and year <= @sim_end"
                    ).reset_index(drop=True)

                    # ersem paths

                    ersem_paths = list(set(df_times.path))
                    ersem_paths.sort()
                    # write to the report

                    write_report("### Matchup summary for observational point data")
                    min_year = df_times.year.min()
                    write_report(f"Model output start year: {min_year}")
                    max_year = df_times.year.max()
                    write_report(f"Model output end year: {max_year}")
                    write_report(
                        f"Number of years in model output: {max_year - min_year + 1}"
                    )
                    write_report(f"Number of paths: {len(ersem_paths)}")
                    # list of files
                    write_report("List of files:")

                    for ff in ersem_paths:
                        write_report(ff)

                    if depths != "surface":
                        with warnings.catch_warnings(record=True) as w:
                            # extract the thickness dataset
                            if thickness is not None:
                                ds_thickness = nc.open_data(thickness, checks=False)
                                if len(ds_thickness.variables) != 1:
                                    raise ValueError(
                                        "The thickness file has more than one variable. Please provide a single variable!"
                                    )
                                ds_thickness.rename({ds_thickness.variables[0]: "e3t"})
                            else:
                                ds_thickness = nc.open_data(ensemble[0], checks=False)

                            ds_thickness.subset(time=0, variables="e3t")
                            ds_thickness.as_missing(0)
                            #####
                            # now output the bathymetry if it does not exists
                            if not os.path.exists("matched/model_bathymetry.nc"):
                                ds_bath = ds_thickness.copy()
                                ds_bath.vertical_sum()
                                ds_bath.to_nc("matched/model_bathymetry.nc", zip=True)

                            # thickness needs to be inverted if the sea surface is at the bottom

                            if surface_level == "bottom":
                                ds_thickness.cdo_command("invertlev")
                            ds_thickness.run()
                            ds_depths = ds_thickness.copy()

                            ds_depths.vertical_cumsum()
                            ds_thickness / 2
                            ds_depths - ds_thickness
                            ds_depths.run()
                            ds_depths.rename({ds_depths.variables[0]: "depth"})
                            if surface_level == "bottom":
                                ds_depths.cdo_command("invertlev")
                            ds_depths.run()

                        for ww in w:
                            if str(ww.message) not in session_warnings:
                                session_warnings.append(str(ww.message))

                    def point_match(variable, layer="all", ds_depths=None):
                        with warnings.catch_warnings(record=True) as w:
                            point_variable = variable
                            if variable == "pft":
                                point_variable = "chlorophyll"
                            ersem_variable = list(
                                all_df.query(
                                    "variable == @point_variable"
                                ).model_variable
                            )[0]
                            paths = glob.glob(
                                f"{data_dir}/point/nws/**/{variable}/**{variable}**.feather"
                            )
                            if variable == "pft":
                                point_variable = "pft"
                            source = os.path.basename(paths[0]).split("_")[0]
                            if depths == "surface":
                                paths = [x for x in paths if "all" in x]
                            else:
                                paths = [x for x in paths if depths in x]

                            if variable == "pft":
                                paths = [x for x in paths if "pft" in x]
                            else:
                                paths = [x for x in paths if f"{point_variable}/" in x]
                            for exc in exclude:
                                paths = [
                                    x
                                    for x in paths
                                    if f"{exc}" not in os.path.basename(x)
                                ]

                            df = pd.concat([pd.read_feather(x) for x in paths])
                            # if it exists, coerce year to int
                            if "year" in df.columns:
                                df = df.assign(year=lambda x: x.year.astype(int))
                            if "month" in df.columns:
                                df = df.assign(month=lambda x: x.month.astype(int))
                            if "day" in df.columns:
                                df = df.assign(day=lambda x: x.day.astype(int))

                            if variable == "doc":
                                # go from mole to g of C
                                df = df.assign(
                                    observation=lambda x: x.observation * 12.011
                                )

                            for x in [
                                x
                                for x in ["year", "month", "day"]
                                if x not in point_time_res
                            ]:
                                if x in df.columns:
                                    df = df.drop(columns=x)
                            if depths == "surface":
                                if "depth" in df.columns:
                                    df = df.query("depth < 5").reset_index(drop=True)
                                    # grouping
                                    # drop depth
                                    df = df.drop(columns="depth")
                                    grouping = [
                                        x
                                        for x in df.columns
                                        if x
                                        in [
                                            "lon",
                                            "lat",
                                            "year",
                                            "month",
                                            "day",
                                            "source",
                                        ]
                                    ]
                                    df = df.groupby(grouping).mean().reset_index()
                                # add in a nominal depth
                            # restrict the lon_lat
                            lon_min = lons[0]
                            lon_max = lons[1]
                            lat_min = lats[0]
                            lat_max = lats[1]
                            df = df.query(
                                "lon >= @lon_min and lon <= @lon_max and lat >= @lat_min and lat <= @lat_max"
                            ).reset_index(drop=True)

                            if variable == "temperature" and mld:
                                df_include = pd.read_feather(
                                    f"{data_dir}/point/nws/mld_profiles.feather"
                                )
                                df = df.merge(df_include).reset_index(drop=True)
                            sel_these = point_time_res
                            sel_these = [x for x in df.columns if x in sel_these]
                            if variable not in ["carbon", "benbio"]:
                                paths = list(
                                    set(
                                        df.loc[:, sel_these]
                                        .drop_duplicates()
                                        .merge(df_times)
                                        .path
                                    )
                                )
                            paths = list(set(df_times.path))
                            if len(paths) == 0:
                                print(f"No matching times for {variable}")
                                return None

                            manager = Manager()
                            # time to subset the df to the lon/lat ranges

                            with warnings.catch_warnings(record=True) as w:
                                ds_grid = nc.open_data(paths[0], checks=False)
                                ds_grid.subset(variables=ds_grid.variables[0])
                                ds_grid.top()
                                ds_grid.subset(time=0)
                                amm7 = False
                                if max(ds_grid.contents.npoints) == 111375:
                                    amm7 = True
                                    ds_grid.fix_amm7_grid()
                                ds_xr = ds_grid.to_xarray()
                            for ww in w:
                                if str(ww.message) not in session_warnings:
                                    session_warnings.append(str(ww.message))
                            # extract the minimum latitude and longitude
                            lon_name = [x for x in list(ds_xr.coords) if "lon" in x][0]
                            lon_min = ds_xr[lon_name].values.min()
                            lon_max = ds_xr[lon_name].values.max()
                            lat_name = [x for x in list(ds_xr.coords) if "lat" in x][0]
                            lat_min = ds_xr[lat_name].values.min()
                            lat_max = ds_xr[lat_name].values.max()
                            df = df.query(
                                "lon >= @lon_min and lon <= @lon_max and lat >= @lat_min and lat <= @lat_max"
                            ).reset_index(drop=True)
                        for ww in w:
                            if str(ww.message) not in session_warnings:
                                session_warnings.append(str(ww.message))

                        valid_cols = [
                            "lon",
                            "lat",
                            "day",
                            "month",
                            "year",
                            "depth",
                            "observation",
                        ]
                        select_these = [x for x in df.columns if x in valid_cols]
                        if variable != "pft":
                            df = df.loc[:, select_these]

                        df_all = manager.list()

                        grid_setup = False
                        pool = multiprocessing.Pool(cores)

                        pbar = tqdm(total=len(paths), position=0, leave=True)
                        results = dict()
                        for ff in paths:
                            if grid_setup is False:
                                if True:
                                    with warnings.catch_warnings(record=True) as w:

                                        ds_grid = nc.open_data(ff, checks=False)
                                        var = ds_grid.variables[0]
                                        ds_grid.subset(variables=var)
                                        if surface_level == "top":
                                            ds_grid.top()
                                        else:
                                            ds_grid.bottom()
                                        ds_grid.as_missing(0)
                                        if max(ds_grid.contents.npoints) == 111375:
                                            ds_grid.fix_amm7_grid()
                                        df_grid = (
                                            ds_grid.to_dataframe()
                                            .reset_index()
                                            .dropna()
                                        )
                                        columns = [
                                            x
                                            for x in df_grid.columns
                                            if "lon" in x or "lat" in x
                                        ]
                                        df_grid = df_grid.loc[
                                            :, columns
                                        ].drop_duplicates()
                                        if not os.path.exists("matched"):
                                            os.makedirs("matched")
                                        df_grid.to_csv(
                                            "matched/model_grid.csv", index=False
                                        )
                                    for ww in w:
                                        if str(ww.message) not in session_warnings:
                                            session_warnings.append(str(ww.message))

                            grid_setup = True
                            if layer == "surface":
                                top_layer = True
                            else:
                                top_layer = False
                            if depths == "surface":
                                ds_depths = None
                            # raise ValueError("stoping")
                            bottom_layer = False
                            if surface_level == "bottom":
                                if layer == "surface":
                                    bottom_layer = True
                                    top_layer = False
                            if vv == "benbio":
                                bottom_layer = False
                                top_layer = False

                            # ff, ersem_variable, df, df_times, ds_depths, ices_variable, df_all, top_layer=False
                            temp = pool.apply_async(
                                # point_match,
                                mm_match,
                                [
                                    ff,
                                    ersem_variable,
                                    df,
                                    df_times,
                                    ds_depths,
                                    point_variable,
                                    df_all,
                                    top_layer,
                                    bottom_layer,
                                ],
                            )

                            results[ff] = temp

                        for k, v in results.items():
                            value = v.get()
                            pbar.update(1)

                        df_all = list(df_all)
                        df_all = [x for x in df_all if x is not None]
                        # do nothing when there is no data
                        if len(df_all) == 0:
                            print(f"No data for {variable}")
                            return None

                        df_all = pd.concat(df_all)
                        if amm7:
                            df_all = (
                                df_all.query("lon > -19")
                                .query("lon < 9")
                                .query("lat > 41")
                                .query("lat < 64.3")
                            )
                        change_this = [
                            x
                            for x in df_all.columns
                            if x
                            not in [
                                "lon",
                                "lat",
                                "year",
                                "month",
                                "day",
                                "depth",
                                "observation",
                            ]
                        ][0]
                        #
                        if vv != "pft":
                            df_all = df_all.rename(
                                columns={change_this: "model"}
                            ).merge(df)
                            # add model to name column names with frac in them
                        df_all = df_all.dropna().reset_index(drop=True)

                        grouping = copy.deepcopy(point_time_res)
                        grouping.append("lon")
                        grouping.append("lat")
                        grouping.append("depth")
                        grouping = [x for x in grouping if x in df_all.columns]
                        grouping = list(set(grouping))
                        df_all = df_all.dropna().reset_index(drop=True)
                        df_all = df_all.groupby(grouping).mean().reset_index()

                        out = f"matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}.csv"
                        # create directory for out if it does not exists
                        if not os.path.exists(os.path.dirname(out)):
                            os.makedirs(os.path.dirname(out))
                        out1 = out.replace(os.path.basename(out), "paths.csv")
                        pd.DataFrame({"path": paths}).to_csv(out1, index=False)
                        if variable == "doc":
                            df_all = df_all.assign(
                                model=lambda x: x.model + (40 * 12.011)
                            )
                        if lon_lim is not None:
                            df_all = df_all.query(
                                f"lon > {lon_lim[0]} and lon < {lon_lim[1]}"
                            )
                        if lat_lim is not None:
                            df_all = df_all.query(
                                f"lat > {lat_lim[0]} and lat < {lat_lim[1]}"
                            )

                        if vv == "pft":
                            print("Fixing pft")
                            # We now need to convert Chl to PFTs
                            ds = nc.open_data(ff, checks=False)
                            ds_contents = ds.contents
                            nano = [
                                x
                                for x in ds_contents.long_name
                                if "chloroph" in x and "nano" in x
                            ]
                            nano = ds.contents.query("long_name in @nano").variable
                            pico = [
                                x
                                for x in ds_contents.long_name
                                if "chloroph" in x and "pico" in x
                            ]
                            pico = ds.contents.query("long_name in @pico").variable
                            micro = [
                                x
                                for x in ds_contents.long_name
                                if "chloroph" in x and ("micro" in x or "diatom" in x)
                            ]
                            micro = ds.contents.query("long_name in @micro").variable
                            # convert to lists
                            nano = nano.tolist()
                            pico = pico.tolist()
                            micro = micro.tolist()
                            # do a row sum
                            df_all["nano_frac"] = df_all.loc[:, nano].sum(axis=1)
                            df_all["pico_frac"] = df_all.loc[:, pico].sum(axis=1)
                            df_all["micro_frac"] = df_all.loc[:, micro].sum(axis=1)
                            valid_vars = [
                                "lon",
                                "lat",
                                "year",
                                "month",
                                "day",
                                "nano_frac",
                                "pico_frac",
                                "micro_frac",
                            ]
                            valid_vars = [x for x in valid_vars if x in df_all.columns]
                            df_all = df_all.loc[:, valid_vars]
                            df_all.rename(
                                columns={
                                    "nano_frac": "nano_frac_model",
                                    "pico_frac": "pico_frac_model",
                                    "micro_frac": "micro_frac_model",
                                },
                                inplace=True,
                            )

                            df = df.rename(
                                columns={
                                    "nano_frac": "nano_frac_obs",
                                    "pico_frac": "pico_frac_obs",
                                    "micro_frac": "micro_frac_obs",
                                }
                            )

                            df_all = df_all.merge(df)

                        if len(df_all) > 0:
                            df_all.to_csv(out, index=False)
                            out_unit = f"matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}_unit.csv"
                            ds = nc.open_data(paths[0], checks=False)
                            ds_contents = ds.contents
                            ersem_variable = ersem_variable.split("+")[0]
                            ds_contents = ds_contents.query(
                                "variable == @ersem_variable"
                            )
                            ds_contents.to_csv(out_unit, index=False)
                        else:
                            print(f"No data for {variable}")

                    vv_variable = vv
                    if vv == "ph":
                        vv_variable = "pH"
                    if vv in ["poc", "doc"]:
                        # upper case
                        vv_variable = vv.upper()
                    if depths == "all":
                        print(
                            f"Matching up model {vv_variable} with vertically resolved bottle and CDT {vv_variable}"
                        )
                    else:
                        if depths == "surface":
                            print(
                                f"Matching up model {vv_variable} with observational surface point {vv_variable} data"
                            )
                        if depths == "bottom":
                            print(
                                f"Matching up model {vv_variable} with near-bottom point {vv_variable} data"
                            )
                    if depths == "benthic":
                        print(
                            f"Matching up model {vv_variable} with benthic point data"
                        )
                    print("**********************")
                    if depths == "surface":
                        point_match(vv, layer="surface")
                    else:
                        point_match(vv, ds_depths=ds_depths)

                    output_warnings = []
                    for ww in session_warnings:
                        if ww is not None:
                            if ww in output_warnings:
                                continue
                            output_warnings.append(str(ww))

                    if len(output_warnings) > 0:
                        output_warnings = list(set(output_warnings))
                        print(f"Warnings for {vv_variable}")
                        for ww in output_warnings:
                            warnings.warn(message=ww)
                    # empty session warnings
        while len(session_warnings) > 0:
            session_warnings.pop()

    strict = True

    gridded_matchup(
        df_mapping=df_mapping,
        folder=folder,
        var_choice=surface,
        exclude=exclude,
        surface=surface_level,
        start=start,
        sim_start=sim_start,
        sim_end=sim_end,
        e3t=thickness,
        domain=model_domain,
        strict=strict,
        lon_lim=lon_lim,
        lat_lim=lat_lim,
        times_dict=times_dict,
    )

    os.system("pandoc matchup_report.md --pdf-engine wkhtmltopdf -o matchup_report.pdf")
