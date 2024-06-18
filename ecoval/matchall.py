import copy

import time
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
from ecoval.utils import extension_of_directory, get_extent
from ecoval.parsers import generate_mapping
from ecoval.gridded import gridded_matchup

# a list of valid variables for validation
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
    "carbon",
    "benbio",
]

session_warnings = Manager().list()

nc.options(parallel=True)
nc.options(progress=False)


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
    df_ff = None

    if ds_depths is not None:
        nc.session.append_safe(ds_depths[0])
    try:
        with warnings.catch_warnings(record=True) as w:
            ds = nc.open_data(ff, checks=False)
            var_match = ersem_variable.split("+")

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

            t_subset = False
            if (
                "year" in df_locs.columns
                or "month" in df_locs.columns
                or "day" in df_locs.columns
            ):
                # idenify if the files have data from multiple days
                if "day" in df_locs.columns:
                    if len(set(df_locs.day)) < 10:
                        df_locs = (
                            df_locs.drop(columns=["month"])
                            .drop_duplicates()
                            .reset_index(drop=True)
                        )
                ff_indices = df_times.query("path == @ff")

                ff_indices = ff_indices.reset_index(drop=True).reset_index()
                ff_indices = ff_indices
                ff_indices = ff_indices.merge(df_locs)
                ff_indices = ff_indices["index"].values
                ff_indices = [int(x) for x in ff_indices]
                ff_indices = list(set(ff_indices))
                t_subset = True
                ds.subset(time=ff_indices)
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
            else:
                return None
        if df_ff is not None:
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
    """
    Get the time resolution of the netCDF files

    Parameters
    -------------
    x : str
        The extension of the file
    folder : str
        The folder containing the netCDF files

    Returns
    -------------
    res : str
        The time resolution of the netCDF files

    """

    final_extension = extension_of_directory(folder)

    if final_extension[0] == "/":
        final_extension = final_extension[1:]

    wild_card = final_extension + x
    wild_card = wild_card.replace("**", "*")
    for x in pathlib.Path(folder).glob(wild_card):
        path = x
        # convert to string
        path = str(path)
        break

    ds = nc.open_data(path, checks=False)
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


random_files = []


def find_paths(folder, exclude=[], n_check=None):
    while True:

        levels = session["levels_down"]

        new_directory = folder + "/"
        if levels > 0:
            for i in range(levels + 1):
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
        if True:
            options = [x for x in options if "part" not in os.path.basename(x)]
            options = [x for x in options if "restart" not in os.path.basename(x)]

        if len([x for x in options if ".nc" in x]) > 0:
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
        random_files.append(ff)
        ds = nc.open_data(ff, checks=False)
        stop = True
        ds_dict = generate_mapping(ds)
        try:
            ds_dict = generate_mapping(ds)
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

    all_df["resolution"] = [get_res(x, folder) for x in all_df.pattern]

    all_df = (
        all_df.sort_values("resolution").groupby("value").head(1).reset_index(drop=True)
    )
    all_df = all_df.rename(columns={"variable": "model_variable"})
    all_df = all_df.rename(columns={"value": "variable"})
    all_df = all_df.drop(columns="resolution")
    all_df = all_df.loc[:, ["variable", "model_variable", "pattern"]]

    return all_df


def matchup(
    sim_dir=None,
    start=None,
    end=None,
    surface_level=None,
    surface="default",
    bottom=["ph", "oxygen"],
    benthic=["carbon", "benbio"],
    pft=False,
    cores=None,
    thickness=None,
    mapping=None,
    mld=False,
    exclude=[],
    levels_down=2,
    point_time_res=["year", "month", "day"],
    lon_lim=None,
    lat_lim=None,
    obs_dir="default",
    n_check=None,
    everything=False,
    overwrite=True,
    point_all=[],
    ask=True,
    out_dir="",
    **kwargs,
):
    """
    Match up model with observational data

    Parameters
    -------------

    sim_dir : str
        Folder containing model output
    start : int
        Start year. First year of the simulations to matchup.
        This must be supplied
    end : int
        End year. Final year of the simulations to matchup.
        This must be supplied
    surface_level : str
        Surface level of the model netCDF files. Either 'top' or 'bottom'. This must be supplied.
    surface : str, list or "dict"
        This defaults to "default", i.e. it will choose all available surface variables, either from gridded or point data.
        When gridded data is available it will use it, otherwise it will use point_data.
        For finer grained control, you can pass a dictionary with two keys, "gridded" and "point".
        The value of each key is a list of variables to matchup.
        Possible variables are: ["temperature", "salinity", "oxygen", "phosphate", "silicate", "nitrate", "ammonium", "alkalinity", "ph", "chlorophyll", "doc", "pco2", "co2flux", "poc"]
        Note: availability of variables depends on the model domain.
        Point data is only available for the North West European Shelf (NWS).
    bottom : list
        List of bottom variables to matchup with observational data.
        This will only work for the north west European shelf currently.
        Full list of options for NWS: ["temperature", "salinity", "oxygen", "phosphate", "silicate", "nitrate", "ammonium", "alkalinity", "ph", "chlorophyll", "doc", "poc"]
    benthic : list
        List of benthic variables.
        This is only available for the North West European Shelf (NWS).
        Full list of options for NWS: ["carbon", "benbio"]
    pft : bool
        Matchup phytoplankton functional types. Default is False.
    cores : int
        Number of cores to use for parallel extraction and matchups of data.
        Default is None, which means all cores are used.
        If you use a large number of cores you may run into RAM issues, so keep an eye on things.
    thickness : str
        Path to a thickness file, i.e. cell vertical thickness. This only needs to be supplied if the variable is missing from the raw data.
        If the e3t variable is in the raw data, it will be used, and thickness does not need to be supplied.
    mapping : str
        Path to mapping file. This is a csv. A starting point can be generated by running `matchup` and saying you are not happy with the matchups.
    levels_down : int
        Number of levels down to look for netCDF files. Default is 2, ie. the files are of the format */*/*.nc.
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
    obs_dir : str
        Path to validation data directory. Default is 'default'. If 'default', the data directory is taken from the session_info dictionary.
    n_check : int
        Number of files when identifying mapping. Default is None, which means all files are checked.
        The mapping is identified by looking at a random output file directory on the assumption on all directories have the same structure.
        In general each directory will only have a small number of files. Only set n_check if there are many files.
    everything : bool
        If True, all possible variables at the surface and near-bottom are matched up. Default is False.
        In most cases this is overkill because point data may not tell you much gridded does not.
    point_all : list
        List of all point variables to matchup for all depths. Default is [].
    kwargs: dict
        Additional arguments
    ask : bool
        If True, the user will be asked if they are happy with the matchups. Default is True.
    out_dir : str
        Path to output directory. Default is "", so the output will be saved in the current directory.

    Returns
    -------------
    None
    Data will be stored in the matched directory.


    Examples
    -------------

    If you wanted to matchup temperature, salinity and oxygen at the surface for gridded data and the near-bottom, for the year 2002 you would run:

    >>> matchup(folder = "path/to/folder", start = 2002, end = 2002, surface = {"gridded": ["temperature", "salinity", "oxygen"], "point": None}, surface_level = "top")



    """

    # check everything is valid

    if start is None:
        raise ValueError("Please provide a start year")

    if end is None:
        raise ValueError("Please provide an end year")

    if isinstance(start, int) is False:
        raise TypeError("Start must be an integer")

    if isinstance(end, int) is False:
        raise TypeError("End must be an integer")

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

    if lon_lim is not None or lat_lim is not None:
        # check both are lists
        if not isinstance(lon_lim, list) or not isinstance(lat_lim, list):
            raise TypeError("lon_lim and lat_lim must be lists")

    # check if the sim_dir exists
    if sim_dir is None:
        raise ValueError("Please provide a sim_dir directory")

    if not os.path.exists(sim_dir):
        raise ValueError(f"{sim_dir} does not exist")

    if not isinstance(cores, int):
        raise ValueError("Please set cores to int")
    nc.options(cores=cores)

    if surface_level is None:
        raise ValueError(
            "You need to specify if the surface is the top or the bottom level"
        )

    if surface_level not in ["top", "bottom"]:
        raise ValueError("surface_level must be top or bottom")

    # set up session info, which will be needed by gridded_matchup
    session_info["overwrite"] = overwrite

    # add out_dir to session_info
    if out_dir != "":
        session_info["out_dir"] = out_dir + "/"
    else:
        session_info["out_dir"] = ""

    if levels_down is not None:
        session["levels_down"] = levels_down
    else:
        session["levels_down"] = 2

    if obs_dir != "default":
        session_info["user_dir"] = True

    surface_default = False
    if surface == "default":
        surface_default = True
    if isinstance(surface, list):
        surface_default = True

    if surface == "default":
        surface = copy.deepcopy(valid_vars)
    # add overwrite to session_info

    if everything:
        pft = True
        mld = True

    if isinstance(surface, dict):
        # make sure there are no more than 2 keys
        if len(surface.keys()) > 2:
            raise ValueError("surface dictionary can only have two keys")
        # loop through the keys

    # coerce bottom to list
    if isinstance(bottom, str):
        bottom = [bottom]

    if obs_dir != "default":
        if not os.path.exists(obs_dir):
            raise ValueError(f"{obs_dir} does not exist")
        session_info["obs_dir"] = obs_dir
    else:
        obs_dir = session_info["obs_dir"]

    # check that lon_lim and lat_lim and valid when either is not None

    if isinstance(benthic, str):
        benthic = [benthic]

    point_surface = []
    surf_dict = False

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
        if point_surface is None:
            point_surface = []
        if surface is None:
            surface = []

    if surface is None:
        surface = []

    # fix benthic if something like "benthic biomass" is an element
    if benthic is None:
        benthic = []
    if isinstance(benthic, str):
        benthic = [benthic]

    for pp in benthic:
        if "ben" in pp and "bio" in pp:
            benthic.remove(pp)
            benthic.append("benbio")

    if bottom is None:
        bottom = []
    if isinstance(bottom, str):
        bottom = [bottom]

    surface_req = copy.deepcopy(surface)
    bottom_req = copy.deepcopy(bottom)
    point_surface_req = copy.deepcopy(point_surface)
    benthic_req = copy.deepcopy(benthic)

    if isinstance(exclude, str):
        exclude = [exclude]

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

    sim_start = -1000
    sim_end = 10000
    for key in kwargs:
        key_failed = True
        if key[:3] == "fol":
            if sim_dir is None:
                sim_dir = kwargs[key]
                key_failed = False
        if key_failed:
            raise ValueError(f"{key} is not a valid argument")

    if end is not None:
        sim_end = end

    if start is not None:
        sim_start = start

    # check validity of variables chosen

    all_df = None
    if isinstance(mapping, pd.DataFrame):
        all_df = mapping
    if isinstance(mapping, str):
        all_df = pd.read_csv(mapping)

    # create lists for working out which variables are needed for point matchups
    # change point_all to list if str
    if isinstance(point_all, str):
        point_all = [point_all]
    # check point_all is a list
    if not isinstance(point_all, list):
        raise ValueError("point_all must be a list")
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

    # restrict surface to valids

    # surface = [x for x in surface if x in valid_surface]

    if all_df is None:
        all_df = find_paths(sim_dir, exclude=exclude, n_check=n_check)

        # add in anything that is missing
        all_vars = valid_vars

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

    # pattern = all_df.iloc[0, :].pattern

    # df = all_df
    # df = df.dropna()
    # df = df.iloc[0:1, :]
    # pattern = list(df.pattern)[0]
    # pattern = pattern.replace("//", "/")
    pattern = all_df.reset_index(drop=True).iloc[0, :].pattern

    final_extension = extension_of_directory(sim_dir)

    if final_extension[0] == "/":
        final_extension = final_extension[1:]

    wild_card = final_extension + pattern
    wild_card = wild_card.replace("**", "*")
    for x in pathlib.Path(sim_dir).glob(wild_card):
        path = x
        # convert to string
        path = str(path)
        break

    ds = nc.open_data(path, checks=False)
    ds_extent = get_extent(ds[0])
    lon_max = ds_extent[1]
    lon_min = ds_extent[0]
    lat_max = ds_extent[3]
    lat_min = ds_extent[2]

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

    if session_info["user_dir"]:
        if global_grid:
            valid_points = list(
                set([x for x in glob.glob(obs_dir + "/point/user/all/*")])
            )
        else:
            valid_points = list(
                set([x for x in glob.glob(obs_dir + "/point/**/all/*")])
            )
    # extract directory base name
    valid_points = [os.path.basename(x) for x in valid_points]
    for pp in point_surface:
        if pp not in valid_points:
            raise ValueError(f"{pp} is not a valid point dataset")

    if global_grid:
        if session_info["user_dir"]:
            valid_surface = [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/user/*")
            ]
            valid_surface += [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/global/*")
            ]
        else:
            valid_surface = [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/global/*")
            ]
    else:
        if session_info["user_dir"]:
            valid_surface = [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/nws/*")
            ]
            valid_surface += [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/user/*")
            ]
        else:
            valid_surface = [
                os.path.basename(x) for x in glob.glob(obs_dir + "/gridded/nws/*")
            ]

    valid_benthic = [
        os.path.basename(x) for x in glob.glob(obs_dir + "/point/nws/benthic/*")
    ]

    if session_info["user_dir"]:
        valid_bottom = [
            os.path.basename(x) for x in glob.glob(obs_dir + "/point/user/bottom/*")
        ]
        if len(valid_bottom) == 0:
            valid_bottom = [
                os.path.basename(x) for x in glob.glob(obs_dir + "/point/nws/bottom/*")
            ]

    else:
        valid_bottom = [
            os.path.basename(x) for x in glob.glob(obs_dir + "/point/nws/bottom/*")
        ]

    if everything:
        surface = valid_surface
        # only valid variables
        surface = [x for x in surface if x in valid_vars]
        point_surface = valid_points
        point_surface = [x for x in point_surface if x in valid_vars]
        benthic = valid_benthic
        benthic = [x for x in benthic if x in valid_vars]
        bottom = valid_bottom
        bottom = [x for x in bottom if x in valid_vars]

    if global_grid:
        point_surface = []
        point_benthic = []
        point_bottom = []
        surface = [x for x in surface if x in valid_surface]
    else:
        point_surface = [x for x in point_surface if x in valid_points]
        point_benthic = [x for x in point_benthic if x in valid_benthic]
        point_bottom = [x for x in point_bottom if x in valid_bottom]
        surface = [x for x in surface if x in valid_surface]

    vars_available = list(
        all_df
        # drop rows where pattern is None
        .dropna()
        # get all variables
        .variable
    )
    # check variables chosen are valid

    for vv in surface_req:
        if vv not in valid_surface:
            raise ValueError(f"{vv} is not a valid surface dataset")
    for vv in bottom_req:
        if vv not in valid_bottom:
            raise ValueError(f"{vv} is not a valid bottom dataset")
    for vv in point_surface_req:
        if vv not in valid_points:
            raise ValueError(f"{vv} is not a valid point dataset")
    for vv in benthic_req:
        if vv not in valid_benthic:
            raise ValueError(f"{vv} is not a valid benthic dataset")

    surface = [x for x in surface if x in vars_available]
    point_surface = [x for x in point_surface if x in vars_available]
    point_bottom = [x for x in point_bottom if x in vars_available]
    point_benthic = [x for x in point_benthic if x in vars_available]
    var_chosen = surface + bottom + point_benthic + point_bottom + point_surface
    var_chosen = list(set(var_chosen))

    if len(point_bottom) > 0 or mld or len(point_all) > 0:
        ds_depths = False
        if True:
            if True:
                if True:
                    try:
                        if True:
                            with warnings.catch_warnings(record=True) as w:
                                # extract the thickness dataset
                                e3t_found = False
                                if thickness is not None:
                                    ds_thickness = nc.open_data(thickness, checks=False)
                                    if len(ds_thickness.variables) != 1:
                                        raise ValueError(
                                            "The thickness file has more than one variable. Please provide a single variable!"
                                        )
                                    ds_thickness.rename(
                                        {ds_thickness.variables[0]: "e3t"}
                                    )
                                    e3t_found = True
                                else:
                                    print(
                                        "Vertical thickness is required for your matchups, but they are not supplied"
                                    )
                                    print(
                                        "Searching through simulation output to find it"
                                    )
                                    for ff in random_files:
                                        ds_thickness = nc.open_data(ff, checks=False)
                                        if "e3t" in ds_thickness.variables:
                                            break
                                            e3t_found
                                if not e3t_found:
                                    raise ValueError("Unable to find e3t")

                                ds_thickness.subset(time=0, variables="e3t")
                                ds_thickness.as_missing(0)
                                #####
                                # now output the bathymetry if it does not exists
                                if not os.path.exists("matched/model_bathymetry.nc"):
                                    ds_bath = ds_thickness.copy()
                                    ds_bath.vertical_sum()
                                    ds_bath.to_nc(
                                        "matched/model_bathymetry.nc", zip=True
                                    )

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
                    except:
                        pass
        if ds_depths is False:
            raise ValueError(
                "You have asked for variables that require the specification of thickness"
            )
    if mapping is not None:
        if True:

            # add the global checker here
            # sort all_df alphabetically by variable
            all_df = all_df.sort_values("variable").reset_index(drop=True)
            surface.sort()
            point_surface.sort()
            point_bottom.sort()
            point_benthic.sort()
            # ensure variables are in all_df.variable
            vars_available = list(
                all_df
                # drop rows where pattern is None
                .dropna()
                # get all variables
                .variable
            )
            surface = [x for x in surface if x in vars_available]
            point_surface = [x for x in point_surface if x in vars_available]
            point_bottom = [x for x in point_bottom if x in vars_available]
            point_benthic = [x for x in point_benthic if x in vars_available]

            all_df_print = copy.deepcopy(all_df).reset_index(drop=True)

            # new tidied variable
            new_variable = []
            for i in range(len(all_df_print)):
                if all_df.variable[i] in var_chosen:
                    if all_df.pattern[i] is not None:
                        new_variable.append(all_df.variable[i] + "**")
                    else:
                        new_variable.append(all_df.variable[i])
                else:
                    new_variable.append(all_df.variable[i])
            all_df_print["variable"] = new_variable
            print(all_df_print)
            print(
                "Note: all possible variables are listed, not just those requested. Variables that will be matched up are starred."
            )
            print("Variables that will be matched up")
            print("******************************")

            if len(surface) > 0:
                print(
                    f"The following variables will be matched up with gridded surface data: {','.join(surface)}"
                )
                missing_surface = [x for x in valid_surface if x not in surface]
                if len(missing_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_surface)}"
                    )
            else:
                print("No variables will be matched up with gridded surface data")
                missing_surface = [x for x in valid_surface if x not in surface]
                if len(missing_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_surface)}"
                    )

            if len(point_surface) > 0:
                print(
                    f"The following variables will be matched up with in-situ near-bottom data: {','.join(point_surface)}"
                )
                missing_point_surface = [
                    x for x in valid_points if x not in point_surface
                ]
                if len(missing_point_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_point_surface)}"
                    )
            else:
                print("No variables will be matched up with in-situ surface data")
                missing_point_surface = [
                    x for x in valid_points if x not in point_surface
                ]
                if len(missing_point_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_point_surface)}"
                    )

            if len(point_bottom) > 0:
                print(
                    f"The following variables will be matched up with in-situ near-bottom data: {','.join(point_bottom)}"
                )
                missing_bottom = [x for x in valid_bottom if x not in point_bottom]
                if len(missing_bottom) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_bottom)}"
                    )
            else:
                print("No variables will be matched up with in-situ bottom data")
                missing_bottom = [x for x in valid_bottom if x not in point_bottom]
                if len(missing_bottom) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_bottom)}"
                    )

            if len(point_benthic) > 0:
                print(
                    f"The following variables will be matched up with in-situ benthic data: {','.join(point_benthic)}"
                )
                missing_benthic = [x for x in valid_benthic if x not in point_benthic]
                if len(missing_benthic) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_benthic)}"
                    )

            else:
                print("No variables will be matched up with in-situ benthic data")
                missing_benthic = [x for x in valid_benthic if x not in point_benthic]
                if len(missing_benthic) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_benthic)}"
                    )

            print("Are you happy with this? Y/N")

            if ask:
                x = input()
            else:
                x = "y"
            if x == "n":
                return None

    if mapping is None:
        if True:
            # add the global checker here
            # sort all_df alphabetically by variable
            all_df = all_df.sort_values("variable").reset_index(drop=True)
            surface.sort()
            point_surface.sort()
            point_bottom.sort()
            point_benthic.sort()
            print("Variables that will be matched up")
            print("******************************")
            if len(surface) > 0:
                print(
                    f"The following variables will be matched up with gridded surface data: {','.join(surface)}"
                )
                missing_surface = [x for x in valid_surface if x not in surface]
                if len(missing_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_surface)}"
                    )
            else:
                print("No variables will be matched up with gridded surface data")
                missing_surface = [x for x in valid_surface if x not in surface]
                if len(missing_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_surface)}"
                    )

            if len(point_surface) > 0:
                print(
                    f"The following variables will be matched up with in-situ surface data: {','.join(point_surface)}"
                )
                missing_point_surface = [
                    x for x in valid_points if x not in point_surface
                ]
                if len(missing_point_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_point_surface)}"
                    )
            else:
                print("No variables will be matched up with in-situ surface data")
                missing_point_surface = [
                    x for x in valid_points if x not in point_surface
                ]
                if len(missing_point_surface) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_point_surface)}"
                    )

            if len(point_bottom) > 0:
                print(
                    f"The following variables will be matched up with in-situ near-bottom data: {','.join(point_bottom)}"
                )
                missing_bottom = [x for x in valid_bottom if x not in point_bottom]
                if len(missing_bottom) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_bottom)}"
                    )
            else:
                print("No variables will be matched up with in-situ bottom data")
                missing_bottom = [x for x in valid_bottom if x not in point_bottom]
                if len(missing_bottom) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_bottom)}"
                    )

            if len(point_benthic) > 0:
                print(
                    f"The following variables will be matched up with in-situ benthic data: {','.join(point_benthic)}"
                )
                missing_benthic = [x for x in valid_benthic if x not in point_benthic]
                if len(missing_benthic) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_benthic)}"
                    )

            else:
                print("No variables will be matched up with in-situ benthic data")
                missing_benthic = [x for x in valid_benthic if x not in point_benthic]
                if len(missing_benthic) > 0:
                    print(
                        f"Surface variables that could be validated, but are not requested: {', '.join(missing_benthic)}"
                    )

            print("******************************")
            print(f"** Inferred mapping of model variable names from {sim_dir}")

            all_df_print = copy.deepcopy(all_df).reset_index(drop=True)

            # new tidied variable
            new_variable = []
            for i in range(len(all_df_print)):
                if all_df.variable[i] in var_chosen:
                    if all_df.pattern[i] is not None:
                        new_variable.append(all_df.variable[i] + "**")
                    else:
                        new_variable.append(all_df.variable[i])
                else:
                    new_variable.append(all_df.variable[i])
            all_df_print["variable"] = new_variable
            print(all_df_print)
            print(
                "Note: all possible variables are listed, not just those requested. Variables that will be matched up are starred."
            )

            print("Are you happy with these matchups? Y/N")

            if ask:
                x = input()
            else:
                x = "y"

            if x.lower() not in ["y", "n"]:
                print("Provide Y or N")
                x = input()

            if x.lower() == "n":
                out = get_out()
                print(f"Inferred mapping saved as {out}")
                all_df.to_csv(out, index=False)
                return None

    if session_info["out_dir"] != "":
        out = session_info["out_dir"] + "/matched/mapping.csv"
    else:
        out = "matched/mapping.csv"
    # check directory exists for out
    out_folder = os.path.dirname(out)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    df_out = all_df.dropna().reset_index(drop=True)
    final_extension = extension_of_directory(sim_dir)
    df_out["pattern"] = [sim_dir + final_extension + x for x in df_out.pattern]
    df_out.to_csv(out, index=False)

    if global_grid is None:
        final_extension = extension_of_directory(sim_dir)
        path = glob.glob(sim_dir + final_extension + all_df.pattern[0])[0]
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

    if "ph" in surface and model_domain == "nws":
        surface.remove("ph")
        point_surface.append("ph")
    if "alkalinity" in surface and model_domain == "nws":
        surface.remove("alkalinity")
        point_surface.append("alkalinity")
    # do the same for alkalinity, poc and doc
    if "poc" in surface and model_domain == "nws":
        if surface_default:
            point_surface.append("poc")
    if "doc" in surface and model_domain == "nws":
        if surface_default:
            point_surface.append("doc")

    if pft:
        point_surface.append("pft")

    if type(surface) is str:
        surface = [surface]

    point_surface = list(set(point_surface))

    # combine all variables into a list
    all_vars = surface + bottom + point_surface + benthic + point_all
    all_vars = list(set(all_vars))

    if pft:
        all_vars.append("micro")
        all_vars.append("nano")
        all_vars.append("pico")

    df_variables = all_df.query("variable in @all_vars").reset_index(drop=True)
    # remove rows where model_variable is None
    df_variables = df_variables.dropna().reset_index(drop=True)

    patterns = list(set(df_variables.pattern))

    times_dict = dict()

    print("*************************************")
    thick_found = False
    if thickness is None:
        print("Identifying whether e3t exists in the files")
        for pattern in patterns:
            final_extension = extension_of_directory(sim_dir)
            ensemble = glob.glob(sim_dir + final_extension + pattern)
            for exc in exclude:
                ensemble = [x for x in ensemble if f"{exc}" not in os.path.basename(x)]

            ds = nc.open_data(ensemble[0])
            if "e3t" in ds.variables:
                print(
                    f"Extracting and saving thickness from {ensemble[0]} as matched/e3t.nc"
                )
                ds.subset(variable="e3t")
                ds.subset(time=0)
                ds.as_missing(0)
                if os.path.exists("matched/e3t.nc"):
                    os.remove("matched/e3t.nc")
                ds.to_nc("matched/e3t.nc", zip=True, overwrite=True)
                thickness = "matched/e3t.nc"
                thick_found = True
                break
    if not thick_found:
        if thickness is None:
            print("It was not. Assuming files have z-levels for any vertical matchups.")

    print("*************************************")
    for pattern in patterns:
        print(f"Indexing file time information for {pattern} files")
        final_extension = extension_of_directory(sim_dir)
        ensemble = glob.glob(sim_dir + final_extension + pattern)
        for exc in exclude:
            ensemble = [x for x in ensemble if f"{exc}" not in os.path.basename(x)]

        ds = xr.open_dataset(ensemble[0])
        time_name = [x for x in list(ds.dims) if "time" in x][0]

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
    with open(session_info["out_dir"] + "matched/times_dict.pkl", "wb") as f:
        pickle.dump(times_dict, f)

    print("********************************")

    # figure out the lon/lat extent in the model
    ds_extent = get_extent(ensemble[0])
    lons = [ds_extent[0], ds_extent[1]]
    lats = [ds_extent[2], ds_extent[3]]

    all_df = all_df.dropna().reset_index(drop=True)
    df_mapping = all_df
    good_model_vars = [x for x in all_df.model_variable if x is not None]

    point_surface = list(set(point_surface))

    df_mapping = all_df

    # raise ValueError("here")
    if model_domain == "nws" or session_info["user_dir"]:

        if len(point_all) > 0 or len(point_bottom) > 0:
            print("Matching up with observational point data")
            print("********************************")

        # if model_variable is None remove from all_df

        for depths in ["bottom", "all", "surface", "benthic"]:
            the_vars = list(df_mapping.dropna().variable)
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

            # sort the list
            point_vars.sort()

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
                    final_extension = extension_of_directory(sim_dir)
                    ensemble = glob.glob(sim_dir + final_extension + pattern)
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

                    sim_paths = list(set(df_times.path))
                    sim_paths.sort()
                    # write to the report

                    write_report("### Matchup summary for observational point data")
                    min_year = df_times.year.min()
                    write_report(f"Model output start year: {min_year}")
                    max_year = df_times.year.max()
                    write_report(f"Model output end year: {max_year}")
                    write_report(
                        f"Number of years in model output: {max_year - min_year + 1}"
                    )
                    write_report(f"Number of paths: {len(sim_paths)}")
                    # list of files
                    write_report("List of files:")

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
                            if session_info["user_dir"]:
                                paths = glob.glob(
                                    f"{obs_dir}/point/user/**/{variable}/**{variable}**.feather"
                                )
                                if len(paths) == 0:
                                    paths = glob.glob(
                                        f"{obs_dir}/point/nws/**/{variable}/**{variable}**.feather"
                                    )
                            else:
                                paths = glob.glob(
                                    f"{obs_dir}/point/nws/**/{variable}/**{variable}**.feather"
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
                                    f"{obs_dir}/point/nws/mld_profiles.feather"
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
                            else:
                                paths = list(set(df_times.path))

                            if len(paths) == 0:
                                print(f"No matching times for {variable}")
                                raise ValueError("here")

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
                        if len(df) == 0:
                            print("No data for this variable")
                            return None

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
                                        if not os.path.exists(
                                            session_info["out_dir"] + "matched"
                                        ):
                                            os.makedirs(
                                                session_info["out_dir"] + "matched"
                                            )
                                        df_grid.to_csv(
                                            session_info["out_dir"]
                                            + "matched/model_grid.csv",
                                            index=False,
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
                            bottom_layer = False
                            if surface_level == "bottom":
                                if layer == "surface":
                                    bottom_layer = True
                                    top_layer = False
                            if vv == "benbio":
                                bottom_layer = False
                                top_layer = False

                            temp = pool.apply_async(
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
                            time.sleep(1)
                            return False

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

                        if session_info["out_dir"] != "":
                            out = f"{session_info['out_dir']}/matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}.csv"
                        else:
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
                            # do a row sum
                            nano = (
                                df_mapping.query("variable == 'nano'")
                                .model_variable.values[0]
                                .split("+")
                            )
                            pico = (
                                df_mapping.query("variable == 'pico'")
                                .model_variable.values[0]
                                .split("+")
                            )
                            micro = (
                                df_mapping.query("variable == 'micro'")
                                .model_variable.values[0]
                                .split("+")
                            )
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

                        print("************************")
                        print(out)
                        print("************************")
                        if len(df_all) > 0:
                            df_all.to_csv(out, index=False)
                            if session_info["out_dir"] != "":
                                out_unit = f"{session_info['out_dir']}/matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}_unit.csv"
                            else:
                                out_unit = f"matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}_unit.csv"
                            ds = nc.open_data(paths[0], checks=False)
                            ds_contents = ds.contents
                            ersem_variable = ersem_variable.split("+")[0]
                            ds_contents = ds_contents.query(
                                "variable == @ersem_variable"
                            )
                            ds_contents.to_csv(out_unit, index=False)
                            return None
                        else:
                            print(f"No data for {variable}")
                            time.sleep(1)
                            return False

                    vv_variable = vv
                    if vv == "ph":
                        vv_variable = "pH"
                    if vv in ["poc", "doc"]:
                        # upper case
                        vv_variable = vv.upper()
                    if session_info["out_dir"] != "":
                        out = glob.glob(
                            session_info["out_dir"]
                            + "/"
                            + f"matched/point/{model_domain}/{depths}/{vv}/**_{depths}_{vv}.csv"
                        )

                    else:
                        out = glob.glob(
                            f"matched/point/{model_domain}/{depths}/{vv}/**_{depths}_{vv}.csv"
                        )

                    if len(out) > 0:
                        if session_info["overwrite"] is False:
                            continue

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
                        try:
                            point_match(vv, layer="surface")
                        except:
                            pass
                    else:
                        point_match(vv, ds_depths=ds_depths)
                        try:
                            print(layer)
                            point_match(vv, ds_depths=ds_depths)
                        except:
                            pass

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

    gridded_matchup(
        df_mapping=df_mapping,
        folder=sim_dir,
        var_choice=surface,
        exclude=exclude,
        surface=surface_level,
        sim_start=sim_start,
        sim_end=sim_end,
        domain=model_domain,
        lon_lim=lon_lim,
        lat_lim=lat_lim,
        times_dict=times_dict,
        ds_thickness=thickness,
    )

    os.system("pandoc matchup_report.md --pdf-engine wkhtmltopdf -o matchup_report.pdf")

    if len(session_info["end_messages"]) > 0:
        print("########################################")
        print("########################################")
        print("Important messages about matchups:")
        print("*" * 30)
        for x in session_info["end_messages"]:
            print(x)
        print("########################################")
        print("########################################")

    # currently redundant fvcom code. Add in later...
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
