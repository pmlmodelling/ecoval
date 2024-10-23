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
import xarray as xr
from ecoval.session import session_info
from multiprocessing import Manager
from tqdm import tqdm
from ecoval.utils import extension_of_directory, get_extent, fvcom_regrid
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
        Dataframe of observational data with /erie_0001.nctime information
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
            if variable == "carbon":
                var_match.append("Q7_pen_depth_c")
                var_match.append("Q6_pen_depth_c")

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
                    if variable == "carbon":

                        ds.assign(total1 = lambda x: x.Q6_c * (1 - exp(-0.1 / x.Q6_pen_depth_c)))
                        ds.assign(total2 = lambda x: x.Q7_c * (1 - exp(-0.1 / x.Q7_pen_depth_c)))
                        ds.assign(total =  lambda x: (x.total1 + x.total2)/0.1, drop = True)
                        ds * ds * 1e-6

                        ds.set_units({"total":"kg/m3"})
                    else:
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


def get_time_res(x, folder=None):
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
    # replace double stars with 1
    wild_card = wild_card.replace("**", "*")

    # figure out if the session is fvcom

    if session_info["fvcom"]:
        for x in pathlib.Path(folder).glob(wild_card):
            path = x
            # convert to string
            path = str(path)
            break
    else:
        wild_card = os.path.basename(wild_card)
        for y in pathlib.Path(folder).glob(wild_card):
            path = y
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


def extract_variable_mapping(folder, exclude=[], n_check=None, fvcom = False):
    """
    Find paths to netCDF files
    Parameters
    -------------
    folder : str
        The folder containing the netCDF files
    exclude : list
        List of strings to exclude
    n_check : int
        Number of files to check

    Returns
    -------------
    all_df : pd.DataFrame
        A DataFrame containing the paths to the netCDF files
    """

    # add restart to exclude
    exclude.append("restart")

    while True:

        levels = session_info["levels_down"]

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
        ds_dict = generate_mapping(ds, fvcom = fvcom)
        try:
            ds_dict = generate_mapping(ds, fvcom = fvcom)
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
            # replace strings of the form _12. with _**.
            new_name = re.sub(r"\d{2,}", "**", new_name)
            #new_name = re.sub(r"_\d{2,}\.", "_**.", new_name)

            all_df.append(
                pd.DataFrame.from_dict(new_dict).melt().assign(pattern=new_name)
            )

    all_df = pd.concat(all_df).reset_index(drop=True)

    patterns = set(all_df.pattern)
    resolution_dict = dict()
    for folder in patterns:
        resolution_dict[folder] = get_time_res(folder, new_directory)
    all_df["resolution"] = [resolution_dict[x] for x in all_df.pattern]

    all_df = (
        all_df.sort_values("resolution").groupby("value").head(1).reset_index(drop=True)
    )
    all_df = all_df.rename(columns={"variable": "model_variable"})
    all_df = all_df.rename(columns={"value": "variable"})
    all_df = all_df.drop(columns="resolution")
    all_df = all_df.loc[:, ["variable", "model_variable", "pattern"]]

    return all_df


def differences(
    sim_dir_1 =None,
    sim_dir_2 =None,
    #variables = {"temperature":["votemper", "votemper"]}, # this approach would be smarter....
    variables = "temperature",
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
    n_dirs_down=2,
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
    cores : int
        Number of cores to use for parallel extraction and matchups of data.
        Default is None, which means all cores are used.
        If you use a large number of cores you may run into RAM issues, so keep an eye on things.
    thickness : str
        Path to a thickness file, i.e. cell vertical thickness. This only needs to be supplied if the variable is missing from the raw data.
        If the e3t variable is in the raw data, it will be used, and thickness does not need to be supplied.
    n_dirs_down : int
        Number of levels down to look for netCDF files. Default is 2, ie. the files are of the format */*/*.nc.
    exclude : list
        List of strings to exclude. This is useful if you have files in the directory that you do not want to include in the matchup.
    lon_lim : list
        List of two floats. The first is the minimum longitude, the second is the maximum longitude. Default is None.
    lat_lim : list
        List of two floats. The first is the minimum latitude, the second is the maximum latitude. Default is None.
    n_check : int
        Number of files when identifying mapping. Default is None, which means all files are checked.
        The mapping is identified by looking at a random output file directory on the assumption on all directories have the same structure.
        In general each directory will only have a small number of files. Only set n_check if there are many files.
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
    levels_down = n_dirs_down

    if "fvcom" not in kwargs:
        fvcom = False
    if "erie" not in kwargs:
        erie = False

    # check everything is valid

    if start is None:
        raise ValueError("Please provide a start year")

    if end is None:
        raise ValueError("Please provide an end year")

    if isinstance(start, int) is False:
        raise TypeError("Start must be an integer")

    if isinstance(end, int) is False:
        raise TypeError("End must be an integer")


    if lon_lim is not None or lat_lim is not None:
        # check both are lists
        if not isinstance(lon_lim, list) or not isinstance(lat_lim, list):
            raise TypeError("lon_lim and lat_lim must be lists")

    # check if the sim_dir exists
    if sim_dir_1 is None:
        raise ValueError("Please provide a sim_dir_1 directory")

    if not os.path.exists(sim_dir_1):
        raise ValueError(f"{sim_dir} does not exist")
    
    if sim_dir_2 is None:
        raise ValueError("Please provide a sim_dir_2 directory")

    if sim_dir_2 is not None:
        if not os.path.exists(sim_dir_2):
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
        session_info["levels_down"] = levels_down
    else:
        session_info["levels_down"] = 2

    if obs_dir != "default":
        session_info["user_dir"] = True

    surface_default = False
    if surface == "default":
        surface_default = True
    if isinstance(surface, list):
        surface_default = True

    if obs_dir != "default":
        if not os.path.exists(obs_dir):
            raise ValueError(f"{obs_dir} does not exist")
        session_info["obs_dir"] = obs_dir
    else:
        obs_dir = session_info["obs_dir"]

    # loop through kwargs, if first three characters match arg and arg is None, set arg to value


    sim_start = -1000
    sim_end = 10000
    for key in kwargs:
        key_failed = True
        if key[:3] == "fol":
            if sim_dir is None:
                sim_dir = kwargs[key]
                key_failed = False
        if key == "fvcom":
            fvcom = kwargs[key]
            key_failed = False
        if key == "erie":
            erie = kwargs[key]
            key_failed = False
        if key_failed:
            raise ValueError(f"{key} is not a valid argument")
        
    session_info["fvcom"] = fvcom

    if end is not None:
        sim_end = end

    if start is not None:
        sim_start = start

    # check validity of variables chosen

    all_df = None



    var_choice = [variables]

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

    pattern_list = []
    all_df_list = []
    times_dict_list = []
    pattern_files_list = []
    for sim_dir in [sim_dir_1, sim_dir_2]:
        if all_df is None:
            all_df = extract_variable_mapping(sim_dir, exclude=exclude, n_check=n_check, fvcom = fvcom)

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

            print(all_df)


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
        if fvcom is False:
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
        else:
            if erie is False:
                global_grid = False
                model_domain = "nws"
            else:
                global_grid = True
                model_domain = "global"
        # pattern = pattern.replace("//", "/")
        all_vars = [variables]
        df_variables = all_df.query("variable in @all_vars").reset_index(drop=True)
        # remove rows where model_variable is None
        df_variables = df_variables.dropna().reset_index(drop=True)

        patterns = list(set(df_variables.pattern))

        times_dict = dict()

        pattern_files = {}
        for pattern in patterns:
            print(f"Indexing file time information for {pattern} files")
            final_extension = extension_of_directory(sim_dir)

            ensemble = glob.glob(sim_dir + final_extension + pattern)
            for exc in exclude:
                ensemble = [x for x in ensemble if f"{exc}" not in os.path.basename(x)]
            pattern_files[pattern] = ensemble

            if fvcom is False:
                ds = xr.open_dataset(ensemble[0])
                time_name = [x for x in list(ds.dims) if "time" in x][0]

            for ff in tqdm(ensemble):
                if "restart" in ff:
                    continue
                if fvcom is False:
                    ds = xr.open_dataset(ff)
                    ff_month = [int(x.dt.month) for x in ds[time_name]]
                    ff_year = [int(x.dt.year) for x in ds[time_name]]
                    days = [int(x.dt.day) for x in ds[time_name]]
                else:
                    found_times = False
                    try:
                        ds = xr.open_dataset(ff)
                        time_name = [x for x in list(ds.dims) if "time" in x][0]
                        ff_month = [int(x.dt.month) for x in ds[time_name]]
                        ff_year = [int(x.dt.year) for x in ds[time_name]]
                        days = [int(x.dt.day) for x in ds[time_name]]
                        found_times = True
                    except:
                        found_times = False

                    if not found_times:
                        ds = nc.open_data(ff, checks = False)
                        ds_times = ds.times
                        ff_month = [int(x.month) for x in ds_times]
                        ff_year = [int(x.year) for x in ds_times]
                        days = [int(x.day) for x in ds_times]
                        if len(ds_times) == 0:
                            try:
                                ds = xr.open_dataset(ff, decode_times = False)
                                times = [x for x in ds.Times.values]
                                # decode bytes
                                times = [x.decode() for x in times]
                                # times are of the format YYYY-MM-DDTHH:MM:SSZ
                                times = [x.split('T')[0] for x in times]
                                years = [x.split('-')[0] for x in times]
                                months = [x.split('-')[1] for x in times]
                                days = [x.split('-')[2] for x in times]
                                # convert to int
                                ff_year = [int(x) for x in years]
                                ff_month = [int(x) for x in months]
                                days = [int(x) for x in days]
                            except:
                                raise ValueError("No times found in the file")
                df_ff = pd.DataFrame(
                    {
                        "year": ff_year,
                        "month": ff_month,
                        "day": days,
                    }
                )
                times_dict[ff] = df_ff

        # save this as a pickle
        # make sure directory exist
        if not os.path.exists(session_info["out_dir"] + "matched"):
            os.makedirs(session_info["out_dir"] + "matched")
        with open(session_info["out_dir"] + "matched/times_dict.pkl", "wb") as f:
            pickle.dump(times_dict, f)
        
        # add results to lists
        pattern_list.append(patterns)
        all_df_list.append(all_df)
        times_dict_list.append(times_dict)
        pattern_files_list.append(pattern_files)


    for vv in var_choice:
        # now create a climatology for each variable
        min_year = start
        max_year = end
        for i in range(2):
            print(all_df_list[i])
            pattern = all_df_list[i].query("variable == @vv").reset_index(drop=True).pattern[0]
            pattern_files = pattern_files_list[0][pattern]
            times_dict = times_dict_list[0]
            i_years = []
            for ff in pattern_files:
                ff_times = times_dict[ff]
                i_years += list(set(ff_times.year))
            min_year = max(min(i_years), min_year)
            max_year = min(max(i_years), max_year)
            
        if min_year > max_year:
            raise ValueError("No data found for the specified years")
        
        # now do the climatology
        for i in range(2):
            var_files = []
            # get the files first
            pattern = all_df_list[i].query("variable == @vv").reset_index(drop=True).pattern[0]
            pattern_files = pattern_files_list[i][pattern]
            times_dict = times_dict_list[i]
            for ff in pattern_files:
                ff_times = times_dict[ff]
                ff_times = ff_times.query("year >= @min_year and year <= @max_year")
                if len(ff_times) > 0:
                    var_files.append(ff)
            
            ds = nc.open_data(var_files)
            vv_model = all_df_list[i].query("variable == @vv").reset_index(drop=True).model_variable[0]
            ds.subset(variables = vv_model)
            ds.top()
            ds.merge("time")
            ds.tmean("year")
            ds.tmean()
            sim_name = "sim_" + str(i)
            out_file = f"data/climatologies/{vv}/climatology_{vv}_{sim_name}.nc"
            # ensure directory exists
            if not os.path.exists(f"data/climatologies/{vv}"):
                os.makedirs(f"data/climatologies/{vv}")
            ds.to_nc(out_file)


        