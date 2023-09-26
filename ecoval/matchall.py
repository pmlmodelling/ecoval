import nctoolkit as nc

nc.options(parallel=True)
nc.options(progress=False)
import glob
import multiprocessing
from multiprocessing import Process, Manager
import os
import pandas as pd
import warnings
import time

import numpy as np
import xarray as xr


def bin_value(x, bin_res):
    return np.floor((x + bin_res / 2) / bin_res + 0.5) * bin_res - bin_res / 2


import xarray as xr
from ecoval.ices import generate_mapping
from ecoval.nsbc import nsbc_matchup
from ecoval.fixers import tidy_warnings

import random
import string

ices_variables = [
    "temperature",
    "salinity",
    "oxygen",
    "chlorophyll",
    "phosphate",
    "silicate",
    "nitrate",
    "ph",
    "ammonium",
    "co2flux",
    "pco2",
    "alkalinity",
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
    # valid_keys = ["thread_safe", "lazy", "cores", "precision", "temp_dir"]

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
            # print(key)
        data_path = terms[1].replace(" ", "")

        if os.path.exists(data_path):
            data_dir = data_path
        else:
            raise ValueError(f"{data_path} does not exist")


def mm_match(ff, ersem_variable, df, df_times, ds_depths, ices_variable, df_all):
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

    # return "x"
    nc.session.append_safe(ds_depths[0])
    try:
        with warnings.catch_warnings(record=True) as w:
            ds = nc.open_data(ff, checks=False)
            ds.subset(variables=ersem_variable)
            ds.as_missing(0)
            ds.run()

            df_locs = (
                df_times.query("path == @ff")
                .merge(df)
                .loc[:, ["lon", "lat", "year", "month", "day", "depth"]]
                .drop_duplicates()
                .reset_index(drop=True)
            )
            # idenify if the files have data from multiple days
            if len(set(df_locs.day)) < 10:
                df_locs = (
                    df_locs.drop(columns=["month"])
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

            if len(df_locs) > 0:
                df_ff = ds.match_points(df_locs, depths=ds_depths, quiet=True)
                df_all.append(df_ff)
        tidy_warnings(w)

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

    path = glob.glob(folder + "/**/**/" + x)[0]

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


def find_paths(folder, fvcom=False):
    i = 1
    while True:
        ensemble = glob.glob(folder + "/**")
        ensemble = [x for x in ensemble if os.path.isdir(x)]
        ensemble = [x for x in ensemble if os.path.basename(x).isdigit()]
        ensemble = list(set(ensemble))
        len_ensemble = len(ensemble)
        stop = -1

        if len_ensemble > 10:
            stop = len_ensemble - int(len_ensemble / 2)

        while True:
            bad = True
            for x in ensemble[-i:]:
                if os.path.isdir(x):
                    try:
                        y = int(x.split("/")[-1])
                        bad = False
                        ensemble = glob.glob(x + "/**")
                        final = x
                        break

                    except:
                        blah = "blah"
            if bad:
                if i > stop:
                    break
            i += 1
        options = glob.glob(final + "/**.nc")
        if not fvcom:
            options = [x for x in options if "part" not in os.path.basename(x)]
            options = [x for x in options if "restart" not in os.path.basename(x)]
        #     options = [x for x in options if "1d" in x or "1m" in x and "part" not in os.path.basename(x)]
        # else:
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
    for ff in options:
        ds = nc.open_data(ff, checks=False)
        ds_dict = generate_mapping(ds, fvcom=fvcom)

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
            all_df.append(
                pd.DataFrame.from_dict(new_dict).melt().assign(pattern=new_name)
            )
    all_df = pd.concat(all_df).reset_index(drop=True)

    if fvcom == False:
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


def matchup(
    folder=None,
    spinup=None,
    mapping=None,
    cores=None,
    variables=None,
    fvcom=False,
    e3t=None,
    exclude=[],
    surface=None,
    start=None,
    end=None,
    ices_all=["temperature", "ph", "alkalinity"],
    ices_bottom=[
        "temperature",
        "salinity",
        "oxygen",
        "phosphate",
        "silicate",
        "nitrate",
        "ammonium",
    ],
    **kwargs,
):
    """
    Match up model with observational data

    Parameters
    -------------
    folder: str
        Directory path for model output
    spinup: int
        Number of years to remove from the start of the model output
    mapping: str or pd.DataFrame
        Path to mapping file or dataframe with columns variable, model_variable, pattern
    cores: int
        Number of cores to use for processing. This must be set
    variables: str or list
        Variables to match. If None, all variables will be matched
    fvcom: bool
        If True, the matchups are assumed to be FVCOM
    e3t: str
        Path to e3t file. If None, e3t will be assumed to be in the files
    start: None
        Start year for matchups
    end: None
        End year for matchups
    ices: list
        List of variables to match with ICES point data
        # take them from ices_variables
        These must be from the list:
        ['temperature',
            'salinity',
            'oxygen',
            'chlorophyll',
            'phosphate',
            'silicate',
            'nitrate',
            'ph',
            'ammonium',
            "co2flux",
            "pco2",
            'alkalinity']
    **kwargs:
        Additional arguments to pass to matchup


    """

    # check if the folder exists
    if folder is None:
        raise ValueError("Please provide a folder")
    
    if not os.path.exists(folder):
        raise ValueError(f"{folder} does not exist")

    # loop through kwargs, if first three characters match arg and arg is None, set arg to value

    # check if the ices variables are in ices_variables

    if isinstance(ices_bottom, str):
        ices_bottom = [ices_bottom]
    if isinstance(ices_all, str):
        ices_all = [ices_all]
    if ices_bottom is None:
        ices_bottom = []
    if ices_all is None:
        ices_all = []

    for vv in ices_bottom:
        if vv not in ices_variables:
            raise ValueError(
                f"{vv} is not a valid variable. Please choose from {ices_variables}"
            )

    for vv in ices_all:
        if vv not in ices_variables:
            raise ValueError(
                f"{vv} is not a valid variable. Please choose from {ices_variables}"
            )

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

    if surface is None:
        raise ValueError(
            "You need to specify if the surface is the top or the bottom level"
        )

    if surface not in ["top", "bottom"]:
        raise ValueError("surface must be top or bottom")

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
        if key[:3] == "var":
            if variables is None:
                variables = kwargs[key]

    if start is not None:
        if spinup is not None:
            raise ValueError("You can only provide one of start or spinup")

    if end is not None:
        sim_end = end

    if start is not None:
        sim_start = start

    if variables is None:
        var_choice = [
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
        ]
    else:
        var_choice = variables
        if isinstance(var_choice, str):
            var_choice = [var_choice]
        if isinstance(var_choice, list):
            var_choice = var_choice
        if not isinstance(var_choice, list):
            raise ValueError("Please provide variables as a list or a string")

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
    ]

    for vv in var_choice:
        if vv not in valid_vars:
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

    if all_df is None:
        all_df = find_paths(folder, fvcom=fvcom)

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
    df_out["pattern"] = [folder + "/**/**/" + x for x in df_out.pattern]
    df_out.to_csv(out, index=False)

    if fvcom:
        import xarray as xr

        # matching up when fvcom
        print("Creating gridded data for NSBC matchups")

        vars = [
            "ammonium",
            "chlorophyll",
            "nitrate",
            "phosphate",
            "oxygen",
            "silicate",
            "temperature",
            "salinity",
        ]
        vars = [x for x in vars if x in var_choice]

        ds_total = nc.open_data()

        for vv in vars:
            pattern = all_df.query("variable == @vv").reset_index(drop=True).pattern[0]

            good_to_go = True

            if pattern is not None:
                ersem_paths = glob.glob(folder + "/**/**/" + pattern)
                if len(ersem_paths) > 0:
                    good_to_go = True

            if good_to_go:
                ersem_paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    ersem_paths = [
                        x for x in ersem_paths if f"{exc}" not in os.path.basename(x)
                    ]

                ds_all = nc.open_data()

                for ff in ersem_paths:
                    drop_variables = ["siglay", "siglev"]
                    ds_xr = xr.open_dataset(
                        ff, drop_variables=drop_variables, decode_times=False
                    )
                    model_variables = (
                        all_df.query("variable == @vv")
                        .reset_index(drop=True)
                        .model_variable
                    )
                    ds_xr = ds_xr[model_variables[0].split("+")]
                    ds1 = nc.from_xarray(ds_xr)
                    ds1.nco_command("ncks -d siglay,0,0")
                    if vv == "temp":
                        ds1.nco_command("ncks -O -C -v temp")
                    lon = ds1.to_xarray().lon.values
                    lat = ds1.to_xarray().lat.values
                    grid = pd.DataFrame({"lon": lon, "lat": lat})
                    out_grid = nc.generate_grid.generate_grid(grid)
                    ds1.subset(variables=model_variables[0].split("+"))
                    ds1.run()
                    out_grid = nc.generate_grid.generate_grid(grid)
                    nc.session.append_safe(out_grid)
                    os.path.exists(out_grid)
                    ds2 = ds1.copy()
                    ds2.run()
                    ds2.cdo_command(f"setgrid,{out_grid}")
                    ds2.as_missing(0)

                    if vv == "chlorophyll":
                        command = "-aexpr,chlorophyll=" + model_variables[0]
                        ds2.cdo_command(command)
                        drop_these = model_variables[0].split("+")
                        ds_contents = ds2.contents
                        ds_contents = ds_contents.query("variable in @drop_these")
                        chl_unit = ds_contents.unit[0]
                        ds2.set_units({"chlorophyll": chl_unit})
                        ds2.drop(variables=drop_these)

                    ds_nsbc = nc.open_data(
                        f"{data_dir}/nsbc/level_3/climatological_monthly_mean/NSBC_Level3_phosphate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN_1960_2014.nc",
                        checks=False,
                    )
                    ds2.regrid(ds_nsbc, "nn")
                    # create a netcdf mask for the fvcom grid
                    df_mask = grid.assign(value=1)
                    bin_res = 0.25
                    df_mask["lon"] = bin_value(df_mask["lon"], bin_res)
                    df_mask["lat"] = bin_value(df_mask["lat"], bin_res)
                    df_mask = df_mask.groupby(["lon", "lat"]).sum().reset_index()
                    df_mask = df_mask.set_index(["lat", "lon"])
                    ds_mask = nc.from_xarray(df_mask.to_xarray())
                    os.system(f"cdo griddes {ds_mask[0]} > /tmp/mygrid")
                    # open the text file text.txt and replace the string "generic" with "lonlat"
                    with open("/tmp/mygrid", "r") as f:
                        lines = f.readlines()

                    # write line by line to /tmp/newgrid

                    with open("/tmp/newgrid", "w") as f:
                        for ll in lines:
                            f.write(ll.replace("generic", "lonlat"))
                    ds_mask.cdo_command(f"setgrid,/tmp/newgrid")
                    ds_mask.regrid(ds_nsbc, "bil")
                    ds_mask > 0

                    ds4 = ds2.copy()
                    # rename the variable to the correct name
                    ds4.rename({ds4.variables[0]: vv})
                    ds_all.append(ds4)
                ds_all.merge("time")

                out = "matched/gridded/nsbc/nsbc_" + vv + ".nc"
                if not os.path.exists(os.path.dirname(out)):
                    os.makedirs(os.path.dirname(out))

                ds_total.append(ds_all)

        ds_year = min(ds_total.year)

        ds_total.merge("variables", ["year", "month", "day"])
        ds_total.set_year(ds_year)

        out = "matched/gridded/nsbc/nsbc_model.nc"
        if not os.path.exists(os.path.dirname(out)):
            os.makedirs(os.path.dirname(out))

        ds_total.to_nc(out, zip=True, overwrite=True)

        return None

    pattern = all_df.iloc[0, :].pattern

    variable = "ices_variable"
    # get the units. File inspection could be randomized in case people have put loose files in there...
    import glob

    print("********************************")
    print("Identifying whether grid is global")
    df = pd.read_csv("matched/mapping.csv")
    df = df.dropna()
    df = df.iloc[0:1, :]
    pattern = list(df.pattern)[0]
    pattern = pattern.replace("//", "/")
    while True:
        i = 0
        patterns = pattern.split("/")
        for x in patterns:
            if x == "**":
                break
            i += 1
        new_pattern = glob.glob("/".join(patterns[0:i]) + "/" + "**")[-1].split("/")[-1]
        patterns[i] = new_pattern
        pattern = "/".join(patterns)

        if len([x for x in pattern.split("/") if x == "**"]) == 0:
            break
    paths = glob.glob(pattern)
    ds = nc.open_data(paths[0], checks=False).to_xarray()
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
    if global_grid:
        print("Grid is global")
    else:
        print("Grid is not global")
    print("********************************")

    all_df = all_df.dropna().reset_index(drop = True)
    df_mapping = all_df
    # get rid of any rows where pattern is None

    if not global_grid:

        if len(ices_all) > 0 or len(ices_bottom) > 0:
            print("Matching up with ICES point data")
            print("********************************")

        df_mapping = all_df

        # if model_variable is None remove from all_df
        good_model_vars = [x for x in all_df.model_variable if x is not None]

        if True:
            for depths in ["bottom", "all"]:
                the_vars = list(df_out.dropna().variable)
                var_choice = [x for x in var_choice if x in the_vars]
                if depths == "all":
                    ices_vars = ices_all
                else:
                    ices_vars = ices_bottom
                    if isinstance(ices_bottom, str):
                        ices_vars = [ices_bottom]
                    if ices_bottom is None:
                        ices_bottom = []
                    # do the same for ices_all
                    if isinstance(ices_all, str):
                        ices_all = [ices_all]
                    if ices_all is None:
                        ices_all = []

                for vv in ices_vars:
                    all_df = df_mapping
                    all_df = all_df.query(
                        "model_variable in @good_model_vars"
                    ).reset_index(drop=True)

                    all_df = all_df.dropna()
                    all_df = all_df.query("variable == @vv").reset_index(drop=True)
                    patterns = list(set(all_df.pattern))

                    for pattern in patterns:
                        ensemble = glob.glob(folder + "/**/**/" + pattern)
                        for exc in exclude:
                            ensemble = [
                                x
                                for x in ensemble
                                if f"{exc}" not in os.path.basename(x)
                            ]
                        import xarray as xr

                        ds = xr.open_dataset(ensemble[0])
                        time_name = [x for x in list(ds.dims) if "time" in x][0]
                        from tqdm import tqdm
                        import xarray as xr

                        nc_times = []
                        df_times = []
                        for ff in ensemble:
                            ds = xr.open_dataset(ff)
                            ff_month = [int(x.dt.month) for x in ds[time_name]][0]
                            ff_year = [int(x.dt.year) for x in ds[time_name]][0]
                            df_times.append(
                                pd.DataFrame(
                                    {
                                        "path": [ff],
                                        "month": [ff_month],
                                        "year": [ff_year],
                                    }
                                )
                            )
                        df_times = pd.concat(df_times)

                        df_times = df_times.query(
                            "year >= @sim_start and year <= @sim_end"
                        ).reset_index(drop=True)

                        if spinup is not None:
                            min_year = df_times.year.min() + spinup
                        else:
                            min_year = df_times.year.min()

                        df_times = df_times.query("year >= @min_year").reset_index(
                            drop=True
                        )

                        # ersem paths

                        ersem_paths = list(set(df_times.path))
                        ersem_paths.sort()
                        # write to the report

                        write_report(f"### Matchup summary for ICES point data")
                        min_year = df_times.year.min()
                        write_report(f"Model output start year: {min_year}")
                        max_year = df_times.year.max()
                        write_report(f"Model output end year: {max_year}")
                        write_report(
                            f"Number of years in model output: {max_year - min_year + 1}"
                        )
                        write_report(f"Number of paths: {len(ersem_paths)}")
                        # list of files
                        write_report(f"List of files:")

                        for ff in ersem_paths:
                            write_report(ff)

                        with warnings.catch_warnings(record=True) as w:
                            # extract the thickness dataset
                            if e3t is not None:
                                ds_thickness = nc.open_data(e3t, checks=False)
                                if "e3t" not in ds_thickness.variables:
                                    options = [
                                        x for x in ds_thickness.variables if "e3t" in x
                                    ]
                                    if len(options) != 1:
                                        raise ValueError("e3t not found in e3t file")
                                    ds_thickness.rename({options[0]: "e3t"})
                            else:
                                ds_thickness = nc.open_data(ensemble[0], checks=False)

                            ds_thickness.subset(time=0, variables="e3t")
                            ds_thickness.as_missing(0)
                            #####

                            # thickness needs to be inverted if the sea surface is at the bottom

                            if surface == "bottom":
                                ds_thickness.cdo_command("invertlev")
                            ds_thickness.run()
                            ds_depths = ds_thickness.copy()

                            ds_depths.vertical_cumsum()
                            ds_thickness / 2
                            ds_depths - ds_thickness
                            ds_depths.run()
                            ds_depths.rename({ds_depths.variables[0]: "depth"})
                            if surface == "bottom":
                                ds_depths.cdo_command("invertlev")
                            ds_depths.run()

                        tidy_warnings(w)

                        def ices_match(variable):
                            with warnings.catch_warnings(record=True) as w:
                                ds = xr.open_dataset(ensemble[0])
                                time_name = [x for x in list(ds.dims) if "time" in x][0]
                                ices_variable = variable
                                ersem_variable = list(
                                    all_df.query(
                                        "variable == @ices_variable"
                                    ).model_variable
                                )[0]
                                paths = glob.glob(
                                    f"{data_dir}/ices/**/**/ices_**.csv"
                                )
                                paths = [x for x in paths if depths in x]
                                paths = [x for x in paths if f"{ices_variable}/" in x]
                                for exc in exclude:
                                    paths = [
                                        x
                                        for x in paths
                                        if f"{exc}" not in os.path.basename(x)
                                    ]

                                df = pd.concat([pd.read_csv(x) for x in paths])
                                if variable == "temperature":
                                    df_include = pd.read_csv(
                                        f"{data_dir}/ices/mld_profiles.csv"
                                    )
                                    df = df.merge(df_include).reset_index(drop=True)
                                paths = list(
                                    set(
                                        df.loc[:, ["year", "month", "day"]]
                                        .drop_duplicates()
                                        .merge(df_times)
                                        .path
                                    )
                                )
                                if len(paths) == 0:
                                    print(f"No matching times for {variable}")
                                    return None

                                manager = Manager()
                                # time to subset the df to the lon/lat ranges
                                ds_grid = nc.open_data(paths[0], checks=False)
                                ds_grid.subset(variables=ds_grid.variables[0])
                                ds_grid.top()
                                ds_grid.subset(time=0)
                                if max(ds_grid.contents.npoints) == 111375:
                                    ds_grid.fix_amm7_grid()
                                ds_xr = ds_grid.to_xarray()
                                # extract the minimum latitude and longitude
                                lon_name = [
                                    x for x in list(ds_xr.coords) if "lon" in x
                                ][0]
                                lon_min = ds_xr[lon_name].values.min()
                                lon_max = ds_xr[lon_name].values.max()
                                lat_name = [
                                    x for x in list(ds_xr.coords) if "lat" in x
                                ][0]
                                lat_min = ds_xr[lat_name].values.min()
                                lat_max = ds_xr[lat_name].values.max()
                                df = df.query(
                                    "lon >= @lon_min and lon <= @lon_max and lat >= @lat_min and lat <= @lat_max"
                                ).reset_index(drop=True)
                            tidy_warnings(w)

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
                                            if surface == "top":
                                                ds_grid.top()
                                            else:
                                                ds_grid.bottom()
                                            ds_grid.as_missing(0)
                                            ds = nc.open_data(paths[0], checks=False)
                                            amm7 = False
                                            if max(ds_grid.contents.npoints) == 111375:
                                                ds_grid.fix_amm7_grid()
                                                amm7 = True
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
                                        tidy_warnings(w)

                                grid_setup = True
                                temp = pool.apply_async(
                                    mm_match,
                                    [
                                        ff,
                                        ersem_variable,
                                        df,
                                        df_times,
                                        ds_depths,
                                        ices_variable,
                                        df_all,
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
                            df_all = df_all.rename(
                                columns={df_all.columns[-1]: "model"}
                            ).merge(df)
                            df_all = df_all.dropna().reset_index(drop=True)
                            out = f"matched/ices/{depths}/ices_{depths}_{variable}.csv"
                            # create directory for out if it does not exists
                            if not os.path.exists(os.path.dirname(out)):
                                os.makedirs(os.path.dirname(out))
                            df_all.to_csv(out, index=False)

                        print("**********************")
                        vv_variable = vv
                        if vv == "ph":
                            vv_variable = "pH"
                        if depths == "all":
                            print(
                                f"Matching up model {vv_variable} with vertically resolved ICES bottle and CDT data"
                            )
                        else:
                            print(
                                f"Matching up model {vv_variable} with ICES near-bottom data"
                            )
                        print("**********************")
                        ices_match(vv)

        # def nsbc_matchup(df_mapping = None, folder = None, var_choice = None, exclude = None, surface = None, start = None, spinup = None, sim_start = None, sim_end = None, e3t = None, report = None):
        nsbc_matchup(
            df_mapping=df_mapping,
            folder=folder,
            var_choice=var_choice,
            exclude=exclude,
            surface=surface,
            start=start,
            spinup=spinup,
            sim_start=sim_start,
            sim_end=sim_end,
            e3t=e3t,
        )

        # extracting sst

        if "temperature" in var_choice:
            print("********************************")
            print("Creating gridded data for OSTIA matchups")
            print("********************************")

            vars = ["temperature"]

            out_dir = "matched/gridded/ostia/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                # get the ostia paths
                paths = nc.create_ensemble(
                    f"{data_dir}/ostia"
                )
                years = [
                    int(os.path.basename(x).split("_")[1].replace(".nc", ""))
                    for x in paths
                ]
                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []
                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    # ff_month = ff_times.dt.month.values
                    # ff_year = ff_times.dt.year.values
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                # subset using spinup

                if spinup is not None:
                    min_year = df_times.year.min() + spinup
                else:
                    min_year = df_times.year.min()

                df_times = df_times.query("year >= @min_year")

                df_times = df_times.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                ostia_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    # ff_month = ff_times.dt.month.values
                    # ff_year = ff_times.dt.year.values
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )

                df_times = pd.concat(df_times)

                if spinup is not None:
                    min_year = df_times.year.min() + spinup
                else:
                    min_year = df_times.year.min()

                ersem_paths = (
                    df_times.query("year >= @min_year")
                    .loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                    .merge(ostia_paths.loc[:, ["year"]], on="year", how="inner")
                )

                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)
                ersem_paths = ersem_paths.query("year >= @min_year").reset_index(
                    drop=True
                )

                # write the report
                write_report("### Matchups for temperature")

                # number of paths
                n_paths = len(ersem_paths)
                write_report(f"Number of paths: {n_paths}")
                # minimum year
                write_report(f"Minimum year: {min(ersem_paths.year)}")
                # maximum year
                write_report(f"Maximum year: {max(ersem_paths.year)}")
                # write the list of files
                write_report(f"Files used for temperature:")
                write_report("```")
                for ff in ersem_paths.path:
                    write_report(ff)
                write_report("```")

                # add a pagebreak

                write_report("\\newpage")

                ostia_paths = ostia_paths.merge(
                    ersem_paths.loc[:, ["year"]], on="year", how="inner"
                ).drop_duplicates()
                ersem_paths = ersem_paths.loc[:, ["path"]]
                ostia_paths = ostia_paths.loc[:, ["path"]]

                ersem_paths = list(ersem_paths.path)
                ostia_paths = list(ostia_paths.path)

                out_file = out_dir + "ostia_model.nc"

                # extract the ersem data

                selection = df_mapping.query(
                    "variable == 'temperature'"
                ).model_variable.values[0]
                with warnings.catch_warnings(record=True) as w:

                    ds_model = nc.open_data(ersem_paths, checks=False)
                    ds_model.subset(variables=selection)
                    if surface == "top":
                        ds_model.top()
                    else:
                        ds_model.bottom()
                    ds_model.as_missing(0)
                    ds_model.merge("time")
                    ds_model.tmean(["year", "month"])
                    # get monthly mean for each year in case files represent each day
                    ds_model.tmean(["month", "year"])
                    ds_model.rename({selection: "model"})
                    # check the number of grid cells before fixing nemo ersem grid

                    ds = nc.open_data(paths[0], checks=False)
                    amm7 = False
                    if max(ds_model.contents.npoints) == 111375:
                        ds_model.fix_amm7_grid()
                        amm7 = True
                    ds_model.run()

                tidy_warnings(w)

                # extract the ostia data

                ds_obs = nc.open_data(ostia_paths, checks=False)
                ds_obs.subset(variables="analysed_sst")
                ds_obs.top()
                ds_obs.as_missing(0)
                ds_obs.tmean("month")
                ds_obs.merge("time")
                ds_obs.rename({"analysed_sst": "observation"})
                if amm7:
                    ds_obs.regrid(ds_model)
                else:
                    ds_model.regrid(ds_obs)
                ds_obs.set_precision("F32")
                # change temperature in ds_obs from Kelvin to Celsius
                ds_obs - 273.15

                ds_obs.run()

                years = [
                    x
                    for x in ds_obs.years
                    if x >= min_year and x >= sim_start and x <= sim_end
                ]
                ds_obs.subset(years=years)

                years = [
                    x
                    for x in ds_model.years
                    if x >= min_year and x >= sim_start and x <= sim_end
                ]
                ds_model.subset(years=years)

                ds_model.append(ds_obs)

                ds_model.merge("variable", ["year", "month"])

                ds_model.to_nc(out_file, zip=True, overwrite=True)

    else:
        if "co2flux" in var_choice:
            # extract the CO2 fluxes
            vars = ["co2flux"]

            out_dir = "matched/gridded/co2fluxes/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                print("********************************")
                print("Matching CO2 fluxes")
                print("********************************")
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").reset_index(drop = True).model_variable)[0]

                selection = []
                for vv in vars:
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

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    # ff_times = xr.open_dataset(ff)[time_name]
                    ds = xr.open_dataset(ff)
                    ff_month = [int(x.dt.month) for x in ds[time_name]]
                    ff_year = [int(x.dt.year) for x in ds[time_name]]
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                ersem_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                if spinup is not None:
                    min_year = ersem_paths.year.min() + spinup
                
                if start is not None:
                    if spinup is None:
                        min_year = start
                if end is not None:
                    max_year = end
                    ersem_paths = ersem_paths.query("year <= @max_year")

                ersem_paths = ersem_paths.query("year >= @min_year")
                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                ersem_years = list(set(ersem_paths.year))

                ersem_paths = ersem_paths.loc[:, ["path"]]

                ersem_paths = list(ersem_paths.path)

                out_file = out_dir + "co2fluxes.nc"

                # extract the model data

                selection = df_mapping.query(
                    "variable == 'co2flux'"
                ).reset_index(drop = True).model_variable.values[0]

                with warnings.catch_warnings(record=True) as w:
                    ds_model = nc.open_data(ersem_paths, checks=False)
                    ds_model.subset(variables=selection)
                    if surface == "top":
                        ds_model.top()
                    else:
                        ds_model.bottom()
                    ds_model.as_missing(0)
                    ds_model.tmean("month")
                    ds_model.merge("time")
                    ds_model.rename({selection: "model"})
                    ds_model.run()
                tidy_warnings(w)

                # extract the observation data

                ds_obs = nc.open_data(
                    f"{data_dir}/fluxes/MPI_SOM_FFN_2022_NCEI_OCADS.nc",
                    checks=False,
                )

                ds_obs.subset(variables="fgco2_smoothed")
                ds_obs.top()
                ds_obs.subset(years=ds_model.years)
                ds_obs.rename({"fgco2_smoothed": "observation"})
                # change temperature in ds_obs from Kelvin to Celsius

                ds_obs.run()

                ds_model.regrid(ds_obs)
                ds_model.run()

                years = [x for x in ds_obs.years if x in ds_model.years]
                ds_model.subset(years=years)
                ds_obs.subset(years=years)

                ds_model.append(ds_obs)
                ds_model.merge("variable", "month")
                ds_model.assign(model=lambda x: x.model * -0.365)
                ds_model.set_units({"model": "mol/m2/yr"})
                ds_model.set_units({"observation": "mol/m2/yr"})

                ds_model.to_nc(out_file, zip=True, overwrite=True)

        if "pco2" in var_choice:
            # extract the surface CO2
            vars = ["pco2"]

            out_dir = "matched/gridded/pco2/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                print("******************************")
                print("Matching the surface CO2")
                print("******************************")
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = ff_times.dt.month.values
                    ff_year = ff_times.dt.year.values
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                ersem_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                min_year = ersem_paths.year.min() + spinup

                ersem_paths = ersem_paths.query("year >= @min_year")
                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                ersem_years = list(set(ersem_paths.year))

                ersem_paths = ersem_paths.loc[:, ["path"]]

                ersem_paths = list(ersem_paths.path)

                out_file = out_dir + "pco2.nc"

                # extract the model data

                selection = df_mapping.query(
                    "variable == 'pco2'"
                ).model_variable.values[0]

                with warnings.catch_warnings(record=True) as w:
                    ds_model = nc.open_data(ersem_paths, checks=False)
                    ds_model.subset(variables=selection)
                    if surface == "top":
                        ds_model.top()
                    else:
                        ds_model.bottom()
                    ds_model.as_missing(0)
                    ds_model.tmean("month")
                    ds_model.merge("time")
                    ds_model.rename({selection: "model"})
                    ds_model.run()
                tidy_warnings(w)

                # extract the observation data

                ds_obs = nc.open_data(
                    f"{data_dir}/fluxes/MPI-SOM_FFN_SOCCOMv2018_minus_4uatm_offset.nc",
                    checks=False,
                )

                ds_obs.subset(variables="spco2")
                ds_obs.top()
                ds_obs.subset(years=ds_model.years)
                ds_obs.rename({"spco2": "observation"})
                # change temperature in ds_obs from Kelvin to Celsius

                ds_obs.run()

                ds_model.regrid(ds_obs)
                ds_model.run()

                years = [x for x in ds_obs.years if x in ds_model.years]
                ds_model.subset(years=years)
                ds_obs.subset(years=years)

                ds_model.append(ds_obs)
                ds_model.merge("variable", "month")
                ds_model.to_nc(out_file, zip=True, overwrite=True)

        if True:
            bad_vars = [x for x in df_mapping.model_variable if x is None]
            all_df = df_mapping.query("model_variable not in @bad_vars")

            all_df = all_df.query("variable in @var_choice").reset_index(drop=True)

            variables = list(all_df.variable)

            noaa_paths = glob.glob(f"{data_dir}/woa/*.nc")

            noaa_variables = [
                x.split("/")[-1].split("_")[1].replace(".nc", "") for x in noaa_paths
            ]

            variables = [x for x in variables if x in ["temperature"]]

            if len(variables) > 0:
                print("********************************")
                print("Matching up WOA temperature with model data")
                print("********************************")
                pattern = all_df.query("variable == @variables[0]").pattern.values[0]

                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                ds_obs = nc.open_data(
                    f"{data_dir}/woa/woa_{variables[0]}.nc",
                    checks=False,
                )
                ds_obs.subset(variables="*_an")
                ds_obs.run()
                ds_obs.rename({ds_obs.variables[0]: "observation"})
                model_var = all_df.query(
                    "variable == @variables[0]"
                ).model_variable.values[0]

                years = []
                if len(paths) == 0:
                    raise ValueError("here")

                # time name
                time_name = [
                    x
                    for x in nc.open_data(paths[0], checks=False).to_xarray().dims
                    if "time" in x
                ][0]

                # we need to identify bad paths
                bad_paths = []

                for ff in paths:
                    try:
                        ds = nc.open_data(ff, checks=False).to_xarray()
                        ff_year = [int(x.dt.year) for x in ds[time_name]][0]
                        years.append(ff_year)
                    except:
                        bad_paths.append(ff)

                paths = [x for x in paths if x not in bad_paths]

                df_years = pd.DataFrame({"year": years, "path": paths})

                if spinup is not None:
                    min_year = df_years.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start
                if end is not None:
                    max_year = end
                    df_years = df_years.query("year <= @max_year").reset_index(drop = True)

                df_years = df_years.query("year >= @min_year")
                df_years = df_years.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)
                # limit it to the final 20 years
                min_year = df_years.year.max() - 20
                df_years = df_years.query("year >= @min_year")
                # limit to NOAA WOA years
                df_years = (
                    df_years.query("year >= 2005")
                    .query("year <= 2017")
                    .reset_index(drop=True)
                )
                if len(df_years) == 0:
                    print("No model temperature between 2005 and 2017")
                else:
                    paths = list(df_years.query("year >= @min_year").path.values)

                    with warnings.catch_warnings(record=True) as w:
                        ds_model = nc.open_data(paths, checks=False)
                        # different approach for mission atlantic files
                        if "MissionAtlantic" in ff:
                            nc.options(cores = 1)
                            ds_model.nco_command(f"ncks   -v {model_var}")
                            nc.options(cores = cores)
                        else:
                            ds_model.subset(variables=model_var)
                        

                        ds_model.as_missing(0)
                        ds_model.tmean("month")
                        ds_model.merge("time")
                        ds_model.tmean("month")
                        ds_model.rename({model_var: "model"})
                        ds_model.run()
                    tidy_warnings(w)

                    noaa_levels = ds_obs.levels
                    # replace 0.0 with 2.5 in noaa_levels
                    noaa_levels = [2.5 if x < 2.0 else x for x in noaa_levels]

                    ds_thickness = nc.open_data(paths[0], checks=False)
                    if "e3t" in ds_thickness.variables:
                        ds_thickness.subset(time=0, variables="e3t")
                        ds_thickness.as_missing(0)
                        ds_thickness.run()
                        ds_model.vertical_interp(noaa_levels, ds_thickness=ds_thickness)
                    else:
                        ds_model.vertical_interp(noaa_levels, fixed=True)

                    ds_obs.vertical_interp(noaa_levels, fixed=True)

                    ds_model.regrid(ds_obs)
                    ds_obs.append(ds_model)
                    ds_obs.merge("variable", "month")

                    out_dir = "matched/woa/"
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)

                    out_file = out_dir + f"woa_{variables[0]}.nc"

                    ds_obs.to_nc(out_file, zip=True, overwrite=True)

        if True:

            # extracting sst

            # remove rows where model_variable is None in all_df
            bad_vars = [x for x in df_mapping.model_variable if x is None]
            all_df = df_mapping.query("model_variable not in @bad_vars")

            all_df = all_df.query("variable in @var_choice").reset_index(drop=True)

            variables = list(all_df.variable)

            noaa_paths = glob.glob(f"{data_dir}/woa/*.nc")

            noaa_variables = [
                x.split("/")[-1].split("_")[1].replace(".nc", "") for x in noaa_paths
            ]
            # exclude temperature from this
            noaa_variables = [x for x in noaa_variables if x != "temperature"]

            variables = [x for x in variables if x in noaa_variables]

            for vv in variables:
                print("********************************")
                print(f"Matching up WOA {vv} with model data")
                print("********************************")
                pattern = all_df.query("variable == @vv").pattern.values[0]

                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                ds_obs = nc.open_data(
                    f"{data_dir}/woa/woa_{vv}.nc",
                    checks=False,
                )
                ds_obs.subset(variables="*_an")
                ds_obs.run()
                ds_obs.rename({ds_obs.variables[0]: "observation"})
                model_var = all_df.query("variable == @vv").model_variable.values[0]

                years = []
                if len(paths) == 0:
                    raise ValueError("here")

                # time name
                time_name = [
                    x
                    for x in nc.open_data(paths[0], checks=False).to_xarray().dims
                    if "time" in x
                ][0]

                # we need to identify bad paths
                bad_paths = []

                for ff in paths:
                    try:
                        ds = nc.open_data(ff, checks=False).to_xarray()
                        ff_year = [int(x.dt.year) for x in ds[time_name]][0]
                        years.append(ff_year)
                    except:
                        bad_paths.append(ff)

                paths = [x for x in paths if x not in bad_paths]

                df_years = pd.DataFrame({"year": years, "path": paths})

                if spinup is not None:
                    min_year = df_years.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start
                if end is not None:
                    max_year = end
                    df_years = df_years.query("year <= @max_year").reset_index(drop = True)

                df_years = df_years.query("year >= @min_year")
                df_years = df_years.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)
                # limit it to the final 20 years
                min_year = df_years.year.max() - 20
                df_years = df_years.query("year >= @min_year")
                paths = list(df_years.query("year >= @min_year").path.values)

                with warnings.catch_warnings(record=True) as w:
                    ds_model = nc.open_data(paths, checks=False)
                    # different approach for mission atlantic files
                    if "MissionAtlantic" in ff:
                        ds_model.nco_command(f"ncks -F -d deptht,1 -v {model_var}")
                    else:
                        ds_model.subset(variables=model_var)
                        if surface == "top":
                            ds_model.top()
                        else:
                            ds_model.bottom()
                    ds_model.as_missing(0)
                    ds_model.tmean("month")
                    ds_model.merge("time")
                    ds_model.rename({model_var: "model"})
                    ds_model.tmean("month")
                    ds_model.run()
                tidy_warnings(w)

                ds_model.regrid(ds_obs)
                ds_obs.append(ds_model)
                ds_obs.merge("variable", "month")

                out_dir = "matched/woa/"
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)

                out_file = out_dir + f"woa_{vv}.nc"

                ds_obs.to_nc(out_file, zip=True, overwrite=True)

        if "chlorophyll" in var_choice:
            # match occci chlorophyll
            vars = ["chlorophyll"]

            out_dir = "matched/gridded/occci/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)

            ensemble = nc.create_ensemble(
                "/data/datasets/CCI/v5.0-release/geographic/netcdf/monthly/chlor_a"
            )
            years = [int(os.path.basename(x).split("-")[-2][0:4]) for x in ensemble]
            months = [int(os.path.basename(x).split("-")[-2][4:6]) for x in ensemble]
            paths = ensemble

            for exc in exclude:
                paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

            df_occci = pd.DataFrame({"year": years, "month": months, "path": paths})

            df = df_mapping.query("variable in @vars").reset_index(drop=True)

            if len(df) > 0:
                print("********************************")
                print("Matching OCCI chlorophyll with model data")
                print("********************************")
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                paths = glob.glob(folder + "/**/**/" + pattern)

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                occci_dates = df_occci.loc[:, ["year", "month"]].drop_duplicates()
                all_dates = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .merge(occci_dates)
                    .drop_duplicates()
                )

                ersem_paths = df_times

                if spinup is not None:
                    min_year = ersem_paths.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start
                    
                if end is not None:
                    max_year = end
                    ersem_paths = ersem_paths.query("year <= @max_year")
                ersem_paths = ersem_paths.query("year >= @min_year")

                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                ersem_paths = (
                    ersem_paths.loc[:, ["year", "month"]]
                    .merge(occci_dates)
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(ersem_paths, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                all_dates = all_dates.merge(
                    ersem_paths.loc[:, ["year"]], on="year", how="inner"
                )

                df_occci = df_occci.merge(
                    all_dates, on=["year", "month"], how="inner"
                ).drop_duplicates()

                ersem_years = list(set(ersem_paths.year))

                ersem_paths = ersem_paths.loc[:, ["path"]]
                paths = list(ersem_paths.path)
                # print(paths)
                # print(df_occci)
                # raise ValueError("What is happening?")

                with warnings.catch_warnings(record=True) as w:
                    ds = nc.open_data(paths, checks=False)
                    ds.subset(variables=selection)
                    if surface == "top":
                        ds.top()
                    else:
                        ds.bottom()
                    ds.as_missing(0)
                    ds.tmean("month")

                    if "chlorophyll" in list(df.variable):
                        command = "-aexpr,chlorophyll=" + mapping["chlorophyll"]
                        ds.cdo_command(command)
                        drop_these = mapping["chlorophyll"].split("+")
                        ds_contents = ds.contents
                        ds_contents = ds_contents.query("variable in @drop_these")
                        chl_unit = ds_contents.unit[0]
                        ds.drop(variables=drop_these)
                    ds.run()

                tidy_warnings(w)
                for key in mapping:
                    if key != "chlorophyll":
                        if mapping[key] in ds.variables:
                            ds.rename({mapping[key]: key})
                if "chlorophyll" in list(df.variable):
                    ds.set_units({"chlorophyll": chl_unit})
                    ds.set_longnames({"chlorophyll": "Total chlorophyll concentration"})

                ds.rename({"chlorophyll": "model"})

                ds.run()
                ds.merge("time")

                ds_all = ds.copy()

                ds_obs = nc.open_data(df_occci.path.values, checks=False)
                ds_obs.subset(variables="chlor_a")
                ds_obs.regrid(ds_all)
                ds_obs.rename({"chlor_a": "observation"})
                ds_obs.run()
                ds_obs.merge("time")
                ds_obs.run()
                ds_all.run()
                if len(ds_all.times) != len(ds_obs.times):
                    raise ValueError("Something is wrong with the times")
                ds_all.append(ds_obs)
                ds_all.merge("variable", ["year", "month"])
                ds_all.to_nc(out_dir + "occci_model.nc", zip=True, overwrite=True)

        if "temperature" in var_choice:
            vars = ["temperature"]

            out_dir = "matched/gridded/cobe2/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                print("******************************")
                print("Matching the surface temperature from COBE2 with model data")
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                ersem_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                if spinup is not None:
                    min_year = ersem_paths.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start

                ersem_paths = ersem_paths.query("year >= @min_year")
                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                if end is not None:
                    max_year = end
                    ersem_paths = ersem_paths.query("year <= @max_year").reset_index(drop = True)

                ersem_years = list(set(ersem_paths.year))

                ersem_paths = ersem_paths.loc[:, ["path"]]

                ersem_paths = list(ersem_paths.path)

                out_file = out_dir + "cobe2_model.nc"

                # extract the model data

                selection = all_df.query(
                    "variable == 'temperature'"
                ).model_variable.values[0]

                with warnings.catch_warnings(record=True) as w:
                    ds_model = nc.open_data(ersem_paths, checks=False)
                    # need a hack for mission atlantic data
                    if "MissionAtlantic" in ersem_paths[0]:
                        ds_model.nco_command("ncks -F -d deptht,1 -v thetao_con")
                    else:
                        ds_model.subset(variables=selection)
                        if surface == "top":
                            ds_model.top()
                        else:
                            ds_model.bottom()
                    ds_model.as_missing(0)
                    ds_model.tmean("month")
                    ds_model.merge("time")
                    ds_model.rename({selection: "model"})
                    ds_model.run()
                tidy_warnings(w)

                # extract the observation data

                ds_obs = nc.open_data(
                    f"{data_dir}/cobe2/sst.mon.mean.nc",
                    checks=False,
                )

                ds_obs.subset(variables="sst")
                ds_obs.top()
                ds_obs.subset(years=ds_model.years)
                ds_obs.rename({"sst": "observation"})
                # change temperature in ds_obs from Kelvin to Celsius

                ds_obs.run()

                ds_model.regrid(ds_obs)

                ds_model.append(ds_obs)
                ds_model.merge("variable", "month")

                ds_model.to_nc(out_file, zip=True, overwrite=True)

        if "alkalinity" in var_choice:
            vars = ["alkalinity"]

            out_dir = "matched/gridded/glodap/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                print("******************************")
                print("Matching the surface alkalinity")
                print("******************************")
                mapping = dict()
                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]
                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )
                df_times = pd.concat(df_times)

                ersem_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                # restrict to glodap climatology years

                ersem_paths = ersem_paths.query("year > 1971 and year < 2014")

                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                if spinup is not None:
                    min_year = ersem_paths.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start

                ersem_paths = ersem_paths.query("year >= @min_year")

                if end is not None:
                    max_year = end
                    ersem_paths = ersem_paths.query("year <= @max_year").reset_index(drop = True)

                if len(ersem_paths) > 0:
                    ersem_years = list(set(ersem_paths.year))

                    ersem_paths = ersem_paths.loc[:, ["path"]]

                    ersem_paths = list(ersem_paths.path)

                    out_file = out_dir + "glodap_alkalinity.nc"

                    # extract the model data

                    selection = df_mapping.query(
                        "variable == 'alkalinity'"
                    ).model_variable.values[0]

                    with warnings.catch_warnings(record=True) as w:
                        ds_model = nc.open_data(ersem_paths, checks=False)
                        # need a hack for mission atlantic data
                        if "MissionAtlantic" in ersem_paths[0]:
                            ds_model.nco_command("ncks -F -d deptht,1 -v thetao_con")
                        else:
                            ds_model.subset(variables=selection)
                            if surface == "top":
                                ds_model.top()
                            else:
                                ds_model.bottom()
                        ds_model.as_missing(0)
                        ds_model.tmean("month")
                        ds_model.merge("time")
                        ds_model.rename({selection: "model"})
                        ds_model.run()
                    tidy_warnings(w)

                    ds_model.tmean()
                    ds_model.run()

                    # extract the observation data
                    obs_file = f"{data_dir}/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TAlk.nc"

                    ds_obs = nc.open_data(obs_file, checks=False)

                    ds_obs.subset(variables="TAlk")
                    ds_obs.top()

                    ds_obs.rename({"TAlk": "observation"})
                    ds_obs.run()

                    ds_model.regrid(ds_obs)

                    ds_model.append(ds_obs)
                    ds_model.merge("variable")

                    ds_model.to_nc(out_file, zip=True, overwrite=True)
                else:
                    warnings.warn(
                        "No model files overlap with GLODAP climatology years"
                    )

        # matchup GLODAP pH
        if "ph" in var_choice:
            vars = ["ph"]

            out_dir = "matched/gridded/glodap/"

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            df = df_mapping.query("variable in @vars").reset_index(drop=True)
            if len(df) > 0:
                mapping = dict()

                for vv in df.variable:
                    mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

                selection = []
                for vv in vars:
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

                # extract the ERSEM paths
                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                import xarray as xr

                time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
                    0
                ]

                df_times = []

                for ff in paths:
                    ff_times = xr.open_dataset(ff)[time_name]
                    ff_month = [int(x.dt.month) for x in ff_times]
                    ff_year = [int(x.dt.year) for x in ff_times]
                    df_times.append(
                        pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
                            path=ff
                        )
                    )

                df_times = pd.concat(df_times)

                ersem_paths = (
                    df_times.loc[:, ["year", "month"]]
                    .drop_duplicates()
                    .groupby("year")
                    .count()
                    .query("month == 12")
                    .reset_index()
                    .merge(df_times, on="year", how="left")
                    .loc[:, ["year", "path"]]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )

                # restrict to glodap climatology years

                ersem_paths = ersem_paths.query("year > 1971 and year < 2014")

                if spinup is not None:
                    min_year = ersem_paths.year.min() + spinup
                if start is not None:
                    if spinup is None:
                        min_year = start
                ersem_paths = ersem_paths.query(
                    "year >= @sim_start and year <= @sim_end"
                ).reset_index(drop=True)

                if end is not None:
                    max_year = end
                    ersem_paths = ersem_paths.query("year <= @max_year").reset_index(drop = True)

                ersem_paths = ersem_paths.query("year >= @min_year")

                if len(ersem_paths) > 0:
                    ersem_years = list(set(ersem_paths.year))

                    ersem_paths = ersem_paths.loc[:, ["path"]]

                    ersem_paths = list(ersem_paths.path)

                    out_file = out_dir + "glodap_ph.nc"

                    # extract the model data

                    selection = all_df.query("variable == 'ph'").model_variable.values[
                        0
                    ]

                    with warnings.catch_warnings(record=True) as w:
                        ds_model = nc.open_data(ersem_paths, checks=False)
                        # need a hack for mission atlantic data
                        if "MissionAtlantic" in ersem_paths[0]:
                            ds_model.nco_command("ncks -F -d deptht,1 -v thetao_con")
                        else:
                            ds_model.subset(variables=selection)
                            if surface == "top":
                                ds_model.top()
                            else:
                                ds_model.bottom()
                        ds_model.as_missing(0)
                        ds_model.tmean("month")
                        ds_model.merge("time")
                        ds_model.rename({selection: "model"})
                        ds_model.run()
                    tidy_warnings(w)

                    ds_model.tmean()
                    ds_model.run()

                    # extract the observation data
                    obs_file = f"{data_dir}/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.pHtsinsitutp.nc"

                    ds_obs = nc.open_data(obs_file, checks=False)

                    ds_obs.subset(variables="pHtsinsitutp")
                    ds_obs.top()

                    ds_obs.rename({"pHtsinsitutp": "observation"})
                    ds_obs.run()

                    ds_model.regrid(ds_obs)

                    ds_model.append(ds_obs)
                    ds_model.merge("variable")

                    ds_model.to_nc(out_file, zip=True, overwrite=True)

    os.system("pandoc matchup_report.md -o matchup_report.pdf")
