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
import numpy as np
import xarray as xr

from multiprocessing import Manager
from tqdm import tqdm
from ecoval.utils import get_datadir, session
from ecoval.utils import extension_of_directory
from ecoval.ices import generate_mapping
from ecoval.nsbc import gridded_matchup
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
    data_dir = get_datadir()
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
    tidy_warnings(w)

    df_out.rename(columns={variable: "model"}, inplace=True)
    df_out = df_out.merge(df_wod)

    df_all.append(df_out)


def mm_match(
    ff, ersem_variable, df, df_times, ds_depths, ices_variable, df_all, top_layer=False
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
            ds.as_missing(0)
            ds.run()
            if len(var_match) > 1:
                ds.sum_all()
            valid_locs = ["lon", "lat", "year", "month", "day", "depth"]
            valid_locs = [x for x in valid_locs if x in df.columns]

            valid_times = "year" in df.columns or "month" in df.columns or "day" in df.columns

            if valid_times: 
                df_locs = (
                    df_times.query("path == @ff")
                    .merge(df)
                    .loc[:, valid_locs]
                    .drop_duplicates()
                    .reset_index(drop=True)
                )
            else:
                df_locs = (df
                           .loc[:, valid_locs]
                )

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
                    df_ff = ds.match_points( df_locs, quiet=True, top=top_layer)
                else:
                    df_ff = ds.match_points( df_locs, depths=ds_depths, quiet=True, top=top_layer)
                valid_vars = ["lon", "lat", "year", "month", "day", "depth", ds.variables[0]]
                valid_vars = [x for x in valid_vars if x in df_ff.columns]
                df_ff = df_ff.loc[:, valid_vars]
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
        try:
            ds_dict = generate_mapping(ds, fvcom=fvcom)
        # output error and ff
        except:
            pass

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


def matchup(
    folder=None,
    spinup=None,
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
        "poc",
    ],
    bottom= ["ph", "oxygen"],
    benthic = ["carbon"],
    cores=None,
    e3t=None,
    mapping=None,
    start=None,
    end=None,
    levels=2,
    lon_lim=None,
    lat_lim=None,
    exclude=[],
    fvcom=False,
    point_surface = [],
    strict = True,
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
    fvcom: bool
        If True, the matchups are assumed to be FVCOM
    e3t: str
        Path to e3t file. If None, e3t will be assumed to be in the files
    start: None
        Start year for matchups
    end: None
        End year for matchups
    levels: int
        Number of directories to go down to find the files
        0 means the files are in that directory. 1 means they are in a subdirectory. 2 means they are in a subsubdirectory
    **kwargs:
        Additional arguments to pass to matchup


    """

    add_point_surface = []
    if len(point_surface) > 0:
        if isinstance(point_surface, str):
            point_surface = [point_surface]
        add_point_surface = copy.deepcopy(point_surface)

    # add not implemented error if lon_lim or lat_lim is not None
    if lon_lim is not None:
        raise NotImplementedError("lon_lim not implemented")
    if lat_lim is not None:
        raise NotImplementedError("lat_lim not implemented")

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
    

    if surface is None:
        surface = []
    if isinstance(surface, str):
        surface = [surface]
    if isinstance(bottom, str):
        bottom = [bottom]
    if bottom is None:
        bottom = []
    var_choice = surface + bottom
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
        all_df = find_paths(folder, fvcom=fvcom, exclude=exclude)

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

    if fvcom:

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
                final_extension = extension_of_directory(folder)
                ersem_paths = glob.glob(folder + final_extension + pattern)
                if len(ersem_paths) > 0:
                    good_to_go = True

            if good_to_go:
                final_extension = extension_of_directory(folder)
                ersem_paths = glob.glob(folder + final_extension + pattern)

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

                    if vv == "doc":
                        command = "-aexpr,doc=" + model_variables[0]
                        ds2.cdo_command(command)
                        drop_these = model_variables[0].split("+")
                        ds_contents = ds2.contents
                        ds_contents = ds_contents.query("variable in @drop_these")
                        doc_unit = ds_contents.unit[0]
                        ds2.set_units({"doc": doc_unit})
                        ds2.drop(variables=drop_these)

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

    # get the units. File inspection could be randomized in case people have put loose files in there...
    import glob

    print("********************************")
    print("Identifying whether it is a northwest European shelf domain")
    print("********************************")
    df = pd.read_csv("matched/mapping.csv")
    df = df.dropna()
    df = df.iloc[0:1, :]
    pattern = list(df.pattern)[0]
    pattern = pattern.replace("//", "/")

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
    if lon_max > 50:
        global_grid = True
    if global_grid:
        print("Grid is not NW European shelf")
    else:
        print("Grid is NW European shelf")

    if global_grid:
        model_domain = "global"
    else:
        model_domain = "nws"

    surf_all = False
    if surface == ["all"]:
        surface = all_vars
        surf_all = True

    if "ph" in surface and model_domain == "nws":
        surface.remove("ph")
        point_surface.append("ph")
    # do the same for alkalinity, poc and doc
    if "poc" in surface and model_domain == "nws":
        point_surface.append("poc")
    if "doc" in surface and model_domain == "nws":
        point_surface.append("doc")

    if type(surface) is str:
        surface = [surface]

    for vv in surface:
        if not os.path.exists(f"{data_dir}/gridded/{model_domain}/{vv}"):
            if not os.path.exists(f"{data_dir}/gridded/globl/{vv}"):
                surface.remove(vv)
    
    if surf_all:
        point_surface.append("ph")
        point_surface.append("poc")
        point_surface.append("doc")
        point_surface.append("alkalinity")

    point_surface = list(set(point_surface))

    print("********************************")
    if True:
        # figure out the lon/lat extent in the model
        lons = [lon_min, lon_max]
        lats = [lat_min, lat_max]
        extent = [lons, lats]
        # start of with the raw coords
        # This will not work with nemo, which outputs the grid incorrectly
        # so we will check if the step between the first lon/lat and the second lon/lat is
        # far bigger than the rest. If this is the case, the first should be ignored
        # get the lon/lat values
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

    for x in add_point_surface:
        point_surface.append(x)
    point_surface = list(set(point_surface))

    # get rid of any rows where pattern is None

    # if True:

    # # if "temperature" in var_choice and False:
    #     vars = ["temperature"]

    #     out_dir = "matched/gridded/cobe2/"

    #     if not os.path.exists(out_dir):
    #         os.makedirs(out_dir)

    #     df = df_mapping.query("variable in @vars").reset_index(drop=True)
    #     if len(df) > 0:
    #         print("******************************")
    #         print("Matching the vertically resolved temperature with NOAA World Ocean Database")
    #         mapping = dict()
    #         for vv in df.variable:
    #             mapping[vv] = list(df.query("variable == @vv").model_variable)[0]

    #         selection = []
    #         for vv in vars:
    #             try:
    #                 selection += mapping[vv].split("+")
    #             except:
    #                 selection = selection

    #         patterns = set(df.pattern)
    #         if len(patterns) > 1:
    #             raise ValueError(
    #                 "Something strange going on in the string patterns. Unable to handle this. Bug fix time!"
    #             )
    #         pattern = list(patterns)[0]

    #         # extract the ERSEM paths
    #         final_extension = extension_of_directory(folder)
    #         paths = glob.glob(folder + final_extension + pattern)

    #         for exc in exclude:
    #             paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

    #         import xarray as xr

    #         time_name = [x for x in xr.open_dataset(paths[0]).dims if "time" in x][
    #             0
    #         ]
    #         df_times = []

    #         for ff in paths:
    #             ff_times = xr.open_dataset(ff)[time_name]
    #             ff_month = [int(x.dt.month) for x in ff_times]
    #             ff_year = [int(x.dt.year) for x in ff_times]
    #             df_times.append(
    #                 pd.DataFrame({"year": ff_year, "month": ff_month}).assign(
    #                     path=ff
    #                 )
    #             )
    #         df_times = pd.concat(df_times)

    #         ersem_paths = (
    #             df_times.loc[:, ["year", "month"]]
    #             .drop_duplicates()
    #             .groupby("year")
    #             .count()
    #             .query("month == 12")
    #             .reset_index()
    #             .merge(df_times, on="year", how="left")
    #             .loc[:, ["year", "path"]]
    #             .drop_duplicates()
    #             .reset_index(drop=True)
    #         )

    #         if spinup is not None:
    #             min_year = ersem_paths.year.min() + spinup
    #         if start is not None:
    #             if spinup is None:
    #                 min_year = start

    #         ersem_paths = ersem_paths.query("year >= @min_year")
    #         ersem_paths = ersem_paths.query(
    #             "year >= @sim_start and year <= @sim_end"
    #         ).reset_index(drop=True)

    #         if end is not None:
    #             max_year = end
    #             ersem_paths = ersem_paths.query("year <= @max_year").reset_index(drop = True)

    #         ersem_years = list(set(ersem_paths.year))

    #         ersem_paths = ersem_paths.loc[:, ["path"]]

    #         ersem_paths = list(ersem_paths.path)

    #         wod_output = "matched/wod/wod_matchups.csv"

    #         df_all = manager.list()

    #         pool = multiprocessing.Pool(cores)

    #         pbar = tqdm(total=len(ersem_paths), position=0, leave=True)

    #         var_sel = list(df.model_variable)[0]
    #         results = dict()
    #         for ff in ersem_paths:
    #             # matchup_wod(ff, var_sel, df_all, "fixed")
    #             temp = pool.apply_async(
    #                 matchup_wod,
    #                 [
    #                     ff,
    #                     var_sel,
    #                     df_all,
    #                     "fixed"

    #                 ],
    #             )
    #             results[ff] = temp

    #         for k, v in results.items():
    #             value = v.get()
    #             pbar.update(1)

    #         df_all = list(df_all)

    #         df_all = pd.concat(df_all)
    #         # make sure directory exists for wod_output
    #         if not os.path.exists(os.path.dirname(wod_output)):
    #             os.makedirs(os.path.dirname(wod_output))

    #         df_all.to_csv(wod_output, index=False)

    # raise ValueError("here")

    if len(point_all) > 0 or len(point_bottom) > 0:
        print("Matching up with observational point data")
        print("********************************")

    df_mapping = all_df

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
            all_df = all_df.query("variable == @vv").reset_index(drop=True)
            patterns = list(set(all_df.pattern))

            for pattern in patterns:
                final_extension = extension_of_directory(folder)
                ensemble = glob.glob(folder + final_extension + pattern)
                for exc in exclude:
                    ensemble = [
                        x for x in ensemble if f"{exc}" not in os.path.basename(x)
                    ]

                ds = xr.open_dataset(ensemble[0])
                time_name = [x for x in list(ds.dims) if "time" in x][0]

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

                df_times = df_times.query("year >= @min_year").reset_index(drop=True)

                # ersem paths

                ersem_paths = list(set(df_times.path))
                ersem_paths.sort()
                # write to the report

                write_report(f"### Matchup summary for observational point data")
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

                if depths != "surface":
                    with warnings.catch_warnings(record=True) as w:
                        # extract the thickness dataset
                        if e3t is not None:
                            ds_thickness = nc.open_data(e3t, checks=False)
                            if "e3t" not in ds_thickness.variables:
                                options = [x for x in ds_thickness.variables if "e3t" in x]
                                if len(options) != 1:
                                    raise ValueError("e3t not found in e3t file")
                                ds_thickness.rename({options[0]: "e3t"})
                        else:
                            ds_thickness = nc.open_data(ensemble[0], checks=False)

                        ds_thickness.subset(time=0, variables="e3t")
                        ds_thickness.as_missing(0)
                        #####

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

                    tidy_warnings(w)
                
                def point_match(variable, layer="all", ds_depths = None):
                    with warnings.catch_warnings(record=True) as w:
                        point_variable = variable
                        ersem_variable = list(
                            all_df.query("variable == @point_variable").model_variable
                        )[0]
                        paths = glob.glob(
                            f"{data_dir}/point/nws/**/{variable}/**{variable}**.csv"
                        )
                        source = os.path.basename(paths[0]).split("_")[0]
                        if depths == "surface":
                            paths = [x for x in paths if "all" in x]
                        else:
                            paths = [x for x in paths if depths in x]

                        paths = [x for x in paths if f"{point_variable}/" in x]
                        for exc in exclude:
                            paths = [
                                x for x in paths if f"{exc}" not in os.path.basename(x)
                            ]


                        df = pd.concat([pd.read_csv(x) for x in paths])

                        if variable == "doc":
                            # go from mole to g of C
                            df = df.assign(observation = lambda x: x.observation * 12.011)
                        if not strict:
                            if "year" in df.columns:
                                df = df.drop(columns = "year")
                        if depths == "surface":
                            df = df.query("depth < 5").reset_index(drop=True)
                            # drop depth
                            df = df.drop(columns="depth")
                            # add in a nominal depth
                            # df = df.assign(depth=2.5)
                        # restrict the lon_lat
                        lon_min = lons[0]
                        lon_max = lons[1]
                        lat_min = lats[0]
                        lat_max = lats[1]
                        df = df.query(
                            "lon >= @lon_min and lon <= @lon_max and lat >= @lat_min and lat <= @lat_max"
                        ).reset_index(drop=True)

                        if variable == "temperature":
                            df_include = pd.read_csv(
                                f"{data_dir}/point/nws/mld_profiles.csv"
                            )
                            df = df.merge(df_include).reset_index(drop=True)
                        sel_these = ["year", "month", "day"]
                        if not strict:
                            sel_these = ["month", "day"]
                        if variable != "carbon":
                            paths = list(
                                set(
                                    df.loc[:, sel_these]
                                    .drop_duplicates()
                                    .merge(df_times)
                                    .path
                                )
                            )
                        paths = df_times.path
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
                        lon_name = [x for x in list(ds_xr.coords) if "lon" in x][0]
                        lon_min = ds_xr[lon_name].values.min()
                        lon_max = ds_xr[lon_name].values.max()
                        lat_name = [x for x in list(ds_xr.coords) if "lat" in x][0]
                        lat_min = ds_xr[lat_name].values.min()
                        lat_max = ds_xr[lat_name].values.max()
                        df = df.query(
                            "lon >= @lon_min and lon <= @lon_max and lat >= @lat_min and lat <= @lat_max"
                        ).reset_index(drop=True)
                    tidy_warnings(w)

                    # valid_cols = ["lon", "lat",	"day"	month	year	depth	model	observation	
                    valid_cols = ["lon", "lat", "day", "month", "year", "depth", "observation"]
                    select_these = [x for x in df.columns if x in valid_cols]
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
                                        ds_grid.to_dataframe().reset_index().dropna()
                                    )
                                    columns = [
                                        x
                                        for x in df_grid.columns
                                        if "lon" in x or "lat" in x
                                    ]
                                    df_grid = df_grid.loc[:, columns].drop_duplicates()
                                    if not os.path.exists("matched"):
                                        os.makedirs("matched")
                                    df_grid.to_csv(
                                        "matched/model_grid.csv", index=False
                                    )
                                tidy_warnings(w)

                        grid_setup = True
                        if layer == "surface":
                            top_layer = True
                        else:
                            top_layer = False
                        if depths == "surface":
                            ds_depths = None
                        # raise ValueError("stoping")

    #ff, ersem_variable, df, df_times, ds_depths, ices_variable, df_all, top_layer=False
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
                    change_this = [x for x in df_all.columns if x not in ["lon", "lat", "year", "month", "day", "depth", "observation"]][0]
                    df_all = df_all.rename(columns={change_this: "model"}).merge(
                        df
                    )
                    df_all = df_all.dropna().reset_index(drop=True)
                    
                    out = f"matched/point/{model_domain}/{depths}/{variable}/{source}_{depths}_{variable}.csv"
                    # create directory for out if it does not exists
                    if not os.path.exists(os.path.dirname(out)):
                        os.makedirs(os.path.dirname(out))
                    out1 = out.replace(os.path.basename(out), "paths.csv")
                    pd.DataFrame({"path": paths}).to_csv(out1, index=False)
                    if variable == "doc":
                        df_all = df_all.assign(model = lambda x: x.model + 40*12.011)
                    df_all.to_csv(out, index=False)

                print("**********************")
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
                print("**********************")
                if depths == "surface":
                    point_match(vv, layer="surface")
                else:
                    point_match(vv, ds_depths = ds_depths)

        # def nsbc_matchup(df_mapping = None, folder = None, var_choice = None, exclude = None, surface = None, start = None, spinup = None, sim_start = None, sim_end = None, e3t = None, report = None):

    gridded_matchup(
        df_mapping=df_mapping,
        folder=folder,
        var_choice=surface,
        exclude=exclude,
        surface=surface_level,
        start=start,
        spinup=spinup,
        sim_start=sim_start,
        sim_end=sim_end,
        e3t=e3t,
        domain=model_domain,
        strict= strict
    )

    os.system("pandoc matchup_report.md -o matchup_report.pdf")
