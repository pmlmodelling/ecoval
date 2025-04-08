import nctoolkit as nc
import copy
import re
import glob
import multiprocessing as mp
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
from ecoval.parsers import generate_mapping


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

    path = x

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



def na_tally(ff, var, the_list, the_list_max, the_list_min):
    # catch warnings
    with warnings.catch_warnings(record=True) as w:
        ds = nc.open_data(ff, checks = False)
        ds.subset(variable = var) 
        ds.as_missing(0)
        ds > - 100000
        ds.spatial_sum()
        ds.vertical_sum()
        ds.run()
        df = ds.to_dataframe()
        time_name = [x for x in df.columns if "time" in x][0]
        df = df.rename(columns={time_name: "time"})
        df = df.loc[:,["time", var]]
        # add the file name
        df = df.assign(path = ff)
        the_list.append(df)
        # now to the maximum values
        ds = nc.open_data(ff, checks = False)
        ds.subset(variable = var)
        ds.as_missing(0)
        ds.spatial_max()
        ds.vertical_max()
        df = ds.to_dataframe()
        time_name = [x for x in df.columns if "time" in x][0]
        df = df.rename(columns={time_name: "time"})
        df = df.loc[:,["time", var]]
        df = df.rename(columns={var: "max_value"}) 
        # add the file name
        df = df.assign(path = ff)
        the_list_max.append(df)

        ds = nc.open_data(ff, checks = False)
        ds.subset(variable = var)
        ds.as_missing(0)
        ds.spatial_min()
        ds.vertical_min()
        df = ds.to_dataframe()
        time_name = [x for x in df.columns if "time" in x][0]
        df = df.rename(columns={time_name: "time"})
        df = df.loc[:,["time", var]]
        df = df.rename(columns={var: "min_value"})
        # add the file name
        df = df.assign(path = ff)
        # add the file name
        the_list_min.append(df)


    # return df

def check_simulation(
    folder=None,
    variables = "P1_Chl",
    cores = 1,
    n_dirs_down = 2,
    exclude = []
):
    """
    
    Parameters
    -------------
    folder: str
        Folder containing model output
    variables: str
        Variables to check in the model output
    cores: int
        Number of cores to use for parallel processing
    n_dirs_down: int
        Number of directories to go down to find the files
    exclude: list
        List of strings to exclude from the file search

    """

    if isinstance(exclude, str):
        exclude = [exclude]

    # create an empty markdown document in "simulation_checks.md"
    with open("simulation_checks.md", "w") as f:
        f.write("# Simulation Checks\n")

    if isinstance(variables, str):
        variables = [variables]

    for vv in variables:
        while True:
            levels = n_dirs_down 

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
                #options = [x for x in options if "part" not in os.path.basename(x)]
                options = [x for x in options if "restart" not in os.path.basename(x)]

            if len([x for x in options if ".nc" in x]) > 0:
                break
        df_pattern = []
        options = list(set(options))
        for ff in options:
            # print(ff)
            ds_ff = nc.open_data(ff, checks = False)
            the_vars = ds_ff.variables
            if vv in the_vars:
                ff_pattern = os.path.basename(ff)
                # replace integers with 4 or more digits with **
                ff_pattern = re.sub(r"\d{4,}", "**", ff_pattern)
                # replace strings of the form _12. with _**.
                ff_pattern = re.sub(r"\d{2,}", "**", ff_pattern)
                ff_res = get_time_res(ff, folder)
                pattern = ff_pattern
                df_pattern.append(pd.DataFrame({"pattern": pattern, "time_res": ff_res}, index=[0]))
        df_pattern = pd.concat(df_pattern)
        # sort by time_res
        df_pattern = df_pattern.sort_values(by=["time_res"], ascending=True).reset_index(drop = True)

        search_pattern = folder + "/"
        for i in range(n_dirs_down):
            search_pattern = search_pattern + "/**"
        pattner = df_pattern.reset_index(drop = True).pattern[0]
        search_pattern = search_pattern + "/" + pattern
        ensemble = glob.glob(search_pattern, recursive = True)
        #  remove anything with restart
        ensemble = [x for x in ensemble if "restart" not in x]
        # nothing with spinup
        ensemble = [x for x in ensemble if "spinup" not in x]
        for exc in exclude:
            ensemble = [x for x in ensemble if exc not in x]
        ensemble = list(set(ensemble))
        paths = ensemble

        nc.options(parallel = True)
        with open("simulation_checks.md", "a") as f:
            text = f"## Checks for {vv}\n\n"
            f.write(text + "\n")

        print(f"Identifying whether number of missing values is consistent in all files for {vv}")
        the_list = mp.Manager().list()
        the_list_max = mp.Manager().list()
        the_list_min = mp.Manager().list()
        pool = mp.Pool(cores)
        for ff in paths:
            # na_tally(ff, vv, the_list)
            pool.apply_async(na_tally, [ff, vv, the_list, the_list_max, the_list_min])
        pool.close()
        pool.join()
        the_list = list(the_list)
        df_missing = pd.concat(the_list)
        # sort by time
        df_missing = df_missing.sort_values(by=["time"], ascending=True).reset_index(drop = True)
        min_time = df_missing["time"].min()
        max_time = df_missing["time"].max()
        min_date = pd.to_datetime(min_time).strftime("%d/%m/%Y")
        max_date = pd.to_datetime(max_time).strftime("%d/%m/%Y")
        if df_pattern.time_res[0] == "d":
            # create sequence of dates, 1 day apart
            date_range = pd.date_range(start=min_time, end=max_time, freq="D")
            n_dates = len(date_range)
            # check if any of these dates are missing in df_missing.time
            missing_dates = []
            for i in range(n_dates):
                if date_range[i] not in df_missing.time.values:
                    missing_dates.append(date_range[i])
            if len(missing_dates) > 0:
                text = f"Potential missing dates in {vv}:\n\n"
                for i in range(len(missing_dates)):
                    # human readable
                    missing_dates[i] = pd.to_datetime(missing_dates[i]).strftime("%d/%m/%Y")
                    text = text + str(missing_dates[i]) + "\n"
                with open("simulation_checks.md", "a") as f:
                    f.write(text + "\n")
            else:
                text = f"No missing dates in {vv} from {min_date} to {max_date}\n\n"
                with open("simulation_checks.md", "a") as f:
                    f.write(text + "\n")




        if float(df_missing[vv].max() - df_missing[vv].min()) == 0:
            text = f"All files have the same number of missing values for {vv}\n\n" 
            with open("simulation_checks.md", "a") as f:
                f.write(text + "\n")
        df_max = pd.concat(the_list_max)
        # sort by max_value
        df_max = df_max.sort_values(by=["max_value"], ascending=False).reset_index(drop = True)
        # select the first row
        df_max = df_max.loc[df_max.index[0]]
        # extract the maximum value
        max_value = float(df_max["max_value"])
        # get this to 5 significant figures
        max_value = "{:.5g}".format(max_value)
        max_date = df_max["time"]
        # make this more human readable, i.e. DD/MM/YYYY
        max_date = pd.to_datetime(max_date).strftime("%d/%m/%Y")
        # extract the file name
        file_name = df_max["path"]
        # remove the folder name
        file_name = file_name.replace(folder, "")
        text = f"Maximum value for {vv} is {max_value} on {max_date} in file:\n\n{file_name}"
        with open("simulation_checks.md", "a") as f:
            f.write(text + "\n\n")
        
        df_min = pd.concat(the_list_min)
        # sort by min_value
        df_min = df_min.sort_values(by=["min_value"], ascending=True).reset_index(drop = True)
        # select the first row
        df_min = df_min.loc[df_min.index[0]]
        # extract the minimum value
        min_value = float(df_min["min_value"])
        # get this to 5 significant figures
        min_value = "{:.5g}".format(min_value)
        min_date = df_min["time"]
        # make this more human readable, i.e. DD/MM/YYYY
        min_date = pd.to_datetime(min_date).strftime("%d/%m/%Y")
        # extract the file name
        file_name = df_min["path"]
        # remove the folder name
        file_name = file_name.replace(folder, "")
        text = f"Minimum value for {vv} is {min_value} on {min_date} in file:\n\n{file_name}"
        with open("simulation_checks.md", "a") as f:
            f.write(text + "\n")
        # add a line break
        with open("simulation_checks.md", "a") as f:
            f.write("\n\n")


        
        nc.options(parallel = False)
    
    # convert the markdown file to pdf, ensure words wrap
    # os.system("pandoc simulation_checks.md -o simulation_checks.pdf --pdf-engine=xelatex --variable mainfont='Arial' --variable fontsize=12pt --variable geometry:margin=1in --variable geometry:top=1in --variable geometry:bottom=1in --variable geometry:left=1in --variable geometry:right=1in")
    os.system("pandoc simulation_checks.md -o simulation_checks.pdf")
    # open the file
    os.system("open simulation_checks.pdf")
