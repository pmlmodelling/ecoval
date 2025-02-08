import copy
import time
import nctoolkit as nc
import re
import shutil
import glob
import multiprocessing
import pathlib
import os
import pandas as pd
import string
import molmass
def get_molar_mass(element):
    from molmass import Formula
    f = Formula(element)
    return f.mass
import random
import warnings
import pickle
import xarray as xr
import pkg_resources
from ecoval.session import session_info
from multiprocessing import Manager
from ecoval.utils import extension_of_directory, get_extent, fvcom_regrid
from ecoval.parsers import generate_mapping
from ecoval.gridded import gridded_matchup

def simulation_differences_comparison():
    """
    Compare the model output to observations for a range of variables.
    """

    shutil.copyfile(pkg_resources.resource_filename(__name__, "data/pml_logo.jpg"), f"pml_logo.jpg")

    if os.path.exists("book_comparison"):
        #get user input to decide if it should be removed
        #user_input = input("book directory already exists. This will be emptied and replaced. Do you want to proceed? (y/n): ")
        user_input = "y"
        if user_input.lower() == "y":
            while True:
                files = glob.glob("book_comparison/**/**/**", recursive=True)
                # list all files in book, recursively
                for ff in files:
                    if ff.startswith("book_comparison"):
                        try:
                            os.remove(ff)
                        except:
                            pass
                files = glob.glob("book_comparison/**/**/**", recursive=True)
                # only list files
                files = [x for x in files if os.path.isfile(x)]
                if len(files) == 0:
                    break
        else:
            print("Exiting")
            return None
                # if not os.path.exists(book_dir):
                    # break
    # create the directory
    if not os.path.exists("book_comparison"):
        os.makedirs("book_comparison")
    if not os.path.exists("book_comparison/notebooks"):
        os.makedirs("book_comparison/notebooks")

    # list files with simdiff
    # use pkg_resources to get all available data

    # list paths in "../../data/"
    file_paths = glob.glob("data/climatologies/**/**/**", recursive=True)
    # must be ".nc"
    file_paths = [x for x in file_paths if ".nc" in x]


    if len([x for x in file_paths if "phenology-clim" not in x]) > 0:
        data_path = pkg_resources.resource_filename(__name__, "data/simdiff_spatial_summaries.ipynb")
        if not os.path.exists(f"book_comparison/notebooks/simdiff_spatial_summaries.ipynb"):
            shutil.copyfile(data_path, "book_comparison/notebooks/simdiff_spatial_summaries.ipynb")

        # copy the mapped_differences notebook
        data_path = pkg_resources.resource_filename(__name__, "data/simdiff_maps_notebook.ipynb")
        if not os.path.exists(f"book_comparison/notebooks/simdiff_maps_notebook.ipynb"):
            shutil.copyfile(data_path, "book_comparison/notebooks/simdiff_maps_notebook.ipynb")
        # copy the temporal_differences notebook
        data_path = pkg_resources.resource_filename(__name__, "data/simdiff_temporal_differences.ipynb")
        if not os.path.exists(f"book_comparison/notebooks/simdiff_temporal_differences.ipynb"):
            shutil.copyfile(data_path, "book_comparison/notebooks/simdiff_temporal_differences.ipynb")
    # move the phenology notebook
    if len([x for x in file_paths if "phenology-clim" in x]) > 0:
        data_path = pkg_resources.resource_filename(__name__, "data/simdiff_phenology.ipynb")
        if not os.path.exists(f"book_comparison/notebooks/simdiff_phenology.ipynb"):
            shutil.copyfile(data_path, "book_comparison/notebooks/simdiff_phenology.ipynb")
    # move the depth_profile notebook
    if len([x for x in file_paths if "depth_profile" in x]) > 0:
        data_path = pkg_resources.resource_filename(__name__, "data/simdiff_depth_profile.ipynb")
        if not os.path.exists(f"book_comparison/notebooks/simdiff_depth_profile.ipynb"):
            shutil.copyfile(data_path, "book_comparison/notebooks/simdiff_depth_profile.ipynb")

    data_path = pkg_resources.resource_filename(__name__, "data/_toc.yml")

    out = "book_comparison/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out) 

    shutil.copyfile(data_path, out)

    data_path = pkg_resources.resource_filename(__name__, "data/intro_compare_simulations.md")

    out = "book_comparison/" + "intro.md" 

    shutil.copyfile(data_path, out)

    data_path = pkg_resources.resource_filename(__name__, "data/_config.yml")

    out = "book_comparison/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out)

    os.system("jupyter-book build book_comparison/")
    import webbrowser

    webbrowser.open("file://" + os.path.abspath("book_comparison/_build/html/index.html"))


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


random_files = []


surface_variables = {"benthic_biomass":"Y2_c+Y3_c"}
surface_variables["co2_flux"] = "O3_fair"

surface_variables["chlorophyll"] = "P1_Chl+P2_Chl+P3_Chl+P4_Chl"
bottom_variables = {"oxygen":"O2_o", "pH":"O3_pH"},
integrated_variables = {"phosphate": "N1_p", "nitrate": "N3_n", "silicate": "N5_s", "doc":"R1_c+R2_c+R3_c"}
integrated_variables["poc"] = "P1_c+P2_c+P3_c+P4_c+Z5_c+Z6_c+R4_c+R6_c+R8_c"
integrated_variables["ammonium"] = "N4_n" 
vertmean_variables = {"oxygen":"O2_o"}
vertmean_variables["alkalinity"] = "O3_TA"
vertmean_variables["light_attenuation"] = "light_xEPS"

#  


def compare_simulations(
    sim_1 =None,
    sim_2 =None,
    surface_variables = None,
#     {'chlorophyll': 'P1_Chl+P2_Chl+P3_Chl+P4_Chl',
#  'oxygen': 'O2_o',
#  'phosphate': 'N1_p',
#  'nitrate': 'N3_n',
#  'silicate': 'N5_s',
#  'benthic_biomass': 'Y2_c+Y3_c',
#  'ph': 'O3_pH',
#  'co2_flux': 'O3_fair',
#  'ammonium': 'N4_n',
#  'pCO2': 'O3_pCO2',
#  'nano phytoplankton chlorophyll': 'P2_Chl',
#  'micro phytoplankton chlorophyll': 'P1_Chl+P4_Chl',
#  'pico phytoplankton chlorophyll': 'P3_Chl',
#  'mesozooplankton': 'Z4_c'},
    bottom_variables = None, 
    vertmean_variables = None, 
    integrated_variables = None, 
#  'POC': 'P1_c+P2_c+P3_c+P4_c+Z5_c+Z6_c+R4_c+R6_c+R8_c',
#  'mesozooplankton': 'Z4_c'},
    phenology = None, 
    depth_profile = None,
    start=None,
    end=None,
    surface_level=None,
    cores=None,
    thickness_files=[None, None],
    exclude=[],
    n_dirs_down=2,
    lon_lim=None,
    lat_lim=None,
    n_check=None,
    overwrite=True,
    ask=True,
    out_dir="",
    build = True,
    **kwargs,
):
    """
    Directly compare two simulations 

    Parameters
    -------------

    sim_dir_1 : str
        The directory containing the first model directory  
    sim_dir_2 : str
        The directory containing the second simulation
    surface_variables : dict
        A dictionary of variables to compare at the surface
    bottom_variables : dict
        A dictionary of variables to compare at the bottom
    vertmean_variables : dict
        A dictionary of variables to compare vertically averaged
    integrated_variables : dict
        A dictionary of variables to compare vertically integrated
    phenology : dict
        A dictionary of variables to compare phenology
    depth_profile : dict
        A dictionary of variables to compare depth profiles
    start : int
        The start year
    end : int
        The end year
    surface_level : str
        The level of the surface
    cores : int
        The number of cores to use
    thickness_files : list  
        A list of thickness files
    exclude : list  
        A list of strings to exclude
    n_dirs_down : int   
        The number of directories down
    lon_lim : list  
        A list of longitude limits
    lat_lim : list  
        A list of latitude limits
    n_check : int
        The number of files to check
    overwrite : bool
        Whether to overwrite existing files
    out_dir : str
        The output directory
    build : bool
        Whether to build the book
    kwargs : dict
        Additional arguments


    Returns
    -------------
    None
    Data will be stored in the matched directory.

    Example
    -------------

    If you wanted to compare surface chlorophyll, bottom oxygen, and integrated phosphate, nitrate, and silicate, the phenology of chlorophyll, and the depth profile of oxygen, you would run:  

    ecoval.compare_simulations(
        sim_dir_1 = "/data/foo",
        sim_dir_2 = "/data/bar",
        surface_variables = {"chlorophyll": "P1_Chl+P2_Chl+P3_Chl+P4_Chl"},
        bottom_variables = {"oxygen": "O2_o", "pH": "O3_pH"},
        vertmean_variables = {"oxygen": "O2_o"},
        integrated_variables = {"phosphate": "N1_p", "nitrate": "N3_n", "silicate": "N5_s"},
        phenology = {"chlorophyll": "P1_Chl+P2_Chl+P3_Chl+P4_Chl"},
        depth_profile = {"oxygen": "O2_o"},
        start = 2021,
        end = 2022,
        surface_level = "top",
        cores = 1,
        thickness_files = [None, None],
        exclude = [],
        n_dirs_down = 2,
        lon_lim = None,
        lat_lim = None,
        n_check = None,
        overwrite = True,
        out_dir = "",
        build = True
    )


    """

    # if any dictionarys are None, set them to empty dictionaries
    if surface_variables is None:
        surface_variables = dict()
    if bottom_variables is None:
        bottom_variables = dict()
    if vertmean_variables is None:
        vertmean_variables = dict()
    if integrated_variables is None:
        integrated_variables = dict()
    if phenology is None:
        phenology = dict()
    if depth_profile is None:
        depth_profile = dict()

    mass = None


    # check if sim_1 is a dict
    if isinstance(sim_1, dict):
        sim_dir_1 = list(sim_1.values())[0]
        sim_name_1 = list(sim_1.keys())[0]
    else: 
        raise ValueError("sim_1 must be a dict")
    if isinstance(sim_2, dict):
        sim_dir_2 = list(sim_2.values())[0]
        sim_name_2 = list(sim_2.keys())[0]
    else:
        raise ValueError("sim_2 must be a dict")
    
    
    sim_dict = {"sim0": sim_name_1, "sim1": sim_name_2}
    # save this as a pickle
    with open("sim_dict.pkl", "wb") as f:
        pickle.dump(sim_dict, f) 


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
        raise ValueError(f"{sim_dir_1} does not exist")
    
    if sim_dir_2 is None:
        raise ValueError("Please provide a sim_dir_2 directory")

    if sim_dir_2 is not None:
        if not os.path.exists(sim_dir_2):
            raise ValueError(f"{sim_dir_2} does not exist")

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


    # restrict surface to valids

    # surface = [x for x in surface if x in valid_surface]

    ff = "times_dict.pkl"
    all_times = dict()
    if os.path.exists(ff):
        all_times = pickle.load(open(ff, "rb"))

    for jj in range(6):
        if jj == 0:
            var_dict = copy.deepcopy(surface_variables)
            measure = "top"
        if jj == 1:
            var_dict = copy.deepcopy(bottom_variables)
            measure = "bottom"
        if jj == 2:
            var_dict = copy.deepcopy(vertmean_variables)
            measure = "vertical_mean"
        if jj == 3:
            var_dict = copy.deepcopy(integrated_variables)
            measure = "vertical_integration"
        if jj == 4:
            var_dict = copy.deepcopy(phenology)
            measure = "phenology"
        if jj == 5:
            var_dict = copy.deepcopy(depth_profile)
            measure = "depth_profile"

        pattern_list = []
        all_df_list = []

        for kk in range(len(var_dict)):
            var_choice = list(var_dict.keys())[kk]
            var_choice = var_choice.replace(" ", "_")

            measure_name = measure
            if measure_name == "top":
                measure_name = "surface"
            if measure_name == "bottom":
                measure_name = "bottom"
            if measure_name == "vertical_mean":
                measure_name = "vertically averaged"
            if measure_name == "vertical_integration":
                measure_name = "vertically integrated"
            if measure_name == "phenology":
                measure_name = "phenology of"
            if measure_name == "depth_profile":
                measure_name = "depth_profile of"
            print(f"********* Extracting data for {measure} {var_choice} ************")
            print("Identifying time info and finding relevant files")

            pattern_files_list = []
            for ii in range(2):
                while True:
                    if ii == 0:
                        folder = sim_dir_1
                    else:
                        folder = sim_dir_2
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
                        #options = [x for x in options if "part" not in os.path.basename(x)]
                        options = [x for x in options if "restart" not in os.path.basename(x)]

                    if len([x for x in options if ".nc" in x]) > 0:
                        break

                # remove any files from options if parts of exclude are in them
                for exc in exclude:
                    options = [x for x in options if f"{exc}" not in os.path.basename(x)]

                random.shuffle(options)

                var_choice = list(var_dict.keys())[kk]

                var_choice = var_choice.replace(" ", "_")
                vv = var_choice
                # replace " " with "_"
                model_vars = list(var_dict.values())[kk]

                # if isinstance(model_vars, list):
                    # model_vars = model_vars[ii]

                t_res = None 
                pattern = None

                for ff in options:
                    ds_ff = nc.open_data(ff, checks = False)
                    the_vars = ds_ff.variables
                    if isinstance(model_vars, list):
                        var1 = model_vars[ii].split("+")[0]
                    else:
                        var1 = model_vars.split("+")[0]
                    if var1 in the_vars:
                        ff_pattern = os.path.basename(ff)
                        # replace integers with 4 or more digits with **
                        ff_pattern = re.sub(r"\d{4,}", "**", ff_pattern)
                        # replace strings of the form _12. with _**.
                        ff_pattern = re.sub(r"\d{2,}", "**", ff_pattern)
                        ff_res = get_time_res(ff, folder)
                        if t_res is None:
                            pattern = ff_pattern
                        else:
                            # if ff_res in ["d", "m"]:
                            pattern = ff_pattern
                    i+= 1
                    if n_check is not None:
                        if i > n_check:
                            break 

                if pattern is None:
                    raise ValueError(f"Unable to identify any files with {var1} in them. Double check the variable names supplied!")
                pattern_list.append(pattern)
                search_pattern = folder + "/"
                for i in range(n_dirs_down):
                    search_pattern = search_pattern + "/**"
                search_pattern = search_pattern + "/" + pattern
                ensemble = glob.glob(search_pattern, recursive = True)
                #  remove anything with restart
                ensemble = [x for x in ensemble if "restart" not in x]
                # nothing with spinup
                ensemble = [x for x in ensemble if "spinup" not in x]
                ensemble = list(set(ensemble))
                pattern_files_list.append(ensemble)
                ds = xr.open_dataset(ensemble[0])
                time_name = [x for x in list(ds.dims) if "time" in x][0]
                ensemble = list(set(ensemble))

                # randomize ensemble
                random.shuffle(ensemble)

                for ff in ensemble:
                    if ff in all_times.keys():
                        continue

                    if "restart" in ff:
                        continue
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
                    all_times[ff] = df_ff

            if True:
                # now create a climatology for each variable
                min_year = start
                max_year = end
                for i in range(2):
                    pattern = pattern_list[i] 
                    pattern_files = pattern_files_list[i]
                    i_years = []
                    for ff in pattern_files:
                        ff_times = all_times[ff]
                        i_years += list(set(ff_times.year))
                    min_year = max(min(i_years), min_year)
                    max_year = min(max(i_years), max_year)

                if min_year > max_year:
                    raise ValueError("No data found for the specified years")

                # now do the climatology
                for i in range(2):
                    if i == 0:
                        print(f"Extracting data for {measure} {var_choice} for the first simulation")
                    else:
                        print(f"Extracting data for {measure} {var_choice} for the second simulation")
                    sim_name = "sim_" + str(i)
                    out_file = f"data/climatologies/{vv}/{measure}/{measure}-climatology-{vv}-{sim_name}.nc"
                    if session_info["overwrite"] is False:
                        if os.path.exists(out_file):
                            continue
                    var_files = []
                    # get the files first
                    pattern = pattern_list[i] 
                    pattern_files = pattern_files_list[i]
                    for ff in pattern_files:
                        ff_times = all_times[ff]
                        ff_times = ff_times.query("year >= @min_year and year <= @max_year")
                        if len(ff_times) > 0:
                            var_files.append(ff)

                    ds = nc.open_data(var_files, checks = False)
                    ds_contents = ds.contents
                    if isinstance(model_vars, list):
                        vv_model = model_vars[i].split("+")
                    else:
                        vv_model = model_vars.split("+")

                    try:
                        unit = ds_contents.query("variable in @vv_model").reset_index(drop = True).unit[0]
                    except:
                        unit = ""
                    try:
                        long_name = ds_contents.query("variable in @vv_model").reset_index(drop = True).long_name[0]
                    except:
                        long_name = ""
                    if measure == "vertical_integration":
                        if "mmol" in unit:
                            if "mmol" in unit:
                                mass = None
                                element = unit.split("/")[0].split(" ")[-1]
                                if element == "O_2":
                                    element = "O"
                                try:
                                    mass = get_molar_mass(element)
                                except:
                                    if "alkalinity" in long_name:
                                        pass
                                    else:
                                        if "oxygen" in long_name:
                                            element = "O"
                                        if "phosphate" in long_name:
                                            element = "P"
                                        if "nitrate" in long_name:
                                            element = "N"
                                        if "nitrogen" in long_name:
                                            element = "N"
                                        if "silicate" in long_name:
                                            element = "Si"
                                        if "carbon carbon" in long_name:
                                            element = "C"
                                        mass = get_molar_mass(element)
                                        # figure out if you can modify the variable
                                if mass is not None:
                                    unit = unit.replace("mmol", "mg")
                                    # ds.set_units({vv_model[0]:unit})
                    
                    # capture warning
                    with warnings.catch_warnings(record=True) as w:
                        ds.subset(variables = vv_model)
                        ds.as_missing(0)
                        ds.sum_all()
                        if measure == "top":
                            ds.top()
                        if measure == "phenology":
                            ds.top()
                        if measure == "bottom":
                            if thickness_files[i] is not None:
                                thickness = thickness_files[i]
                            else:
                                thickness = session_info["obs_dir"] +  "/amm7_e3t.nc"
                            n_levels = len(ds.levels) - 1
                            ds.cdo_command(f"sellevidx,1/{n_levels}")
                            ds.cdo_command("bottomvalue")
                            # ds.cdo_command("bottomvalue")
                        if measure == "vertical_mean":
                            if thickness_files[i] is not None:
                                thickness = thickness_files[i]
                            else:
                                thickness = session_info["obs_dir"] +  "/amm7_e3t.nc"
                            n_levels = len(ds.levels) - 1
                            ds_thickness = nc.open_data(thickness, checks = False)
                            ds_thickness.cdo_command(f"sellevidx,1/{n_levels}")
                            ds.cdo_command(f"sellevidx,1/{n_levels}")
                            ds.vertical_mean(thickness = ds_thickness)

                        if measure == "depth_profile":
                            if thickness_files[i] is not None:
                                thickness = thickness_files[i]
                            else:
                                thickness = session_info["obs_dir"] +  "/amm7_e3t.nc"
                            ds_depth = nc.open_data(thickness)
                            ds_depth.subset(variable = "e3t", time = 0)
                            ds_e3t = ds_depth.copy() 
                            ds_e3t / 2
                            ds_depth - ds_e3t
                            ds_depth.vertical_cumsum()
                            ds_depth.run()
                            ds_depth.rename({"e3t": "depth"})
                            ds_e3t = nc.open_data(ff)
                            ds_e3t.subset(variable = "e3t", time = 0)
                            ds_e3t.run()
                            ds_depth.append(ds_e3t)
                            ds_depth.merge("variables")
                            ds.tmean()
                            ds.ensemble_mean()
                            ds.append(ds_depth)
                            ds.merge("variable")
                            ds.assign(total_e3t = lambda x: x.total * x.e3t, total_e3t_depth = lambda x: x.total * x.e3t * x.depth)
                            ds.run()
                            ds.vertical_sum()
                            ds.assign(total_depth = lambda x: x.total_e3t_depth/x.total_e3t, drop = True)

                        if measure == "vertical_integration":
                            if thickness_files[i] is not None:
                                thickness = thickness_files[i]
                            else:
                                thickness = session_info["obs_dir"] +  "/amm7_e3t.nc"
                            ds.vertical_integration(thickness = thickness)
                        ds.merge("time")
                        if measure == "phenology":
                            ds.tmean(["year", "month", "day"])
                        else:
                            ds.tmean(["year", "month"])
                            ds.tmean(["month"])

                        ds.run()
                        try:
                            ds.fix_amm7_grid()
                        except:
                            pass
                        if measure == "phenology":
                            ds.phenology(var = ds.variables[0], metric = "peak")
                            ds.ensemble_mean()
                        else:
                            ds.rename({ds.variables[0]: vv})
                            ds.set_units({vv:unit})
                        if mass is not None:
                            ds * mass
                        sim_name = "sim_" + str(i)
                        out_file = f"data/climatologies/{vv}/{measure}/{measure}-climatology-{vv}-{sim_name}.nc"
                        # if this exists, delte it
                        if os.path.exists(out_file):
                            os.remove(out_file)
                        # ensure directory exists
                        if not os.path.exists(os.path.dirname(out_file)):
                            os.makedirs(os.path.dirname(out_file))
                        ds.set_fill(-9999)
                    ds.to_nc(out_file, zip = True)

    if build:
        simulation_differences_comparison()


        
