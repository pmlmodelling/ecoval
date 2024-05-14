# import what is needed in nsbc_matchup
import os
import glob
import warnings
import time
import numpy as np
import nctoolkit as nc

from ecoval.fixers import tidy_warnings
from ecoval.utils import extension_of_directory, get_extent
from ecoval.session import session_info



def write_report(x):
    # append x to report
    with open("matchup_report.md", "a") as f:
        f.write(x + "\n")
        # add blank line
        f.write("\n")


def gridded_matchup(
    df_mapping=None,
    folder=None,
    var_choice=None,
    exclude=None,
    surface=None,
    start=None,
    spinup=None,
    sim_start=None,
    sim_end=None,
    e3t=None,
    domain="nws",
    strict = True,
    lon_lim = None,
    lat_lim = None
):
    """
    Function to create NSBC matchups for a given set of variables

    Parameters
    ----------
    df_mapping : pandas.DataFrame
        DataFrame containing the mapping between model variables and NSBC variables
    folder : str
        Path to folder containing model data
    var_choice : list
        List of variables to create matchups for
    exclude : list
        List of strings to exclude from the file search
    surface : str
        Surface to use for matchups. Either "top" or "bottom"
    start : int
        Start year for matchups
    spinup : int
        Spinup time for matchups
    sim_start : int
        Start year for model simulations
    sim_end : int
        End year for model simulations
    e3t : xarray.DataArray
        Vertical grid spacing for model data
    domain : str
        Domain to use for matchups. Either "NWS" or "global"
        This indicates whether the matchups use northwest European shelf data or global data

    """
    data_dir = session_info["data_dir"] 

    all_df = df_mapping
    # if model_variable is None remove from all_df
    good_model_vars = [x for x in all_df.model_variable if x is not None]

    all_df = all_df.query("model_variable in @good_model_vars").reset_index(drop=True)

    all_df = all_df.dropna()

    vars = [
        "ammonium",
        "chlorophyll",
        "nitrate",
        "phosphate",
        "oxygen",
        "silicate",
        "salinity",
        "poc",
        "temperature",
        "co2flux",
        "pco2",
        "ph",
        "alkalinity",
        "doc",
    ]
    vars = [x for x in vars if x in var_choice]

    if len(vars) > 0:
        # first up, do the top

        mapping = dict()
        ds_all = nc.open_data()

        for vv in vars:
            # figure out the data source
            dir_var = f"{data_dir}/gridded/{domain}/{vv}"
            # check if it exists
            if not os.path.exists(dir_var):
                dir_var = f"{data_dir}/gridded/global/{vv}"

            vv_source = [
                os.path.basename(x).replace(".txt", "")
                for x in glob.glob(dir_var + "/*")
                if ".txt" in x
            ][0]
            # out_dir = f"matched/gridded/{vv_source}"

            # if not os.path.exists(out_dir):
            #     os.makedirs(out_dir)

            print("**********************")
            #
            vv_name = vv
            if vv == "poc":
                vv_name = "particulate organic carbon"
            if vv == "doc":
                vv_name = "dissolved organic carbon"
            if vv == "co2flux":
                vv_name = "sea-air CO2 flux"
            if vv == "ph":
                vv_name = "pH"

            print(
                f"Matching up surface {vv_name} with {vv_source.upper()} gridded data"
            )
            print("**********************")
            df = df_mapping.query("variable == @vv").reset_index(drop=True)
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
                    if surface == "top":
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
                all_years = ds.years

                if spinup is not None:
                    years = [x for x in all_years if x >= min(all_years) + spinup]

                sim_years = range(sim_start, sim_end + 1)
                if start is not None:
                    years = [x for x in all_years if x in sim_years]

                for ff in paths:
                    try:
                        ds = nc.open_data(ff, checks=False)
                        ds_years = ds.years
                        if len([x for x in ds_years if x in years]) > 0:
                            new_paths.append(ff)
                    except:
                        print(f"Unable to find relevant years in  {ff}")

                paths = list(set(new_paths))
                paths.sort()

                # get the number of paths

                n_paths = len(paths)

                # Write to report
                write_report(f"### Matchups for {vv}")
                write_report(f"Number of paths: {n_paths}")
                # minimum year
                write_report(f"Minimum year: {min(years)}")
                # maximum year
                write_report(f"Maximum year: {max(years)}")
                # write the list of files
                write_report(f"Files used for {vv}:")
                write_report("```")
                paths.sort()
                for ff in paths:
                    write_report(ff)
                write_report("```")

                # add a pagebreak
                write_report("\\newpage")

                # figure out if cdo or nco is faster....

                use_nco = False
                if surface == "top":
                    try:
                        with warnings.catch_warnings(record=True) as w:
                            use_nco = False
                            cdo_time = time.time()

                            ds = nc.open_data(paths[0], checks=False)
                            ds.subset(variables=selection)
                            if surface == "top":
                                ds.top()
                            else:
                                ds.bottom()
                            ds.run()
                            cdo_time = time.time() - cdo_time

                            nco_time = time.time()
                            ds = nc.open_data(paths[0], checks=False)
                            nco_selection = ",".join(selection)
                            ds.nco_command(f"ncks -F -d deptht,1 -v {nco_selection}")
                            ds.run()
                            nco_time = time.time() - nco_time

                            if nco_time < cdo_time:
                                use_nco = True
                        tidy_warnings(w)

                    except:
                        use_nco = False

                with warnings.catch_warnings(record=True) as w:
                    ds = nc.open_data(paths, checks=False)

                    ds.subset(years=years)
                    if vv_source == "glodap":
                        ds_years = ds.years
                        if len([x for x in ds_years if x in range(1971, 2015)]) > 0:
                            ds.subset(years=range(1971, 2015))
                            ds.run()

                    if use_nco:
                        ds.nco_command(f"ncks -F -d deptht,1 -v {nco_selection}")
                    else:
                        ds.subset(variables=selection)
                        if surface == "top":
                            ds.top()
                        else:
                            ds.bottom()

                    ds.as_missing(0)
                    if vv_source == "glodap":
                        ds.merge("time")
                        ds.tmean()

                    # the code below needs to be simplifed
                    # essentially anything with a + in the mapping should be split out
                    # and then the command should be run for each variable

                    ignore_later = []
                    for vv in list(df.variable):
                        if "+" in mapping[vv]:
                            command = f"-aexpr,{vv}=" + mapping[vv]
                            ds.cdo_command(command)
                            drop_these = mapping[vv].split("+")
                            ds_contents = ds.contents
                            ds_contents = ds_contents.query("variable in @drop_these")
                            var_unit = ds_contents.unit[0]
                            ds.drop(variables=drop_these)
                            ignore_later.append(vv)

                    ds.run()
                    for key in mapping:
                        if key not in ignore_later:
                            if mapping[key] in ds.variables:
                                ds.rename({mapping[key]: key})
                    if "chlorophyll" in list(df.variable):
                        ds.set_units({"chlorophyll": var_unit})
                        ds.set_longnames(
                            {"chlorophyll": "Total chlorophyll concentration"}
                        )
                    if "poc" in list(df.variable):
                        ds.set_units({"poc": var_unit})
                        ds.set_longnames(
                            {"poc": "Particulate organic carbon concentration"}
                        )
                    if "doc" in list(df.variable):
                        ds.set_units({"doc": var_unit})
                        ds.set_longnames(
                            {"doc": "Dissolved organic carbon concentration"}
                        )

                    ds.run()
                    ds.merge("time")
                    ds.subset(years = years)
                    ds.run()

                tidy_warnings(w)

                # figure out the start and end year
                start_year = min(ds.years)
                end_year = max(ds.years)

                vv_file = nc.create_ensemble(dir_var)
                # except:
                # vv_file = nc.create_ensemble(dir_var)
                ds_obs = nc.open_data(
                    vv_file,
                    checks=False,
                )
                obs_years = ds_obs.years

                if len(obs_years) == 1:
                    ds.tmean("month")
                else:
                    ds.tmean(["year", "month"])

                amm7 = False
                if max(ds.contents.npoints) == 111375:
                    ds.fix_amm7_grid()
                    amm7 = True
                    ds.subset(lon=[-19, 9], lat=[41, 64.3])
                if lon_lim is not None and lat_lim is not None:
                    ds.subset(lon=lon_lim, lat = lat_lim)

                # dir_var = f"{data_dir}/gridded/{domain}/{vv}"
                # if vv not in ["poc", "temperature"]:
                if vv in ["poc", "doc"]:
                    if strict:
                        ds_obs.subset(years=years)
                    ds_obs.merge("time")
                    ds_obs.tmean("month")
                    ds.tmean("month")

                if vv in ["temperature"]:
                    if strict:
                        ds_obs.subset(years=years)
                    ds_obs.subset(years=years)
                    ds_obs.tmean(["year", "month"])
                    ds_obs.merge("time")
                    ds_obs.tmean(["year", "month"])
                if vv in ["salinity"] and domain != "nws":
                    ds_obs.top()
                    ds_obs.subset(years=years)

                if vv not in ["poc", "temperature"]:
                    if len(ds_obs.times) > 12:
                        ds_obs.subset(years=years)

                if vv_source == "occci" and vv == "chlorophyll":
                    ds_obs.subset(variable="chlor_a")

                ds_xr = ds.to_xarray()
                lon_name = [x for x in ds_xr.coords if "lon" in x]
                lat_name = [x for x in ds_xr.coords if "lat" in x]
                lon = ds_xr[lon_name[0]].values
                lat = ds_xr[lat_name[0]].values
                lon_max = lon.max()
                lon_min = lon.min()
                lat_max = lat.max()
                lat_min = lat.min()

                # figure out the lon/lat extent in the model
                lons = [lon_min, lon_max]
                lats = [lat_min, lat_max]
                # start of with the raw coords
                # This will not work with nemo, which outputs the grid incorrectly
                # so we will check if the step between the first lon/lat and the second lon/lat is
                # far bigger than the rest. If this is the case, the first should be ignored
                # get the lon/lat values
                lon_vals = ds_xr[lon_name[0]].values
                lat_vals = ds_xr[lat_name[0]].values
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

                lon_min_model = lons[0]
                lon_max_model = lons[1]
                lat_min_model = lats[0]
                lat_max_model = lats[1]

                # now do the same for ds_nsbc
                ds_xr = ds_obs.to_xarray()
                lon_name = [x for x in ds_xr.coords if "lon" in x]
                lat_name = [x for x in ds_xr.coords if "lat" in x]
                lon = ds_xr[lon_name[0]].values
                lat = ds_xr[lat_name[0]].values
                lon_max = lon.max()
                lon_min = lon.min()
                lat_max = lat.max()
                lat_min = lat.min()

                lon_min = max(lon_min, lon_min_model)
                lon_max = min(lon_max, lon_max_model)
                lat_min = max(lat_min, lat_min_model)
                lat_max = min(lat_max, lat_max_model)

                lons = [lon_min, lon_max]
                lats = [lat_min, lat_max]


                if domain != "global":
                    ds.subset(lon=lons, lat=lats)
                    ds_obs.subset(lon=lons, lat=lats)

                if domain == "global":
                    model_extent = get_extent(ds[0])
                    obs_extent = get_extent(ds_obs[0])
                    lon_min = max(model_extent[0], obs_extent[0])
                    lon_max = min(model_extent[1], obs_extent[1])
                    lat_min = max(model_extent[2], obs_extent[2])
                    lat_max = min(model_extent[3], obs_extent[3])
                   # make sure lon_min is greater than -180
                    if lon_min < -180:
                        lon_min = -180
                    if lon_max > 180:
                        lon_max = 180
                    if lat_min < -90:
                        lat_min = -90
                    if lat_max > 90:
                        lat_max = 90

                    lons = [lon_min, lon_max]
                    lats = [lat_min, lat_max]
                    print(lons)
                    print(lats)
                    ds.subset(lon=lons, lat=lats)
                    ds_obs.subset(lon=lons, lat=lats)
                

                n1 = ds_obs.contents.npoints[0]
                n2 = ds.contents.npoints[0]

                # if not amm7:
                if n1 >= n2:
                    ds_obs.regrid(ds, method="nn")
                else:
                    ds.regrid(ds_obs, method="nn")

                ds_obs.rename({ds_obs.variables[0]: "observation"})
                ds.rename({ds.variables[0]: "model"})

                # it is possible the years do not overlap, e.g. with satellite Chl
                if len(ds.times) > 12:
                    years1 = ds.years
                    years2 = ds_obs.years
                    all_years = [x for x in years1 if x in years2]
                    if len(all_years) != len(years1):
                        if len(all_years) != len(years2):
                            ds_obs.subset(years=all_years)
                            ds.subset(years=all_years)
                            ds_obs.run()
                            ds.run()
                if len(ds_obs) > 1:
                    ds_obs.merge("time")

                ds_obs.run()
                ds.run()

                if vv == "doc":
                    ds_obs * 12.011
                    ds + (40 * 12.011)

                ds_obs.append(ds)

                if len(ds.times) > 12:
                    ds_obs.merge("variable", match=["year", "month"])
                else:
                    ds_obs.merge("variable", match="month")
                ds_obs.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
                ds_obs.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")
                ds_obs.set_fill(-9999)
                ds_obs.top()
                ds_mask = ds_obs.copy()
                ds_mask.assign( mask_these=lambda x: -1e30 * ((isnan(x.observation) + isnan(x.model)) > 0), drop=True)
                ds_mask.as_missing([-1e40, -1e20])
                ds_obs + ds_mask

                # fix the co2 flux units
                if vv == "co2flux":
                    ds_obs.assign(model=lambda x: x.model * -0.365)
                    ds_obs.set_units({"model": "mol/m2/yr"})
                    ds_obs.set_units({"observation": "mol/m2/yr"})

                # figure out if the temperature is in degrees C
                if vv == "temperature":
                    if ds_obs.to_xarray().model.max() > 100:
                        ds_obs.assign(model=lambda x: x.model - 273.15)
                    if ds_obs.to_xarray().observation.max() > 100:
                        ds_obs.assign(observation=lambda x: x.observation - 273.15)
                    # set the units
                    ds_obs.set_units({"model": "degrees C"})
                    ds_obs.set_units({"observation": "degrees C"})
                # # now, we need to exclude data outside the lon/lat range with data

                out_file = f"matched/gridded/{domain}/{vv}/{vv_source}_{vv}_surface.nc"
                # check directory exists for out_file
                if not os.path.exists(os.path.dirname(out_file)):
                    os.makedirs(os.path.dirname(out_file))
                # remove the file if it exists
                if os.path.exists(out_file):
                    os.remove(out_file)
                ds_obs.set_precision("F32")
                if vv == "salinity" and domain != "nws":
                    ds_obs.tmean("month")
                ds_obs.to_nc(out_file, zip=True, overwrite=True)

                # ds_all.append(ds)

        # if len(ds_all) > 1:
        #     ds_all.merge("variable", "month")
        #     ds_all.run()
        #     # add start year as a global attribute using nco
        #     ds_all.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
        #     ds_all.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")

        #     ds_all.to_nc(f"matched/gridded/{}/nsbc_model.nc", zip=True, overwrite=True)
        # else:
        #     if len(ds_all) == 1:
        #         # add start year as a global attribute using nco
        #         ds_all.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
        #         ds_all.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")
        #         ds_all.to_nc(
        #             "matched/gridded/nsbc/nsbc_model.nc", zip=True, overwrite=True
        #         )
        #     else:
        #         print(f"No NSBC matchups for {vv} for surface")

        # now, do the vertical mean

        return None

        mapping = dict()
        ds_all = nc.open_data()

        for vv in vars:
            df = df_mapping.query("variable == @vv").reset_index(drop=True)
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

                years = range(sim_start, sim_end + 1)

                new_paths = []
                for ff in paths:
                    try:
                        ds = nc.open_data(ff, checks=False)
                        ds_years = ds.years
                        if len([x for x in ds_years if x in years]) > 0:
                            new_paths.append(ff)
                    except:
                        print(f"Unable to find relevant years in  {ff}")

                paths = list(set(new_paths))
                paths.sort()

                # figure out if cdo or nco is faster....

                use_nco = False

                with warnings.catch_warnings(record=True) as w:
                    ds = nc.open_data(paths, checks=False)

                    ds_years = ds.years
                    years = [x for x in ds_years if x >= sim_start and x <= sim_end]

                    if spinup is not None:
                        start_year = min(ds.years) + spinup
                    else:
                        start_year = min(years)
                    years = [x for x in years if x >= start_year]

                    ds.subset(years=years)

                    if use_nco:
                        ds.nco_command(f"ncks -F -d deptht,1 -v {nco_selection}")
                    else:
                        ds.subset(variables=selection)
                    ds.as_missing(0)

                    if "chlorophyll" in list(df.variable):
                        command = "-aexpr,chlorophyll=" + mapping["chlorophyll"]
                        ds.cdo_command(command)
                        drop_these = mapping["chlorophyll"].split("+")
                        ds_contents = ds.contents
                        ds_contents = ds_contents.query("variable in @drop_these")
                        chl_unit = ds_contents.unit[0]
                        ds.drop(variables=drop_these)
                    ds.run()
                    for key in mapping:
                        if key != "chlorophyll":
                            if mapping[key] in ds.variables:
                                ds.rename({mapping[key]: key})
                    if "chlorophyll" in list(df.variable):
                        ds.set_units({"chlorophyll": chl_unit})
                        ds.set_longnames(
                            {"chlorophyll": "Total chlorophyll concentration"}
                        )

                    ds.run()
                tidy_warnings(w)

                ds.merge("time")
                ds.run()
                # figure out the start and end year
                start_year = min(ds.years)
                end_year = max(ds.years)

                ds.tmean("month")

                with warnings.catch_warnings(record=True) as w:
                    if e3t is None:
                        try:
                            ds_thickness = nc.open_data(paths[0], checks=False)
                            ds_thickness.subset(time=0, variables="e3t")
                            ds_thickness.run()
                            e3t = ds_thickness[0]

                        except:
                            pass
                tidy_warnings(w)

                if e3t is not None:
                    ds.vertical_mean(fixed=False, thickness=e3t)
                else:
                    ds.vertical_mean(fixed=False)

                amm7 = False
                if max(ds.contents.npoints) == 111375:
                    ds.fix_amm7_grid()
                    amm7 = True
                ds_nsbc = nc.open_data(
                    f"{data_dir}/nsbc/level_3/climatological_monthly_mean/NSBC_Level3_phosphate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN_1960_2014.nc",
                    checks=False,
                )

                if not amm7:
                    ds.regrid(ds_nsbc)

                ds_all.append(ds)

        if len(ds_all) > 1:
            ds_all.merge("variable", "month")
            ds_all.run()
            # add start year as a global attribute using nco
            ds_all.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
            ds_all.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")

            ds_all.to_nc(
                "matched/gridded/nsbc/nsbc_model_verticalmean.nc",
                zip=True,
                overwrite=True,
            )
        else:
            if len(ds_all) == 1:
                # add start year as a global attribute using nco
                ds_all.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
                ds_all.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")
                ds_all.to_nc(
                    "matched/gridded/nsbc/nsbc_model_verticalmean.nc",
                    zip=True,
                    overwrite=True,
                )
            else:
                print(f"No NSBC matchups for vertical mean for {vv}")

        del ds_all
        del ds
