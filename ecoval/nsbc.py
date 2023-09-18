# import what is needed in nsbc_matchup
import os
import glob
import warnings
import time
import numpy as np
import pandas as pd
import xarray as xr
import nctoolkit as nc

from ecoval.fixers import tidy_warnings

def write_report(x):
    # append x to report
    with open("matchup_report.md", "a") as f:
        f.write(x + "\n")
        # add blank line
        f.write("\n")


def nsbc_matchup(
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

    """

    all_df = df_mapping
    # if model_variable is None remove from all_df
    good_model_vars = [x for x in all_df.model_variable if x is not None]

    all_df = all_df.query("model_variable in @good_model_vars").reset_index(drop=True)

    all_df = all_df.dropna()

    print("Creating gridded data for NSBC matchups")

    vars = [
        "ammonium",
        "chlorophyll",
        "nitrate",
        "phosphate",
        "oxygen",
        "silicate",
        "salinity",
    ]
    vars = [x for x in vars if x in var_choice]

    out_dir = "matched/gridded/nsbc"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if len(vars) > 0:
        # first up, do the top

        mapping = dict()
        ds_all = nc.open_data()

        for vv in vars:
            print("**********************")
            print(f"Matching up {vv} with NSBC gridded data")
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

                paths = glob.glob(folder + "/**/**/" + pattern)

                for exc in exclude:
                    paths = [x for x in paths if f"{exc}" not in os.path.basename(x)]

                new_paths = []
                ds = nc.open_data(paths, checks=False)
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

                    if use_nco:
                        ds.nco_command(f"ncks -F -d deptht,1 -v {nco_selection}")
                    else:
                        ds.subset(variables=selection)
                        if surface == "top":
                            ds.top()
                        else:
                            ds.bottom()

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
                    ds.merge("time")
                    ds.run()

                tidy_warnings(w)

                # figure out the start and end year
                start_year = min(ds.years)
                end_year = max(ds.years)

                ds.tmean("month")

                amm7 = False
                if max(ds.contents.npoints) == 111375:
                    ds.fix_amm7_grid()
                    amm7 = True
                ds_nsbc = nc.open_data(
                    "/data/proteus1/scratch/rwi/evaldata/data/nsbc/level_3/climatological_monthly_mean/NSBC_Level3_phosphate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN_1960_2014.nc",
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

            ds_all.to_nc("matched/gridded/nsbc/nsbc_model.nc", zip=True, overwrite=True)
        else:
            if len(ds_all) == 1:
                # add start year as a global attribute using nco
                ds_all.nco_command(f"ncatted -O -a start_year,global,o,c,{start_year}")
                ds_all.nco_command(f"ncatted -O -a end_year,global,o,c,{end_year}")
                ds_all.to_nc(
                    "matched/gridded/nsbc/nsbc_model.nc", zip=True, overwrite=True
                )
            else:
                print(f"No NSBC matchups for {vv} for surface")

        # now, do the vertical mean

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

                paths = glob.glob(folder + "/**/**/" + pattern)

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
                    "/data/proteus1/scratch/rwi/evaldata/data/nsbc/level_3/climatological_monthly_mean/NSBC_Level3_phosphate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN_1960_2014.nc",
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
