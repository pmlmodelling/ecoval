import nctoolkit as nc
import pandas as pd
import os
import xarray as xr
import warnings
from ecoval.utils import get_obsdir

from pathlib import Path
from netCDF4 import Dataset

obs_dir = get_obsdir()



def generate_mapping(ds):
    """
    Generate mapping of model and ICES data
    """

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
        "doc",
        "carbon",
        "benbio",
        "alkalinity",
        "micro",
        "nano",
        "pico",
    ]
    ds1 = nc.open_data(ds[0], checks=False)
    ds_contents = ds1.contents

    ds_contents["long_name"] = [str(x) for x in ds_contents["long_name"]]

    ds_contents_top = ds_contents.query("nlevels == 1").reset_index(drop=True)
    ds_contents = ds_contents.query("nlevels > 1").reset_index(drop=True)

    ersem_dict = {}
    for vv in ices_variables:
        vv_check = vv
        if vv != "ph":
            the_vars = [
                x
                for x in [str(x) for x in ds_contents.long_name]
                if vv_check.lower() in x
                and "benthic" not in x.lower()
                and "river" not in x.lower()
            ]
        if vv == "doc":
            # doc = [x for x in ds.contents.long_name if "arbon" in x and "iss" in x and " organic" in x and "benthic" not in x]
            the_vars = [
                x
                for x in ds_contents.long_name
                if "arbon" in x
                and "iss" in x
                and " organic" in x
                and "benthic" not in x
            ]
            vars_2 = [
                x for x in ds.contents.long_name if "photolabile" in str(x) and "carbon" in str(x)
            ]
            if len(vars_2) > 0:
                the_vars += vars_2

            the_vars = [x for x in the_vars if " loss " not in x]
            the_vars = [x for x in the_vars if "depth" not in x]

        if vv == "carbon":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "carbon" in x
                and "benthic" in x.lower()
                and ("refractory" in x.lower() or "particul" in x.lower())
                and "penetr" not in x.lower()
            ]

        if vv == "benbio":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "carbon" in x
                and ("feeder" in x.lower() or "predator" in x.lower())
                and "penetr" not in x.lower()
            ]

        if vv == "ph":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "pH" in x and "benthic" not in x.lower() and "river" not in x.lower()
            ]
        if vv == "co2flux":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "co2" in x.lower()
                and "flux" in x.lower()
                and "river" not in x.lower()
            ]

        if vv == "pco2":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "carbonate" in x.lower()
                and "partial" in x.lower()
                and "river" not in x.lower()
            ]

        if vv == "silicate":
            the_vars = [
                x
                for x in ds_contents.long_name
                if ("silicate" in x or "silicic" in x)
                and "benthic" not in x.lower()
                and "river" not in x.lower()
            ]
        if vv == "oxygen":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "oxygen" in x.lower()
                and "benthic" not in x.lower()
                and "river" not in x.lower()
            ]

        if vv == "micro":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x and ("micro" in x or "diatom" in x)
            ]

        if vv == "nano":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x.lower() and "nano" in x
            ]

        if vv == "pico":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x.lower() and "pico" in x
            ]

        if vv in ["carbon", "benbio"]:
            ersem_vars = ds_contents_top.query("long_name in @the_vars").variable
        else:
            if vv != "co2flux" and vv != "pco2":
                ersem_vars = ds_contents.query("long_name in @the_vars").variable
            else:
                if vv == "co2flux":
                    ersem_vars = ds_contents_top.query(
                        "long_name in @the_vars"
                    ).variable
                else:
                    ersem_vars = ds_contents_top.query(
                        "long_name in @the_vars"
                    ).variable

        add = True

        if len(ersem_vars) > 1 and vv not in ["doc", "chlorophyll", "carbon", "benbio", "micro"]:
            add = False

        if add:
            if len(ersem_vars) > 0:
                ersem_dict[vv] = "+".join(ersem_vars)
            else:

                if "nitrate" not in ersem_dict.keys() and vv == "nitrate":
                    if "ammonium" not in ersem_dict.keys():
                        df_nitrogen = (
                            ds_contents.query("long_name.str.contains('nitrogen')")
                            .query("long_name.str.contains('nutrient')")
                            .reset_index(drop=True)
                        )
                        if len(df_nitrogen) == 1:
                            ersem_dict["nitrate"] = df_nitrogen.variable[0]
                            warnings.warn(
                                "No nitrate variable found, using nitrogen nutrient variable"
                            )

                if vv not in ersem_dict.keys():
                    ersem_dict[vv] = None
    return ersem_dict


# redundant fvcom contents function. User later
# def fvcom_contents(ds):
#     import xarray as xr

#     drop_variables = ["siglay", "siglev"]

#     ff = ds[0]

#     ds_xr = xr.open_dataset(ff, drop_variables=drop_variables, decode_times=False)
#     # ds_nc = nc.open_data(ff, checks = False)
#     nc_ds = Dataset(ff)
#     variables = ds_xr.variables
#     # variables = list(ds_xr.variables)
#     good_vars = []
#     longs = []
#     for vv in variables:
#         try:
#             x = nc_ds[vv].long_name
#             good_vars.append(vv)
#             longs.append(x)
#         except:
#             print(vv)
#     contents = pd.DataFrame({"variable": good_vars, "long_name": longs})
#     return contents
