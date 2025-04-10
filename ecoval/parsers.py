import nctoolkit as nc
import re
import pandas as pd
import os
import xarray as xr
import warnings
from ecoval.session import session_info

from pathlib import Path
from netCDF4 import Dataset

bad_conc_vars = ["medium", "pod", "size"]

def standard_mapping(model = None):
    if model is None:
        return None
    if model == "ersem":
        mapping = {}
        mapping["alkalinity"] = "O3_TA"
        mapping["ammonium"] = "N4_n"
        mapping["benbio"] = "Y2_c+Y3_c"
        mapping["benthic_carbon_flux"] = "net_PelBen_POC_result"
        mapping["carbon"] = "Q6_c+Q7_c"
        mapping["chlorophyll"] = "P1_Chl+P2_Chl+P3_Chl+P4_Chl"
        mapping["co2flux"] = "O3_fair"
        mapping["doc"] = "R1_c+R2_c+R3_c"
        mapping["mesozoo"] = "Z4_c"
        mapping["micro"] = "P1_Chl+P4_Chl"
        mapping["nano"] = "P2_Chl"
        mapping["nitrate"] = "N3_n"
        mapping["oxycons"] = "ben_nit_nrate+Ymacro_fYG3c_result+Y4_fYG3c+H1_fHG3c+H2_fHG3c"
        mapping["oxygen"] = "O2_o"
        mapping["pco2"] = "O3_pCO2"
        mapping["ph"] = "O3_pH"
        mapping["phosphate"] = "N1_p"
        mapping["pico"] = "P3_Chl"
        mapping["poc"] = None
        mapping["salinity"] = "vosaline"
        mapping["silicate"] = "N5_s"
        mapping["spm"] = "P4_Chl"
        mapping["temperature"] = "votemper"
        mapping["kd"] = "light_xEPS"
        return mapping
    # throw an error
    raise ValueError("Model not supported")




def fvcom_contents(ds):
    import xarray as xr
    drop_variables = ['siglay','siglev']

    ff = ds[0]

    ds_xr = xr.open_dataset(ff, drop_variables = drop_variables, decode_times=False)
    #ds_nc = nc.open_data(ff, checks = False)
    nc_ds = Dataset(ff)
    variables = ds_xr.variables
    # variables = list(ds_xr.variables)
    good_vars = []
    longs = []
    for vv in variables:
        try:
            x = nc_ds[vv].long_name
            good_vars.append(vv)
            longs.append(x)
        except:
            pass
    contents = pd.DataFrame({"variable":good_vars,"long_name":longs})
    contents = contents.assign(nlevels = 50)
    return contents

def generate_mapping(ds, fvcom = False):
    """
    Generate mapping of model and observational variables
    """

    candidate_variables = [
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
        "benthic_carbon_flux",
        "mesozoo",
        "spm",
        "oxycons",
        "kd"
    ]
    if fvcom is False:
        ds1 = nc.open_data(ds[0], checks=False)
        ds_contents = ds1.contents
    else:
        ds_contents = fvcom_contents(ds)

    ds_contents["long_name"] = [str(x) for x in ds_contents["long_name"]]

    ds_contents_top = ds_contents.query("nlevels == 1").reset_index(drop=True)
    ds_contents = ds_contents.query("nlevels > 1").reset_index(drop=True)
    ds_variables = ds_contents.variable

    model_dict = {}
    for vv in candidate_variables:
        vv_check = vv
        if vv != "ph":
            the_vars = [
                x
                for x in [str(x) for x in ds_contents.long_name]
                if vv_check.lower() in x.lower()
                and "benthic" not in x.lower()
                and "river" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None:
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)
        if vv == "spm":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x and ("micro" in x)
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
                x
                for x in ds_contents.long_name
                if "photolabile" in str(x) and "carbon" in str(x)
            ]
            if len(vars_2) > 0:
                the_vars += vars_2

            the_vars = [x for x in the_vars if " loss " not in x]
            the_vars = [x for x in the_vars if "depth" not in x]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)
        
        if vv == "benthic_carbon_flux":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "carbon" in x.lower()
                and "benth" in x.lower()
                and "penetr" not in x.lower()
                and "flux" in x.lower()
                and "organic" in x.lower()
                and "inorganic" not in x.lower()
            ]
            the_vars = list(set(the_vars))
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "carbon":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "carbon" in x
                and "benthic" in x.lower()
                and ("refractory" in x.lower() or "particul" in x.lower())
                and "penetr" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "benbio":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "carbon" in x
                and ("feeder" in x.lower() or "predator" in x.lower())
                and "penetr" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "ph":
            the_vars = [
                x
                for x in ds_contents.long_name
                if ("pH" in x) or (" ph " in x) and "benthic" not in x.lower() and "river" not in x.lower()
            ]
            if len(the_vars) == 0:
                if len(ds_contents.query("variable == 'ph'")) > 0:
                    the_vars = list(ds_contents.query("variable == 'ph'").long_name)
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "co2flux":
            the_vars = [
                x
                for x in ds_contents_top.long_name
                if "co2" in x.lower()
                and "flux" in x.lower()
                and "river" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None:
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "pco2":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "carbonate" in x.lower()
                and "partial" in x.lower()
                and "river" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None:
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "silicate":
            the_vars = [
                x
                for x in ds_contents.long_name
                if ("silicate" in x.lower() or "silicic" in x.lower())
                and "benthic" not in x.lower()
                and "river" not in x.lower()
            ]

            if len(the_vars) > 1:
                the_vars_2 = [
                    x
                    for x in ds_contents.long_name
                    if "benthic" not in x.lower()
                    and "river" not in x.lower()
                ]
                the_vars_2 = [x for x in the_vars_2 if re.match(r"silicate.silicate", x)]
                if len(the_vars_2) == 1:
                    the_vars = the_vars_2
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "oxygen":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "oxygen" in x.lower()
                and "benthic" not in x.lower()
                and "saturation" not in x.lower()
                and "util" not in x.lower()
                and "flux" not in x.lower()
                and "river" not in x.lower()
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "micro":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x and ("micro" in x or "diatom" in x)
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "mesozoo":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "mesozoo" in x and
                "ingestion" not in x.lower()  and
                "mortality" not in x.lower()  and
                "loss" not in x.lower()  

            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None:
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)
        if vv == "kd":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "atten" in x.lower() and "coeff" in x.lower()
            ]


        if vv == "nano":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x.lower() and "nano" in x
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv == "pico":
            the_vars = [
                x
                for x in ds_contents.long_name
                if "chloroph" in x.lower() and "pico" in x
            ]
            if len(the_vars) > 1:
                if session_info["model"] is not None:
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)
        
        if vv == "oxycons":
            oxy_con_vars = list(set(["Ymacro_fYG3c_result", "Y4_fYG3c", "H1_fHG3c", "H2_fHG3c", "ben_nit_nrate"]))
            the_vars = [x for x in ds_contents_top.variable if x in oxy_con_vars]
            if len(the_vars) != len(oxy_con_vars):
                the_vars = []
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)


        if vv == "chlorophyll":
            if len(the_vars) > 1:
                the_vars = [x for x in the_vars if "total" not in x.lower()]
            if len(the_vars) > 1:
                if session_info["model"] is not None: 
                    standard_var = standard_mapping(session_info["model"])[vv]
                    if len(ds_contents.query("variable == @standard_var")) > 0:
                        the_vars = list(ds_contents.query("variable == @standard_var").long_name)

        if vv in ["carbon", "benbio", "benthic_carbon_flux"]:
            model_vars = ds_contents_top.query("long_name in @the_vars").variable
            if len(the_vars) > 0:
                if vv == "benthic_carbon_flux":
                    model_vars = [x for x in model_vars if "net" in x]
        else:
            if vv != "co2flux" and vv != "pco2":
                model_vars = ds_contents.query("long_name in @the_vars").variable
            else:
                if vv == "co2flux":
                    model_vars = ds_contents_top.query(
                        "long_name in @the_vars"
                    ).variable
                else:
                    model_vars = ds_contents.query(
                        "long_name in @the_vars"
                    ).variable

        add = True
        if vv == "oxycons":
            model_vars = the_vars

        if len(model_vars) > 1 and vv not in [
            "doc",
            "chlorophyll",
            "carbon",
            "benbio",
            "micro",
            "oxycons"
        ]:
            add = False

        if add:
            if len(model_vars) > 0:
                model_dict[vv] = "+".join(model_vars)
            else:

                if "nitrate" not in model_dict.keys() and vv == "nitrate":
                    if "ammonium" not in model_dict.keys():
                        df_nitrogen = (
                            ds_contents.query("long_name.str.contains('nitrogen')")
                            .query("long_name.str.contains('nutrient')")
                            .reset_index(drop=True)
                        )
                        if len(df_nitrogen) == 1:
                            model_dict["nitrate"] = df_nitrogen.variable[0]
                            warnings.warn(
                                "No nitrate variable found, using nitrogen nutrient variable"
                            )

                if vv not in model_dict.keys():
                    model_dict[vv] = None
    return model_dict


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
