
import nctoolkit as nc
import pandas as pd
import os
import xarray as xr
from ecoval.utils import get_datadir
data_dir = get_datadir()
import warnings
from pathlib import Path
from netCDF4 import Dataset

def split_path(path):
    x = path.split("/")
    for i in range(len(x)):
        if "*" in x[i]:
            break
    part1 = "/".join(x[:i])
    part2 = "/".join(x[i:])
    part2 = part2.replace("**", "*")
    return [part1, part2]

def get_units(ds, mapping):
    """
    """
    base_units = pd.read_csv(f"{data_dir}/ices/ices_units.csv")
    ds1 = nc.open_data(ds[0], checks = False)
    the_contents = ds1.contents
    keys = []
    m_units = []
    o_units = []
    for key in mapping.keys():
        if mapping[key] is not None:
            vv = mapping[key].split("+")[0]
            m_units.append(str(the_contents.query("variable == @vv").unit.values[0]))
            o_units.append(str(base_units.query("variable == @key").unit.values[0]))
            keys.append(key)


    df_out = pd.DataFrame({"key":keys, "model_unit":m_units, "observation_unit":o_units})

    return df_out

def get_dates(ds):
    print("Getting dates in dataset!")
    years = []
    months = []
    time_name = None
    for ff in ds:
        ds_xr = xr.open_dataset(ff)
        if time_name is None:
            time_name = [x for x in list(ds_xr.dims) if "time" in x][0]
        years += list((ds_xr[time_name].dt.year.values))
        months += list((ds_xr[time_name].dt.month.values))
    return pd.DataFrame({"year":years, "month":months}).drop_duplicates()

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
            print(vv)
    contents = pd.DataFrame({"variable":good_vars,"long_name":longs})
    return contents

def generate_mapping(ds, fvcom = False):
    """
    Generate mapping of model and ICES data
    """

    ices_variables = ['temperature',
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
     "doc",
     "carbon",
     'alkalinity']
    ds1 = nc.open_data(ds[0], checks = False)
    try:
        ds_contents = ds1.contents
    except:
        try:
            ds_contents = fvcom_contents(ds1)
        except:
            raise ValueError("Could not read contents of dataset!")

    ds_contents["long_name"] = [str(x) for x in ds_contents["long_name"]]
    if fvcom == False:
        ds_contents_top = ds_contents.query("nlevels == 1").reset_index(drop = True)
        ds_contents = ds_contents.query("nlevels > 1").reset_index(drop = True)

    ersem_dict = {}
    for vv in ices_variables:
        vv_check = vv
        if vv != "ph":
            the_vars = [x for x in [str(x) for x in ds_contents.long_name] if vv_check.lower() in x 
                        and "benthic" not in x.lower()
                        and "river" not in x.lower()
                       ]
        if vv == "doc":
        # doc = [x for x in ds.contents.long_name if "arbon" in x and "iss" in x and " organic" in x and "benthic" not in x]
            the_vars = [x for x in ds_contents.long_name if "arbon" in x and "iss" in x and " organic" in x and "benthic" not in x]
            vars_2 = [x for x in ds.contents.long_name if "photolabile" in x and "carbon" in x]
            if len(vars_2) > 0:
                the_vars += vars_2

            the_vars = [x for x in the_vars if " loss " not in x] 
            the_vars = [x for x in the_vars if "depth" not in x] 


        if vv == "carbon":
            the_vars = [x for x in ds_contents_top.long_name if "carbon" in x 
                        and "benthic" in x.lower()
                        and ("refractory" in x.lower() or "particul" in x.lower())
                        and "penetr" not in x.lower()
                       ]

        if vv == "ph":
            the_vars = [x for x in ds_contents.long_name if "pH" in x 
                        and "benthic" not in x.lower()
                        and "river" not in x.lower()
                       ]
        if vv == "co2flux":
            if fvcom == False:
                the_vars = [x for x in ds_contents_top.long_name if "co2" in x.lower() 
                            and "flux"  in x.lower()
                            and "river" not in x.lower()
                           ]
            else:
                the_vars = [x for x in ds_contents_top.long_name if "co2" in x.lower() 
                            and "flux"  in x.lower()
                            and "river" not in x.lower()
                           ]
                

        if vv == "pco2":
            if fvcom == False:
                the_vars = [x for x in ds_contents.long_name if "pco2" in x.lower() 
                            and "carbonate"  in x.lower()
                            and "river" not in x.lower()
                           ]
            else:
                the_vars = [x for x in ds_contents_top.long_name if "pco2" in x.lower() 
                            and "carbonate"  in x.lower()
                            and "river" not in x.lower()
                           ]


        if vv == "silicate":
            if fvcom == False:
                the_vars = [x for x in ds_contents.long_name if ("silicate" in x  or "silicic" in x)
                            and "benthic" not in x.lower()
                            and "river" not in x.lower()
                           ]
            else:
                the_vars = [x for x in ds_contents_top.long_name if ("silicate silicate" in x  or "silicic" in x)
                            and "benthic" not in x.lower()
                            and "river" not in x.lower()
                           ]
        if vv == "oxygen":
            if fvcom == False:
                the_vars = [x for x in ds_contents.long_name if "oxygen" in x.lower() 
                            and "benthic" not in x.lower()
                            and "river" not in x.lower()
                           ]
            else:
                the_vars = [x for x in ds_contents_top.long_name if "oxygen oxygen" in x.lower() 
                            and "benthic" not in x.lower()
                            and "river" not in x.lower()
                           ]

        
        if vv == "carbon":
            ersem_vars = ds_contents_top.query("long_name in @the_vars").variable
        else:
            if vv != "co2flux" and vv != "pco2":
                ersem_vars = ds_contents.query("long_name in @the_vars").variable
            else:
                if vv == "co2flux":
                    if fvcom == False:
                        ersem_vars = ds_contents_top.query("long_name in @the_vars").variable
                    else:
                        ersem_vars = ds_contents.query("long_name in @the_vars").variable
                else:
                    if fvcom == False:
                        ersem_vars = ds_contents_top.query("long_name in @the_vars").variable
                    else:
                        ersem_vars = ds_contents.query("long_name in @the_vars").variable

        add = True

        if len(ersem_vars) > 1 and vv not in  ["doc","chlorophyll", "carbon"]:
            add = False

        if add:
            if len(ersem_vars) > 0:
                ersem_dict[vv] = "+".join(ersem_vars)
            else:
    
                if "nitrate" not in ersem_dict.keys() and vv == "nitrate":
                    if "ammonium" not in ersem_dict.keys():
                        df_nitrogen = (
                            ds_contents
                            .query("long_name.str.contains('nitrogen')")
                            .query("long_name.str.contains('nutrient')")
                            .reset_index(drop = True)
                        )
                        if len(df_nitrogen) == 1:
                            ersem_dict["nitrate"] = df_nitrogen.variable[0]
                            warnings.warn("No nitrate variable found, using nitrogen nutrient variable")

                if vv not in ersem_dict.keys():
                    ersem_dict[vv] = None 

    return ersem_dict



def match_data(ds, mapping, nan = None):
    """
    Match data
    """

    if not isinstance(mapping, dict):
        return TypeError("mapping must be a dictionary")

    the_dates = get_dates(ds)
    the_years = set(the_dates.year) 

    all_out = dict()

    # add ability to handle +

    for key in mapping.keys():
        if mapping[key] is not None:
            all_df = []
            for yy in the_years:
                infile = f"{data_dir}/ices/{key}/ices_{key}_{yy}.csv"
                if os.path.exists(infile):
                    i_df = pd.read_csv(infile)
                    all_df.append(i_df)

            
            if len(all_df) > 0:
                all_df = pd.concat(all_df).reset_index(drop = True)
                all_locs = all_df.drop(columns = key).merge(the_dates)

                if len(all_locs) == 0:
                    print(f"There are no matching times for {key}")
                else:
                    try:
                        if "+" in mapping[key]:
                            vars = mapping[key].split("+")
                            df_match = ds.match_points(all_locs, variables = vars, nan = nan)
                            df_match = eval(f"df_match.assign({key} " + "= lambda x: " + "+".join(["x." + x for x in  mapping[key].split("+")]) + ")")

                        else:
                            df_match = ds.match_points(all_locs, variables = mapping[key], nan = nan)
                        df_match = df_match.rename(columns = {df_match.columns[-1] : "model"})
                        df_match = df_match.merge(all_df)
                        df_match = df_match.rename(columns = {key:"observation"})
                        all_out[key] = df_match
                    except Exception as e:
                        warnings.warn(message = f"{e} for {key}")

    return all_out



def match_ices(path, out_dir = None, mapping = None):

    if not os.path.exists(out_dir):
        raise ValueError("out does not exist")

    if not os.path.isdir(out_dir):
        raise ValueError("out is not a directory")

    if mapping is None:
        if isinstance(path, list):
            ds = nc.open_data(path[0])
            the_map = generate_mapping(ds)
        else:

            for file_path in Path(split_path(path)[0]).glob(split_path(path)[1]):
                ds = nc.open_data(str(file_path), checks = False)
                the_map = generate_mapping(ds)
                break
        print(the_map)

        x = input("Is this mapping correct? Y/N. If not, please provide a correct one.")

        if x.lower() == "n":
            return None
    else:
        the_map = mapping


    ds = nc.open_data(path, checks = False)

    print("Matching up")
    new = match_data(ds, the_map, nan = 0)
    for key in new.keys():
        df_kk = new[key].dropna()
        if len(df_kk) > 0:
            out = f"{out_dir}/ersemobs_{key}.csv"
            df_kk.to_csv(out)

    


    #new = mm.ices.match_data(ds, the_map, nan = 0)




