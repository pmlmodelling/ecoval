
import nctoolkit as nc
import pandas as pd
import os
import xarray as xr
import warnings

def get_units(ds, mapping):
    """
    """
    base_units = pd.read_csv("/data/proteus1/scratch/rwi/ersemval/data/noaa/noaa_units.csv")
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
    print("Getting years in dataset!")
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

def generate_mapping(ds):
    """
    Generate mapping of model and NOAA data
    """

    noaa_variables = [
            "doc",
     ]
    ds1 = nc.open_data(ds[0], checks = False)
    ds_contents = ds1.contents
    ersem_dict = {}
    for vv in noaa_variables:
        print(vv)
        vv_check = vv
        if vv == "doc":
            the_vars = [x for x in ds.contents.long_name if " organic" in x and "carbon" in x and "benthic" not in x and "dissolved" in x]
        
        ersem_vars = ds_contents.query("long_name in @the_vars").variable
        if len(ersem_vars) > 0:
            ersem_dict[vv] = "+".join(ersem_vars)
        else:
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
            infile = f"data/noaa/noaa_{key}.csv"
            all_df = pd.read_csv(infile)
            if len(all_df) > 0:
                all_locs = all_df.drop(columns = key).merge(the_dates)

                if len(all_locs) == 0:
                    print(f"There are no matching times for {key}")
                else:
                    #return all_locs
                    #return all_locs
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


