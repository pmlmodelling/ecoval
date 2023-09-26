
from IPython.display import Markdown as md, display
from plotnine import *
import nctoolkit as nc
nc.options(parallel = True)
import pandas as pd
import numpy as np
import pkg_resources
import os

def global_regionals(ds_model, ds_obs, variable, i_figure):
    """
    This function is used to generate the regional plots for the global model.

    Parameters
    ----------
    ds_model : nctoolkit.Dataset
        The model dataset.
    ds_obs : nctoolkit.Dataset
        The observation dataset.
    variable : str
        The variable to plot.
    i_figure : int
        The figure number.


    """
    display(md(f"## Performance of the model in oceans"))
    Variable = ds_model.contents.long_name.values[0].title()

    df_regions = []

    for rr in ["atlantic", "pacific"]:
        ocean = rr.title() + " Ocean" 
        # extract the region from the package data

        ff = pkg_resources.resource_filename(__name__, f"data/region_{rr}.nc")

        ds_rr = nc.open_data(ff, checks = False)
        ds_rr.regrid(ds_model)
        ds_rr.run()

        ds_obs_rr = ds_obs.copy()
        ds_obs_rr * ds_rr

        ds_model_rr = ds_model.copy()
        ds_model_rr * ds_rr
        ds_model_rr.as_missing(0)
        ds_obs_rr.as_missing(0)

        try:
            ds_model_rr_zonal = ds_model_rr.copy()
            ds_obs_rr_zonal = ds_obs_rr.copy()
            ds_model_rr_zonal.zonal_mean()
            ds_obs_rr_zonal.zonal_mean()
        except:
            ds_model_rr_zonal = ds_model_rr.copy()
            ds_obs_rr_zonal = ds_obs_rr.copy()
            ds_model_rr_zonal.to_latlon(lon = [-179.5, 179.5], lat =  [-89.5, 89.5], res = 1)
            ds_obs_rr_zonal.to_latlon(lon = [-179.5, 179.5], lat =  [-89.5, 89.5], res = 1)

            ds_model_rr_zonal.zonal_mean()
            ds_obs_rr_zonal.zonal_mean()


        # extract the name of model time

        model_time = [x for x in ds_model.to_xarray().dims if "time" in x][0]
        obs_time = [x for x in ds_obs.to_xarray().dims if "time" in x][0]

        df_hover = pd.concat([(
            ds_model_rr_zonal
            .to_dataframe()
            .reset_index()
            .dropna()
            .rename(columns = {model_time: "time"})
            .assign(month = lambda x: x.time.dt.month)
            .loc[:,["month", "lat", "model"]]
            .rename(columns = {"model": "value"})
            .assign(region = rr)
            .assign(source = "model")
        ), (
            ds_obs_rr_zonal
            .to_dataframe()
            .reset_index()
            .dropna()
            .rename(columns = {obs_time: "time"})
            .assign(month = lambda x: x.time.dt.month)
            .loc[:,["month", "lat", "observation"]]
            .rename(columns = {"observation": "value"})
            .assign(region = rr)
            .assign(source = "observation")
        )]
        )

        df_regions.append(df_hover)

        del ds_model_rr, ds_obs_rr, ds_model_rr_zonal, ds_obs_rr_zonal, ds_rr

    df_regions = pd.concat(df_regions)

    label = Variable + " (" +  ds_model.contents.unit.values[0] + ")"

    df_regions = (
        df_regions
        # change region to title
        .assign(region = lambda x: x.region.str.title())
        # change source to title
        .assign(source = lambda x: x.source.str.title())
    )

    gg = (
        ggplot(df_regions)+
        geom_raster(aes(x = "month", y = "lat", fill = "value"))+
        facet_grid("region~source")+
        theme_bw()+
        scale_y_continuous(breaks = np.arange(-90, 100, 30), labels = ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])+
        scale_x_continuous(breaks = [1, 4, 7, 10], labels = ["Jan", "Apr", "Jul", "Oct"])+ 
        labs(y = None, fill = label, x= "Month")+
        theme(legend_position = "top")

    )


    # check max and min of value in df_regions

    max_value = df_regions.value.max()
    min_value = df_regions.value.min()

    # check if max_value is larger than 0 and min_value is smaller than 0

    if max_value > 0 and min_value < 0:
        # check max_value and min_value have similar magnitude
        if (max_value /-min_value) < 2 and (max_value /-min_value) > 0.5:
            gg = gg + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)



    gg = gg.draw()
    return gg
    gg.show()
    #display(gg)
    md_result = md(f"**Figure {i_figure}**: Zonal mean {variable} for model and observations.")
    i_figure += 1

    del ds_model, ds_obs


    nc.cleanup()



    return  md_result, i_figure




