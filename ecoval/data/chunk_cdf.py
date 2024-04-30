# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Cumulative distribution function 

# %% tags=["remove-input"]
md(f"The ability of the model to reproduce the broad-scale statistical distribution of {layer} {vv_name} is assessed by comparing the cumulative distribution function (CDF) of the modelled and observed {vv_name}. The CDF is a function that maps the probability that a random variable is less than or equal to a given value. The CDF is calculated by counting the number of values less than or equal to a given value and dividing by the total number of values. The CDF is a non-parametric measure of the statistical distribution of a random variable. It is a more robust measure of the statistical distribution than the mean and standard deviation, which are sensitive to outliers.")
ds = nc.open_data(f"../../results/monthly_mean/monthlymean_{variable}.nc")
ds_regions = nc.open_data(f"{data_dir}/amm7_val_subdomains.nc")
ds_regions.subset(variable = "Shelf")
ds_regions.as_missing(0)
ds_regions.regrid(ds)
ds * ds_regions
df = (
    ds
    .to_dataframe()
    .dropna()
    .reset_index()
    )
time_name = [x for x in df.columns if "time" in x][0]
df.rename(columns = {time_name: "time"}, inplace = True)
lon_name = [x for x in df.columns if "lon" in x][0]
lat_name = [x for x in df.columns if "lat" in x][0]
df.rename(columns = {lon_name: "lon", lat_name: "lat"}, inplace = True)
df = (
    df
    .assign(month = lambda x: x.time.dt.month)
    .loc[:,["lon", "lat", "model", "observation", "month"]]
    .drop_duplicates()
    )
lon_name = [x for x in df.columns if "lon" in x][0]
lat_name = [x for x in df.columns if "lat" in x][0]
# rename
df = df.melt(["lon", "lat", "month"])
units = ds.contents.unit[0]


# %% tags=["remove-input"]

%%capture --no-display
%%R -i df -i units -i variable -w 10 -h 10 --units in -r 100
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggthemes)

df <- df %>%
    group_by(variable, month) %>%
    # calculate 98th percentile of value
    summarize(limit = quantile(value, 0.98)) %>%
    inner_join(df) %>%
    # filter to values below 98th percentile
    filter(value < limit) %>%
    # remove limit column
    select(-limit)


# ggplot(iris, aes(x = Sepal.Length, y = Species, group = Species)) + 
#   geom_density_ridges(fill = "#00AFBB")

# df <- df %>%
#     mutate(month = as.factor(month))
# convert month to month name
df$month <- month.name[df$month]
# convert month to factor
df$month <- as.factor(df$month)
# ensure month is ordered
df$month <- factor(df$month, levels = month.name)


# edf plot

ggplot(df, aes(x = value, colour = variable)) +
    stat_ecdf()+
    facet_wrap(~month)+
    labs(y = "Cumulative probability", x = str_glue("{str_to_title(variable)} ({units})"))+
    theme_bw(base_size = 14)+
    theme(legend.position = "top")+
    labs(colour = NULL)+
    scale_color_fivethirtyeight()



# %% tags=["remove-input"]
md(f"**Figure {i_figure}**: Empirical distribution function of {layer} {vv_name} for model and observations for each month. This compares the distributions on the shelf across the entire domain using grid cells with model-observation matchups.")


