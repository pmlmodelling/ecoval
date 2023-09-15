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

# %% tags=["remove-cell"]
ds1 = ds_model.copy()
ds1.cdo_command("setname,model")
ds1.run()
ds2 = ds_obs.copy()
ds2.cdo_command("setname,observation")
ds2.run()
ds_cor = nc.open_data([ds1.current[0], ds2.current[0]])
ds_cor.merge(match=["month"])
ds_cor.run()
ds_ts = ds_cor.copy()
ds_cor.cor_time("model", "observation")
title = f"Seasonal temporal correlation between {variable} for model and observations"
ds_cor.run()

# output to nc

out = f"../../results/temporals/{variable}_cor.nc"
if not os.path.exists(os.path.dirname(out)):
    os.makedirs(os.path.dirname(out))
ds_cor.to_nc(out, zip = True, overwrite = True)

# output to csv

df_cor = ds_cor.to_dataframe().reset_index()
df_cor = df_cor.dropna()
out = f"../../results/temporals/{variable}_cor.csv"
if not os.path.exists(os.path.dirname(out)):
    os.makedirs(os.path.dirname(out))
df_cor.to_csv(out, index = False)


# %% tags=["remove-input"]

df_cor = ds_cor.to_dataframe().reset_index()
# get range of lon and lat without missing values of cor in df_cor
lon_min = df_cor.dropna().lon.min()
lon_max = df_cor.dropna().lon.max()
lat_min = df_cor.dropna().lat.min()
lat_max = df_cor.dropna().lat.max()
lon_range = [lon_min, lon_max]
lat_range = [lat_min, lat_max]
ds_cor.subset(lon = lon_range, lat = lat_range)
ds_cor.pub_plot()



# %% tags=["remove-input"]
md(f"**Figure {i_figure}**: Seasonal temporal correlation between model and observations for {variable}. This is the Pearson correlation coefficient between climatology monthly mean values in the model and observations.")
i_figure += 1

# %% tags=["remove-cell"]
ds_regions = nc.open_data("/data/proteus1/scratch/rwi/evaldata/data/amm7_val_subdomains.nc")
ds_regions.as_missing(0)
ds_regions.set_fill(-9999)
ds_regions.run()
ds_regions.regrid(ds_model)
regions_contents = ds_regions.contents

# %% tags=["remove-cell"]

lon_size = len(ds_regions.to_xarray().lon)
lat_size = len(ds_regions.to_xarray().lat)
region_size = lon_size * lat_size

# get the model grid size

lon_size = len(ds_model.to_xarray().lon)
lat_size = len(ds_model.to_xarray().lat)
model_size = lon_size * lat_size

if model_size == region_size:
    regional = True
else:
    regional = False




# %% tags=["remove-cell"]
shape = gpd.read_file("/data/proteus1/scratch/rwi/evaldata//data/mapping/TM_WORLD_BORDERS-0.3.shp")
mod_var = ds_model.variables[0]
obs_var = ds_obs.variables[0]
# create xlim using mod_var
time_name = [x for x in list(ds_model.to_xarray().coords) if "time" in x][0]
time_name

df_model = (
    ds_model
    .to_dataframe()
    .reset_index()
    .rename(columns = {time_name: "time"})
    .assign(month = lambda x: x.time.dt.month)
    .loc[:,["lon", "lat", "month", mod_var ]]
    .rename(columns = {mod_var: "model"})
    .dropna()
)

xlim = np.array([df_model.lon.min(), df_model.lon.max()])
ylim = np.array([df_model.lat.min(), df_model.lat.max()])

df_obs = (
    ds_obs
    .to_dataframe()
    .reset_index()
)

time_name = [x for x in list(ds_obs.to_xarray().coords) if "time" in x][0]


df_obs = (
    df_obs
    .rename(columns = {time_name: "time"}  )
    .assign(month = lambda x: x.time.dt.month)
    .loc[:,["lon", "lat", "month", obs_var ]]
    .rename(columns = {obs_var: "observation"})
    .dropna()

)

df_diff = (
    df_model
    .merge(df_obs, on = ["lon", "lat", "month"])
    .assign(diff = lambda x: x.model - x.observation)
)
model_unit = ds_model.contents.unit[0]

# %% tags=["remove-input"]
# %%capture --no-display
# %%R -i model_unit -i df_model -i df_obs -i df_diff -w 800 -h 600 -i fast_plot
library(ggplot2, warn.conflicts = FALSE)
library(cowplot, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

df_model_raw <- drop_na(df_model)
df_obs_raw <- drop_na(df_obs)
df_diff_raw <- drop_na(df_diff)

if (model_unit == "degC")
    model_unit = "°C"

df_model <- drop_na(df_model)
df_model <- df_model %>%
    filter(month %in% c(1,2,3,4,5,6))
df_obs <- df_obs %>%
    filter(month %in% c(1,2,3,4,5,6))
df_diff <- df_diff %>%
    filter(month %in% c(1,2,3,4,5,6))
# 

df_model_raw %>%
    group_by(month) %>%
    summarize(model_98 = quantile(model, probs = 0.98)) %>%
    ungroup() %>%
    summarize(model_98 = mean(model_98)) %>%
    ungroup() %>%
    pull(model_98) -> model_98

obs_98 <- df_obs_raw %>%
    group_by(month) %>%
    summarize(obs_98 = quantile(observation, probs = 0.98)) %>%
    ungroup() %>%
    summarize(obs_98 = mean(obs_98)) %>%
    ungroup() %>%
    pull(obs_98)


df_model <- df_model %>%
    mutate(model = ifelse(model > model_98, model_98, model))

df_obs <- df_obs %>%
    mutate(observation = ifelse(observation > obs_98, obs_98, observation))

diff_02 <- df_diff_raw %>%
    group_by(month) %>%
    summarize(diff_02 = quantile(diff, probs = 0.02)) %>% 
    ungroup() %>%
    summarize(diff_02 = mean(diff_02)) %>%
    ungroup() %>%
    pull(diff_02)

diff_98 <- df_diff_raw %>%
    group_by(month) %>%
    summarize(diff_98 = quantile(diff, probs = 0.98)) %>%
    ungroup() %>%
    summarize(diff_98 = mean(diff_98)) %>%
    ungroup() %>%
    pull(diff_98)

df_diff <- df_diff %>% 
    mutate(diff = ifelse(diff < diff_02, diff_02, diff)) %>%
    mutate(diff = ifelse(diff > diff_98, diff_98, diff))



xlim = c(min(df_model$lon), max(df_model$lon))
ylim = c(min(df_model$lat), max(df_model$lat))

world_map <- map_data("world") 

# get 98th percentile of df_model$

# function to convert month number in int to month name
month_name <- function(x){
    month.abb[x]
}

# vectorize month_name function
# 
month_name <- Vectorize(month_name) 



# convert month number to month name in dataframes

# 

df_model <- df_model %>%
    mutate(month = month_name(month))

df_obs <- df_obs %>%
    mutate(month = month_name(month))

df_diff <- df_diff %>%
    mutate(month = month_name(month))

# first six months of the year in dataframes

df_model <- df_model %>%
    filter(month %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))

df_obs <- df_obs %>%
    filter(month %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))

df_diff <- df_diff %>%  
    filter(month %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))


# change month to suitable factor in dataframes

df_model$month <- factor(df_model$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))

df_obs$month <- factor(df_obs$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))

df_diff$month <- factor(df_diff$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))







gg1 = ggplot(df_model)+
geom_tile(aes(x = lon, y = lat, fill = model))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
theme_bw(base_size = 12)+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
labs(x = NULL, y = NULL, title = "Model", fill = model_unit)

if (fast_plot == FALSE)
    gg1 <- gg1 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

gg2 = ggplot(df_obs)+
geom_tile(aes(x = lon, y = lat, fill = observation))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
theme_bw(base_size = 12)+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
labs(x = NULL, y = NULL, title = "Observation", fill = model_unit)

if (fast_plot == FALSE)
    gg2 <-  gg2 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")


gg3 = ggplot(df_diff)+
geom_tile(aes(x = lon, y = lat, fill = diff))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
theme_bw(base_size = 12)+
scale_fill_gradient2(guide  = guide_colourbar(title.position = "right"), low = "blue", high = "red")+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Model - Observation", fill = model_unit)

if (fast_plot == FALSE)
    gg3 <- gg3 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")


cowplot::plot_grid(gg1, gg2, gg3, ncol = 1)



# %% tags=["remove-input"]
# %%capture --no-display
# %%R -i model_unit -i df_model -i df_obs -i df_diff -w 800 -h 600 -i fast_plot
library(ggplot2, warn.conflicts = FALSE)
library(cowplot, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

df_model_raw <- drop_na(df_model)
df_obs_raw <- drop_na(df_obs)
df_diff_raw <- drop_na(df_diff)


df_model_raw %>%
    group_by(month) %>%
    summarize(model_98 = quantile(model, probs = 0.98)) %>%
    ungroup() %>%
    summarize(model_98 = mean(model_98)) %>%
    ungroup() %>%
    pull(model_98) -> model_98

obs_98 <- df_obs_raw %>%
    group_by(month) %>%
    summarize(obs_98 = quantile(observation, probs = 0.98)) %>%
    ungroup() %>%
    summarize(obs_98 = mean(obs_98)) %>%
    ungroup() %>%
    pull(obs_98)



diff_02 <- df_diff_raw %>%
    group_by(month) %>%
    summarize(diff_02 = quantile(diff, probs = 0.02)) %>% 
    ungroup() %>%
    summarize(diff_02 = mean(diff_02)) %>%
    ungroup() %>%
    pull(diff_02)

diff_98 <- df_diff_raw %>%
    group_by(month) %>%
    summarize(diff_98 = quantile(diff, probs = 0.98)) %>%
    ungroup() %>%
    summarize(diff_98 = mean(diff_98)) %>%
    ungroup() %>%
    pull(diff_98)


world_map <- map_data("world") %>% 
    fortify()

df_model <- df_model_raw
df_obs <- df_obs_raw
df_diff <- df_diff_raw 

df_model <- df_model %>%
    filter(month > 6) 
df_obs <- df_obs %>%
    filter(month > 6) 
df_diff <- df_diff %>%
    filter(month > 6) 
# 


df_model <- df_model %>%
    mutate(model = ifelse(model > model_98, model_98, model))

df_obs <- df_obs %>%
    mutate(observation = ifelse(observation > obs_98, obs_98, observation))


df_diff <- df_diff %>% 
    mutate(diff = ifelse(diff < diff_02, diff_02, diff)) %>%
    mutate(diff = ifelse(diff > diff_98, diff_98, diff))



xlim = c(min(df_model$lon), max(df_model$lon))
ylim = c(min(df_model$lat), max(df_model$lat))

world_map <- map_data("world") %>% 
    fortify()

# function to convert month number in int to month name
month_name <- function(x){
    month.abb[x]
}

# vectorize month_name function
# 
month_name <- Vectorize(month_name) 



# convert month number to month name in dataframes

# 

df_model <- df_model %>%
    mutate(month = month_name(month))

df_obs <- df_obs %>%
    mutate(month = month_name(month))

df_diff <- df_diff %>%
    mutate(month = month_name(month))

    # select the final six months of the year in dataframes

df_model <- df_model %>%    
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_obs <- df_obs %>%
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_diff <- df_diff %>%
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

    # convert month name to factor

df_model$month <- factor(df_model$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_obs$month <- factor(df_obs$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_diff$month <- factor(df_diff$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))



# get 98th percentile of df_model$

gg4 = ggplot(df_model)+
geom_tile(aes(x = lon, y = lat, fill = model))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
theme_bw(base_size = 12)+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
labs(x = NULL, y = NULL, title = "Model", fill = model_unit)


if( fast_plot == FALSE)
    gg4 <- gg4 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")



gg5 = ggplot(df_obs)+
geom_tile(aes(x = lon, y = lat, fill = observation))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
theme_bw(base_size = 12)+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
labs(x = NULL, y = NULL, title = "Observation", fill = model_unit)

if (fast_plot == FALSE)
    gg5 <- gg5 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")   


gg6 = ggplot(df_diff )+
geom_tile(aes(x = lon, y = lat, fill = diff))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim,    ylim = ylim)+ 
scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))+ 
theme_bw(base_size = 12)+
scale_fill_gradient2(guide  = guide_colourbar(title.position = "right"), low = "blue", high = "red")+
theme(legend.title = element_text(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Model - Observation", fill = model_unit)


if (fast_plot == FALSE)
    gg6 <- gg6 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

cowplot::plot_grid( gg4, gg5, gg6, ncol = 1)



# %% tags=["remove-input"]
md(f"**Figure {i_figure}**: Monthly mean {variable} for the model, observation and the difference between model and observations. For clarity, the maximum values are capped to the 98th percentiles") 
i_figure += 1
