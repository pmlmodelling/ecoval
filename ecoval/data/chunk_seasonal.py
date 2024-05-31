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
if fix_grid:
    ds_cor_plot = ds_cor.copy()
    ds_cor_plot.to_latlon(lon = [lon_min , lon_max], lat = [lat_min, lat_max], res = [lon_res, lat_res])
    ds_cor_plot.pub_plot("cor")
else:
    ds_cor.pub_plot("cor")



# %% tags=["remove-input"]
md(f"**Figure {chapter}{i_figure}**: Seasonal temporal correlation between model and observations for {layer} {vv_name}. This is the Pearson correlation coefficient between climatology monthly mean values in the model and observations.")
i_figure += 1

# %% tags=["remove-cell"]
lon_range = lon_max - lon_min
lat_range = lat_max - lat_min

global_grid = False
if lon_range > 340:
    if lat_range > 160:
        global_grid = True

if not global_grid:
    data_path = pkg_resources.resource_filename("ecoval", "data/amm7_val_subdomains.nc")
else:
    data_path = pkg_resources.resource_filename("ecoval", "data/global_subdomains.nc")
ds_regions = nc.open_data(data_path, checks = False)
# pull this in from the package data

ds_regions.as_missing(0)
ds_regions.set_fill(-9999)
ds_regions.run()
ds_regions.regrid(ds_model)
regions_contents = ds_regions.contents

# figure out if you can sensibly do a regional analysis for nws
grid = pd.read_csv("../../matched/model_grid.csv")
lon = grid.loc[:,[x for x in grid.columns if "lon" in x]].values
lon = np.unique(lon)
lon.sort()
lat = grid.loc[:,[x for x in grid.columns if "lat" in x]].values
lat = np.unique(lat)
lat.sort()
# get unique values in grid and sort them
lon = np.unique(lon)
lon.sort()
lon_min = lon.min()
lon_max = lon.max()
lat_min = lat.min()
lat_max = lat.max()
regional = False
if lon_min > -30:
    if lon_max < 15:
        if lat_min > 35:
            if lat_max < 70:
                regional = True
if global_grid:
    regional = True


# %% tags=["remove-cell"]





# %% tags=["remove-cell"]
#shape = gpd.read_file(f"{data_dir}/mapping/TM_WORLD_BORDERS-0.3.shp")
mod_var = ds_model.variables[0]
obs_var = ds_obs.variables[0]
# create xlim using mod_var
time_name = [x for x in list(ds_model.to_xarray().coords) if "time" in x][0]
lon_name = [x for x in list(ds_obs.to_xarray().coords) if "lon" in x][0]
lat_name = [x for x in list(ds_obs.to_xarray().coords) if "lat" in x][0]
time_name

df_model = (
    ds_model
    .to_dataframe()
    .reset_index()
    .rename(columns = {time_name: "time"})
    .rename(columns = {lon_name: "lon"})
    .rename(columns = {lat_name: "lat"})
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
lon_name = [x for x in list(ds_obs.to_xarray().coords) if "lon" in x][0]
lat_name = [x for x in list(ds_obs.to_xarray().coords) if "lat" in x][0]


df_obs = (
    df_obs
    .rename(columns = {time_name: "time"}  )
    .rename(columns = {lon_name: "lon"}  )
    .rename(columns = {lat_name: "lat"}  )
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
from ecoval.utils import get_extent
raw_extent = get_extent(ds_annual[0])
if np.abs(raw_extent[0] - df_model.lon.min()) > 3:
    # convert longitude to -180-180
    df_model["lon" ] = [x if x < 180 else x -360 for x in df_model.lon]
    df_obs["lon" ] = [x if x < 180 else x -360 for x in df_obs.lon]
    df_diff["lon" ] = [x if x < 180 else x -360 for x in df_diff.lon]

# generate a temporary csv file name in /tmp
# create adhoc dir if not
if not os.path.exists("adhoc/tmp"):
    os.makedirs("adhoc/tmp")
tmp_csv = f"adhoc/df_obs_model.feather"
df_obs.to_feather(tmp_csv)
tmp_csv = f"adhoc/df_model_model.feather"
df_model.to_feather(tmp_csv)
tmp_csv = f"adhoc/df_diff_model.feather"
df_diff.to_feather(tmp_csv)



# %% tags=["remove-input"]
# %%capture --no-display
# %%R -i model_unit  -w 800 -h 600 -i fast_plot -i variable -i raw_extent -r 100


library(tidyverse, warn.conflicts = FALSE)
library(cowplot, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

df_model <- arrow::read_feather("adhoc/df_model_model.feather")


if("month" %in% colnames(df_model)){

    model_unit = str_replace(model_unit, "/m\\^3", "m<sup>-3</sup>")

    df_model_raw <- drop_na(df_model)
    df_obs <- arrow::read_feather("adhoc/df_obs_model.feather")
    df_diff <- arrow::read_feather("adhoc/df_diff_model.feather")


#df_model_raw <- drop_na(df_model)
df_obs_raw <- drop_na(df_obs)
df_diff_raw <- drop_na(df_diff)

if (variable == "temperature" )
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




# lon_label = c("20°W", "10°W", "0°", "10°E")
# lat_label = c("45°N", "50°N", "55°N", "60°N", "65°N")
# Create sensible lon/lat labels based on the min/min lon and lat
# This needs to work on any data, including global data
lon_breaks = c(xlim[1], xlim[1] + 10, 0, xlim[2] - 10)
lat_breaks = c(ylim[1], ylim[1] + 5, ylim[1] + 10, ylim[1] + 15, ylim[2])
lon_labels = c(paste0(round(xlim[1]), "°W"), paste0(round(xlim[1] + 10), "°W"), "0°", paste0(round(xlim[2] - 10), "°E"))
lat_labels = c(paste0(round(ylim[1]), "°N"), paste0(round(ylim[1] + 5), "°N"), paste0(round(ylim[1] + 10), "°N"), paste0(round(ylim[1] + 15), "°N"), paste0(round(ylim[2]), "°N"))



gg1 = ggplot(df_model)+
geom_tile(aes(x = lon, y = lat, fill = model))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
theme_bw(base_size = 12)+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
labs(x = NULL, y = NULL, title = "Model", fill = model_unit)

y_labels <-  as.numeric(na.omit(layer_scales(gg1)$y$break_positions()))
x_labels <- as.numeric(na.omit(layer_scales(gg1)$x$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate
# i.e. 10 should be 10 °N, -10 should be 10 °S

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg1 <- gg1 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)

if (fast_plot == FALSE)
    gg1 <- gg1 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

gg2 = ggplot(df_obs)+
geom_tile(aes(x = lon, y = lat, fill = observation))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
theme_bw(base_size = 12)+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Observation", fill = model_unit)



y_labels <-  as.numeric(na.omit(layer_scales(gg2)$y$break_positions()))
x_labels <- as.numeric(na.omit(layer_scales(gg2)$x$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate
# i.e. 10 should be 10 °N, -10 should be 10 °S

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg2 <- gg2 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)

if (fast_plot == FALSE)
    gg2 <-  gg2 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

gg3 = ggplot(df_diff)+
geom_tile(aes(x = lon, y = lat, fill = diff))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
theme_bw(base_size = 12)+
scale_fill_gradient2(guide  = guide_colourbar(title.position = "right"), low = "blue", high = "red")+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Model - Observation", fill = model_unit)



y_labels <-  as.numeric(na.omit(layer_scales(gg3)$y$break_positions()))
x_labels <- as.numeric(na.omit(layer_scales(gg3)$x$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate
# i.e. 10 should be 10 °N, -10 should be 10 °S

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg3 <- gg3 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)

if (fast_plot == FALSE)
    gg3 <- gg3 + 
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")


gg1 <- gg1 + 
   theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) # Left ggmargin 


# gg1
# reduce the size of the plot
# options(repr.plot.width = 10, repr.plot.height = 3)
# gg1
cowplot::plot_grid(gg1, gg2, gg3, ncol = 1)

}



# %% tags=["remove-input"]
# %%capture --no-display
# %%R -i model_unit -w 800 -h 600 -i fast_plot -i variable -i raw_extent -r 100

library(tidyverse, warn.conflicts = FALSE)
library(cowplot, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

df_model <- arrow::read_feather("adhoc/df_model_model.feather")

if ("month" %in% colnames(df_model)){

    df_model_raw <- drop_na(df_model)
    df_obs <- arrow::read_feather("adhoc/df_obs_model.feather")
    df_diff <- arrow::read_feather("adhoc/df_diff_model.feather")


df_obs_raw <- drop_na(df_obs)
df_diff_raw <- drop_na(df_diff)

if (variable == "temperature" )
    model_unit = "°C"

    model_unit = str_replace(model_unit, "/m\\^3", "m<sup>-3</sup>")

df_model <- drop_na(df_model)
df_model <- df_model %>%
    filter(month > 6) 
df_obs <- df_obs %>%
    filter(month > 6) 
df_diff <- df_diff %>%
    filter(month > 6) 

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
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_obs <- df_obs %>%
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_diff <- df_diff %>%
    filter(month %in% c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

# change month to suitable factor in dataframes

df_model$month <- factor(df_model$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_obs$month <- factor(df_obs$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

df_diff$month <- factor(df_diff$month, levels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

# Create sensible lon/lat labels based on the min/min lon and lat
# This needs to work on any data, including global data
lon_breaks = c(xlim[1], xlim[1] + 10, 0, xlim[2] - 10)
lat_breaks = c(ylim[1], ylim[1] + 5, ylim[1] + 10, ylim[1] + 15, ylim[2])
lon_labels = c(paste0(round(xlim[1]), "°W"), paste0(round(xlim[1] + 10), "°W"), "0°", paste0(round(xlim[2] - 10), "°E"))
lat_labels = c(paste0(round(ylim[1]), "°N"), paste0(round(ylim[1] + 5), "°N"), paste0(round(ylim[1] + 10), "°N"), paste0(round(ylim[1] + 15), "°N"), paste0(round(ylim[2]), "°N"))



gg1 = ggplot(df_model)+
geom_tile(aes(x = lon, y = lat, fill = model))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
theme_bw(base_size = 12)+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
labs(x = NULL, y = NULL, title = "Model", fill = model_unit)


x_labels <- as.numeric(na.omit(layer_scales(gg1)$x$break_positions()))
y_labels <- as.numeric(na.omit(layer_scales(gg1)$y$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg1 <- gg1 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")


gg2 = ggplot(df_obs)+
geom_tile(aes(x = lon, y = lat, fill = observation))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
scale_fill_viridis_c(guide  = guide_colourbar(title.position = "right"))+
theme_bw(base_size = 12)+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Observation", fill = model_unit)


x_labels <- as.numeric(na.omit(layer_scales(gg2)$x$break_positions()))
y_labels <- as.numeric(na.omit(layer_scales(gg2)$y$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg2 <- gg2 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

gg3 = ggplot(df_diff)+
geom_tile(aes(x = lon, y = lat, fill = diff))+
facet_wrap(~month, ncol = 6)+
coord_cartesian(xlim = xlim, ylim = ylim)+ 
theme_bw(base_size = 12)+
scale_fill_gradient2(guide  = guide_colourbar(title.position = "right"), low = "blue", high = "red")+
theme(legend.title = ggtext::element_markdown(angle = -90), legend.title.align = 0.5)+
labs(x = NULL, y = NULL, title = "Model - Observation", fill = model_unit)


x_labels <- as.numeric(na.omit(layer_scales(gg3)$x$break_positions()))
y_labels <- as.numeric(na.omit(layer_scales(gg3)$y$break_positions()))
x_breaks <- x_labels
y_breaks <- y_labels

# y labels are north-south coordinates. Make them more appropriate

y_labels <- ifelse(y_labels >= 0, paste0(y_labels, "°N"), paste0(abs(y_labels), "°S"))
x_labels <- ifelse(x_labels >= 0, paste0(x_labels, "°E"), paste0(abs(x_labels), "°W"))

gg3 <- gg3 + scale_y_continuous(breaks = y_breaks, labels = y_labels)+
    scale_x_continuous(breaks = x_breaks, labels = x_labels)+
        geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey")

gg1 <- gg1 + 
   theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) # Left ggmargin 

# figure out if it's a global file
if(abs(raw_extent[1] - raw_extent[2]) > 350){
    gg1 <- gg1 + 
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), labels = c("180°W", "90°W", "0°", "90°E", "180°E"))+
    scale_y_continuous(breaks = c(-90, -45, 0, 45, 90), labels = c("90°S", "45°S", "0°", "45°N", "90°N"))

    gg2 <- gg2 +
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), labels = c("180°W", "90°W", "0°", "90°E", "180°E"))+
    scale_y_continuous(breaks = c(-90, -45, 0, 45, 90), labels = c("90°S", "45°S", "0°", "45°N", "90°N"))

    gg3 <- gg3 +
    scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), labels = c("180°W", "90°W", "0°", "90°E", "180°E"))+
    scale_y_continuous(breaks = c(-90, -45, 0, 45, 90), labels = c("90°S", "45°S", "0°", "45°N", "90°N"))

}

# appropriate plotting for nws 
if((raw_extent[1] > -30) & (raw_extent[2] < 20)){
    gg1 + 
    scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))
    gg2 + 
    scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))
    gg3 +
    scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous(breaks = c(45, 50, 55, 60, 65), labels = c("45°N", "50°N", "55°N", "60°N", "65°N"))
}

# gg1
# reduce the size of the plot
# options(repr.plot.width = 10, repr.plot.height = 3)
# gg1
cowplot::plot_grid(gg1, gg2, gg3, ncol = 1)

}


# %% tags=["remove-input"]
md(f"**Figure {chapter}{i_figure}**: Monthly mean {layer} {vv_name} for the model, observation and the difference between model and observations. For clarity, the maximum values are capped to the 98th percentiles") 
i_figure += 1
