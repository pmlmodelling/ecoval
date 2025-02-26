# %% tags=["remove-input"]
md(f"## Assessing model bias for {layer} {vv_name}")
#
# A critical metric for model performance is the bias between model and observed values. Here the bias is calculated as the mean difference between model and observed values. A positive bias indicates that the model overestimates the observed values, while a negative bias indicates that the model underestimates the observed values. 



# %% tags=["remove-input"]
ds_bias = ds_model.copy()
ds_bias - ds_obs
ds_bias.tmean()
ds_bias.run()
# model plot

df_bias = ds_bias.to_dataframe().reset_index()
if "model" == ds_bias.variables[0]:
    ds_bias.set_longnames({"model": Variable + " bias"})
else:
    ds_bias.set_longnames({variable: Variable + " bias"})

# get th
ds_ave = ds_bias.copy()
ds_ave.spatial_mean()
ds_ave.rename({ds_ave.variables[0]: "bias"})
ave_bias = ds_ave.to_dataframe().bias.values[0]
ds_summary = ds_bias.copy()
ds_summary > 0 
ds_summary.spatial_mean()
ds_summary.rename({ds_summary.variables[0]: "bias"})
positive_bias = ds_summary.to_dataframe().bias.values[0] * 100
# as a percentage to 1 dp
positive_bias = round(positive_bias, 1)
the_unit = ds_summary.contents.unit.values[0]

md(f"Figure {chapter}{i_figure} shows the average bias of {layer} {vv_name} simulated by the model. A positive bias indicates that the model overestimates the observation, while a negative bias indicates that the model overpredicts the observation.") 
if positive_bias > 50:
    positive_bias = str(positive_bias) + "%"
    md(f"The spatial average bias of {layer} {vv_name} is {ave_bias:.2f} {the_unit}. Overall, the model overestimates the observations in {positive_bias} of the model domain.")
else:
    negative_bias = str(100-positive_bias) + "%"
    md(f"The spatial average bias of {layer} {vv_name} is {ave_bias:.2f} {the_unit}. Overall, the model underestimates the observations in {negative_bias} of the model domain.")



# %% tags=["remove-input"]

if fix_grid:
    ds_bias_1 = ds_bias.copy()
    ds_bias_1.to_latlon(lon = [lon_min , lon_max], lat = [lat_min, lat_max], res = [lon_res, lat_res])
    # change the name
    ds_bias_1.set_longnames({ds_bias_1.variables[0]: "Model bias"})
    ds_bias_1.pub_plot(robust = True)
else:
    # change the name
    ds_bias.set_longnames({ds_bias.variables[0]: "Model bias"})
    ds_bias.pub_plot(robust = True)


# %% tags=["remove-input"]
md(f"**Figure {chapter}{i_figure}**: Bias of {layer} {vv_name} from the model. A positive bias indicates that the model overestimates the observation.  For clarity, the colorbar is limited to the 2nd and 98th percentile of the data.")
i_figure += 1

