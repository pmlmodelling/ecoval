# %% [markdown]
# ## Assessing model bias 
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

ds_bias.pub_plot(robust = True)


# %% tags=["remove-input"]
md(f"**Figure {i_figure}**: Bias of {layer} {vv_name} from the model. A positive bias indicates that the model overestimates the observation.  For clarity, the colorbar is limited to the 2nd and 98th percentile of the data.")
i_figure += 1
