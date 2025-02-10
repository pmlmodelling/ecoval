# %% [markdown]

# %% tags=["remove-cell"]

ds_annual = ds_model.copy()
ds_annual.rename({ds_annual.variables[0]: "model"})
ds_annual.append(ds_obs)
ds_annual.tmean()
ds_annual.merge("variables")
ds_annual.rename({ds_obs.variables[0]: "observation"})
out_dir = "../../results/annual_mean/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = out_dir + f"annualmean_{variable}_{source}.nc"
ds_annual.to_nc(out_file, zip = True, overwrite = True)
# Calculate the monthly mean and output it

ds_monthly = ds_model.copy()
ds_monthly.rename({ds_monthly.variables[0]: "model"})
ds_monthly.append(ds_obs)
ds_monthly.tmean("month")
ds_monthly.merge("variables", ["month"])
ds_monthly.rename({ds_obs.variables[0]: "observation"})
out_dir = "../../results/monthly_mean/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = out_dir + f"monthlymean_{variable}_{source}.nc"
ds_monthly.to_nc(out_file, zip = True, overwrite = True)

 
