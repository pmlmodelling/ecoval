Examples
============


Matching up data
---------------------------

**Basic requirements**

Validating a simulation is a two-step process. First, you need to match up the simulation data with observational data. Then, you can validate the simulation against the matched-up data.

If you want to matchup simulation and observational data with ecoval, you need to use the `matchup` function. This has some essential requirements.

- `sim_dir`: The directory containing the simulation data.
- `surface_level`: The surface level of the simulation data. This tells ecoval if the top of bottom level of the files are the sea surface. This is normally "top".
- `cores`: The number of cores to use for the matchup. Set to a high value to speed things up.
- `start`: The first year to use for the matchups.
- `end`: The last year to use for the matchups. 

For example, if you had data in the /foo/bar directory, you could run the following code to match up the data for the years 2000 to 2010, using the top surface level and 6 cores:

.. code:: ipython3

   import ecoval
   ecoval.matchup("/foo/bar", surface_level = "top",  cores = 6, start = 2000, end = 2010)


Once you have done this, you can run `validate` to validate the simulation data against the matched-up data. 


.. code:: ipython3

   ecoval.validate()

By default, this will produce an html document. To get a pdf do the following:

.. code:: ipython3

   ecoval.validate("pdf")

Note: **validate must be run from the directory with the matched up data.**



Specifying which variables to validate
--------------------------------------

Three `matchup` arguments specify which variables to validate: `surface`, `bottom` and `benthic`.

If you do not want any of these three to be matched up, set them to None. For example, if you do not want to validate any benthic variables, you would do the following:

.. code:: ipython3

   ecoval.matchup(.., benthic = None)

Matchups for the surface come in two flavours, either using gridded or point data.
If you want to specifify what you want, you need to provide a dictionary to matchup.
So, if you only wanted to validate gridded temperature data and point chlorophyll data, you would do the following:

.. code:: ipython3

   ecoval.matchup(.., surface = {"gridded": ["temperature"], "point": ["chlorophyll"]})


Getting a list of available observational data
------------------------------------------------

A useful helper function is `list_observational_data`, which will give you a list of all the available observational data that can be used for matchups.

.. code:: ipython3

   ecoval.list_observational_data()

This will provide a complete list of observational data available that can pruned down to what you need.





Carrying out an ultra-extensive validation
------------------------------------------------

If you want to carry out a validation which is as extensive as possible, you can do the following:

.. code:: ipython3

   ecoval.matchup(.., everything = True)

This will matchup all the simulation with all available gridded and observational data across all depth layers. If you have a big simulation, this could take a weekend.


Validating plankton functional types
------------------------------------

If you want to validate plankton functional types, you need to use the pft argument in `matchup` as follows:

.. code:: ipython3

   ecoval.matchup(.., pft = True)

Note: PFT data is only available for 2011 and 2012.

Specifying the cell thickness
--------------------------------------

To carry out 3D interpolation of the simulation data, ecoval needs to know the thickness of the cells in the model.
For NEMO simulations it will assume it is the e3t variable and by default it will try to find this in the simulation files.

However, if this is not available in the files you will need to provide it for any 3D matchups with observational data.

Just create a file with the vertical thicknesses of each cell and point `matchup` towards it as follows:

.. code:: ipython3

   ecoval.matchup(.., thickenss = "/foo/bar/thickness.nc")

Controlling temporal accuracy of point matchups
------------------------------------------------

By default, ecoval will do precise matchups for point data. In other words it will work out the matchup for the precise day the observation comes from.

Sometimes this is not ideal, as you might only have one year of simulation and you want to matchup with all of the observational data to get an (imperfect) picture of how things look against observations.

In this case, you can relax the accuracy of the matchups using the point_time_res argument. This will tell you how precise things will be.

By default, this will be a list off ["year", "month", "day"], i.e. precise to the day in each year. To make things less fine-grained, change this to ["month", "day"], or "month" to get matchups to the precise day or month of the year, while ignoring the actual year in the observational data.

So, for example if you did the following:

.. code:: ipython3

   ecoval.matchup("/foo/bar", point_time_res = ["month"])

It will matchup all observational values for January, February and so on and come it with the average monthly value for the simulation.

You can also have some additional control with the `point_years` argument. This will tell ecoval which years to use from the point data.
For example, if you only wanted to use the years 2000 to 2010, you would do the following:

.. code:: ipython3

   ecoval.matchup("/foo/bar", point_years = [2000, 2010])

This can be useful if you have 1 year of simulation, but only want to compare it with data for the decade prior, not the decades long observational dataset.



Handling the location of simulation files
------------------------------------------------

The `matchup` function will look for simulation files in the directory you specify. However, it needs to make some assumptions about where the files are located.

By default, it will assume that simulation files are stored 2 directories down from the base directory, i.e. files look like ../2000/01/foo_bar.nc.

If the structure is different, you can specify the `n_dirs_down` argument. For example, if your files are all in the `sim_dir` directory and not a subdirectory, do the matchups as follows:

.. code:: ipython3

   ecoval.matchup(.., n_dirs_down = 0)

Handing dubious files in the simulation directory
------------------------------------------------

ecoval will automatically scan through the simulation directory and figure out which files are which, identify where variables are stored and so on.

In general, this works fairly well. However, it is possible you will have files stored that cause confusion. For example, you might have some post-processed files in among raw simulation output.

If you want to ignore certain files, use the `exclude` argument. This will take a list of strings and any files that partially match the string will be ignored.

So for example, if you want to ignore all files with "initial_conditions" in them, you would do the following:

.. code:: ipython3

   ecoval.matchup(.., exclude = ["initial_conditions"])


Validating a spatial subset
------------------------------------------------

Sometimes you only want to validate a spatial subset of the simulation data. For example, you might want to ignore regions close to the model boundary.

In this case you can specifify `lon_lim` and `lat_lim`, which will tell you the minimum and maximum latitudes to consider.

This would work as follows:

.. code:: ipython3

   ecoval.matchup(.., lon_lim = [-10, 10], lat_lim = [40, 50])

if you wanted to validate a region between 10 degrees west and 10 degrees east, and between 40 and 50 degrees north.


Speeding up file identification
------------------------------------------------

To identify files in the simulation directory, ecoval will look at the files in a random subdirectory and identify a mapping from variables to file names, e.g. ***foo**bar**.nc.
This is normally fast enough, as there are typically only a few files in a subdirectory. 
However, occasionally you could have all of the simulation files in a single directory. In this case you want to specify `n_check`, which tells ecoval how many randomly selected files to check.

For example, if you had 1000s of files in a directory, you might want to set `n_check` to 20 to identify things quickly.

.. code:: ipython3

   ecoval.matchup(.., n_check = 20)

Not asking for user input when running matchup
------------------------------------------------

By default, the `matchup` function will ask you for confirmation that you are happy with the proposed matchups. This is to avoid you accidentally running a matchup that will take a long time and not give you the results you want. 
If you want to skip this step, you can set the `ask` argument to False. This will run the matchups without asking for confirmation.
.. code:: ipython3

   ecoval.matchup(.., ask = False)

Validating mixed layer depth and stratifcation
------------------------------------------------

If you want to validate mixed layer depth and stratification, you can do so by specifying the `mld` argument in `matchup` as follows:

.. code:: ipython3

   ecoval.matchup(.., mld = True)

Note: this has the potential to matchup a lot of data, and maybe require **a lot** of interpolation. 


Specifying where ecoval stores matched up data
------------------------------------------------
By default, ecoval will store matched up data the directory it is run from.
If you want to specify a different directory, you can do so using the `out_dir` argument in `matchup`. For example, if you wanted to store the matched up data in "foo_bar", you would do the following: 
.. code:: ipython3

   ecoval.matchup(.., out_dir = "foo_bar")


Modifying ecoval's variable mapping
------------------------------------------------

ecoval will automatically map variables, e.g. temperature, to file patterns based on the metadata in netCDF files.
Sometimes this will fail because ecoval can either not figure things out correctly or the metadata info could be ambiguous.
If you are not happy with the way ecoval is matching up variables, type "no" when it asks you to confirm the matchups and then edit the resulting csv file.
Once that is done, re-run the `matchup` function with the `mapping` argument set to the path of the edited csv file as follows:

.. code:: ipython3

   ecoval.matchup(.., mapping = "/foo/bar/mapping.csv")


Validating using a local copy of the observational data
------------------------------------------------

By default, ecoval will use observational data stored on PML's servers. If you have a copy of this data on your local machine, you can use it instead by specifying the `obs_dir` argument in `matchup` as follows:

.. code:: ipython3

   ecoval.matchup(.., obs_dir = "/foo/bar/obs_data") 


Validating a subset of matched up variables
------------------------------------------------
If you only want to validate some of the variables you have matched up, specifify the `variables` argument in `validate` as follows:

.. code:: ipython3

   ecoval.validate(variables = ["temperature", "chlorophyll"])

Validating a spatial subset of the matched up data
------------------------------------------------

If you want to validate a spatial subset of the matched up data, you can specify `lon_lim` and `lat_lim` in `validate`. This will limit the validation to the specified longitude and latitude ranges.
For example, if you wanted to validate a region between 10 degrees west and 10 degrees east, and between 40 and 50 degrees north, you would do the following:
.. code:: ipython3

   ecoval.validate(lon_lim = [-10, 10], lat_lim = [40, 50])