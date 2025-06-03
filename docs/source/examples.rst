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


Specifying which variables to validate
--------------------------------------

Three `matchup` arguments specify which variables to validate: `surface`, `bottom` and `benthic`.


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

