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



