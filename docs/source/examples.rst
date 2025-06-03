Examples
============


Matching up data
---------------------------

**Basic requirements**

If you want to matchup simulation and observational data with ecoval, you need to use the `matchup` function. This has some essential requirements.

- sim_dir: The directory containing the simulation data.
- surface_level: The surface level of the simulation data. This tells ecoval if the top of bottom level of the files are the sea surface. This is normally "top".
- cores: The number of cores to use for the matchup. Set to a high value to speed things up.
- start: The first year to use for the matchups.
- end: The last year to use for the matchups. 

For example, if you had data in the /foo/bar directory, you could run the following code to match up the data for the years 2000 to 2010, using the top surface level and 6 cores:

.. code:: ipython3
   import ecoval
   ecoval.matchup("/foo/bar", surface_level = "top",  cores = 6, start = 2000, end = 2010)

