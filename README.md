# ecoval
An ERSEM validation tool


## Steps to carry out automated validation of model output 


First, clone this directory:

```sh
git clone https://github.com/pmlmodelling/ecoval.git
```

Then move to this directory.

```sh
cd ecoval
```


Second, set up a conda environment.

```sh
mamba env create --name ecoval -f ecoval.yml
```

Activate this environment. 


Now, install the package.

```sh
pip install .

```


You can now build the docs in two steps. First, matchup the data in Python. You need to specify a spinup period, in this case 5 years. This process might take a couple of hours to run, depending on the size of the simulation. Increase the number of cores to get faster matchups.


```sh
import ecoval
ecoval.matchup("/data/sthenno1/scratch/hpo/LOCATE/new_production_run/data/", cores = 6, spinup = 5)

```
This will put all relevant matchup data into a folder called matched. Note: this could take a couple of hours if you have a large simulation.

Ideally, the data directory specified will only have model simulation output in it, and it should have a consistent structure. The matchup function will infer the folder structure and read in all the relevant data. But if things are inconsistent, or you have stray files, things could go wrong.

Once this is done you can build the docs. This should take 10-15 minutes.


```sh
import ecoval
ecoval.validate()
```

## Corrupt files

Please ensure that there are no corrupt files in the data directory. If there are, the matchup function will probably fail. You can check for corrupt files using the following command:

```sh
    import nctoolkit as nc
    ds = nc.open_data("foo.nc")
    ds.is_corrupt()
```


## Minimal simulation requirements

Simulations should have at least one year of complete data. Matchups for gridded data will require the model to have at least monthly resolution; if it is daily gridded model output will be averaged in each month to matchup with gridded observations. 

Point observation matchups will do a strict day/month/year matchup. If you have monthly output only, any point matchups will be relatively ineffective.



