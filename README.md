[![Documentation Status](https://readthedocs.org/projects/ecoval-pml/badge/?version=latest)](https://ecoval-pml.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/pmlmodelling/ecoval/branch/main/graph/badge.svg)](https://codecov.io/gh/pmlmodelling/ecoval)
# ecoVal
Marine ecosystem model validation made easy in Python

ecoVal is a fully automated marine ecosystem model tool that reduces ecosystem model validation to a simple to follow step-by-step producedure. It is designed to validate ecosystem models stored in NetCDF format against a range of in-situ and gridded observations. Currently, it uses curated data sets for validations of models on the northwest European shelf, and global datasets elsewhere. 

If you are interested in a demo and discussions of how it can be used for your purposes, please contact [Robert Wilson](mailto:rwi@pml.ac.uk).


**Note: At present, if you want to use ecoVal, you will have to request access to the validation data. This will change in future once various rights issues have been reviewed.**

## Website
A website for ecoVal is available at [https://ecoval-pml.readthedocs.io/](https://ecoval-pml.readthedocs.io/). This contains a full user guide, and information on how to use the package.


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

```sh
mamba activate ecoval 
```


Now, install the package.

```sh
pip install .

```


You can now build the docs in two steps. First, matchup the data in Python. You can specify a start and end years for the comparisons. In this only the years from 2000-2005 are validated. This process might take a couple of hours to run, depending on the size of the simulation. Increase the number of cores to get faster matchups.


```sh
import ecoval
ecoval.matchup(sim_dir = "/foo/bar", cores = 6, start = 2000, end = 2005, surface_level = "top")

```
This will put all relevant matchup data into a folder called matched. Note: this could take a couple of hours if you have a large simulation. Note: you will have to specify whether the surface is the top or bottom level in the file structure. This is almost always the top level.

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


## Modifying jupyter notebooks produced

If you want to tweak the analysis produced by ecoval, you can do so by changing the notebooks ecoval uses to produce the validation summary.

Internally, ecoval will create a number of notebook, run them and then generate an html file. If you would like to work with one of the notebooks, you can open them in the book/notebooks directory. So long as you are using the conda environment created using the commands above, you should able to run the notebooks problem free. 

Once you have modified notebooks, you can then rebuild the validation docs. Do this by running the following commands from the same directory as `validate` was run before::



```sh
import ecoval
ecoval.rebuild()
```

## Validation data sources for the northwest European Shelf

In-situ  and gridded historical observations are used for the following variables for the northwest European shelf.

| Variable | In-situ data | Gridded data source |
| --- | --- | --- | 
| Alkalinity | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | |
| Ammonium | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Chlorophyll | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Dissolved Inorganic Carbon | PANGAEA | |
| Nitrate | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Oxygen | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) |[NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| pCO2 | [SOCAT]( https://www.socat.info/) | |
| Plankton Functional Types | [Cefas](https://www.cefas.co.uk/data-and-publications/dois/north-sea-phytoplankton-pigments-2010-to-2011/) |  |
| pH | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | |
| Phosphate | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Particulate Organic Carbon | PANGAEA | [NCEO](https://catalogue.ceda.ac.uk/uuid/480394782006415897d7715dd3f5ceda) |
| Salinity | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Silicate | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx)| [NSBC](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html) |
| Temperature | [ICES](https://www.ices.dk/data/data-portals/Pages/ocean.aspx) | [CMEMS](https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/description) |










