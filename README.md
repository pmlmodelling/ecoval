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

Once this is done you can build the docs. This should take 10-15 minutes.


```sh
import ecoval
ecoval.validate()
```

