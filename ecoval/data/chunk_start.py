import nctoolkit as nc
nc.options(parallel=True)
import os
import pkg_resources
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from IPython.display import Markdown as md
%load_ext rpy2.ipython
import jellyfish
import pickle
import ecoval
import copy
import calendar
data_dir = ecoval.get_obsdir(level = -2)
import nctoolkit as nc
import hvplot.xarray
import geopandas as gpd
import xarray as xr
import holoviews as hv
from holteandtalley import HolteAndTalley
import numpy as np
from mask import mask_all, mask_shelf
import pandas as pd
build = "book_build"
if build == "pdf":
    pd.set_option('styler.render.repr', 'latex')
chapter = "book_chapter"
if build == "pdf":
    chapter = f"{chapter}."
else:
    chapter = ""
import glob
import cmocean as cm
from tqdm import tqdm
import hvplot.pandas
from plotnine import *
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
i_figure = 1
i_table = 1
stamp = nc.session_info["stamp"]
out = ".trackers/" + stamp
fast_plot = fast_plot_value
if not os.path.exists(".trackers"):
    os.makedirs(".trackers")
# save out as empty file
with open(out, 'w') as f:
    f.write("")
