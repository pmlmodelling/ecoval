
ff = "data/sst.mon.mean.nc"

import ecoval
import nctoolkit as nc
import numpy as np
import pandas as pd
import glob
import os
import shutil


class TestFinal:
    def test_point(self):

        x = ecoval.matchup_starting_point("data/example/2000/01/amm7_1d_20000101_20000131_grid_T.nc", obs_dir = "data/evaldata")

        x_split = x.split("\n")
        y = "ecoval.matchup(folder = 'data/example', cores = 6, start = -10000, end = 10000, surface_level = 'top/bottom', surface = {'gridded': ['nitrate','temperature'], 'point': ['temperature','nitrate']}, bottom = [], benthic = [])"
        y_split = "ecoval.matchup(folder = 'data/example', cores = 6, start = -10000, end = 10000, surface_level = 'top/bottom', surface = {'gridded': ['nitrate','temperature'], 'point': ['temperature','nitrate']}, bottom = [], benthic = [])".split("\n")

        # compare each element in x_split and y_split
        for i in range(len(x_split)):
            assert x_split[i] == y_split[i]


        assert len(x) == len(y)