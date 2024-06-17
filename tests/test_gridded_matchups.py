
ff = "data/sst.mon.mean.nc"

import ecoval
import nctoolkit as nc
import numpy as np
import pandas as pd
import glob
import os
import shutil


class TestFinal:

    def test_gridded(self):

        surface = {"gridded": ["nitrate", "temperature"], "point":[]}
        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False, out_dir = "/tmp")

        assert os.path.exists("/tmp/matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        assert os.path.exists("/tmp/matched/gridded/nws/temperature/model_temperature_surface.nc")

        ds = nc.open_data("/tmp/matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        df = ds.to_dataframe().reset_index()
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        ds = nc.open_data("/tmp/matched/gridded/nws/temperature/model_temperature_surface.nc")
        df = ds.to_dataframe().reset_index()
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        if os.path.exists("/tmp/matched/gridded/nws/nitrate/model_nitrate_surface.nc"):
            os.remove("/tmp/matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        if os.path.exists("/tmp/matched/gridded/nws/temperature/model_temperature_surface.nc"):
            os.remove("/tmp/matched/gridded/nws/temperature/model_temperature_surface.nc")

        shutil.rmtree("/tmp/matched")




