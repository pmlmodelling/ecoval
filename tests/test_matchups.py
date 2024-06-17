
ff = "data/sst.mon.mean.nc"

import ecoval
import numpy as np
import pandas as pd
import os

surface = {"gridded": None, "point":["nitrate", "temperature"]}

class TestFinal:
    def test_matchups(self):
        if os.path.exists("matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("matched/point/nws/surface/temperature/model_surface_temperature.csv")

        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False)

        assert os.path.exists("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        assert os.path.exists("matched/point/nws/surface/temperature/model_surface_temperature.csv")

        df = pd.read_csv("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        df = pd.read_csv("matched/point/nws/surface/temperature/model_surface_temperature.csv")
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        if os.path.exists("matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("matched/point/nws/surface/temperature/model_surface_temperature.csv")
