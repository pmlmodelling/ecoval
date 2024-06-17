
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
        import tempfile

        print(tempfile.gettempdir())
        assert tempfile.gettempdir() == "/tmp"
        surface = {"gridded": None, "point":["nitrate", "temperature"]}
        # get the name of the temporary directory


        if os.path.exists("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv")

        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False, out_dir = "/tmp")
        # list files in /tmp recursively
        paths = glob.glob("/tmp/matched/*", recursive = True)
        print(paths)

        assert os.path.exists("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        assert os.path.exists("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv")

        df = pd.read_csv("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        df = pd.read_csv("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv")
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        if os.path.exists("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("/tmp/matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("/tmp/matched/point/nws/surface/temperature/model_surface_temperature.csv")
       
        # list all files in "/tmp/matched"
        paths = glob.glob("/tmp/matched/*")

        for ff in paths:
            if "/tmp/matched" in ff:
                try:
                    os.remove(ff)
                except:
                    pass


        paths = glob.glob("/tmp/matched/*")

        for ff in paths:
            if "/tmp/matched" in ff:
                try:
                    os.removedirs(ff, recursive = True)
                except:
                    pass

        shutil.rmtree("/tmp/matched")
        #os.removedirs("/tmp/matched")

