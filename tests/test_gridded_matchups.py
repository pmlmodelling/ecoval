
import ecoval
import PyPDF2
import nctoolkit as nc
import numpy as np
import pandas as pd
import glob
import os
import shutil


class TestFinal:

    def test_gridded(self):

        surface = {"gridded": ["nitrate", "temperature"], "point":[]}
        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False)

        assert os.path.exists("matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        assert os.path.exists("matched/gridded/nws/temperature/model_temperature_surface.nc")

        ds = nc.open_data("matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        df = ds.to_dataframe().reset_index()
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01

        ds = nc.open_data("matched/gridded/nws/temperature/model_temperature_surface.nc")
        df = ds.to_dataframe().reset_index()
        assert np.corrcoef(df["observation"], df["model"])[0,1] > 0.999
        # ensure average abs difference < 0.01
        assert np.mean(np.abs(df["observation"] - df["model"])) < 0.01
        if os.path.exists("book_html"):
            shutil.rmtree("book_html")
        # results directory
        if os.path.exists("results"):
            shutil.rmtree("results")

        ecoval.validate(build = "html", test = True)

        ff = [x for x in glob.glob("book_html/_build/html/notebooks/*") if "nitrate" in x]
        assert len(ff) == 1
        ff = ff[0]
        # ff = "book_html/_build/html/notebooks/001_model_nitrate.html"
        line = "This is getting to the end!"
        #read in and  identify if line is in the file
        with open(ff, 'r') as file:
            filedata = file.read()
            assert line in filedata

        # ff = "book_html/_build/html/notebooks/002_model_temperature.html"
        ff = [x for x in glob.glob("book_html/_build/html/notebooks/*") if "temperature" in x]
        assert len(ff) == 1
        ff = ff[0]
        line = "This is getting to the end!"
        #read in and  identify if line is in the file
        with open(ff, 'r') as file:
            filedata = file.read()
            assert line in filedata
        # remove book_pdf directory
        shutil.rmtree("book_html")
        # results directory
        shutil.rmtree("results")

        #if os.path.exists("matched/gridded/nws/nitrate/model_nitrate_surface.nc"):
        #    os.remove("matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        #if os.path.exists("matched/gridded/nws/temperature/model_temperature_surface.nc"):
        #    os.remove("matched/gridded/nws/temperature/model_temperature_surface.nc")
        shutil.rmtree("matched")




