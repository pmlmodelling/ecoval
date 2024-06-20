
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

        directory = "matched/point/nws/surface/nitrate/"
        # create the directory, recursively
        os.makedirs(directory, exist_ok = True)
        directory = "matched/point/nws/surface/temperature/"
        # create the directory, recursively
        os.makedirs(directory, exist_ok = True)

        if os.path.exists("matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("matched/point/nws/surface/temperature/model_surface_temperature.csv")

        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False)
        # list files in /tmp recursively
        # paths = glob.glob("/tmp/matched/*", recursive = True)
        ff =  "matched/point/nws/surface/nitrate/model_surface_nitrate.csv"
        # when was this file modified
        ff_time = os.path.getmtime(ff)
        ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = [], benthic = None, cores = 1, ask = False,  overwrite = False)

        ff = "matched/point/nws/surface/nitrate/model_surface_nitrate.csv"
        # when was this file modified
        new_time = os.path.getmtime(ff)
        assert ff_time == new_time

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

        # remove book_pdf directory if it exists
        if os.path.exists("book_html"):
            shutil.rmtree("book_html")
        # remove results directory if it exists
        if os.path.exists("results"):
            shutil.rmtree("results")

        if os.path.exists("validation_report.pdf"):
            os.remove("validation_report.pdf")
        ecoval.validate(build = "html", test = True)

        if os.path.exists("matched/point/nws/surface/nitrate/model_surface_nitrate.csv"):
            os.remove("matched/point/nws/surface/nitrate/model_surface_nitrate.csv")
        if os.path.exists("matched/point/nws/surface/temperature/model_surface_temperature.csv"):
            os.remove("matched/point/nws/surface/temperature/model_surface_temperature.csv")

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
        ff = [x for x in glob.glob("book_html/_build/html/notebooks/*") if "temperature" in x][0]
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

        shutil.rmtree("matched")
        #os.removedirs("/tmp/matched")

