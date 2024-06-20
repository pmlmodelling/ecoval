
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

        # remove book_pdf directory if it exists
        if os.path.exists("book_pdf"):
            shutil.rmtree("book_pdf")
        # remove results directory if it exists
        if os.path.exists("results"):
            shutil.rmtree("results")

        if os.path.exists("validation_report.pdf"):
            os.remove("validation_report.pdf")
        ecoval.validate(build = "pdf")


        assert os.path.exists("validation_report.pdf")
        pdf = PyPDF2.PdfFileReader(open("validation_report.pdf", "rb"))
        # make sure this is 20 pages long
        assert pdf.numPages == 25
        os.remove("validation_report.pdf")
        # remove book_pdf directory
        shutil.rmtree("book_pdf")
        # results directory
        shutil.rmtree("results")

        if os.path.exists("matched/gridded/nws/nitrate/model_nitrate_surface.nc"):
            os.remove("matched/gridded/nws/nitrate/model_nitrate_surface.nc")
        if os.path.exists("matched/gridded/nws/temperature/model_temperature_surface.nc"):
            os.remove("matched/gridded/nws/temperature/model_temperature_surface.nc")
        shutil.rmtree("matched")




