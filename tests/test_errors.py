
ff = "data/sst.mon.mean.nc"

import ecoval
import pytest
import nctoolkit as nc
import numpy as np
import pandas as pd
import glob
import os
import shutil


class TestFinal:
    def test_errors(self):
        import tempfile

        print(tempfile.gettempdir())
        assert tempfile.gettempdir() == "/tmp"
        surface = {"gridded": None, "point":["nitrate", "temperature"]}
        # get the name of the temporary directory
        sim_dir = "data/example"
        obs_dir = "data/evaldata"


        # test if something throws a ValueError

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = surface, bottom = ["pco2"], benthic = None, cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top", surface = "benbio", bottom = ["pco2"], benthic = None, cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = 2000, surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata",  end = 2000, surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(TypeError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = "2000", end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(TypeError):
            ecoval.matchup(sim_dir = "data/example", obs_dir = "data/evaldata", start = 2000, end = "2000",  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, point_time_res = 1, start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, point_time_res = [1], start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(TypeError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, lon_lim = 1, start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(TypeError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, lat_lim = 1, start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = "foo/bar", obs_dir = obs_dir,  start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")

        with pytest.raises(ValueError):
            ecoval.matchup( obs_dir = obs_dir,  start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature", cores = 1, ask = False, out_dir = "/tmp")
        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature",  ask = False, out_dir = "/tmp")

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, start = 2000, end = 2000,   benthic = "temperature",  ask = False, out_dir = "/tmp", cores = 1)

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, start = 2000, end = 2000,  surface_level = "foobar",  benthic = "temperature",  ask = False, out_dir = "/tmp", cores = 1)

        with pytest.raises(ValueError):
            ecoval.matchup(sim_dir = sim_dir, obs_dir = obs_dir, start = 2000, end = 2000,  surface_level = "top",  benthic = "temperature",  ask = False, out_dir = "/tmp", cores = 0)



        # list files in /tmp recursively