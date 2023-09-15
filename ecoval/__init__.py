import ecoval.ices as ices
import ecoval.noaa as noaa
from ecoval.ices import match_ices
from ecoval.matchall import matchup
from ecoval.fixers import tidy_name
from ecoval.regionals import global_regionals
import pandas as pd
import glob

import nctoolkit as nc
import webbrowser
from ecoval.chunkers import add_chunks
import glob
import os
import pandas as pd

def fix_toc():
    paths = glob.glob("book/notebooks/*.ipynb")
    variables = list(pd.read_csv("matched/mapping.csv").variable)
    variables.sort()

    vv_dict = dict()
    for vv in variables:
        if vv != "ph":
            vv_paths = [os.path.basename(x) for x in paths if vv in x]
            if len(vv_paths) > 0: 
                vv_dict[vv] = vv_paths
        else:
            vv_paths = [os.path.basename(x) for x in paths if vv in x and "phos" not in x and "chlo" not in x]
            if len(vv_paths) > 0: 
                vv_dict[vv] = vv_paths
    # get summary docs
    ss_paths = [os.path.basename(x) for x in paths if "summary" in x]

    out = "book/_toc.yml"
    # write line by line to out
    with open(out, "w") as f:
        #"format: jb-book"
        x = f.write("format: jb-book\n")
        x = f.write("root: intro\n")
        x = f.write("parts:\n")
        x = f.write(f"- caption: Introduction\n")
        x = f.write("  chapters:\n")
        x = f.write(f"  - file: notebooks/000_info.ipynb\n")

        # loop over variables in each vv_dict
        # value is the file in the chapter section
        # key is the variable name, so is the section
        for vv in vv_dict.keys():
            # capitalize if not ph
            if vv != "ph":
                vv_out = vv.capitalize()
            # correct ph
            if vv == "ph":
                vv_out = "pH"
            x = f.write(f"- caption: {vv_out}\n")
            x = f.write("  chapters:\n")
            for file in vv_dict[vv]:
                x = f.write(f"  - file: notebooks/{file}\n")

        x = f.write(f"- caption: Summaries\n")
        x = f.write("  chapters:\n")
        for file in ss_paths:
            x = f.write(f"  - file: notebooks/{file}\n")




def compare(model_dict = None):
    """
    This function will compare the validition output from multiple models. 
    """

    # make a folder called book/compare
    import shutil

    import pkg_resources
    import os
    data_path = pkg_resources.resource_filename(__name__, "data/mask.py")
    if not os.path.exists("book/compare/notebooks/mask.py"):
        if not os.path.exists("book/compare/notebooks"):
            # create directory recursively
            os.makedirs("book/compare/notebooks")
        shutil.copyfile(data_path, "book/compare/notebooks/mask.py")


    if not os.path.exists("book/compare"):
        # create directory recursively
        os.makedirs("book/compare")

    # move toc etc to book/compare

    data_path = pkg_resources.resource_filename(__name__, "data/_toc.yml")

    out = "book/compare/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out)

    data_path = pkg_resources.resource_filename(__name__, "data/requirements.txt")

    out = "book/compare/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out)

    data_path = pkg_resources.resource_filename(__name__, "data/intro.md")

    out = "book/compare/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out)

    # copy config

    data_path = pkg_resources.resource_filename(__name__, "data/_config.yml")

    out = "book/compare/" + os.path.basename(data_path)

    shutil.copyfile(data_path, out)


    # copy the comparison_seasonal notebook

    # make sure the directory exists

    if not os.path.exists("book/compare/notebooks"):
        # create directory recursively
        os.makedirs("book/compare/notebooks")

    file1 = pkg_resources.resource_filename(__name__, "data/comparison_seasonal.ipynb")
    if len(glob.glob("book/compare/notebooks/*comparison_seasonal.ipynb")) == 0:
        shutil.copyfile(file1, "book/compare/notebooks/comparison_seasonal.ipynb")

    model_dict_str = str(model_dict)

    with open("book/compare/notebooks/comparison_seasonal.ipynb", "r") as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("model_dict_str", model_dict_str)


    # Write the file out again

    with open("book/compare/notebooks/comparison_seasonal.ipynb", "w") as file:
        file.write(filedata)


    # now sort out the comparison_spatial notebook

    file1 = pkg_resources.resource_filename(__name__, "data/comparison_spatial.ipynb")
    if len(glob.glob("book/compare/notebooks/*comparison_spatial.ipynb")) == 0:
        shutil.copyfile(file1, "book/compare/notebooks/comparison_spatial.ipynb")

    model_dict_str = str(model_dict)

    with open("book/compare/notebooks/comparison_spatial.ipynb", "r") as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("model_dict_str", model_dict_str)

    # Write the file out again

    with open("book/compare/notebooks/comparison_spatial.ipynb", "w") as file:
        file.write(filedata)


    # now to comparison_bias

    file1 = pkg_resources.resource_filename(__name__, "data/comparison_bias.ipynb")

    if len(glob.glob("book/compare/notebooks/*comparison_bias.ipynb")) == 0:
        shutil.copyfile(file1, "book/compare/notebooks/comparison_bias.ipynb")

    model_dict_str = str(model_dict)

    with open("book/compare/notebooks/comparison_bias.ipynb", "r") as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("model_dict_str", model_dict_str)

    # Write the file out again

    with open("book/compare/notebooks/comparison_bias.ipynb", "w") as file:
        file.write(filedata)

    # add the chunks


    # sync the notebooks

    os.system("jupytext --set-formats ipynb,py:percent book/compare/notebooks/*.ipynb")

    add_chunks()

    # fix the chunks
    os.system("jupytext --sync book/compare/notebooks/*.ipynb")

    # loop through notebooks and change fast_plot_value to fast_plot

    for ff in glob.glob("book/compare/notebooks/*.ipynb"):
        with open(ff, "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("fast_plot_value", "False")

        # Write the file out again
        with open(ff, "w") as file:
            file.write(filedata)


    os.system("jupyter-book build book/compare/")
    import webbrowser
    webbrowser.open("file://" + os.path.abspath("book/compare/_build/html/index.html"))






def validate(title="Automated model evaluation", author=None):
    # docstring
    """
    This function will run the model evaluation for all of the available datasets.

    Parameters
    ----------
    title : str
        The title of the book. Default is "Automated model evaluation"
    author : str
        The author of the book. Default is None

    Returns
    -------
    None
    """

    fast_plot = False


    # figure out of "book/notebooks" has ipynb files

    empty = True

    if len(glob.glob("book/notebooks/*.ipynb")) > 0:
        empty = False

    import os

    if empty:
        from shutil import copyfile

        if not os.path.exists("book"):
            os.mkdir("book")
        if not os.path.exists("book/notebooks"):
            os.mkdir("book/notebooks")

        # from importlib.resources import files
        import pkg_resources

        data_path = pkg_resources.resource_filename(__name__, "data/mask.py")
        if not os.path.exists("book/notebooks/mask.py"):
            copyfile(data_path, "book/notebooks/mask.py")

        data_path = pkg_resources.resource_filename(__name__, "data/000_info.ipynb")
        if not os.path.exists("book/notebooks/000_info.ipynb"):
            copyfile(data_path, "book/notebooks/000_info.ipynb")

        # data_path = files("ecoval.data").joinpath("toc.yml")
        data_path = pkg_resources.resource_filename(__name__, "data/_toc.yml")

        out = "book/" + os.path.basename(data_path)
        copyfile(data_path, out)

        data_path = pkg_resources.resource_filename(__name__, "data/requirements.txt")
        out = "book/" + os.path.basename(data_path)
        copyfile(data_path, out)

        data_path = pkg_resources.resource_filename(__name__, "data/intro.md")
        out = "book/" + os.path.basename(data_path)
        copyfile(data_path, out)

        # copy config

        data_path = pkg_resources.resource_filename(__name__, "data/_config.yml")
        out = "book/" + os.path.basename(data_path)
        # change project_title in _config.yml to title

        with open(data_path, "r") as file:
            filedata = file.read()

        # Replace the target string
        if author is not None:
            filedata = filedata.replace("author:", f"author: {author}")
        # replace title
        filedata = filedata.replace("title:", f"title: {title}")

        # Write the file out again
        with open(out, "w") as file:
            file.write(filedata)

        # copyfile(data_path, out)

        books = [
            "co2fluxes.ipynb",
            "ices_template.ipynb",
            "nsbc_template.ipynb",
            "nsbc_template_verticals.ipynb",
            "pco2.ipynb",
            "pco2fluxes.ipynb",
            "woa_template.ipynb",
        ]

        path_df = []

        if os.path.exists("matched/gridded/occci/occci_model.nc"):
            file1 = pkg_resources.resource_filename(__name__, "data/occci_chl.ipynb")
            if len(glob.glob("book/notebooks/*occci_chl.ipynb")) == 0:
                copyfile(file1, "book/notebooks/occci_chl.ipynb")
                variable = "chl"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/occci_chl.ipynb"],
                        }
                    )
                )

        if os.path.exists("matched/gridded/ostia/ostia_model.nc"):
            file1 = pkg_resources.resource_filename(
                __name__, "data/ostia_temperature.ipynb"
            )
            if len(glob.glob("book/notebooks/*ostia_temperature.ipynb")) == 0:
                copyfile(file1, "book/notebooks/ostia_temperature.ipynb")
                variable = "temperature"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/ostia_temperature.ipynb"],
                        }
                    )
                )

        if os.path.exists("matched/gridded/pco2/pco2.nc"):
            file1 = pkg_resources.resource_filename(__name__, "data/pco2.ipynb")
            if len(glob.glob("book/notebooks/*pco2.ipynb")) == 0:
                copyfile(file1, "book/notebooks/pco2.ipynb")
                variable = "pco2"
                path_df.append(
                    pd.DataFrame(
                        {"variable": [variable], "path": [f"book/notebooks/pco2.ipynb"]}
                    )
                )

        # co2 fluxes

        if os.path.exists("matched/gridded/co2fluxes/co2fluxes.nc"):
            file1 = pkg_resources.resource_filename(__name__, "data/co2fluxes.ipynb")
            if len(glob.glob("book/notebooks/*co2fluxes.ipynb")) == 0:
                copyfile(file1, "book/notebooks/co2fluxes.ipynb")
                variable = "co2fluxes"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/co2fluxes.ipynb"],
                        }
                    )
                )

        # cobe2 sst

        if os.path.exists("matched/gridded/cobe2/cobe2_model.nc"):
            file1 = pkg_resources.resource_filename(
                __name__, "data/cobe2_temperature.ipynb"
            )
            if len(glob.glob("book/notebooks/*cobe2_temperature.ipynb")) == 0:
                copyfile(file1, "book/notebooks/cobe2_temperature.ipynb")
                variable = "temperature"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/cobe2_temperature.ipynb"],
                        }
                    )
                )

        for vv in [
            "salinity",
            "silicate",
            "nitrate",
            "oxygen",
            "phosphate",
            "temperature",
        ]:
            variable = vv
            Variable = variable.title()
            if os.path.exists(f"matched/woa/woa_{vv}.nc"):
                file1 = pkg_resources.resource_filename(
                    __name__, "data/woa_template.ipynb"
                )
                if len(glob.glob(f"book/notebooks/*woa_{variable}.ipynb")) == 0:
                    count = 0
                    with open(file1, "r") as file:
                        filedata = file.read()

                    # Replace the target string
                    filedata = filedata.replace("template_variable", variable)
                    filedata = filedata.replace("template_title", Variable)

                    # Write the file out again
                    with open(f"book/notebooks/woa_{variable}.ipynb", "w") as file:
                        file.write(filedata)

                    variable = vv
                    path_df.append(
                        pd.DataFrame(
                            {
                                "variable": [variable],
                                "path": [f"book/notebooks/woa_{variable}.ipynb"],
                            }
                        )
                    )

        # ices notebooks

        for vv in ["ph", "alkalinity"]:
            variable = vv
            if vv != "ph":
                Variable = variable.title()
            else:
                Variable = "pH"
            vv_file = "../../" + f"matched/ices/all/ices_all_{vv}.csv"
            vv_file_find = f"matched/ices/all/ices_all_{vv}.csv"
            if not os.path.exists(vv_file_find):
                print(f"Not finding {vv_file}")
            if os.path.exists(vv_file_find):
                if len(glob.glob(f"book/notebooks/*ices_{variable}.ipynb")) == 0:
                    file1 = pkg_resources.resource_filename(
                        __name__, "data/ices_template.ipynb"
                    )
                    count = 0
                    with open(file1, "r") as file:
                        filedata = file.read()

                    # Replace the target string
                    out = f"book/notebooks/ices_{variable}.ipynb"
                    filedata = filedata.replace("template_file_name", vv_file)
                    filedata = filedata.replace("template_title", Variable)

                    # Write the file out again
                    with open(f"book/notebooks/ices_{variable}.ipynb", "w") as file:
                        file.write(filedata)
                    variable = vv
                    path_df.append(
                        pd.DataFrame(
                            {
                                "variable": [variable],
                                "path": [f"book/notebooks/ices_{variable}.ipynb"],
                            }
                        )
                    )
        
        # set up the ices_bottom notebooks

        #ices_vars = ['temperature', 'salinity', 'oxygen', 'phosphate', 'silicate', 'nitrate', 'ammonium']   
        for vv in ["temperature", "salinity", "oxygen", "phosphate", "silicate", "nitrate", "ammonium"]:
        # for vv in [ "nitrate"]:
            variable = vv
            if vv != "ph":
                Variable = variable.title()
            else:
                Variable = "pH"
            vv_file = "../../" + f"matched/ices/bottom/ices_bottom_{vv}.csv"
            vv_file_find = f"matched/ices/bottom/ices_bottom_{vv}.csv"
            if not os.path.exists(vv_file_find):
                print(f"Not finding {vv_file}")
            if os.path.exists(vv_file_find):
                if len(glob.glob(f"book/notebooks/*ices_bottom_{variable}.ipynb")) == 0:
                    file1 = pkg_resources.resource_filename(
                        __name__, "data/ices_bottom.ipynb"
                    )
                    count = 0
                    with open(file1, "r") as file:
                        filedata = file.read()

                    # Replace the target string
                    out = f"book/notebooks/ices_bottom_{variable}.ipynb"
                    filedata = filedata.replace("ices_variable", vv)
                    filedata = filedata.replace("template_title", Variable)

                    # Write the file out again
                    with open(f"book/notebooks/ices_bottom_{variable}.ipynb", "w") as file:
                        file.write(filedata)
                    variable = vv
                    path_df.append(
                        pd.DataFrame(
                            {
                                "variable": [variable],
                                "path": [f"book/notebooks/ices_bottom_{variable}.ipynb"],
                            }
                        )
                    )


        # see if the ices temperature matchups exist.
        if os.path.exists(f"matched/ices/all/ices_all_temperature.csv"):
            if len(glob.glob(f"book/notebooks/temperature_vertical.ipynb")) == 0:
                file1 = pkg_resources.resource_filename(
                    __name__, "data/temperature_vertical.ipynb"
                )
                count = 0
                with open(file1, "r") as file:
                    filedata = file.read()

                # Replace the target string
                vv_file = "../../" + f"matched/ices_temperature.csv"
                filedata = filedata.replace("template_file_name", vv_file)
                filedata = filedata.replace("template_title", "Temperature")

                # Write the file out again
                with open(
                    f"book/notebooks/temperature_vertical.ipynb", "w"
                ) as file:
                    file.write(filedata)
                variable = vv
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": ["temperature"],
                            "path": [f"book/notebooks/temperature_vertical.ipynb"],
                        }
                    )
                )

        # nsbc notebooks

        import nctoolkit as nc

        if os.path.exists("matched/gridded/nsbc/nsbc_model.nc"):
            ds = nc.open_data("matched/gridded/nsbc/nsbc_model.nc")

            for vv in ds.variables:
                variable = vv
                if not os.path.exists(f"book/notebooks/nsbc_{variable}.ipynb"):
                    if variable == "ph":
                        Variable = "pH"
                    else:
                        Variable = variable.title()
                    file1 = pkg_resources.resource_filename(
                        __name__, "data/nsbc_template.ipynb"
                    )
                    if len(glob.glob(f"book/notebooks/*nsbc_{variable}.ipynb")) == 0:
                        count = 0
                        with open(file1, "r") as file:
                            filedata = file.read()

                        # Replace the target string
                        filedata = filedata.replace("template_variable", variable)
                        filedata = filedata.replace("template_title", Variable)

                        # Write the file out again
                        with open(f"book/notebooks/nsbc_{variable}.ipynb", "w") as file:
                            file.write(filedata)

                        variable = vv
                        path_df.append(
                            pd.DataFrame(
                                {
                                    "variable": [variable],
                                    "path": [f"book/notebooks/nsbc_{variable}.ipynb"],
                                }
                            )
                        )

            # repeate for the verticals version
        if os.path.exists("matched/gridded/nsbc/nsbc_model_verticalmean.nc"):
            ds = nc.open_data("matched/gridded/nsbc/nsbc_model_verticalmean.nc")

            for vv in ds.variables:
                variable = vv
                if not os.path.exists(
                    f"book/notebooks/nsbc_{variable}_verticals.ipynb"
                ):
                    if variable == "ph":
                        Variable = "pH"
                    else:
                        Variable = variable.title()
                    file1 = pkg_resources.resource_filename(
                        __name__, "data/nsbc_template_verticals.ipynb"
                    )
                    if (
                        len(
                            glob.glob(
                                f"book/notebooks/*nsbc_{variable}_verticals.ipynb"
                            )
                        )
                        == 0
                    ):
                        count = 0
                        with open(file1, "r") as file:
                            filedata = file.read()

                        # Replace the target string
                        filedata = filedata.replace("template_variable", variable)
                        filedata = filedata.replace("template_title", Variable)

                        # Write the file out again
                        with open(
                            f"book/notebooks/nsbc_{variable}_verticals.ipynb", "w"
                        ) as file:
                            file.write(filedata)

                        variable = vv
                        path_df.append(
                            pd.DataFrame(
                                {
                                    "variable": [variable],
                                    "path": [
                                        f"book/notebooks/nsbc_{variable}_verticals.ipynb"
                                    ],
                                }
                            )
                        )

        # glodap notebooks

        # glodap alkalinity

        if os.path.exists("matched/gridded/glodap/glodap_alkalinity.nc"):
            variable = "alkalinity"
            Variable = variable.title()
            file1 = pkg_resources.resource_filename(
                __name__, "data/glodap_template.ipynb"
            )
            if len(glob.glob(f"book/notebooks/*glodap_{variable}.ipynb")) == 0:
                count = 0
                with open(file1, "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("template_variable", variable)
                filedata = filedata.replace("template_title", Variable)

                # Write the file out again
                with open(f"book/notebooks/glodap_{variable}.ipynb", "w") as file:
                    file.write(filedata)

                variable = "alkalinity"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/glodap_{variable}.ipynb"],
                        }
                    )
                )

        # glodap ph

        if os.path.exists("matched/gridded/glodap/glodap_ph.nc"):
            variable = "pH"
            if variable == "ph":
                Variable = variable.title()
            else:
                Variable = "pH"
            file1 = pkg_resources.resource_filename(
                __name__, "data/glodap_template.ipynb"
            )
            if len(glob.glob(f"book/notebooks/*glodap_{variable}.ipynb")) == 0:
                count = 0
                with open(file1, "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("template_variable", variable)
                filedata = filedata.replace("template_title", Variable)

                # Write the file out again
                with open(f"book/notebooks/glodap_{variable}.ipynb", "w") as file:
                    file.write(filedata)

                variable = "pH"
                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": [variable],
                            "path": [f"book/notebooks/glodap_{variable}.ipynb"],
                        }
                    )
                )

        # loop through all notebooks and replace paths

        i = 0

        # need to start by figuring out whether anything has already been run...

        i = 0
        tried = False

        for ff in [x for x in glob.glob("book/notebooks/*.ipynb") if "info" not in x]:
            try:
                i_ff = int(os.path.basename(ff).split("_")[0])
                if i_ff > i:
                    i = i_ff
                    tried = True
            except:
                pass

        i_orig = i

        if len(path_df) > 0:
            path_df = pd.concat(path_df)
            path_df = path_df.sort_values("variable").reset_index(drop=True)

        for i in range(len(path_df)):
            file1 = path_df.path.values[i]
            print(file1)
            # pad i with zeros using zfill
            i_pad = str(i + 1).zfill(3)
            new_file = (
                os.path.dirname(file1) + "/" + i_pad + "_" + os.path.basename(file1)
            )
            os.rename(file1, new_file)
            # print(key, value)

        # copy the summary.ipynb notebook and add i_pad to the name

        i = i + 2
        i_pad = str(i).zfill(3)

        shelf = False
        # figure out if we need the shelf
        try:
            ds_regions = nc.open_data(
                "/data/proteus1/scratch/rwi/evaldata/data/amm7_val_subdomains.nc"
            )
            ds_xr = ds_regions.to_xarray()
            lon_size = len(ds_xr.lon)
            lat_size = len(ds_xr.lat)

            ensemble = nc.create_ensemble("matched")
            ff = ensemble[0]
            ds = nc.open_data(ff)
            ds_xr = ds.to_xarray()
            lon_size_ds = len(ds_xr.lon)
            lat_size_ds = len(ds_xr.lat)

            if lon_size_ds == lon_size:
                if lat_size_ds == lat_size:
                    shelf = True
                else:
                    shelf = False
        except:
            shelf = False

        # copy the summary notebook for the full domain

        if shelf:
            file1 = pkg_resources.resource_filename(__name__, "data/summary.ipynb")
            if len(glob.glob(f"book/notebooks/*summary.ipynb")) == 0:
                copyfile(file1, f"book/notebooks/{i_pad}_summary_shelf.ipynb")
                # change domain_title to "Shelf"
                with open(f"book/notebooks/{i_pad}_summary_shelf.ipynb", "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("domain_title", "Shelf")
                # replace shelf_mask with True
                filedata = filedata.replace("shelf_mask", "True")

                # Write the file out again
                with open(f"book/notebooks/{i_pad}_summary_shelf.ipynb", "w") as file:
                    file.write(filedata)

                i += 1

                i_pad = str(i).zfill(3)

        file1 = pkg_resources.resource_filename(__name__, "data/summary.ipynb")
        if len(glob.glob(f"book/notebooks/*summary.ipynb")) == 0:
            copyfile(file1, f"book/notebooks/{i_pad}_summary.ipynb")
        else:
            if i > (i_orig + 1):
                initial = glob.glob(f"book/notebooks/*summary.ipynb")[0]
                copyfile(initial, f"book/notebooks/{i_pad}_summary.ipynb")
                i += 1

        # change domain_title to "Full domain"

        with open(f"book/notebooks/{i_pad}_summary.ipynb", "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("domain_title", "Full domain")
        # replace shelf_mask with False
        filedata = filedata.replace("shelf_mask", "False")

        # Write the file out again
        with open(f"book/notebooks/{i_pad}_summary.ipynb", "w") as file:
            file.write(filedata)

        # pair the notebooks using jupyter text

        os.system("jupytext --set-formats ipynb,py:percent book/notebooks/*.ipynb")

        # add the chunks
        add_chunks()

        # sync the notebooks
        #
        os.system("jupytext --sync book/notebooks/*.ipynb")

    #        return None

    # return None

    # loop through notebooks and change fast_plot_value to fast_plot

    for ff in glob.glob("book/notebooks/*.ipynb"):
        with open(ff, "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("fast_plot_value", str(fast_plot))

        # Write the file out again
        with open(ff, "w") as file:
            file.write(filedata)

    # fix the toc using the function

    fix_toc()

    os.system("jupyter-book build book/")

    import os

    stamps = [os.path.basename(x) for x in glob.glob("book/notebooks/.trackers/*")]
    stamps.append("nctoolkit_rwi_uhosarcenctoolkittmp")

    delete = []
    for x in stamps:
        delete += glob.glob("/tmp/*" + x + "*")

    for ff in delete:
        if os.path.exists(ff):
            if "nctoolkit" in x:
                os.remove(ff)
    webbrowser.open("file://" + os.path.abspath("book/_build/html/index.html"))


try:
    from importlib.metadata import version as _version
except ImportError:
    from importlib_metadata import version as _version

try:
    __version__ = _version("ecoval")
except Exception:
    __version__ = "999"
