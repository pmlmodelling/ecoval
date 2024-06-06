from ecoval.utils import get_datadir
import pkg_resources
import pandas as pd
import shutil
import glob
import subprocess
import nctoolkit as nc
from ecoval.matchall import matchup
from ecoval.trends import trends
from ecoval.fixers import tidy_name
from ecoval.regionals import global_regionals
from ecoval.helpers import matchup_starting_point 
from ecoval.session import session_info
import webbrowser
from ecoval.chunkers import add_chunks
import os
import re


def fix_toc(book_dir):
    paths = glob.glob(f"{book_dir}/notebooks/*.ipynb")
    variables = list(pd.read_csv("matched/mapping.csv").variable)
    if len([x for x in paths if "pft" in x]) > 0:
        variables.append("pft")
    variables.sort()

    vv_dict = dict()
    for vv in variables:
        if vv != "ph":
            vv_paths = [os.path.basename(x) for x in paths if vv in x]
            if len(vv_paths) > 0:
                vv_dict[vv] = vv_paths
        else:
            vv_paths = [
                os.path.basename(x)
                for x in paths
                if vv in x and "phos" not in x and "chlo" not in x and "occci" not in x
            ]
            if len(vv_paths) > 0:
                vv_dict[vv] = vv_paths
    # get summary docs
    ss_paths = [os.path.basename(x) for x in paths if "summary" in x]

    out = f"{book_dir}/_toc.yml"

    # write line by line to out
    i_chapter = 1
    with open(out, "w") as f:
        # "format: jb-book"
        x = f.write("format: jb-book\n")
        x = f.write("root: intro\n")
        x = f.write("parts:\n")
        x = f.write(f"- caption: Introduction\n")
        x = f.write("  chapters:\n")
        x = f.write(f"  - file: notebooks/000_info.ipynb\n")
        # open notebook and replace book_chapter with i_chapter
        with open(f"{book_dir}/notebooks/000_info.ipynb", "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("book_chapter", str(i_chapter))

        # Write the file out again
        with open(f"{book_dir}/notebooks/000_info.ipynb", "w") as file:
            file.write(filedata)
        i_chapter += 1




        x = f.write(f"- caption: Summaries\n")
        x = f.write("  chapters:\n")
        for ff in ss_paths:
            x = f.write(f"  - file: notebooks/{ff}\n")
            # open notebook and replace book_chapter with i_chapter
            with open(f"{book_dir}/notebooks/{ff}", "r") as file:
                filedata = file.read()

            filedata = filedata.replace("book_chapter", str(i_chapter))

            # Replace the target string

            # Write the file out again
            with open(f"{book_dir}/notebooks/{ff}", "w") as file:
                file.write(filedata)
            i_chapter += 1


        # loop over variables in each vv_dict
        # value is the file in the chapter section
        # key is the variable name, so is the section
        for vv in vv_dict.keys():
            # capitalize if not ph
            if vv != "ph":
                vv_out = vv.capitalize()
            if vv.lower() in ["poc", "doc"]:
                if vv.lower() == "poc":
                    vv_out = "Particulate Organic Carbon"
                if vv.lower() == "doc":
                    vv_out = "Dissolved Organic Carbon"
            if vv == "pco2":
                vv_out = "pCO2"
            # correct ph
            if vv == "ph":
                vv_out = "pH"
            if vv.lower() == "benbio":
                vv_out = "Benthos"
            if vv.lower() == "pft":
                vv_out = "Plankton Functional Types"
            x = f.write(f"- caption: {vv_out}\n")
            x = f.write("  chapters:\n")
            for ff in vv_dict[vv]:
                x = f.write(f"  - file: notebooks/{ff}\n")

                # open notebook and replace book_chapter with i_chapter
                with open(f"{book_dir}/notebooks/{ff}", "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("book_chapter", str(i_chapter))

                # Write the file out again

                with open(f"{book_dir}/notebooks/{ff}", "w") as file:
                    file.write(filedata)
                i_chapter += 1




def fix_toc_comparison():
    book_dir = "book"

    out = f"{book_dir}/compare/_toc.yml"
    # write line by line to out
    with open(out, "w") as f:
        # "format: jb-book"
        x = f.write("format: jb-book\n")
        x = f.write("root: intro\n")
        x = f.write("parts:\n")
        x = f.write(f"- caption: Comparisons with gridded surface observations\n")
        x = f.write("  chapters:\n")
        x = f.write(f"  - file: notebooks/comparison_bias.ipynb\n")
        x = f.write(f"  - file: notebooks/comparison_spatial.ipynb\n")
        x = f.write(f"  - file: notebooks/comparison_seasonal.ipynb\n")
        x = f.write(f"  - file: notebooks/comparison_regional.ipynb\n")
        x = f.write(f"- caption: Comparisons with point observations\n")
        x = f.write("  chapters:\n")
        x = f.write(f"  - file: notebooks/comparison_point_surface.ipynb\n")
        x = f.write(f"  - file: notebooks/comparison_point_benthic.ipynb\n")
        x = f.write(f"  - file: notebooks/comparison_point_bottom.ipynb\n")


def compare(model_dict=None):
    """
    This function will compare the validition output from multiple models.

    Parameters
    ----------
    model_dict : dict
        A dictionary of model names and the paths to the validation output. Default is None.
        Example: {"model1": "/path/to/model1", "model2": "/path/to/model2"}
        If the models have different grids, put the model with the smallest grid first.


    """
    if os.path.exists("book"):
        #get user input to decide if it should be removed
        user_input = input("book directory already exists. This will be emptied and replaced. Do you want to proceed? (y/n): ")
        if user_input.lower() == "y":
            while True:
                files = glob.glob("book/**/**/**", recursive=True)
                # list all files in book, recursively
                for ff in files:
                    if ff.startswith("book"):
                        try:
                            os.remove(ff)
                        except:
                            pass
                files = glob.glob("book/**/**/**", recursive=True)
                # only list files
                files = [x for x in files if os.path.isfile(x)]
                if len(files) == 0:
                    break
        else:
            print("Exiting")
            return None
                # if not os.path.exists(book_dir):
                    # break

    # make a folder called book/compare

    data_path = pkg_resources.resource_filename(__name__, "data/mask.py")
    if not os.path.exists(f"book/compare/notebooks/mask.py"):
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

    fix_toc_comparison()

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

    # move the regional book

    file1 = pkg_resources.resource_filename(__name__, "data/comparison_regional.ipynb")
    if len(glob.glob("book/compare/notebooks/*comparison_regional.ipynb")) == 0:
        shutil.copyfile(file1, "book/compare/notebooks/comparison_regional.ipynb")

    model_dict_str = str(model_dict)

    with open("book/compare/notebooks/comparison_regional.ipynb", "r") as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("model_dict_str", model_dict_str)

    # Write the file out again

    with open("book/compare/notebooks/comparison_regional.ipynb", "w") as file:
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

    # figure out if both simulations have point data

    i = 0    
    
    if i == 0:
        for ss in ["surface", "bottom", "benthic"]:
            file1 = pkg_resources.resource_filename(__name__, "data/comparison_point.ipynb")

            if len(glob.glob(f"book/compare/notebooks/*comparison_point_{ss}.ipynb")) == 0:
                shutil.copyfile(file1, f"book/compare/notebooks/comparison_point_{ss}.ipynb")

            model_dict_str = str(model_dict)

            with open(f"book/compare/notebooks/comparison_point_{ss}.ipynb", "r") as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace("model_dict_str", model_dict_str)

            # Write the file out again

            with open(f"book/compare/notebooks/comparison_point_{ss}.ipynb", "w") as file:
                file.write(filedata)
            # replace layer in the notebook with ss
            with open(f"book/compare/notebooks/comparison_point_{ss}.ipynb", "r") as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace("layer", ss)

            # Write the file out again

            with open(f"book/compare/notebooks/comparison_point_{ss}.ipynb", "w") as file:
                file.write(filedata)


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


def validate(title="Automated model evaluation", author=None, variables = "all", r_warnings = False, build = "html", model = None):
    # docstring
    """
    This function will run the model evaluation for all of the available datasets.

    Parameters
    ----------
    title : str
        The title of the book. Default is "Automated model evaluation"
    author : str
        The author of the book. Default is None
    variables : str or list
        The variables to run the model evaluation for. Default is "all"
    r_warnings : bool
        Whether to suppress R warnings. Default is False
    build : str
        Whether to build the book as "html" or "pdf". Default is "html"
    model : str
        The name of the model. This is only for providing model info and a schematic.
        The only option right now is "ersem". 


    Returns
    -------
    None
    """
    import os
    path_df = []

    fast_plot = False

    # figure out of "book/notebooks" has ipynb files

    empty = True

    # book directory is book, book1, book2, book10 etc.

    # create a new name if one already exists
    i = 0
    book_dir = "book"

    if os.path.exists("book"):
        #get user input to decide if it should be removed
        user_input = input("book directory already exists. This will be emptied and replaced. Do you want to proceed? (y/n): ")
        if user_input.lower() == "y":
            while True:
                files = glob.glob("book/**/**/**", recursive=True)
                # list all files in book, recursively
                for ff in files:
                    if ff.startswith("book"):
                        try:
                            os.remove(ff)
                        except:
                            pass
                files = glob.glob("book/**/**/**", recursive=True)
                # only list files
                files = [x for x in files if os.path.isfile(x)]
                if len(files) == 0:
                    break
        else:
            print("Exiting")
            return None
                # if not os.path.exists(book_dir):
                    # break

    import os

    if variables != "all":
        if isinstance(variables, str):
            variables = [variables]

    if empty:
        from shutil import copyfile

        if not os.path.exists(book_dir):
            os.mkdir(book_dir)
        if not os.path.exists(f"{book_dir}/notebooks"):
            os.mkdir(f"{book_dir}/notebooks")

        # from importlib.resources import files
        import pkg_resources

        data_path = pkg_resources.resource_filename(__name__, "data/mask.py")
        if not os.path.exists(f"{book_dir}/notebooks/mask.py"):
            copyfile(data_path, f"{book_dir}/notebooks/mask.py")

        data_path = pkg_resources.resource_filename(__name__, "data/000_info.ipynb")
        if not os.path.exists(f"{book_dir}/notebooks/000_info.ipynb"):
            copyfile(data_path, f"{book_dir}/notebooks/000_info.ipynb")

        # open this file and replace model_name with model


        # Replace the target string
        if model is not None:
            with open(f"{book_dir}/notebooks/000_info.ipynb", "r") as file:
                filedata = file.read()
            filedata = filedata.replace("model_name", model)

            # Write the file out again
            with open(f"{book_dir}/notebooks/000_info.ipynb", "w") as file:
                file.write(filedata)

        # data_path = files("ecoval.data").joinpath("toc.yml")
        data_path = pkg_resources.resource_filename(__name__, "data/_toc.yml")

        out = f"{book_dir}/" + os.path.basename(data_path)
        copyfile(data_path, out)

        data_path = pkg_resources.resource_filename(__name__, "data/requirements.txt")
        out = f"{book_dir}/" + os.path.basename(data_path)
        copyfile(data_path, out)

        data_path = pkg_resources.resource_filename(__name__, "data/intro.md")
        out = f"{book_dir}/" + os.path.basename(data_path)
        copyfile(data_path, out)

        # copy config

        data_path = pkg_resources.resource_filename(__name__, "data/_config.yml")
        out = f"{book_dir}/" + os.path.basename(data_path)
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

        path_df = []

        # loop through the point matchups and generate notebooks

        point_paths = glob.glob("matched/point/**/**/**/**.csv")
        point_paths = [x for x in point_paths if "paths.csv" not in x]
        point_paths = [x for x in point_paths if "pft" not in x]
        point_paths = [x for x in point_paths if "unit" not in x]
        # loop through the paths
        for pp in point_paths:
            vv = os.path.basename(pp).split("_")[2].replace(".csv", "")
            if True:
                if vv not in variables:
                    continue
            source = os.path.basename(pp).split("_")[0]
            variable = vv
            layer = os.path.basename(pp).split("_")[1].replace(".csv", "")
            if True:
                if vv != "ph":
                    Variable = variable
                else:
                    Variable = "pH"
                if vv == "co2flux":
                    Variable = "Air-sea CO2 fluxes"
                if vv in ["poc", "doc"]:
                    Variable = Variable.upper()
                if vv == "pco2":
                    Variable = "pCO2"
                if vv == "benbio":
                    Variable = "macrobenthos biomass"
                if vv == "pft":
                    Variable = "Plankton Functional Types"
                vv_file = pp
                vv_file_find = pp.replace("../../", "")

                if os.path.exists(vv_file_find):
                    if (
                        len(glob.glob(f"{book_dir}/notebooks/*point_{layer}_{variable}.ipynb"))
                        == 0
                    ):
                        file1 = pkg_resources.resource_filename(
                            __name__, "data/point_template.ipynb"
                        )
                        with open(file1, "r") as file:
                            filedata = file.read()

                        # Replace the target string
                        out = f"{book_dir}/notebooks/{source}_{layer}_{variable}.ipynb"
                        filedata = filedata.replace("point_variable", variable)
                        filedata = filedata.replace("point_layer", layer)
                        filedata = filedata.replace("template_title", Variable)

                        # Write the file out again
                        with open(out, "w") as file:
                            file.write(filedata)
                        variable = vv
                        if variable in ["poc", "doc"]:
                            variable = variable.upper()
                        if variable == "pco2":
                            variable = "pCO2"
                        path_df.append(
                            pd.DataFrame(
                                {
                                    "variable": [variable],
                                    "path": out,
                                }
                            )
                        )


        # Loop through the gridded matchups and generate notebooks
        # identify gridded variables in matched data
        gridded_paths = glob.glob("matched/gridded/**/**/**.nc")

        if len(gridded_paths) > 0:
            for vv in [
                os.path.basename(x).split("_")[1].replace(".nc", "") for x in gridded_paths
            ]:
                source = os.path.basename(glob.glob(f"matched/gridded/**/**/**_{vv}_**.nc")[0]).split("_")[0]

                variable = vv
                if variables != "all":
                    if vv not in variables:
                        continue
                if not os.path.exists(f"{book_dir}/notebooks/{source}_{variable}.ipynb"):
                    if variable == "ph":
                        Variable = "pH"
                    else:
                        Variable = variable
                    if variable == "co2flux":
                        Variable = "air-sea CO2 flux"
                    if variable in ["poc", "doc"]:
                        Variable = Variable.upper()
                    if variable == "benbio":
                        Variable = "macrobenthos biomass"
                    if variable == "pco2":
                        Variable = "pCO2"
                    if variable == "pft":
                        Variable = "PFT"
                    file1 = pkg_resources.resource_filename(
                        __name__, "data/gridded_template.ipynb"
                    )
                    if len(glob.glob(f"{book_dir}/notebooks/*{source}_{variable}.ipynb")) == 0:
                        with open(file1, "r") as file:
                            filedata = file.read()

                        # Replace the target string
                        filedata = filedata.replace("template_variable", variable)
                        filedata = filedata.replace("template_title", Variable)

                        # Write the file out again
                        with open(f"{book_dir}/notebooks/{source}_{variable}.ipynb", "w") as file:
                            file.write(filedata)

                        variable = vv
                        path_df.append(
                            pd.DataFrame(
                                {
                                    "variable": [variable],
                                    "path": [f"{book_dir}/notebooks/{source}_{variable}.ipynb"],
                                }
                            )
                        )

        # figure out if mld needs to be calculated...

        mld_paths = glob.glob("matched/point/**/all/temperature/**temperature**.csv")

        run = True
        if variables != "all":
            if "mld" not in variables:
                run = False
        if len(mld_paths) > 0 and run:

            if not os.path.exists(f"{book_dir}/notebooks/temperature_mld.ipynb"):
                file1 = pkg_resources.resource_filename(__name__, "data/mld_template.ipynb")
                with open(file1, "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("template_variable", "temperature")
                filedata = filedata.replace("template_title", "Mixed layer depth")

                # Write the file out again
                with open(f"{book_dir}/notebooks/temperature_mld.ipynb", "w") as file:
                    file.write(filedata)

                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": ["temperature"],
                            "path": [f"{book_dir}/notebooks/temperature_mld.ipynb"],
                        }
                    )
                )
        

        # copy pft ntoebook over

        #/data/proteus1/scratch/rwi/pygetm/pft/matched/point/nws/surface/pft/cefas_surface_pft.csv
        if os.path.exists("matched/point/nws/surface/pft/cefas_surface_pft.csv"):
            if not os.path.exists(f"{book_dir}/notebooks/surface_pft.ipynb"):
                file1 = pkg_resources.resource_filename(__name__, "data/pft_template.ipynb")
                with open(file1, "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("template_variable", "pft")
                filedata = filedata.replace("template_title", "Primary production")

                # Write the file out again
                with open(f"{book_dir}/notebooks/surface_pft.ipynb", "w") as file:
                    file.write(filedata)

                path_df.append(
                    pd.DataFrame(
                        {
                            "variable": ["pft"],
                            "path": [f"{book_dir}/notebooks/surface_pft.ipynb"],
                        }
                    )
                )



        # loop through all notebooks and replace paths

        i = 0

        # need to start by figuring out whether anything has already been run...

        i = 0

        for ff in [x for x in glob.glob("{book_dir}/notebooks/*.ipynb") if "info" not in x]:
            try:
                i_ff = int(os.path.basename(ff).split("_")[0])
                if i_ff > i:
                    i = i_ff
            except:
                pass

        i_orig = i

        if len(path_df) > 0:
            path_df = pd.concat(path_df)
            path_df = path_df.sort_values("variable").reset_index(drop=True)

        for i in range(len(path_df)):
            file1 = path_df.path.values[i]
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
            ds_regions = nc.open_data(f"{data_dir}/amm7_val_subdomains.nc")
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
            if len(glob.glob(f"{book_dir}/notebooks/*summary.ipynb")) == 0:
                copyfile(file1, f"{book_dir}/notebooks/{i_pad}_summary_shelf.ipynb")
                # change domain_title to "Shelf"
                with open(f"{book_dir}/notebooks/{i_pad}_summary_shelf.ipynb", "r") as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace("domain_title", "Shelf")
                # replace shelf_mask with True
                filedata = filedata.replace("shelf_mask", "True")

                # Write the file out again
                with open(f"{book_dir}/notebooks/{i_pad}_summary_shelf.ipynb", "w") as file:
                    file.write(filedata)

                i += 1

                i_pad = str(i).zfill(3)

        file1 = pkg_resources.resource_filename(__name__, "data/summary.ipynb")
        if len(glob.glob(f"{book_dir}/notebooks/*summary.ipynb")) == 0:
            copyfile(file1, f"{book_dir}/notebooks/{i_pad}_summary.ipynb")
        else:
            if i > (i_orig + 1):
                initial = glob.glob(f"{book_dir}/notebooks/*summary.ipynb")[0]
                copyfile(initial, f"{book_dir}/notebooks/{i_pad}_summary.ipynb")
                i += 1

        # change domain_title to "Full domain"

        with open(f"{book_dir}/notebooks/{i_pad}_summary.ipynb", "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("domain_title", "Full domain")
        # replace shelf_mask with False
        filedata = filedata.replace("shelf_mask", "False")

        # Write the file out again
        with open(f"{book_dir}/notebooks/{i_pad}_summary.ipynb", "w") as file:
            file.write(filedata)

        # pair the notebooks using jupyter text

        os.system(f"jupytext --set-formats ipynb,py:percent {book_dir}/notebooks/*.ipynb")

        # add the chunks
        add_chunks(build)

        # loop through the notebooks and set r warnings options
        for ff in glob.glob(f"{book_dir}/notebooks/*.py"):
            with open(ff, "r") as file:
                filedata = file.read()
        
            # loop through line by line, and rewrite the original file
            lines = filedata.split("\n")
            new_lines = []
            for line in lines:
                if "%%R" in line:
                    new_lines.append(line)
                    new_lines.append("options(warn=-1)")
                else:
                    new_lines.append(line)
            
            # write the new lines to the file
            with open(ff, "w") as file:
                for line in new_lines:
                    file.write(line + "\n")







        # sync the notebooks
        #
        os.system(f"jupytext --sync {book_dir}/notebooks/*.ipynb")

    #        return None

    # return None

    # loop through notebooks and change fast_plot_value to fast_plot

    for ff in glob.glob(f"{book_dir}/notebooks/*.ipynb"):
        with open(ff, "r") as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace("fast_plot_value", str(fast_plot))

        # fix linees using the above
        def fix_r_magic(x):
            if "%%R" in x:
                # x = x.sp
                x =  re.sub(r" -r\s+\d+", "", x)
                x = x.replace("%%R", "%%R -r 120 ") 
            return x

        # lines = filedata.split("\n")
        # new_lines = []
        # for line in lines:
        #     new_lines.append(fix_r_magic(line))

        # filedata = "\n".join(new_lines)

        # Write the file out again
        with open(ff, "w") as file:
            file.write(filedata)

    # fix the toc using the function

    # raise ValueError("here")
    fix_toc(book_dir)

    for ff in glob.glob(f"{book_dir}/notebooks/*.ipynb"):
        ff_clean = ff.replace(".ipynb", ".py")
        if os.path.exists(ff_clean):
            os.remove(ff_clean)

    if build == "html":
        os.system(f"jupyter-book build {book_dir}/")
    else:
        os.system(f"jupyter-book build {book_dir}/ --builder pdflatex")
        # Now the latex file needs to be modified and re-compiled
        # open the latex file and ensure tables have an hline after the header 
        with open(f"{book_dir}/_build/latex/python.tex", "r") as file:
            filedata = file.read()

        # Replace the target string
            filedata = filedata.replace("\\begin{tabular}{ ", "\\begin{tabular}{\n \\hline\n")
        
        # read line by line and modify

        lines = filedata.split("\n")
        new_lines = []
        i = 3
        for line in lines:
            if "\\begin{tabular}" in line:
                new_lines.append(line)
                i = 0
            else:
                if i == 0:
                    # ensure everything is bold in the latex in this line
                    def make_header_bold(x):
                         # header is of format
                         #Variable & Model mean & Observed mean & Model bias & Percentage bias \\
                         # stuff in between & needs to be bold
                         # split the string by "&"
                         y = x.split("&")
                         y = ["\\textbf{" + x.replace("\\", "").strip() + "}" for x in y]
                         # make sure the 2 in R2 is a superscript
                         y = [x.replace("R2", "R$^2$") for x in y]
                         y = " & ".join(y) + " \\\\"
                         # remove the last element

                         return y
                    if "textbf" not in line:
                        new_lines.append(make_header_bold(line))
                    else:
                        new_lines.append(line)
                else:
                    new_lines.append(line)
                if i == 0:
                    new_lines.append("\\hline")
                    i = 1
        # now replace noindent with center in lines with sphinxincludegraphics
        for i in range(len(new_lines)):
            if "sphinxincludegraphics" in new_lines[i]:
                new_lines[i] = new_lines[i].replace("noindent", "center")
        
        # we do not want lines that start with\part{
        new_lines = [x for x in new_lines if not x.startswith("\\part{")]

        filedata = "\n".join(new_lines)

        # Write the file out again
        with open(f"{book_dir}/_build/latex/python.tex", "w") as file:
            file.write(filedata)
        makecmd = os.environ.get("MAKE", "make")
        # remove the pdf file
        if os.path.exists(f"{book_dir}/_build/latex/python.pdf"):
            os.remove(f"{book_dir}/_build/latex/python.pdf")
        output_path = f"{book_dir}/_build/latex"
        output = subprocess.run([makecmd, "all-pdf"], cwd=output_path)
        if output.returncode != 0:
            _error("Error: Failed to build pdf")
            return output.returncode
        print(f"A PDF of your validation can be found at: {output_path} ")







    import os

    stamps = [os.path.basename(x) for x in glob.glob(f"{book_dir}/notebooks/.trackers/*")]
    stamps.append("nctoolkit_rwi_uhosarcenctoolkittmp")

    delete = []
    for x in stamps:
        delete += glob.glob("/tmp/*" + x + "*")

    for ff in delete:
        if os.path.exists(ff):
            if "nctoolkit" in x:
                os.remove(ff)
    # create a symlink  
    out_ff = f"{book_dir}/_build/html/index.html"
    #if os.path.exists("validation.html"):
    #    os.remove("validation.html")
    #os.symlink(out_ff, "validation.html")
    # clean up with the python files that are now useless

    if build == "html":
        webbrowser.open("file://" + os.path.abspath(f"{book_dir}/_build/html/index.html"))
    else:
        out_ff = f"validation_report.pdf"
        # make a copy of the pdf
        shutil.copyfile(f"{book_dir}/_build/latex/python.pdf", out_ff)

        webbrowser.open("file://" + os.path.abspath(f"{book_dir}/_build/latex/python.pdf"))


def rebuild():
    os.system(f"jupyter-book build book/")
    webbrowser.open("file://" + os.path.abspath("book/_build/html/index.html"))


try:
    from importlib.metadata import version as _version
except ImportError:
    from importlib_metadata import version as _version

try:
    __version__ = _version("ecoval")
except Exception:
    __version__ = "999"
