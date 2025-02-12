
from ecoval.matchall import matchup
import pandas as pd
import glob
import os

import nctoolkit as nc
import webbrowser

import pkg_resources

import re
def is_chunk(x):
    x = x.replace("\n", "" )
    try:     
        return re.findall(r'chunk_.[a-z]*', x)[0] == x
    except:
        return False

def add_chunks(build = "html", cdf = False, dir = None):

    nws = False
    if len(glob.glob("matched/gridded/nws/**/*.nc")) > 0:
        nws = True

    paths = glob.glob(f"book_{build}/notebooks/*.py")
    paths += glob.glob("book/compare/notebooks/*.py")
    if dir is not None:
        paths += glob.glob(f"{dir}/*.py")

    for path in paths:
        # read file line by line
        # Using readlines()
        file1 = open(path, 'r')
        Lines = file1.readlines()

        count = 0
        # Strips the newline character
        # generate new lines
        new_lines = []
        for line in Lines:
            if is_chunk(line):
                # get the file names
                chunk_file = line.replace("\n", "") + ".py"
                # check if globals is in chunk file
                if cdf == False:
                    if "chunk_cdf" in chunk_file:
                        continue
                if "globals" in chunk_file:
                #     # now, we need to figure out if it's a global grid
                    try:
                        file_paths = nc.create_ensemble("matched/gridded/")
                        if len(file_paths) == 0:
                            chunk_file = "chunk_empty.py" 
                        else:
                            file_path = file_paths[0]
                            ds = nc.open_data(file_path)
                            ds = ds.to_xarray()
                            lon_name = [x for x in ds.coords if "lon" in x][0]
                            lat_name = [x for x in ds.coords if "lat" in x][0]
                            # now, get the lon and lat min/max
                            lon_min = float(ds[lon_name].min())
                            lon_max = float(ds[lon_name].max())
                            lat_min = float(ds[lat_name].min())
                            lat_max = float(ds[lat_name].max())
                            lon_range = lon_max - lon_min
                            lat_range = lat_max - lat_min
                            if lon_range > 340 and lat_range > 160:
                                chunk_file = "chunk_globals.py"
                            else:
                                chunk_file = "chunk_empty.py"
                            raise ValueError(chunk_file)
                    except:
                        chunk_file = chunk_file


                if not nws:
                    if "cdf" in chunk_file:
                        continue
                data_path = pkg_resources.resource_filename(__name__, f"data/{chunk_file}")

                # read the chunk file in line by line

                # Using readlines()
                file2 = open(data_path, 'r')
                chunk_lines = file2.readlines()

                # add the chunk lines to the new lines
                for chunk_line in chunk_lines:

                    if build == "pdf":
                        new_lines.append(chunk_line.replace("book_build", "pdf"))
                    else:
                        new_lines.append(chunk_line)

                # close the chunk file
                file2.close()

                count += 1
            else:
                new_lines.append(line)

        # close the original file
        file1.close()

        # write the new lines to the file
        file1 = open(path, 'w')
        file1.writelines(new_lines)
        file1.close()

        

