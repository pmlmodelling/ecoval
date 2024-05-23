
import ecoval.ices as ices
import ecoval.noaa as noaa
from ecoval.ices import match_ices
from ecoval.matchall import matchup
from ecoval.fixers import tidy_name
from ecoval.regionals import global_regionals
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

def add_chunks(build = "html"):

    nws = False
    if len(glob.glob("matched/gridded/nws/**/*.nc")) > 0:
        nws = True

    paths = glob.glob("book/notebooks/*.py")
    paths += glob.glob("book/compare/notebooks/*.py")

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

