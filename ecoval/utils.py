

import os
import glob





def find_config(level = 0):
    # first look in the working directory
    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(ff):
            return ff

    if level == -2:
        for ff in ["../../.ecovalrc", "../../ecovalrc"]:
            if os.path.exists(ff):
                return ff

    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(ff):
            return ff

    # now look in the home directory....
    from os.path import expanduser

    home = expanduser("~")
    for ff in [".ecovalrc", "ecovalrc"]:
        if os.path.exists(home + "/" + ff):
            return home + "/" + ff

    return None

def get_datadir(level = 0):

    data_dir = "/data/proteus1/scratch/rwi/evaldata/data/"
    config_file = find_config(level = level)

    if config_file is not None:
        # valid_keys = ["thread_safe", "lazy", "cores", "precision", "temp_dir"]

        file1 = open(config_file, "r")
        Lines = file1.readlines()

        count = 0
        # Strips the newline character
        for line in Lines:
            text = line.replace(" ", "").strip()
            if text.count(":") != 1:
                if len(text) > 0:
                    raise ValueError(f"Line in {config_file} is invalid: {line}")

        for line in Lines:
            text = line.replace(" ", "").strip()
            if len(text) > 0:
                terms = text.split(":")
                key = terms[0]
                value = None
                # print(key)
            data_path = terms[1].replace(" ", "")

            if os.path.exists(data_path):
                data_dir = data_path
            else:
                raise ValueError(f"{data_path} does not exist")
    return data_dir
    

def extension_of_directory(starting_directory, exclude = []):
    n_remove = len(starting_directory)

    n_max = 0
    for root, directories, files in os.walk(starting_directory):
        r_n = (len(root[n_remove:].split("/")))
        paths = root[n_remove+1:].split("/")
        paths = [x for x in paths if len(x) > 0]
        r_n = len(paths)
        if r_n > n_max:
            n_max = r_n 
        checks =  glob.glob(root +   "/*.nc")
        # remove anything in exclude from checks
        for exc in exclude:
            checks = [x for x in checks if exc not in x]
        checks = [x for x in checks if "part" not in os.path.basename(x)]
        checks = [x for x in checks if "restart" not in os.path.basename(x)]
        if len(checks) > 0:
            n_max = r_n
            break
    ## repeat "**" n_max
    return "/" + "/".join(["**" for i in range(n_max)]) + "/"