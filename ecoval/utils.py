

import os




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
    