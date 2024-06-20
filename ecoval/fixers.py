import warnings


# function to convert list to string with , separator, with an "and" at the end
def list_to_string(lst):
    if len(lst) == 1:
        return lst[0]
    elif len(lst) == 2:
        return " and ".join(lst)
    else:
        return ", ".join(lst[:-1]) + ", and " + lst[-1]


def ignore_warning(x):
    """
    Parameters
    -------------
    x: str
        Warning message

        Returns
        -------------
        True if the warning should be ignored
        False if the warning should not be ignored
    """
    if "Adding a time series with the same number of time steps" in x:
        return True
    # did not have valid years
    if "did not have valid years" in x:
        return True
    if "There is only file in the dataset. No need to merge" in x:
        return True
    if "time bounds unsupporte" in x:
        return True
    if "deflate" in x:
        return True
    if "None of the points are contained" in x:
        return True
    if "0 as the fill value" in x:
        return True
    if "found more than one time variabl" in x:
        return True
    if "coordinates variable time" in x and "be assigned" in x:
        return True
    return False


def tidy_name_1(x, lower=False):
    """
    A function to create a better name for tables etc.

    """

    if "flux" in x.lower():
        if lower:
            return "air-sea CO2 flux"
        else:
            return "Air-sea CO2 flux"

    if "pco2" in x.lower():
        return "pCO2"
    if "sst" in x.lower():
        return "SST"
    if "nitrate" in x.lower():
        if lower:
            return "nitrate"
        else:
            return "Nitrate"
    if "poc" in x.lower():
        return "POC"
    if x.lower() == "doc":
        return "DOC"
    if x.lower() == "ph":
        return "pH"
    if lower:
        return x.lower()
    else:
        return x.title()


# vectorize the function

def tidy_name(x, lower=False):
    if isinstance(x, str):
        x = [x]
    out = [tidy_name_1(xx, lower=lower) for xx in x]
    # sort
    out = sorted(out)
    out = list_to_string(out)
    return out


def tidy_warnings(w):
    # A function to tidy up the warnings
    out_warnings = []
    for x in w:
        x_message = str(x.message)
        bad = True
        bad = ignore_warning(x_message) is False
        if bad:
            out_warnings.append(x_message)
    for ww in out_warnings:
        warnings.warn(ww)
