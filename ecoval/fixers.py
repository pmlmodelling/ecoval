import warnings

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
    if "0 as the fill value" in x:
        return True
    if "found more than one time variabl" in x:
        return True
    if "coordinates variable time" in x and "be assigned" in x:
        return True
    return False

def tidy_name(x):
    """
    A function to create a better name for tables etc.

    """

    if "flux" in x.lower():
        return "Air-sea CO2 flux"
    if "pco2" in x.lower():
        return "pCO2"
    if "sst" in x.lower():
        return "SST"
    if "nitrate" in x.lower():
        return "Nitrate"
    if x.lower() == "ph":
        return "pH"
    
    return x.title()

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




