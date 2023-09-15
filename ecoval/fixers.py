

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




