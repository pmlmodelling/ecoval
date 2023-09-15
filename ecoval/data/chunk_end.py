# delete all objects starting with ds
for name in dir():
    if name.startswith("ds"):
        del globals()[name]
for name in dir():
    if name.endswith("mask"):
        del globals()[name]
nc.cleanup()
for ff in nc.session.get_safe():
    nc.session.remove_safe(ff)
nc.cleanup()