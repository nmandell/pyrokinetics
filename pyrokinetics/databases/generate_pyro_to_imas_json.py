from pyrokinetics import Pyro
from pyrokinetics.databases.IMAS import GKODS
from cleverdict import CleverDict
import json
from operator import getitem
from functools import reduce


pyro = Pyro(gk_file='input.cgyro')

imas = GKODS(pyro)

ods = imas.ods

mapping_dict = CleverDict({})

# Settable parameters
# Code
ods["code"]["name"] = ["gk_code"]


# Flux surface
ods["flux_surface"]["r_minor_norm"] = ["local_geometry", "rho"]

with open("imas_to_pyro.json", 'w+') as f:
    json.dump(ods, f)


