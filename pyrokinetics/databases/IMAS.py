from omas.omas_utils import imas_structure
from omas import latest_imas_version
from cleverdict import CleverDict
from typing import Optional
from ..pyro import Pyro
import json
from pathlib import Path
import operator
from functools import reduce


class GKODS:
    """
    Class representing the OMAS gyrokinetics structure

    """
    def __init__(self, 
                 pyro: Optional[Pyro] = None,
                 ods: Optional[dict] = None,
             ):
        
        self.ods = imas_structure(latest_imas_version, "gyrokinetics")
        
        self.json_map_file =  Path(__file__).parent / "imas_to_pyro.json"

        with open(self.json_map_file) as f:
            self.pyro_to_ods_dict = json.load(f)

        if pyro:
            # FIXME store pyro or pass as argument?
            #self.pyro = pyro
            self.map_pyro_to_ods(pyro)

    def map_pyro_to_ods(self, 
                        pyro: Pyro,
                    ):
        for key_1, val_1 in self.pyro_to_ods_dict.items():
            if key_1 in ['ids_properties']:
                continue

            for key_2, val_2 in val_1.items():
                if key_2 in ['library', ':']:
                    continue

                if val_2:
                    self.ods[key_1][key_2] = reduce(getattr, val_2, pyro)


    def map_ods_to_pyro(self,
                    ):

        raise NotImplementedError


    def _pyro_to_ods_key(self,):
        
        data_dict = {"code" : "test"}

        return data_dict
        


def get_from_dict(data_dict, map_list):
    """
    Gets item in dict given location as a list of string
    """
    return reduce(operator.getitem, map_list, data_dict)


def set_in_dict(data_dict, map_list, value):
    """
    Sets item in dict given location as a list of string
    """
    get_from_dict(data_dict, map_list[:-1])[map_list[-1]] = copy.deepcopy(value)
