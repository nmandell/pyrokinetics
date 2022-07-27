import numpy as np
from cleverdict import CleverDict
from copy import copy
from typing import Dict, Any, Optional
from ..typing import PathLike
from ..constants import pi, sqrt2, electron_charge
from ..local_species import LocalSpecies
from ..local_geometry import (
    LocalGeometry,
    LocalGeometryMiller,
    default_miller_inputs,
)
from ..numerics import Numerics
from ..templates import gk_templates
from .GKInput import GKInput


class GKInputSTELLA(GKInput):
    """
    Class that can read stella input files, and produce
    Numerics, LocalSpecies, and LocalGeometry objects
    """

    code_name = "stella"
    default_file_name = "input.in"

    pyro_stella_miller = {
        "rho": ["millergeo_parameters", "rhoc"],
        "Rmaj": ["millergeo_parameters", "rmaj"],
        "q": ["millergeo_parameters", "qinp"],
        "kappa": ["millergeo_parameters", "kappa"],
        "shat": ["millergeo_parameters", "shat"],
        "shift": ["millergeo_parameters", "shift"],
        "beta_prime": ["millergeo_parameters", "betaprim"],
    }

    pyro_stella_species = {
        "mass": "mass",
        "z": "z",
        "dens": "dens",
        "temp": "temp",
#        "nu": "vnewk",
        "a_lt": "tprim",
        "a_ln": "fprim",
#        "a_lv": "uprim",
    }

    def read(self, filename: PathLike) -> Dict[str, Any]:
        """
        Reads stella input file into a dictionary
        """
        result = super().read(filename)
        return result

    def read_str(self, input_string: str) -> Dict[str, Any]:
        """
        Reads stella input file given as string
        Uses default read_str, which assumes input_string is a Fortran90 namelist
        """
        result = super().read_str(input_string)
        return result

    def verify(self, filename: PathLike):
        """
        Ensure this file is a valid stella input file, and that it contains sufficient
        info for Pyrokinetics to work with
        """
        # The following keys are not strictly needed for a stella input file,
        # but they are needed by Pyrokinetics
        expected_keys = [
            "knobs",
            "parameters",
            "geo_knobs",
            "millergeo_parameters",
            "zgrid_parameters",
            "species_knobs",
            "kt_grids_knobs",
            "physics_flags",
        ]
        if not self.verify_expected_keys(filename, expected_keys):
            raise ValueError(f"Unable to verify {filename} as stella file")

    def write(self, filename: PathLike, float_format: str = ""):
        super().write(filename, float_format=float_format)

    def is_nonlinear(self) -> bool:
        try:
            is_box = self.data["kt_grids_knobs"]["grid_option"] == "box"
            is_nonlinear = self.data["physics_flags"]["nonlinear"] == True
            return is_box and is_nonlinear
        except KeyError:
            return False

    def add_flags(self, flags) -> None:
        """
        Add extra flags to stella input file
        """
        super().add_flags(flags)

    def get_local_geometry(self) -> LocalGeometry:
        """
        Returns local geometry. Delegates to more specific functions
        """
        stella_eq = self.data["geo_knobs"]["geo_option"]

        local_eq = True
        if stella_eq not in ["miller", "local", "default"]:
            local_eq = False
            raise NotImplementedError(
                f"stella equilibrium option {stella_eq} not implemented"
            )

#        geotype = self.data["theta_grid_parameters"].get("geotype", 0)
#        if geotype != 0:
#            raise NotImplementedError("GS2 Fourier options are not implemented")

        return self.get_local_geometry_miller()

    def get_local_geometry_miller(self) -> LocalGeometryMiller:
        """
        Load Miller object from stella file
        """
        miller_data = default_miller_inputs()

        for pyro_key, (stella_param, stella_key) in self.pyro_stella_miller.items():
            miller_data[stella_key] = self.data[stella_param][stella_key]

        rho = miller_data["rho"]
        kappa = miller_data["kappa"]
        miller_data["delta"] = np.sin(self.data["millergeo_parameters"]["tri"])
        miller_data["s_kappa"] = (
            self.data["millergeo_parameters"]["kapprim"] * rho / kappa
        )
        miller_data["s_delta"] = self.data["millergeo_parameters"]["triprim"] * rho

        # Get beta and beta_prime normalised to R_major(in case R_geo != R_major)
        r_geo = self.data["millergeo_parameters"].get("rgeo", miller_data["rmaj"])

        beta = self.data["parameters"]["beta"] * (miller_data["rmaj"] / r_geo) ** 2
        miller_data["beta_prime"] *= -0.5*(miller_data["rmaj"] / r_geo) ** 2

        # Assume pref*8pi*1e-7 = 1.0
        # FIXME Is this assumption general enough? Can't we get pref from local_species?
        # FIXME B0 = None can cause problems when writing
        miller_data["B0"] = np.sqrt(1.0 / beta) if beta != 0.0 else None

        # must construct using from_gk_data as we cannot determine bunit_over_b0 here
        return LocalGeometryMiller.from_gk_data(miller_data)

    def get_local_species(self):
        """
        Load LocalSpecies object from stella file
        """
        # Dictionary of local species parameters
        local_species = LocalSpecies()

        ion_count = 0

        # Load each species into a dictionary
        for i_sp in range(self.data["species_knobs"]["nspec"]):

            species_data = CleverDict()

            stella_key = f"species_parameters_{i_sp + 1}"

            stella_data = self.data[stella_key]

            for pyro_key, stella_key in self.pyro_stella_species.items():
                species_data[pyro_key] = stella_data[stella_key]

            species_data.vel = 0.0
            species_data.a_lv = 0.0

            if species_data.z == -1:
                name = "electron"
            else:
                ion_count += 1
                name = f"ion{ion_count}"

            dens = stella_data["dens"]
            temp = stella_data["temp"]
            charge = stella_data["z"]
            mass = stella_data["mass"]
            # Account for sqrt(2) in vth
            species_data.nu = stella_data["vnew_ref"] * sqrt2 * dens * charge**4 / sqrt(mass) / temp**1.5
                
            species_data.name = name

            # Add individual species data to dictionary of species
            local_species.add_species(name=name, species_data=species_data)

        local_species.normalise()
        return local_species

    def get_numerics(self) -> Numerics:
        """Gather numerical info (grid spacing, time steps, etc)"""

        numerics_data = {}

        # Set no. of fields
        numerics_data["phi"] = self.data["knobs"].get("fphi", 0.0) > 0.0
        numerics_data["apar"] = self.data["knobs"].get("fapar", 0.0) > 0.0
        numerics_data["bpar"] = self.data["knobs"].get("fbpar", 0.0) > 0.0

        # Set time stepping
        delta_time = self.data["knobs"].get("delt", 0.005) / sqrt2
        numerics_data["delta_time"] = delta_time
        numerics_data["max_time"] = self.data["knobs"].get("nstep", 50000) * delta_time

        # Fourier space grid
        # Linear simulation
        if self.is_linear():
            numerics_data["nky"] = 1
            numerics_data["nkx"] = 1
            numerics_data["ky"] = self.data["kt_grids_range_parameters"]["aky_min"] / sqrt2
            numerics_data["kx"] = 0.0
            numerics_data["theta0"] = self.data["kt_grids_range_parameters"].get(
                "theta0_min", 0.0
            )
            numerics_data["nonlinear"] = False
        # Nonlinear/multiple modes in box
        # kt_grids_knobs.grid_option == "box"
        else:
            box = self.data["kt_grids_box_parameters"]
            keys = box.keys()

            # Set up ky grid
            if "ny" in keys:
                numerics_data["nky"] = int((box["ny"] - 1) / 3 + 1)
            elif "n0" in keys:
                numerics_data["nky"] = box["n0"]
            elif "nky" in keys:
                numerics_data["nky"] = box["naky"]
            else:
                raise RuntimeError(f"ky grid details not found in {keys}")

            if "y0" in keys:
                if box["y0"] < 0.0:
                    numerics_data["ky"] = -box["y0"] / sqrt2
                else:
                    numerics_data["ky"] = 1 / box["y0"] / sqrt2
            else:
                raise RuntimeError(f"Min ky details not found in {keys}")

            if "nx" in keys:
                numerics_data["nkx"] = int((2 * box["nx"] - 1) / 3 + 1)
            elif "ntheta0" in keys():
                numerics_data["nkx"] = int((2 * box["ntheta0"] - 1) / 3 + 1)
            else:
                raise RuntimeError("kx grid details not found in {keys}")

            shat_params = self.pyro_stella_miller["shat"]
            shat = self.data[shat_params[0]][shat_params[1]]
            if abs(shat) > 1e-6:
                numerics_data["kx"] = (
                    numerics_data["ky"] * shat * 2 * pi / box["jtwist"]
                )
            else:
                numerics_data["kx"] = 2 * pi / (box["x0"] * sqrt2)

            try:
                numerics_data["nonlinear"] = (
                    self.data["physics_flags"]["nonlinear"] == True
                )
            except KeyError:
                numerics_data["nonlinear"] = False

        # Theta grid
        numerics_data["ntheta"] = self.data["zgrid_parameters"]["nzed"]
        numerics_data["nperiod"] = self.data["zgrid_parameters"]["nperiod"]

        # Parallel velocity grid (nvpa = 2 * nvgrid)
        try:
            numerics_data["nvpa"] = (
                self.data["vpamu_grids_parameters"]["nvgrid"] * 2
            )
        except KeyError:
            numerics_data["nvpa"] = self.data["vpamu_grids_parameters"]["nvgrid"]

        numerics_data["nmu"] = self.data["vpamu_grids_parameters"]["nmu"]

        return Numerics(numerics_data)

    def set(
        self,
        local_geometry: LocalGeometry,
        local_species: LocalSpecies,
        numerics: Numerics,
        template_file: Optional[PathLike] = None,
        **kwargs,
    ):
        """
        Set self.data using LocalGeometry, LocalSpecies, and Numerics.
        These may be obtained via another GKInput file, or from Equilibrium/Kinetics
        objects.
        """
        # If self.data is not already populated, fill in defaults from a given
        # template file. If this is not provided by the user, fall back to the
        # default.
        if self.data is None:
            if template_file is None:
                template_file = gk_templates["stella"]
            self.read(template_file)

        # Set Miller Geometry bits
        if not isinstance(local_geometry, LocalGeometryMiller):
            raise NotImplementedError(
                f"LocalGeometry type {local_geometry.__class__.__name__} for stella not supported yet"
            )

        # Ensure Miller settings
        self.data["geo_knobs"]["geo_option"] = 'miller'

        # Assign Miller values to input file
        for key, val in self.pyro_stella_miller.items():
            self.data[val[0]][val[1]] = local_geometry[key]

        self.data["millergeo_parameters"]["kapprim"] = (
            local_geometry.s_kappa * local_geometry.kappa / local_geometry.rho
        )
        self.data["miilergeo_parameters"]["tri"] = np.arcsin(local_geometry.delta)
        self.data["millergeo_parameters"]["triprim"] = (
            local_geometry["s_delta"] / local_geometry.rho
        )
        self.data["millergeo_parameters"]["rgeo"] = local_geometry.Rmaj

        # Set local species bits
        self.data["species_knobs"]["nspec"] = local_species.nspec
        for iSp, name in enumerate(local_species.names):

            # add new outer params for each species
            species_key = f"species_parameters_{iSp + 1}"

            if name == "electron":
                self.data[species_key]["type"] = "electron"
            else:
                try:
                    self.data[species_key]["type"] = "ion"
                except KeyError:
                    self.data[species_key] = copy(self.data["species_parameters_1"])
                    self.data[species_key]["type"] = "ion"

            for key, val in self.pyro_stella_species.items():
                self.data[species_key][val] = local_species[name][key]

            # set the reference collision frequency (electron-electron collisions)
            # NEEDS SORTING OUT
            self.data["parameters"]["vnew_ref"] = local_species[name]["nu"] / sqrt2
            
        # If species are defined calculate beta
        if local_species.nref is not None:

            pref = local_species.nref * local_species.tref * electron_charge
            # FIXME local_geometry.B0 may be set to None
            bref = local_geometry.B0

            beta = pref / bref**2 * 8 * pi * 1e-7

        # Calculate from reference  at centre of flux surface
        else:
            if local_geometry.B0 is not None:
                beta = 1 / local_geometry.B0**2
            else:
                beta = 0.0

        self.data["parameters"]["beta"] = beta

        # Set numerics bits
        # Set no. of fields
        self.data["knobs"]["fphi"] = 1.0 if numerics.phi else 0.0
        self.data["knobs"]["fapar"] = 1.0 if numerics.apar else 0.0
        self.data["knobs"]["fbpar"] = 1.0 if numerics.bpar else 0.0

        # Set time stepping
        self.data["knobs"]["delt"] = numerics.delta_time * sqrt2
        self.data["knobs"]["nstep"] = int(numerics.max_time / numerics.delta_time)

        if numerics.nky == 1:
            self.data["kt_grids_knobs"]["grid_option"] = "range"

            if "kt_grids_range_parameters" not in self.data.keys():
                self.data["kt_grids_range_parameters"] = {}

            self.data["kt_grids_range_parameters"]["aky_min"] = numerics.ky * sqrt2
            self.data["kt_grids_range_parameters"]["aky_max"] = numerics.ky * sqrt2
            self.data["kt_grids_range_parameters"]["theta0_min"] = numerics.theta0
            self.data["kt_grids_range_parameters"]["theta0_max"] = numerics.theta0
            self.data["zgrid_parameters"]["nperiod"] = numerics.nperiod

        else:
            self.data["kt_grids_knobs"]["grid_option"] = "box"

            if "kt_grids_box_parameters" not in self.data.keys():
                self.data["kt_grids_box_parameters"] = {}

            self.data["kt_grids_box_parameters"]["nx"] = int(
                ((numerics.nkx - 1) * 3 / 2) + 1
            )
            self.data["kt_grids_box_parameters"]["ny"] = int(
                ((numerics.nky - 1) * 3) + 1
            )

            self.data["kt_grids_box_parameters"]["y0"] = -numerics.ky * sqrt2

            # Currently forces NL sims to have nperiod = 1
            self.data["theta_grid_parameters"]["nperiod"] = 1

            shat = local_geometry.shat
            if abs(shat) < 1e-6:
                self.data["kt_grids_box_parameters"]["x0"] = (
                    2 * pi / numerics.kx / sqrt2
                )
            else:
                self.data["kt_grids_box_parameters"]["jtwist"] = int(
                    (numerics.ky * shat * 2 * pi / numerics.kx) + 0.1
                )

        self.data["zgrid_parameters"]["nzed"] = numerics.ntheta

        self.data["vpamu_grids_parameters"]["nvgrid"] = numerics.nvpa // 2
        self.data["vpamu_grids_parameters"]["nmu"] = numerics.nmu

        if numerics.nonlinear:
            if "physics_flags" not in self.data.keys():
                self.data["physics_flags"] = {}

            self.data["physics_flags"]["nonlinear"] = True
        else:
            try:
                self.data["physics_flags"]["nonlinear"] = False
            except KeyError:
                pass
