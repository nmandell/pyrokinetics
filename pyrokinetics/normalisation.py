"""Proxy objects for pint units

The problem we're trying to solve is that each unique simulation run
(i.e. `Pyro` object) will have its own set of reference values to
normalise that run to. Separately, each code has its own convention
for normalising its quantities, including pyrokinetics, as well as
other external conventions, such as GKDB/IMAS. When we have multiple
simulations, we (potentially) need the full Cartesian product of
``(simulation reference values) x (normalisation conventions)`` in
order to, for example, get all of the simulations into a single,
comparable normalisation.

To do this, we use the `pint` library, which allows us to create
unique units for each combination of reference value and normalisation
convention for a given simulation. We then wrap this up in a set of
proxy objects which gives us nice names for these units.


Normalisation for single simulation =>
Conventions =>
Unique units

"""


from dataclasses import dataclass
from typing import Optional, Dict
import numpy as np

import pint

from pyrokinetics.constants import electron_charge, mu0
from pyrokinetics.kinetics import Kinetics
from pyrokinetics.local_geometry import LocalGeometry


@dataclass
class Convention:
    """A description of a normalisation convention, including what
    species reference values use, whether the velocity includes a
    ``sqrt(2)`` factor, what length scales are normalised to, and so
    on.

    TODO: Do we need to specifiy "kind" of deuterium mass? Actual mass vs 2*m_p?

    Attributes
    ----------
    tref_species:
        The species to normalise temperatures to
    nref_species:
        The species to normalise densities to
    mref_species:
        The species to normalise masses to
    vref_multiplier:
        Velocity multiplier
    lref_type:
        What to normalise length scales to
    bref_type:
        Magnetic field normalisation. Must be either ``B0`` or ``Bunit``

    """

    name: str
    tref_species: str = "electron"
    nref_species: str = "electron"
    mref_species: str = "deuterium"
    vref_multiplier: float = 1.0
    lref_type: str = "minor_radius"
    bref_type: str = "B0"


NORMALISATION_CONVENTIONS = {
    "pyrokinetics": Convention("pyrokinetics"),
    "cgyro": Convention("cgyro", bref_type="Bunit"),
    "gs2": Convention("gs2", vref_multiplier=np.sqrt(2)),
    "gene": Convention("gene", lref_type="major_radius"),
    "gkdb": Convention("gkdb", vref_multiplier=np.sqrt(2)),
    "imas": Convention("imas", vref_multiplier=np.sqrt(2)),
}
"""Particular normalisation conventions"""


def _create_unit_registry() -> pint.UnitRegistry:
    """Create a default pint.UnitRegistry with some common features we need"""
    ureg = pint.UnitRegistry()
    ureg.enable_contexts("boltzmann")
    ureg.define("deuterium_mass = 3.3435837724e-27 kg")

    return ureg


ureg = _create_unit_registry()
"""Default unit registry"""


class SimulationNormalisation:
    """Holds the normalisations for a given simulation for all the
    known conventions.

    Has a current convention which sets what the short names refer to.

    Examples
    --------

        >>> norm = SimulationNormalisation("run001")
        >>> norm.set_lref_bref(local_geometry)
        >>> norm.set_kinetic_references(kinetics)
        >>> norm.lref      # Current convention's lref
        1 <Unit('lref_pyrokinetics_run001')>
        >>> norm.gs2.lref  # Specific convention's lref
        1 <Unit('lref_gs2_run001')>
        # Change the current default convention
        >>> norm.default_convention = "gkdb"
        >>> norm.gkdb.lref
        1 <Unit('lref_gkdb_run001')>

    """

    def __init__(
        self,
        name: str,
        default_convention: str = "pyrokinetics",
        registry: pint.UnitRegistry = ureg,
        geometry: Optional[LocalGeometry] = None,
        kinetics: Optional[Kinetics] = None,
        psi_n: Optional[float] = None,
    ):
        self.units = ureg
        self.name = name

        self._conventions: Dict[str, ConventionNormalisation] = {
            name: ConventionNormalisation(self.name, convention, self.units)
            for name, convention in NORMALISATION_CONVENTIONS.items()
        }

        self.default_convention = default_convention

        if geometry:
            self.set_lref_bref(geometry)
        if kinetics:
            self.set_kinetic_references(kinetics, psi_n)

    def __getattr__(self, item):
        try:
            return self._conventions[item]
        except KeyError:
            raise AttributeError(name=item, obj=self)

    @property
    def default_convention(self):
        """Change the current convention that the short names refer to"""
        return self._current_convention

    @default_convention.setter
    def default_convention(self, convention):
        self._current_convention = self._conventions[convention]
        self._update_references()

    def _update_references(self):
        """Update all the short names to the current convention's
        actual units"""

        # Note that this relies on private details of the unit registry
        for key in list(self.units._cache.root_units.keys()):
            if self.name in key:
                del self.units._cache.root_units[key]

        self.bref = self._current_convention.bref
        self.lref = self._current_convention.lref
        self.mref = self._current_convention.mref
        self.nref = self._current_convention.nref
        self.qref = self._current_convention.qref
        self.tref = self._current_convention.tref
        self.vref = self._current_convention.vref
        self.beta = self._current_convention.beta
        self.rhoref = self._current_convention.rhoref

    def set_lref_bref(self, local_geometry: LocalGeometry):
        """Set the length and magnetic field reference values for all
        the conventions from the local geometry"""
        for convention in self._conventions.values():
            convention.set_lref_bref(local_geometry)
        self._update_references()

    def set_kinetic_references(self, kinetics: Kinetics, psi_n: float):
        """Set the temperature, density, and mass reference values for
        all the conventions"""

        for convention in self._conventions.values():
            convention.set_kinetic_references(kinetics, psi_n)
        self._update_references()


class ConventionNormalisation:
    """A concrete set of reference values/normalisations.

    You should call `ConventionNormalistion.set_lref_bref` and then
    `ConventionNormalistion.set_kinetic_references` (in that order)
    before attempting to use most of these units

    Parameters
    ----------
    run_name:
        Name of the specific simulation run
    convention:
        Object describing how particular reference values should be set
    registry:
        The pint registry to add these units to
    definitions:
        Dictionary of definitions for each reference value. If not
        given, the default set will be used

    """

    REF_DEFS = {
        "deuterium_mass": {"def": "3.3435837724e-27 kg"},
        "bref": {"def": "nan tesla", "base": "tesla"},
        "lref": {"def": "nan metres", "base": "meter"},
        "mref": {"def": "deuterium_mass", "base": "gram"},
        "nref": {"def": "nan m**-3"},
        "qref": {"def": "elementary_charge"},
        "tref": {"def": "nan eV", "base": "kelvin"},
        "vref": {"def": "(tref / mref)**(0.5)"},
        "beta": {"def": "2 * mu0 * nref * tref / bref**2"},
        "rhoref": {"def": "mref * vref / qref / bref"},
    }

    def __init__(
        self,
        run_name: str,
        convention: Convention,
        registry: pint.UnitRegistry,
        definitions: Optional[Dict[str, Dict[str, str]]] = None,
    ):
        self.convention = convention
        self.name = f"{self.convention.name}_{run_name}"

        self._registry = registry
        self._system = registry.get_system(self.name)

        self.definitions = definitions or self.REF_DEFS

        for unit, definition in self.definitions.items():
            convention_unit = f"{unit}_{self.name}"

            unit_def = definition["def"]
            for unit_name in list(self.REF_DEFS.keys()):
                unit_def = unit_def.replace(unit_name, f"{unit_name}_{self.name}")

            if unit == "vref":
                unit_def = f"{self.convention.vref_multiplier} * {unit_def}"

            self._registry.define(f"{convention_unit} = {unit_def}")

            if "base" in definition:
                self._system.base_units[definition["base"]] = {convention_unit: 1.0}

        self._system.base_units["second"] = {
            f"lref_{self.name}": 1.0,
            f"vref_{self.name}": -1.0,
        }

        # getattr rather than []-indexing as latter returns a quantity
        # rather than a unit (??)
        self.bref = getattr(self._registry, f"bref_{self.name}")
        self.lref = getattr(self._registry, f"lref_{self.name}")
        self.mref = getattr(self._registry, f"mref_{self.name}")
        self.nref = getattr(self._registry, f"nref_{self.name}")
        self.qref = getattr(self._registry, f"qref_{self.name}")
        self.tref = getattr(self._registry, f"tref_{self.name}")
        self.vref = getattr(self._registry, f"vref_{self.name}")
        self.beta = getattr(self._registry, f"beta_{self.name}")
        self.rhoref = getattr(self._registry, f"rhoref_{self.name}")

    def set_lref_bref(self, local_geometry: LocalGeometry):
        BREF_TYPES = {
            "B0": local_geometry.B0,
            "Bunit": local_geometry.B0 * local_geometry.bunit_over_b0,
        }
        LREF_TYPES = {
            "minor_radius": local_geometry.a_minor,
            "major_radius": local_geometry.Rmaj,
        }

        bref_type = self.convention.bref_type
        if bref_type not in BREF_TYPES:
            raise ValueError(
                f"Unrecognised bref_type: got '{bref_type}', expected one of {list(BREF_TYPES.keys())}"
            )

        lref_type = self.convention.lref_type
        if lref_type not in LREF_TYPES:
            raise ValueError(
                f"Unrecognised lref_type: got '{lref_type}', expected one of {list(LREF_TYPES.keys())}"
            )

        bref = BREF_TYPES[bref_type]
        self._registry.define(f"bref_{self.name} = {bref} tesla")

        lref = LREF_TYPES[lref_type]
        self._registry.define(f"lref_{self.name} = {lref} metres")

    def set_kinetic_references(self, kinetics: Kinetics, psi_n: float):
        tref = kinetics.species_data[self.convention.tref_species].get_temp(psi_n)
        nref = kinetics.species_data[self.convention.nref_species].get_dens(psi_n)
        mref = kinetics.species_data[self.convention.mref_species].get_mass()

        self._registry.define(f"tref_{self.name} = {tref} eV")
        self._registry.define(f"nref_{self.name} = {nref} m**-3")
        self._registry.define(f"mref_{self.name} = {mref} kg")


class Normalisation:
    """A concrete set of normalisation parameters following a given convention

    Attributes
    ----------
    tref:
        Reference temperature
    nref:
        Reference density
    mref:
        Reference mass
    vref:
        Reference velocity
    lref:
        Reference length scale
    bref:
        Reference magnetic field
    """

    def __init__(
        self,
        convention: Optional[int] = "pyrokinetics",
        tref: Optional[float] = None,
        nref: Optional[float] = None,
        mref: Optional[float] = None,
        vref: Optional[float] = None,
        lref: Optional[float] = None,
        bref: Optional[float] = None,
        beta: Optional[float] = None,
        rhoref: Optional[float] = None,
    ):

        self.nocos = self.choose_convention(convention)
        self.tref = tref
        self.nref = nref
        self.mref = mref
        self.vref = vref
        self.lref = lref
        self.bref = bref
        self.beta = self._calculate_beta()
        self.rhoref = self._calculate_rhoref()

    def __repr__(self):
        return (
            f"Normalisation(nocos='{self.nocos.name}', "
            f"tref={self.tref}, "
            f"nref={self.nref}, "
            f"mref={self.mref}, "
            f"vref={self.vref}, "
            f"lref={self.lref}, "
            f"bref={self.bref}, "
            f"beta={self.beta}, "
            f"rhoref={self.rhoref}"
            ")"
        )

    @staticmethod
    def choose_convention(convention: str = "pyrokinetics"):
        """Set normalisation convention"""
        try:
            return NORMALISATION_CONVENTIONS[convention]
        except KeyError:
            raise NotImplementedError(f"NOCOS value {convention} not yet supported")

    def _calculate_beta(self):
        """Return beta from normalised value"""

        if self.bref is None:
            return None

        if self.nref is None:
            return 1.0 / self.bref**2

        return self.nref * self.tref * electron_charge / (self.bref**2 / (2 * mu0))

    def _calculate_rhoref(self):
        """Return reference Larmor radius"""

        if self.vref is None or self.bref is None:
            return None

        return self.mref * self.vref / electron_charge / self.bref

    @classmethod
    def from_kinetics(
        cls,
        kinetics: Kinetics,
        psi_n: float,
        convention: str = "pyrokinetics",
        lref: Optional[float] = None,
        bref: Optional[float] = None,
    ):
        """Create a `Normalisation` using local normalising species data from kinetics object"""

        nocos = cls.choose_convention(convention)

        tref = kinetics.species_data[nocos.tref_species].get_temp(psi_n)
        nref = kinetics.species_data[nocos.nref_species].get_dens(psi_n)
        mref = kinetics.species_data[nocos.mref_species].get_mass()
        vref = np.sqrt(electron_charge * tref / mref) * nocos.vref_multiplier

        return cls(
            convention, tref=tref, nref=nref, mref=mref, vref=vref, lref=lref, bref=bref
        )

    @classmethod
    def from_local_geometry(
        cls, local_geometry: LocalGeometry, convention: str = "pyrokinetics", **kwargs
    ):
        """Create a `Normalisation` using local normalising field from `LocalGeometry` Object.

        This really only sets `bref`, and you'll likely need to pass that into `from_kinetics`
        """

        nocos = cls.choose_convention(convention)

        if nocos.bref_type == "B0":
            bref = local_geometry.B0
        elif nocos.bref_type == "Bunit":
            bref = local_geometry.B0 * local_geometry.bunit_over_b0
        else:
            raise ValueError(f"bref_type : {nocos.bref_type} is not recognised")

        return cls(convention, bref=bref, **kwargs)