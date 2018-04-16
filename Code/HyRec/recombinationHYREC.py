from .baseconfig import CAMB_Structure, dll_import
from ctypes import c_bool, c_int, c_double

# logical
Do21cm = dll_import(c_double, "recombination", "do21cm")
# Do21cm.value = False

# logical
doTmatTspin = dll_import(c_bool, "recombination", "dotmattspin")
# doTmatTspin.value = False

recombination_saha_z = dll_import(c_double, "recombination", "recombination_saha_z")

recombination_saha_tau = dll_import(c_double, "recombination", "recombination_saha_tau")


# ---Derived Types in recombination.f90

#PDM 2016 for HYREC
class RecombinationParams(CAMB_Structure):
    """
    Hold parameters for the HYREC recombination model.
    """
    _fields_ = [
        ("Pann", c_double), #Annihilation
        ("FineSRat", c_double), #fine-structure constant
        ("EMassRat", c_double), #Electron mass
    ]
