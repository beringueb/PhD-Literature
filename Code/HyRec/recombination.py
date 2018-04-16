from .baseconfig import CAMB_Structure, dll_import
from ctypes import c_bool, c_int, c_double

# logical
#Do21cm = dll_import(c_double, "recombination", "do21cm")
#Do21cm.value = False

# logical
#doTmatTspin = dll_import(c_bool, "recombination", "dotmattspin")
#doTmatTspin.value = False

#recombination_saha_z = dll_import(c_double, "recombination", "recombination_saha_z")

#recombination_saha_tau = dll_import(c_double, "recombination", "recombination_saha_tau")


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
    
    def set_params(self, DM_Pann=1.0, FineS=1.0, EMass=1.0):
        """
        Set parameters for modified recombination
        :param DM_Pann: Dark matter annihilation parameter P_ann = f<v sigma>/m_dm in cm^3/s/GeV /10^-27 (w HyRec v jan 2015 only)
        :param FineS: Fine structure constant relative to measured value = 1.0
        :param EMass: Electron Mass relative to measured value = 1.0
        """

        self.Pann = DM_Pann*1.0e-27
        self.FineSRat = FineS
        self.EMassRat = EMass

        return self
