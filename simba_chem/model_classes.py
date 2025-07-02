"""
SIMBA Model Classes Module

This module defines the core data structures used in SIMBA using Python's dataclass 
framework. It provides structured containers for chemical species, reactions, and 
environmental conditions with built-in validation and information display.

Key Components:
- Elements: Container for basic chemical element data
- Species: Container for chemical species properties and abundances
- Gas: Container for gas phase properties
- Dust: Container for dust grain properties
- Reactions: Container for chemical reaction network data
- Environment: Container for physical environmental conditions
- Parameters: Container for simulation parameters and constants

Features:
- Type checking and validation for all attributes
- Default values for required parameters
- Info methods for displaying current state
- Error checking and validation

Dependencies:
   numpy: Required for numerical arrays
   dataclasses: Used for class definitions
   typing: Used for type hints
   numpy.typing: Used for array type hints
"""


import numpy as np
from dataclasses import dataclass, field
from typing import List
from numpy.typing import NDArray


############
# ELEMENTS #
############

class Elements:
    """
    Container for chemical elements data.
    Currently stores only element names.

    Attributes:
        name (List[str]): List of chemical element names.

    Notes:
        - Element names should follow standard chemical notation (e.g., 'H', 'He', 'C').
        - Default value is empty list.
        - Values are validated when set through input file.
    """
    
    name: List[str] = field(default_factory=list)
    
    def validate(self):
        """Validate chemical elements data."""
        # Type checks
        if not isinstance(self.name, list):
            raise TypeError("Element names must be in a list")
        if not all(isinstance(x, str) for x in self.name):
            raise TypeError("Element names must be strings")
            
        # Value checks
        if any(not x for x in self.name):
            raise ValueError("Element names cannot be empty")
        if len(self.name) != len(set(self.name)):
            raise ValueError("Duplicate element names found")
            
    def info(self):
        """Print basic information about elements."""
        print(f"Number of elements: {len(self.name)}")
        print(f"Elements included: {', '.join(self.name)}")


###########
# SPECIES #
###########

@dataclass
class Species:
    """
    Container for chemical species data in an astrochemical network.
    
    Attributes:
        name (List[str]): List of chemical species names.
        abundance (NDArray): Array of species abundances relative to total H.
        number (NDArray): Array of number densities in cm^-3.
        mass (NDArray): Array of masses in atomic mass units (amu).
        charge (NDArray): Array of charges in elementary charge units (e).
    
    Notes:
        - All arrays must have same length as name list.
        - Abundances must be between 0 and 1.
        - Number densities and masses must be positive.
        - Charges must be whole numbers, typically -4 to +4.
        - Values are validated when set through input file.
    """
    
    name: List[str] = field(default_factory=list)
    abundance: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    number: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    mass: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    charge: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    
    def validate(self):
        """Validate chemical species data."""
        # Type checks
        if not isinstance(self.name, list):
            raise TypeError("Species names must be in a list")
        if not all(isinstance(x, str) for x in self.name):
            raise TypeError("Species names must be strings")
            
        if not all(isinstance(x, np.ndarray) for x in [self.abundance, self.number, self.mass, self.charge]):
            raise TypeError("Numerical data must be NumPy arrays")
            
        # Check charges are whole numbers
        if not np.allclose(self.charge, np.round(self.charge)):
            raise ValueError("Charges must be whole numbers")
            
        # Length checks
        length = len(self.name)
        if not all(len(x) == length for x in [self.abundance, self.number, self.mass, self.charge]):
            raise ValueError("All species arrays must have the same length")
            
        # Value checks
        if np.any((self.abundance < 0) | (self.abundance > 1)):
            raise ValueError("Abundances must be between 0 and 1")
        if np.any(self.number < 0):
            raise ValueError("Number densities cannot be negative")
        if np.any(self.mass < 0):
            raise ValueError("Masses must be positive")
        if np.any(np.abs(self.charge) > 4):
            raise ValueError("Charges outside expected range (-4 to +4)")
    
    def info(self):
        """Print basic information about species."""
        n_species = len(self.name)
        n_neutral = np.sum(np.abs(self.charge) < 1e-10)  # Allow for floating point
        n_ions = n_species - n_neutral
        print("\nSpecies Information")
        print(f"Total species: {n_species}")
        print(f"Neutral species: {n_neutral}")
        print(f"Ionic species: {n_ions}")
        print(f"Mass range: {self.mass.min():.1f} to {self.mass.max():.1f} amu")
        print(f"Charge range: {int(self.charge.min()):+d} to {int(self.charge.max()):+d}")
        

#######
# GAS #
#######

@dataclass
class Gas:
    """
    Container for gas properties.
    
    Attributes:
        n_gas (float): Number density of the gas in cm^-3.
        t_gas (float): Gas temperature in K.
        h2_col (float): H2 column density in cm^-2.
    
    Notes:
        - n_gas and temperature must be positive.
        - h2_col must be non-negative.
        - Values are validated when set through input file.
    """
    
    n_gas: float = 1.0
    t_gas: float = 1.0
    h2_col: float = 1.0
    
    def validate(self):
        """Validate gas properties."""
        # Type checks
        if not all(isinstance(x, (int, float)) for x in [self.n_gas, self.t_gas, self.h2_col]):
            raise TypeError("All gas properties must be numerical")
            
        # Value checks
        if self.n_gas <= 0:
            raise ValueError("Gas number density must be positive")
        if self.t_gas <= 0:
            raise ValueError("Gas temperature must be positive")
        if self.h2_col < 0:
            raise ValueError("H2 column density cannot be negative")
    
    def info(self):
        """Print basic information about gas properties."""
        print("\nGas Properties")
        print(f"Number density: {self.n_gas:.2e} cm^-3")
        print(f"t_gas: {self.t_gas:.1f} K")
        print(f"H2 column density: {self.h2_col:.2e} cm^-2")


########
# DUST #
########

@dataclass
class Dust:
    """
    Container for dust grain properties.
    
    Attributes:
        n_dust (float): Number density of dust grains in cm^-3.
        t_dust (float): Dust temperature in K.
        radius (float): Grain radius in cm.
        binding_sites (float): Number of molecular binding sites per grain.
    
    Notes:
        - All values must be positive.
        - Default radius is 1e-5 cm (0.1 micron), typical for ISM grains.
        - Default binding_sites is 1e6, typical for 0.1 micron grains.
        - Values are validated when set through input file.
    """
    
    n_dust: float = 1.0
    t_dust: float = 1.0
    radius: float = 1e-5
    binding_sites: float = 1e6
    
    def validate(self):
        """Validate dust properties."""
        # Type checks
        if not all(isinstance(x, (int, float)) for x in 
                  [self.n_dust, self.t_dust, self.radius, 
                   self.binding_sites]):
            raise TypeError("All dust properties must be numerical")
            
        # Value checks
        if self.n_dust <= 0:
            raise ValueError("Dust number density must be positive")
        if self.t_dust <= 0:
            raise ValueError("Temperature must be positive")
        if self.radius <= 0:
            raise ValueError("Grain radius must be positive")
        if self.binding_sites <= 0:
            raise ValueError("Number of binding sites must be positive")
    
    def info(self):
        """Print basic information about dust properties."""
        print("\nDust Properties")
        print(f"Dust number density: {self.n_dust:.2e} cm^-3")
        print(f"Dust temperature: {self.t_dust:.1f} K")
        print(f"Grain radius: {self.radius:.2e} cm")
        print(f"Binding sites per grain: {self.binding_sites:.2e}")


#############
# REACTIONS #
#############

@dataclass
class Reactions:
    """
    Container for chemical reaction data.
    
    Attributes:
        educts (List[List[str]]): Reactants for each reaction.
        products (List[List[str]]): Products for each reaction.
        reaction_id (NDArray): Unique identifier for each reaction.
        itype (NDArray): Reaction type identifiers.
        a (NDArray): Rate coefficient parameter 'a'.
        b (NDArray): Rate coefficient parameter 'b'.
        c (NDArray): Rate coefficient parameter 'c'.
        temp_min (NDArray): Minimum valid temperature in K.
        temp_max (NDArray): Maximum valid temperature in K.
        pd_data (NDArray): Photodissociation data.
        k (NDArray): Calculated rate coefficients.
        labels (List[str]): Human-readable reaction descriptors.
    
    Notes:
        - All arrays must have the same length (number of reactions).
        - Reaction IDs must be unique positive integers.
        - Temperature limits must be physically reasonable (temp_max > temp_min ≥ 0).
        - Rate coefficients k must be non-negative.
        - Values are validated when set through input file.
    """
    
    educts: List[List[str]] = field(default_factory=list)
    products: List[List[str]] = field(default_factory=list)
    reaction_id: NDArray = field(default_factory=lambda: np.array([], dtype=np.int32))
    itype: NDArray = field(default_factory=lambda: np.array([], dtype=np.int32))
    a: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    b: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    c: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp_min: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp_max: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    pd_data: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    k: NDArray = field(default_factory=lambda: np.array([], dtype=np.float64))
    labels: List[str] = field(default_factory=list)
    
    def validate(self):
        """Validate reaction data."""
        # Skip validation if no reactions loaded
        if not self.educts:
            return
            
        # Type checks
        if not isinstance(self.educts, list) or not all(isinstance(x, list) for x in self.educts):
            raise TypeError("Educts must be a list of lists")
        if not isinstance(self.products, list) or not all(isinstance(x, list) for x in self.products):
            raise TypeError("Products must be a list of lists")
        if not all(isinstance(x, str) for sublist in self.educts + self.products for x in sublist):
            raise TypeError("Species names must be strings")
            
        arrays = [self.reaction_id, self.itype, self.a, self.b, self.c,
                 self.temp_min, self.temp_max]
        if not all(isinstance(x, np.ndarray) for x in arrays):
            raise TypeError("Numerical data must be NumPy arrays")
            
        # Length checks
        length = len(self.educts)
        check_lists = [self.products, self.labels]
        check_arrays = [self.reaction_id, self.itype, self.a, self.b, self.c,
                       self.temp_min, self.temp_max, self.pd_data, self.k]
        
        if not all(len(x) == length for x in check_lists + check_arrays):
            raise ValueError("All reaction properties must have the same length")
            
        # Value checks
        if len(np.unique(self.reaction_id)) != length:
            raise ValueError("Reaction IDs must be unique")
        if np.any(self.reaction_id < 0):
            raise ValueError("Reaction IDs must be positive")
            
    
    def info(self):
        """Print basic information about reactions."""
        if not self.educts:
            print("\nNo reactions loaded")
            return
            
        n_reactions = len(self.educts)
        reaction_types = np.unique(self.itype)
        
        print("\nReaction Network Information")
        print(f"Total reactions: {n_reactions}")
        print(f"Unique reaction types: {len(reaction_types)}")
        print("\nReaction Types:")
        for rt in reaction_types:
            n_type = np.sum(self.itype == rt)
            print(f"Type {rt:2d}: {n_type:4d} reactions")
        

###############
# ENVIRONMENT #
###############

@dataclass
class Environment:
    """
    Container for physical conditions.
    
    Attributes:
        gtd (float): Gas-to-dust mass ratio.
        Av (float): Visual extinction in magnitudes.
        G_0 (float): FUV radiation field strength in Habing units.
        G_0_unatt (float): Unattenuated FUV field in Habing units.
        Zeta_X (float): X-ray ionization rate in s^-1.
        Zeta_CR (float): Cosmic ray ionization rate in s^-1.
        pah_ism (float): PAH abundance relative to ISM value.
        dg100 (float): Normalized gas-to-dust ratio (gtd/100).
    
    Notes:
        - All values must be non-negative.
        - G_0_unatt is calculated from G_0 and Av.
        - dg100 is derived from gtd and includes He correction factor.
        - Values are validated when set through input file.
    """
    
    gtd: float = 1.0
    Av: float = 1.0
    G_0: float = 1.0
    G_0_unatt: float = 1.0
    Zeta_X: float = 1.0
    Zeta_CR: float = 1.0
    pah_ism: float = 1.0
    dg100: float = 1.0
    
    def validate(self):
        """Validate physical conditions."""
        # Type checks
        if not all(isinstance(x, (int, float)) for x in 
                  [self.gtd, self.Av, self.G_0, self.G_0_unatt,
                   self.Zeta_X, self.Zeta_CR, self.pah_ism, self.dg100]):
            raise TypeError("All conditions must be numerical")
            
        # Value checks
        if self.gtd <= 0:
            raise ValueError("Gas-to-dust ratio must be positive")
        if self.Av < 0:
            raise ValueError("Visual extinction cannot be negative")
        if self.G_0 < 0:
            raise ValueError("Radiation field strength cannot be negative")
        if self.G_0_unatt < 0:
            raise ValueError("Unattenuated radiation field cannot be negative")
        if self.Zeta_X < 0:
            raise ValueError("X-ray ionization rate cannot be negative")
        if self.Zeta_CR < 0:
            raise ValueError("Cosmic ray ionization rate cannot be negative")
        if self.pah_ism < 0:
            raise ValueError("PAH abundance cannot be negative")
        if self.dg100 <= 0:
            raise ValueError("Normalized gas-to-dust ratio must be positive")
    
    def info(self):
        """Print basic information about physical conditions."""
        print("\nPhysical Conditions")
        print(f"Gas-to-dust ratio: {self.gtd:.1f}")
        print(f"Visual extinction: {self.Av:.2f} mag")
        print(f"FUV field (G0): {self.G_0:.1e} Draine")
        print(f"Unattenuated FUV: {self.G_0_unatt:.1e} Draine")
        print(f"X-ray ionization rate: {self.Zeta_X:.2e} s^-1")
        print(f"Cosmic ray ionization rate: {self.Zeta_CR:.2e} s^-1")
        print(f"PAH abundance (rel. to ISM): {self.pah_ism:.2f}")
        print(f"Normalized gas-to-dust: {self.dg100:.2f}")
        
        
##############
# PARAMETERS #
##############

@dataclass
class Parameters:
    """
    Container for simulation parameters and physical constants.
    
    Attributes:
        n_elements (int): Number of chemical elements in network.
        n_species (int): Number of chemical species in network.
        n_reactions (int): Number of reactions in network.
        time_initial (float): Initial time in seconds.
        time_final (float): Final time in seconds.
        timeout_duration (float): Timeout solver in seconds
        delta_v (float): Velocity dispersion in km/s for self-shielding.
        av_nH (float): Column to visual extinction conversion factor.
        self_shielding (bool): Flag for self-shielding.
        column (bool): Flag for using column-based self-shielding.
        k_B (float): Boltzmann constant in erg/K.
        yr_sec (float): Number of seconds in a year.
        m_p (float): Proton mass in g.
        save_logs (bool): Flag to control file logging.
        verbose (bool): Flag to control console output.
        verbosity_level (str): Logging level - "ERROR", "WARNING", "INFO", or "DEBUG".
    
    Notes:
        - Network size parameters must be non-negative integers.
        - Time values must be positive, with final > initial.
        - Physical constants are in cgs units and must be positive.
        - Values are validated when set through input file.
        - A timeout is used to prevent the solver from running too long, and will 
          automatically move on to the next integration method (default = 10 minutes).
    """
    
    # Network size parameters
    n_elements: int = 0
    n_species: int = 0
    n_reactions: int = 0
    
    # Time parameters
    time_initial: float = 1.0
    time_final: float = 1.0
    timeout_duration: int = 600
    
    # Physical parameters
    delta_v: float = 0.2
    av_nH: float = 1/5.34e-22
    self_shielding: bool = True
    column: bool = True
    
    # Physical constants in cgs
    k_B: float = 1.3806504e-16
    yr_sec: float = 3.1556926e7
    m_p: float = 1.660538e-24
    
    # Logging
    save_logs: bool = True
    verbose: bool = True
    verbosity_level: str = "INFO"


    
    def validate(self):
        """Validate parameters."""
        # Type checks
        for name in ['n_elements', 'n_species', 'n_reactions']:
            if not isinstance(getattr(self, name), int):
                raise TypeError(f"{name} must be an integer")
                
        for name in ['time_initial', 'time_final', 'delta_v', 'av_nH',
                    'k_B', 'yr_sec', 'm_p']:
            if not isinstance(getattr(self, name), (int, float)):
                raise TypeError(f"{name} must be numerical")
    
        if not isinstance(self.self_shielding, bool):
            raise TypeError("self_shielding must be boolean")
        
        if not isinstance(self.column, bool):
            raise TypeError("column must be boolean")
        
        if not isinstance(self.save_logs, bool):
            raise TypeError("save_logs must be boolean")
            
        # Value checks
        if any(x < 0 for x in [self.n_elements, self.n_species, self.n_reactions]):
            raise ValueError("Network size parameters cannot be negative")
        if self.time_initial < 0:
            raise ValueError("Initial time cannot be negative")
        if self.time_final <= 0:
            raise ValueError("Final time must be positive")
        if self.time_final < self.time_initial:
            raise ValueError("Final time must be greater than initial time")
        if self.delta_v <= 0:
            raise ValueError("Velocity dispersion must be positive")
        if self.av_nH <= 0:
            raise ValueError("Column conversion factor must be positive")
        if any(x <= 0 for x in [self.k_B, self.yr_sec, self.m_p]):
            raise ValueError("Physical constants must be positive")
    
    def info(self):
        """Print basic information about parameters."""
        print("\nNetwork Parameters")
        print(f"Elements: {self.n_elements}")
        print(f"Species: {self.n_species}")
        print(f"Reactions: {self.n_reactions}")
        print("\nTime Setup")
        print(f"Initial time: {self.time_initial:.2e} s")
        print(f"Final time: {self.time_final:.2e} s")
        print(f"Duration: {(self.time_final - self.time_initial)/self.yr_sec:.2e} yr")
        print("\nPhysical Setup")
        print(f"Velocity dispersion: {self.delta_v:.2f} km/s")
        print(f"Column density mode: {'On' if self.column else 'Off'}")