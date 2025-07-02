"""
SIMBA: Solver for Inferring Molecular aBundances in Astrophysical environments

This code implements a comprehensive solver for chemical reaction networks
in astrophysical environments. It is designed to model and simulate complex
chemical processes in various cosmic settings such as the ISM, molecular clouds,
and protoplanetary disks.

Key Features:
- Initialization of chemical species, reactions, and environmental parameters
- Efficient solving of stiff ODEs representing chemical reactions
- Support for various reaction types including gas-phase, grain-surface, and photochemistry
- Integration of self-shielding factors for specific molecules (H2, CO, N2, C)
- Optimization using Numba JIT compilation for performance-critical functions
- Comprehensive logging and progress tracking

Main Components:
- Classes: Elements, Species, Gas, Dust, Environment, Reactions, Parameters, Simba
- Key Functions: calculate_derivatives, calculate_jacobian, compute_rate_coefficients, solve_network

Usage:
Initialize a Simba object, set up the network parameters, and call solve_network()
to integrate the chemical evolution over time.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time
from tqdm.auto import tqdm
import signal
import time
from contextlib import contextmanager
import platform
import logging
import sys
import os
from . import helpers
from . import model_classes
from . import calculus
from . import self_shielding as ss
from .data import CO_SELFSHIELDING_FILE, N2_SELFSHIELDING_FILE
from .helpers import safe_exp

def setup_windows_console():
    """
    Simple function to enable UTF-8 support on Windows console.
    Call this once at the start of your program.
    """
    if platform.system() == 'Windows':
        try:
            # Try to enable UTF-8 mode for Windows console
            os.system('chcp 65001 > nul')
            # Set stdout/stderr encoding to utf-8
            if hasattr(sys.stdout, 'reconfigure'):
                sys.stdout.reconfigure(encoding='utf-8')
                sys.stderr.reconfigure(encoding='utf-8')
        except:
            # If anything fails, just continue - the error handler will catch issues
            pass

class SafeConsoleHandler(logging.StreamHandler):
    """
    Custom logging handler that gracefully handles Unicode encoding errors.
    Falls back to ASCII representation when Unicode fails.
    """
    def emit(self, record):
        try:
            # Try normal logging first
            super().emit(record)
        except UnicodeEncodeError:
            # If Unicode fails, convert the message to ASCII and try again
            try:
                # Get the formatted message
                msg = self.format(record)
                # Convert Unicode to ASCII equivalents
                ascii_msg = msg.encode('ascii', 'replace').decode('ascii')
                # Create a new record with the ASCII message
                ascii_record = logging.LogRecord(
                    record.name, record.levelno, record.pathname, record.lineno,
                    ascii_msg, record.args, record.exc_info, record.funcName, record.stack_info
                )
                # Emit the ASCII version
                logging.StreamHandler.emit(self, ascii_record)
            except:
                # Ultimate fallback - just print a simple message
                print(f"[{record.levelname}] Unicode encoding error in log message")

class Simba:
    def __init__(self):

        setup_windows_console()
 
        self.elements = model_classes.Elements()
        self.species = model_classes.Species()
        self.gas = model_classes.Gas()
        self.dust = model_classes.Dust()
        self.environment = model_classes.Environment()
        self.reactions = model_classes.Reactions()
        self.parameters = model_classes.Parameters()
        self.abundance_history = []  
        self._progress_bar = None
        
        
        ##########
        # Logger #
        ##########
        
        logger = logging.getLogger()
        
        # Clear any existing handlers
        if logger.hasHandlers():
            logger.handlers.clear()
        
        # Set root logger to lowest level needed (DEBUG)
        logger.setLevel(logging.DEBUG)
        
        # Create console handler - only active if verbose is True
        if self.parameters.verbose:
            console_handler = SafeConsoleHandler()
            # Set level based on verbosity_level parameter
            level = getattr(logging, self.parameters.verbosity_level)
            console_handler.setLevel(level)
            console_format = logging.Formatter('%(message)s')
            console_handler.setFormatter(console_format)
            logger.addHandler(console_handler)
        
        # Optionally create file handler if user wants logs saved
        if self.parameters.save_logs:
            file_handler = logging.FileHandler('simba.log')
            file_handler.setLevel(logging.INFO)  # Could make this configurable too
            file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_format)
            logger.addHandler(file_handler)


    def set_verbosity(self, verbose=True, level="INFO"):

        logger = logging.getLogger()
        
        # Update parameters
        self.parameters.verbose = verbose
        self.parameters.verbosity_level = level
        
        # Remove existing console handlers
        for handler in logger.handlers[:]:
            if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                logger.removeHandler(handler)
        
        # Add new console handler if verbose
        if verbose:
            console_handler = logging.StreamHandler()
            log_level = getattr(logging, level)
            console_handler.setLevel(log_level)
            console_format = logging.Formatter('%(message)s')
            console_handler.setFormatter(console_format)
            logger.addHandler(console_handler)


    ##############
    # Initialise #
    ##############
    
    def init_simba(self, input_file):
        """Initialize the SIMBA chemical network solver"""

        logging.info("\n")
        logging.info("┏" + ("━" * 68) + "┓")
        logging.info("┃" + "".center(68) + "┃")
        logging.info("┃" + "┏┓┳┳┳┓┳┓┏┓".center(68) + "┃")
        logging.info("┃" + "┗┓┃┃┃┃┣┫┣┫".center(68) + "┃")
        logging.info("┃" + "┗┛┻┛ ┗┻┛┛┗".center(68) + "┃")
        logging.info("┃" + "Astrochemical Network Solver".center(68) + "┃")
        logging.info("┃" + "by Luke Keyte".center(68) + "┃")
        logging.info("┃" + "".center(68) + "┃")
        logging.info("┃" + "Version 1.0.2 | 2025".center(68) + "┃")
        logging.info("┃" + "".center(68) + "┃")
        logging.info("┗" + ("━" * 68) + "┛")
        logging.info("\n")
                                
        
        input_data = helpers.read_input_file(input_file)
        

        # Input parameters
        helpers.log_section("INITIALIZATION")
        try:   
            logging.info(" ◆ Loading input parameters")
            self.gas.n_gas = input_data['n_gas']
            self.gas.t_gas = input_data['t_gas']
            self.dust.n_dust = input_data['n_dust']
            self.dust.t_dust = input_data['t_dust']
            self.environment.gtd = input_data['gtd']
            self.environment.Av = input_data['Av']
            self.environment.G_0 = input_data['G_0']
            self.environment.G_0_unatt = input_data['G_0'] * np.exp(np.minimum(3.02 * input_data['Av'], 700))
            self.environment.Zeta_X = input_data['Zeta_X']
            self.environment.Zeta_CR = input_data['Zeta_CR']
            self.environment.pah_ism = input_data['pah_ism']
            self.environment.dg100 = ((1/self.environment.gtd) * 100) / 1.4
            self.parameters.time_final = input_data['t_chem'] * self.parameters.yr_sec
            self.parameters.self_shielding = input_data['self_shielding']
            self.parameters.column = input_data['column']
            if self.parameters.column:
                self.gas.h2_col = input_data['h2_col']
            logging.info("    ► Parameters loaded successfully\n")
        except AttributeError as e:
            logging.error(f"Missing required input parameter: {str(e)}")
            raise
        except Exception as e:
            logging.error(f"Failed to load input parameters: {str(e)}")
            raise


        # Chemical network
        try:
            logging.info(" ◆ Loading chemical network")
            network_data = helpers.read_chemnet(input_data['network'])
            
            self.parameters.n_elements = network_data[0]
            self.elements.name = network_data[1]
            self.parameters.n_species = network_data[2]
            self.species.name = network_data[3]
            self.species.abundance = network_data[4]
            self.species.mass = network_data[5]
            self.species.charge = network_data[6]
            self.parameters.n_reactions = network_data[7]
            self.reactions.educts = network_data[8]
            self.reactions.products = network_data[9]
            self.reactions.reaction_id = network_data[10]
            self.reactions.itype = network_data[11]
            self.reactions.a = network_data[12]
            self.reactions.b = network_data[13]
            self.reactions.c = network_data[14]
            self.reactions.temp_min = network_data[15]
            self.reactions.temp_max = network_data[16]
            self.reactions.pd_data = network_data[17]
            self.reactions.k = np.zeros(self.parameters.n_reactions)
            self.reactions.labels = helpers.create_reaction_labels(self.parameters.n_reactions, self.reactions.educts, self.reactions.products)  
            
            # PAH abundance
            i_PAH0 = self.species.name.index('PAH0')
            i_PAHp = self.species.name.index('PAH+')
            i_PAHm = self.species.name.index('PAH-')
            i_PAHH = self.species.name.index('PAH_H')
            if self.environment.G_0 > 1e6:
                # PAHS are photodistroyed when G0 > 1e6 (Visser et a. 2007)
                self.species.abundance[i_PAH0] = 1e-10
                self.species.abundance[i_PAHp] = 0.0
                self.species.abundance[i_PAHm] = 0.0
                self.species.abundance[i_PAHH] = 0.0
            else:
                self.species.abundance[i_PAH0] = 4.17e-07 * self.environment.pah_ism
                self.species.abundance[i_PAHp] = 0.0
                self.species.abundance[i_PAHm] = 0.0
                self.species.abundance[i_PAHH] = 0.0
            
            # Absolute abundances
            self.species.number = self.species.abundance * self.gas.n_gas
            
            logging.info("    ► Network successfully loaded")
            logging.info(f"      • {self.parameters.n_elements} elements")
            logging.info(f"      • {self.parameters.n_species} atomic/molecular species")
            logging.info(f"      • {self.parameters.n_reactions} reactions\n")
            
        except Exception as e:
            logging.error("Failed to load chemical network")
            logging.error(f"Error details: {str(e)}")
            
        
        # Self-shielding factors
        logging.info(" ◆ Loading self-shielding data")
        try:
            self.ss_co = ss.read_selfshielding_co(CO_SELFSHIELDING_FILE)
            logging.info("    ► CO self-shielding: Loaded")
        except FileNotFoundError:
            logging.error("CO self-shielding data file not found: data_selfshielding_c.dat")
            raise
        except Exception as e:
            logging.error(f"Error reading CO self-shielding data: {str(e)}")
            raise
        try:
            self.ss_n2 = ss.read_selfshielding_n2(N2_SELFSHIELDING_FILE)
            logging.info("    ► N₂ self-shielding: Loaded\n")
        except FileNotFoundError:
            logging.error("N2 self-shielding data file not found: data_selfshielding_n2.dat")
            raise
        except Exception as e:
            logging.error(f"Error reading N2 self-shielding data: {str(e)}")
            raise
        
        # Validation
        self.elements.validate()
        self.species.validate()
        self.gas.validate()
        self.dust.validate()
        self.environment.validate()
        self.parameters.validate()
        self.reactions.validate()
        
        # Setup reaction indices
        self.setup_reaction_indices()

        # Print system parameters
        helpers.log_section("SYSTEM PARAMETERS")
        helpers.log_table_header(" ◆ PHYSICAL CONDITIONS")
        helpers.log_table_row("Gas Density", self.gas.n_gas, "cm^-3")
        helpers.log_table_row("Dust Density", self.dust.n_dust, "cm^-3")
        helpers.log_table_row("Gas Temperature", self.gas.t_gas, "K")
        helpers.log_table_row("Dust Temperature", self.dust.t_dust, "K")
        helpers.log_table_row("Gas/Dust Ratio", self.environment.gtd)
        helpers.log_table_row("Visual Extinction", self.environment.Av)
        helpers.log_table_footer()
        helpers.log_table_header(" ◆ RADIATION FIELD")
        helpers.log_table_row("UV Field", self.environment.G_0, "G0")
        helpers.log_table_row("X-ray Rate", self.environment.Zeta_X, "s^-1")
        helpers.log_table_row("Cosmic Ray Rate", self.environment.Zeta_CR, "s^-1")
        helpers.log_table_row("PAH/ISM Ratio", self.environment.pah_ism)
        helpers.log_table_footer()
        helpers.log_table_header(" ◆ CALCULATION PARAMETERS")
        helpers.log_table_row("Chemical Time", self.parameters.time_final/self.parameters.yr_sec, "years")
        helpers.log_table_row("Self-shielding", self.parameters.self_shielding)
        helpers.log_table_row("Use H Column?", self.parameters.column)
        if self.parameters.column:
            helpers.log_table_row("H Column", self.gas.h2_col, "cm^-2")
        helpers.log_table_footer()

        logging.info("")


    ##########################
    # Setup reaction indices #
    ##########################
    
    def setup_reaction_indices(self):
        """
        Create index arrays for reactions to speed up calculations
        """
        self.idx = np.zeros((self.parameters.n_reactions, 8), dtype=np.int32)
        
        for i in range(self.parameters.n_reactions):
            # Convert educts and products to indices in species array
            for j, educt in enumerate(self.reactions.educts[i]):
                if educt != '' and educt != 'M' and educt != 'PHOTON' and educt != 'CRP' and educt != 'CRPHOT' and educt != 'XELECTRON':
                    self.idx[i,j] = self.species.name.index(educt)
                else:
                    self.idx[i,j] = -1
                    
            for j, product in enumerate(self.reactions.products[i]):
                if product != '' and product != 'M' and product != 'PHOTON' and product != 'CRP' and product != 'CRPHOT' and product != 'XELECTRON':
                    self.idx[i,j+3] = self.species.name.index(product)
                else:
                    self.idx[i,j+3] = -1

    
    #############################
    # Compute rate coefficients #
    #############################
    
    def compute_rate_coefficients(self, y):
        """
        Computes reaction rate coefficients for all chemical reactions in the network.

        This method calculates the rate coefficients for various types of reactions including:
        gas-phase reactions, surface reactions, photochemistry, cosmic ray interactions,
        and grain-surface chemistry. It handles different reaction types through their
        itype identifiers, matching those used in DALI (Bruderer et al. 2012, 2013).

        Included:
            - H2 formation on dust grains (Cazaux & Tielens 2002/04, Bosman+22)
            - Hydrogenation reactions (Visser+11)
            - Photodesorption processes (Visser+11)
            - Gas-phase reactions with temperature dependencies
            - Photodissociation with self-shielding 
                (SS factors from Draine & Bertoldi 1996, Kamp & Bertoldi 2000, Visser+09, Visser+18)
            - Cosmic ray induced reactions (Stauber+05, Heays+14, Visser+18)
            - X-ray induced reactions (Stauber+05, Bruderer+09b)
            - PAH/grain charge exchange
            - Thermal desorption and freezeout (Visser+09a, Visser+11)
            - H2 excitation processes (Tielens & Hollenbach 1985, Visser+18)

        Parameters:
            y (numpy.ndarray): Current abundances of all species in the network,
                            in units of cm^-3.

        Updates:
            self.reactions.k (numpy.ndarray): Array of rate coefficients for all reactions.
            self.disso_H2 (float): H2 dissociation rate, used for reaction types 90 & 91.

        Raises:
            ValueError: If an unknown reaction type is encountered or if an invalid
                    rate coefficient is calculated (NaN or Inf).

        Notes:
            - Self-shielding factors are applied for H2, CO, N2, and C when enabled
            - PAH and grain surface chemistry are included with appropriate scaling
            - Temperature and UV field limits are enforced for relevant reaction types
            - All rates are validated to prevent invalid values
        """
                        
        # Index of some key species
        i_H   = self.species.name.index('H')
        i_H2  = self.species.name.index('H2')
        i_He  = self.species.name.index('He')
        i_CO  = self.species.name.index('CO')
        i_H2O = self.species.name.index('H2O')
        i_C   = self.species.name.index('C')
        i_N2  = self.species.name.index('N2')
        
        # Initialise H2 dissociation rate
        self.disso_H2 = 0.0
        
        # Dust properties
        ngrndust100 = 1e-12                                                     # gtd number density for gtd=100
        n_gr        = self.gas.n_gas * ngrndust100 * self.environment.dg100     # n_grains used for chemistry (not self.dust.ndust!!)
        n_gr        = np.maximum(n_gr, 1e-30)                                   # safeguard for n_gr to prevent underflow
        
        # Inputs to prefactors for hydrogenation/freezeout
        n_hydro = 0.0
        n_ice   = 0.0
        for i in range(0, self.parameters.n_reactions):
            ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5 = self.idx[i,:]
            if self.reactions.itype[i] == 11:
                n_hydro += y[ir1]
            if self.reactions.itype[i] == 80:
                n_ice += y[ir1]
        
        # Safeguard for n_hydro and n_ice to prevent division by zero
        n_hydro = np.maximum(n_hydro, 1e-30)
        n_ice = np.maximum(n_ice, 1e-30)

        # Loop through all reactions
        for i in range(0, self.parameters.n_reactions):
            ir1, ir2, ir3, ip1, ip2, ip3, ip4, ip5 = self.idx[i,:]             # indices for this reaction

            # 10: H2 formation on grains
            if self.reactions.itype[i] == 10:
                if self.dust.t_dust < 10.0:
                    eta   = 1.0
                    stick = 1.0/(1.0 + 0.04 * np.sqrt(self.gas.t_gas + self.dust.t_dust) + 2e-3 * self.gas.t_gas + 8e-6 * self.gas.t_gas**2)
                else:
                    # Calculate mean velocity for monolayers per second
                    v_mean_h2  = np.sqrt(8.0 * self.parameters.k_B * self.gas.t_gas / (np.pi * self.parameters.m_p))
                    stick      = 1.0/(1.0 + 0.04 * np.sqrt(self.gas.t_gas + self.dust.t_dust) + 2e-3 * self.gas.t_gas + 8e-6 * self.gas.t_gas**2)
                    f_mlps     = v_mean_h2 * y[i_H] * np.pi * self.dust.radius**2 / self.dust.binding_sites * stick                    
                    f_mlps     = np.maximum(f_mlps, 1e-30)
                    sqterm     = (1.0 + np.sqrt((30000.0-200.0)/(600.0-200.0)))**2
                    beta_H2    = 3e12 * safe_exp(-320.0/self.dust.t_dust)
                    beta_alpha = 0.25 * sqterm * safe_exp(-200.0/self.dust.t_dust)
                    xi         = 1.0/(1.0 + 1.3e13/(2.0 * f_mlps) * safe_exp(-1.5 * 30000.0/self.dust.t_dust) * sqterm)
                    eta        = xi/(1.0 + 0.005*f_mlps/(2.0*beta_H2) + beta_alpha)
    
                s_eta = stick * eta
                k = s_eta * self.reactions.a[i] * (self.gas.t_gas**self.reactions.b[i]) / (1e-10 + y[i_H]) * self.gas.n_gas * self.environment.dg100
            
            # 11: Hydrogenation
            elif self.reactions.itype[i] == 11:
                if self.gas.n_gas * self.environment.dg100 < 1e3:
                    # prehydro = 1e-99
                    prehydro = np.pi * self.dust.radius**2 * n_gr / np.maximum(n_hydro, self.dust.binding_sites*n_gr)
                else:
                    prehydro = np.pi * self.dust.radius**2 * n_gr / np.maximum(n_hydro, self.dust.binding_sites*n_gr)
                k = prehydro * np.sqrt(8.0*self.parameters.k_B/(np.pi*self.parameters.m_p))*np.sqrt(self.gas.t_gas)
                  
            # 12: Photodesorption
            elif self.reactions.itype[i] == 12:
                if self.gas.n_gas * self.environment.dg100 < 1e3:
                    # pregrain = 1e-99
                    pregrain = np.pi * self.dust.radius**2 * n_gr / np.maximum(n_ice, self.dust.binding_sites*n_gr)
                else:
                    pregrain = np.pi * self.dust.radius**2 * n_gr / np.maximum(n_ice, self.dust.binding_sites*n_gr)
                k = pregrain * self.reactions.a[i] * 1e8 * (self.environment.G_0+1e-4*((self.environment.Zeta_CR+self.environment.Zeta_X)/5e-17))

            # 20: Normal gas-phase reaction        
            elif self.reactions.itype[i] == 20:
                k = self.reactions.a[i] * (self.gas.t_gas/300.0)**self.reactions.b[i] * safe_exp(-self.reactions.c[i]/self.gas.t_gas)
                # Not the highest temperature reaction available
                if self.gas.t_gas >= np.abs(self.reactions.temp_max[i]) and self.reactions.temp_max[i] < 0.0:
                   k = 0.0
                # Keep constant above maximum temperature
                if self.gas.t_gas > np.abs(self.reactions.temp_max[i]) and self.reactions.temp_max[i] > 0.0:
                    k = self.reactions.a[i] * (self.reactions.temp_max[i] / 300.0) ** self.reactions.b[i] * safe_exp(-self.reactions.c[i] / self.reactions.temp_max[i])
                # Switch off below minimum temperature
                if self.gas.t_gas < self.reactions.temp_min[i]:
                    k = 0.0
                    
            # 21: Normal gas-phase reaction (do not extrapolate in temperature)
            elif self.reactions.itype[i] == 21:
                k = self.reactions.a[i] * (self.gas.t_gas/300.0)**self.reactions.b[i] * safe_exp(-self.reactions.c[i]/self.gas.t_gas)
                if self.gas.t_gas > self.reactions.temp_max[i]:
                    k = self.reactions.a[i] * (self.reactions.temp_max[i]/300.0)**self.reactions.b[i] * safe_exp(-self.reactions.c[i]/self.reactions.temp_max[i])
                if self.gas.t_gas < self.reactions.temp_min[i]:
                    k = self.reactions.a[i] * (self.reactions.temp_min[i]/300.0)**self.reactions.b[i] * safe_exp(-self.reactions.c[i]/self.reactions.temp_min[i])
                    
            # 22: Normal gas-phase reaction (switch off outside temperature range)
            elif self.reactions.itype[i] == 22:
                k = self.reactions.a[i] * (self.gas.t_gas/300.0)**self.reactions.b[i] * safe_exp(-self.reactions.c[i]/self.gas.t_gas)
                if self.gas.t_gas > self.reactions.temp_max[i] or self.gas.t_gas < self.reactions.temp_min[i]:
                    k = 0.0

            # 30: Photodissociation
            elif self.reactions.itype[i] in (30, 34, 35, 36, 37, 39):
                safe_G0_unatt = np.clip(self.environment.G_0_unatt, 0.0, 1e50)
                k = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av)                

            # 31: H2 dissociation inc. self-shielding
            elif self.reactions.itype[i] == 31:
                safe_G0_unatt = np.clip(self.environment.G_0_unatt, 0.0, 1e50)
                if self.parameters.self_shielding:
                    if self.parameters.column:
                        col_h2 = self.gas.h2_col
                    else:
                        col_h2 = self.environment.Av * self.parameters.av_nH * (np.maximum(y[i_H2], 0.0)/self.gas.n_gas)
                    ssfact_H2  = ss.calc_selfshielding_h2(col_h2, self.parameters.delta_v) 
                    k = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av) * ssfact_H2
                else:
                    k = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av)
                self.disso_H2 = k   # needed for reaction types 90 & 91

            # 32: CO dissociation inc. self-shielding
            elif self.reactions.itype[i] == 32:
                safe_G0_unatt = np.clip(self.environment.G_0_unatt, 0.0, 1e50)
                if self.parameters.self_shielding:
                    if self.parameters.column:
                        col_h2 = self.gas.h2_col
                    else:
                        col_h2 = self.environment.Av * self.parameters.av_nH * (np.maximum(y[i_H2], 0.0)/self.gas.n_gas)
                    col_co = (col_h2/np.maximum(1e-20, (np.maximum(y[i_H2], 0.0)/self.gas.n_gas))) * (np.maximum(y[i_CO], 0.0)/self.gas.n_gas)
                    ssfactor_co = ss.calc_selfshielding_co(self.ss_co[0], self.ss_co[1], self.ss_co[2], col_h2, col_co)
                    k  = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av) * ssfactor_co
                else:
                    k  = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av)

            # 33: C ionization inc. self-shielding
            elif self.reactions.itype[i] == 33:
                safe_G0_unatt = np.clip(self.environment.G_0_unatt, 0.0, 1e50)
                if self.parameters.self_shielding:
                    if self.parameters.column:
                        col_h2 = self.gas.h2_col
                    else:
                        col_h2 = self.environment.Av * self.parameters.av_nH * (np.maximum(y[i_H2], 0.0)/self.gas.n_gas)  
                    col_c = (col_h2/np.maximum(1e-20, (np.maximum(y[i_H2], 0.0)/self.gas.n_gas))) * (np.maximum(y[i_C], 0.0)/self.gas.n_gas)
                    ssfactor_c = ss.calc_selfshielding_c(col_h2, col_c, self.gas.t_gas)
                    k  = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av) * ssfactor_c
                else:
                    k  = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av)

            # 38: N2 photodissociation inc. self-shielding
            elif self.reactions.itype[i] == 38:
                safe_G0_unatt = np.clip(self.environment.G_0_unatt, 0.0, 1e50)
                if self.parameters.self_shielding:
                    if self.parameters.column:
                        col_h2 = self.gas.h2_col
                    else:
                        col_h2 = self.environment.Av * self.parameters.av_nH * (y[i_H2]/self.gas.n_gas)
                    col_h  = (col_h2/np.maximum(1e-20, (np.maximum(y[i_H2], 0.0)/self.gas.n_gas))) * (np.maximum(y[i_H], 0.0)/self.gas.n_gas)
                    col_n2 = (col_h2/np.maximum(1e-20, (np.maximum(y[i_H2], 0.0)/self.gas.n_gas))) * (np.maximum(y[i_N2], 0.0)/self.gas.n_gas)
                    ssfactor_n2 = ss.calc_selfshielding_n2(self.ss_n2[0], self.ss_n2[1], self.ss_n2[2], self.ss_n2[3], col_h2, col_h, col_n2)
                    k = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av) * ssfactor_n2
                else:
                    k = self.reactions.a[i] * self.reactions.b[i] * safe_G0_unatt * safe_exp(-self.reactions.c[i]*self.environment.Av)

            # 40: Direct cosmic ray ionization
            elif self.reactions.itype[i] == 40:
                k = self.reactions.a[i] * (self.environment.Zeta_CR/1.35e-17)     # scaled value (see eg. UMIST22 paper Millar+22 §2.2)
            
            # 41: Cosmic ray / X-ray induced FUV reaction
            elif self.reactions.itype[i] == 41:
                k = self.reactions.a[i] * ((self.environment.Zeta_CR+self.environment.Zeta_X)/1.35e-17) * self.reactions.c[i]/(1.0-0.5) * (self.gas.t_gas/300.0)**self.reactions.b[i]

            # 42: Cosmic ray induced FUV reaction CO dissociation with self-shielding (Visser et al. 2018 & Heays et al. 2014)
            elif self.reactions.itype[i] == 42:
                co_abun = np.maximum(y[i_CO] / self.gas.n_gas, 1e-12)
                pd_eff  = 56.14 / (5.11e4 * (co_abun**0.792) + 1.0) + 4.3         # Visser et al. (2018) Eq. 2
                k = pd_eff * (self.environment.Zeta_CR+self.environment.Zeta_X)

            # 43: Cosmic ray induced FUV reaction He decay from 2.1P
            elif self.reactions.itype[i] == 43:
                k = self.reactions.a[i] * ((self.environment.Zeta_CR+self.environment.Zeta_X)*0.0107) * np.maximum(y[i_He], 0.0) / np.maximum(self.reactions.b[i] * np.maximum(y[i_H2], 0.0) + self.reactions.c[i] * np.maximum(y[i_H], 0.0), 1e-20)

            # 60: X-ray secondary ionization of H
            elif self.reactions.itype[i] == 60:
                k = self.reactions.a[i] * self.environment.Zeta_X * 0.56          # 0.56 from Staeuber/Doty code Eixion(1)/Eixion(2)
                
            # 61: X-ray secondary ionization of H2
            elif self.reactions.itype[i] == 61:
                k = self.reactions.a[i] * self.environment.Zeta_X
            
            # 62: X-ray secondary ionization of other molecules
            elif self.reactions.itype[i] == 62:
                k = self.reactions.a[i] * self.environment.Zeta_X
  
            # 70: Photoelectron production from PAHs/grains
            elif self.reactions.itype[i] == 70:
                k = self.reactions.a[i] * (self.environment.G_0+1e-4)             # +1e-4 simulates CR ionization
          
            # 71: Charge exchange/recombination with PAHs/grains
            elif self.reactions.itype[i] == 71:
                phi_pah = 0.5
                k = self.reactions.a[i] * phi_pah * ((self.gas.t_gas/100.0)**self.reactions.b[i])
            
            # 72: Charge exchange/recombination with PAHs/grains (species heavier than H)
            elif self.reactions.itype[i] == 72:
                phi_pah = 0.5
                k = self.reactions.a[i] * phi_pah * ((self.gas.t_gas/100.0)**self.reactions.b[i]) * 1/np.sqrt(self.species.mass[ir1])  
              
            # 80: Thermal desorption
            elif self.reactions.itype[i] == 80:
                if self.gas.n_gas * self.environment.dg100 < 1e3:
                    # pregrain = 1e-99
                    pregrain = (np.pi * self.dust.radius**2 * n_gr) / np.maximum(n_ice, self.dust.binding_sites*n_gr)
                else:
                    pregrain = (np.pi * self.dust.radius**2 * n_gr) / np.maximum(n_ice, self.dust.binding_sites*n_gr)
                k = 4.0 * pregrain * self.reactions.a[i] * safe_exp(-self.reactions.b[i]/self.dust.t_dust)

            # 81: Freezeout
            elif self.reactions.itype[i] == 81:
                if self.gas.n_gas * self.environment.dg100 < 1e3:
                    # prefreeze = 1e-99
                    prefreeze = np.pi * self.dust.radius**2 * n_gr * np.sqrt( 8.0 * self.parameters.k_B / (np.pi * self.parameters.m_p) )
                else:
                    prefreeze = np.pi * self.dust.radius**2 * n_gr * np.sqrt( 8.0 * self.parameters.k_B / (np.pi * self.parameters.m_p) )
                k = self.reactions.a[i] * prefreeze * np.sqrt(self.gas.t_gas/self.species.mass[ir1])

            # 90: Pumping of H2 to H2*
            elif self.reactions.itype[i] == 90:
                tinv    = (self.gas.t_gas / 1000.0) + 1.0  # downward collision rate from Le Bourlot et al. (1999)
                col_HH2 = 10**(-11.058+0.05554/tinv-2.3900/(tinv*tinv))*(np.maximum(y[i_H], 0.0)) + 10**(-11.084-3.6706/tinv-2.0230 /(tinv*tinv))*np.maximum(y[i_H2], 0.0)
                k = self.disso_H2 * 10.0 + col_HH2 * safe_exp(-30163.0 / self.gas.t_gas)      

            # 91: Radiative and collisional de-excitation of H2* to H2
            elif self.reactions.itype[i] == 91:
                tinv    = (self.gas.t_gas / 1000.0) + 1.0  # downward collision rate from Le Bourlot et al. (1999)
                col_HH2 = 10**(-11.058+0.05554/tinv-2.3900/(tinv*tinv))*(np.maximum(y[i_H], 0.0)) + 10**(-11.084-3.6706/tinv-2.0230 /(tinv*tinv))*np.maximum(y[i_H2], 0.0)
                k = self.disso_H2 * 10.0 + 2e-7 + col_HH2

            # 92: Further reactions with H2* + XX
            elif self.reactions.itype[i] == 92:
                k = self.reactions.a[i] * (self.gas.t_gas/300.0)**self.reactions.b[i] * safe_exp(-np.maximum(0.0, self.reactions.c[i]-30163.0)/self.gas.t_gas)
                
            # ERROR: Reaction ID not found
            else:
                logging.error(f'Unknown reaction type {int(self.reactions.itype[i])} for reaction {int(self.reactions.reaction_id[i])}')
                raise ValueError('ABORTED. Invalid reaction type')
                
            # ERROR: NaN or Inf rate coefficient
            if np.isnan(k) or np.isinf(k):
                logging.error(f'Invalid rate coefficient (k={k}) for reaction {self.reactions.reaction_id[i]}')
                # Instead of aborting, set a small value and continue
                k = 1e-30
                logging.warning(f'Setting reaction {self.reactions.reaction_id[i]} rate to {k}')

            self.reactions.k[i] = k

        # Final safety check for NaN or Inf values
        if np.any(np.isnan(self.reactions.k)) or np.any(np.isinf(self.reactions.k)):
            # Replace any remaining invalid values with a small number
            self.reactions.k = np.nan_to_num(self.reactions.k, nan=1e-30, posinf=1e20, neginf=1e-30)
            logging.warning("Corrected invalid rate coefficients in final check")


    ##################
    # RUN THE SOLVER #
    ##################

    def solve_network(self):
        """
        Solves the chemical reaction network using a stiff ODE solver.

        This method integrates the system of ordinary differential equations
        that govern the chemical reaction network using SciPy's solve_ivp` function 
        with the backward differentiation formula (BDF) method.

        The function tracks species abundances and reaction rates over time, updates 
        species properties at the final time, and stores the integration results for analysis.

        Includes fallback incase of integration failure, switching to Radau method, then 
        relaxed tolerances.

        Integration is aborted if not succesful after time set in self.parameters.timeout_duration
        (default = 600 seconds (10 minutes))

        Returns:
            dict: A dictionary containing:
                - `time` (numpy.ndarray): Time points of the solution.
                - `abundances` (numpy.ndarray): Species abundances over time.
                - `rates` (numpy.ndarray): Reaction rates over the evaluated time points.
                - `success` (bool): Whether the integration was successful.
                - `message` (str): Solver message (success or error).
                - `species` (list): Names of the chemical species.
                - `reaction_labels` (list): List of formatted strings describing each reaction.
        """
        
        logging.info("\n")
        helpers.log_section("SOLVING NETWORK")
        logging.info(" ◆ Starting chemical network integration...")

        # Define solver configurations to try
        solver_configs = [
            {
                "name": "BDF method",
                "method": "BDF",
                "rtol": 1e-3,
                "atol": 1e-8,
                "max_step": self.parameters.time_final/10,
                "first_step": self.parameters.time_initial/1e4
            },
            {
                "name": "Radau method",
                "method": "Radau",
                "rtol": 1e-3,
                "atol": 1e-8,
                "max_step": self.parameters.time_final/10,
                "first_step": self.parameters.time_initial/1e4
            },
            {
                "name": "Radau method (relaxed tolerance)",
                "method": "Radau",
                "rtol": 1e-2,
                "atol": 1e-8,
                "max_step": self.parameters.time_final/10,
                "first_step": self.parameters.time_initial/1e4
            }
        ]
        
        # Try each configuration until one succeeds
        for i, config in enumerate(solver_configs):
            logging.info(f"   (attempt {i+1}/{len(solver_configs)}: integration with {config['name']})\n")
            result = self._try_integration(config)
            
            if result.get('success', False):
                result['solver_config'] = config['name']    # record which configuration worked
                return result
            else:
                logging.warning(f"Integration failed with {config['name']}: {result.get('message', 'Unknown error')}")
        
        # If all attempts failed
        logging.error("All integration attempts failed.")
        return {
            'success': False,
            'message': "Failed after trying multiple solver configurations"
        }

    def _try_integration(self, config):
        """
        Helper method that attempts to solve the network with a specific configuration.
        Uses a timeout to prevent solver from running too long.
        """
        
        @contextmanager
        def timeout_handler(seconds):
            if platform.system() == 'Windows':
                # On Windows, skip timeout entirely
                yield
            else:
                # On Mac/Linux, use original signal-based timeout
                def handle_timeout(signum, frame):
                    raise TimeoutError(f"Integration timed out after {seconds} seconds")
                    
                # Set the timeout handler
                original_handler = signal.signal(signal.SIGALRM, handle_timeout)
                signal.alarm(seconds)
                
                try:
                    yield
                finally:
                    # Restore original handler and cancel alarm
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, original_handler)
        
        # Setup initial conditions with explicit float64 precision
        t0  = np.float64(self.parameters.time_initial)
        tf  = np.float64(self.parameters.time_final)
        y0  = np.array(self.species.number, dtype=np.float64)
        eps = 1e-10 * tf                                         # add small buffer to prevent floating point issues
        
        # Create logarithmically spaced time points
        t_eval = np.unique(np.concatenate([
            np.logspace(np.log10(t0), np.log10(t0 + (tf-t0)*0.1), 50),
            np.logspace(np.log10(t0 + (tf-t0)*0.1), np.log10(tf-eps), 50)
        ]))
        
        logging.debug(f"Integration timespan: {t0:.2e} to {tf:.2e} seconds")
        logging.debug(f"Number of evaluation points: {len(t_eval)}")

        # Initialise storage for rates
        self.rate_history = np.zeros((len(t_eval), self.parameters.n_reactions))
        
        # Initialize lists to store ALL rates and times during integration
        self._all_rates = []
        self._all_times = []
        
        # Initialize progress bar
        start_time = time.time()
        self._progress_bar = tqdm(ascii=" ■", ncols=70, total=100, bar_format="   [{bar}] {percentage:3.0f}% [{elapsed}<{remaining}]", colour='#6495ED')   # cornflowerblue :)

        # Wrapper for derivatives()
        def wrapper_derivatives(t, y):
            # Update progress bar
            progress = min(100, int((t - t0) / (tf - t0) * 100))
            if self._progress_bar is not None:
                self._progress_bar.n = progress
                self._progress_bar.refresh()        
                    
            # Compute rate coefficients
            self.compute_rate_coefficients(y) 

            # Calculate derivatives
            ydot = np.zeros(self.parameters.n_species, dtype=np.float64)
            ydot, rates = calculus.calculate_derivatives(y, self.reactions.k, self.idx, ydot, self.parameters.n_reactions)

            # Store time and rates
            self._all_times.append(t)
            self._all_rates.append(rates)
            return ydot

        # Wrapper for jacobian()
        def wrapper_jacobian(t, y):
            return calculus.calculate_jacobian(y, self.reactions.k, self.idx, self.parameters.n_reactions)


        ########################
        # Solve the ODE system #
        ########################

        try: 
            with timeout_handler(self.parameters.timeout_duration):
                solution = solve_ivp(
                    fun=wrapper_derivatives,
                    t_span=(t0, tf),
                    y0=y0,
                    method=config['method'],
                    t_eval=t_eval,
                    jac=wrapper_jacobian,
                    rtol=config['rtol'],
                    atol=config['atol'],
                    max_step=config['max_step'],
                    first_step=config['first_step']               
                )
            
            # Ensure progress bar reaches 100% if successful
            if solution.success and self._progress_bar is not None:
                self._progress_bar.n = 100
                self._progress_bar.refresh()
            
            end_time = time.time()
            
            if solution.success:
                if self._progress_bar is not None:
                    self._progress_bar.close()
                    self._progress_bar = None  # Set to None to prevent any future references

                logging.info(f"\n")
                logging.info(f"   ► Integration successful!")
                logging.info(f"     • Runtime:   {end_time - start_time:.2f} seconds")
                logging.info(f"     • Timesteps: {len(solution.t)}")
                logging.info(f"\n")          
                    
                # Convert lists to arrays
                all_times = np.array(self._all_times)
                all_rates = np.array(self._all_rates)
                
                # For each evaluation time, find the closest actual computed time
                for i, eval_time in enumerate(t_eval):
                    idx = np.argmin(np.abs(all_times - eval_time))
                    self.rate_history[i] = all_rates[idx]
                
                # Store results
                self.abundance_history = solution.y.T
                self.time_points = solution.t
                self.species.number = solution.y[:,-1]
                
                
                ######################
                # Summary Statistics #
                ######################

                helpers.log_section("SOLUTION ANALYSIS")

                # Most abundant species
                helpers.log_table_header_analysis(' ◆ SPECIES ABUNDANCE [X/H] (top 10)', 'Species', 'Abundance')
                top_idx = np.argsort(solution.y[:,-1])[-10:][::-1]
                for idx in top_idx:
                    helpers.log_table_row_species(self.species.name[idx], solution.y[:,-1][idx]/self.gas.n_gas)
                helpers.log_table_footer_analysis()
                
                # Most efficient reactions
                helpers.log_table_header_analysis(' ◆ DOMINANT REACTIONS', 'Reaction', 'Rate (cm⁻³ s⁻¹)')
                top_idx = np.argsort(self.rate_history[-1, :])[-10:][::-1]
                for idx in top_idx:
                    # logging.info(f"  * {self.reactions.labels[idx]:}: {self.rate_history[:,idx][-1]:.2e} cm^-3 s^-1")
                    helpers.log_table_row_reactions(self.reactions.labels[idx], self.rate_history[:,idx][-1])
                helpers.log_table_footer_analysis()
                
                # Mass conservation
                mass_i         = np.sum(solution.y[:, 0] * self.species.mass)
                mass_f         = np.sum(solution.y[:, -1] * self.species.mass)
                mass_change    = abs((mass_f-mass_i)/ mass_i) * 100
                mass_tolerance = 0.1  

                logging.info(" ◆ CONSERVATION DIAGNOSTICS")
                logging.info(f"   ► Mass conservation")
                logging.info(f"      • Initial total mass : {mass_i:.3e} amu cm^-3")
                logging.info(f"      • Final total mass   : {mass_f:.3e} amu cm^-3")
                logging.info(f"      • Difference         : {mass_change:.1e} %")
                
                if mass_change < mass_tolerance:
                    logging.info(f"      [OK] Mass conservation satisfied")
                else:
                    logging.warning(f"      [WARNING] Mass conservation violated")
                    logging.warning(f"      [WARNING] Difference exceeds tolerance by {mass_change/mass_tolerance:.1e}x")
                
                # Charge conservation
                charge_i         = np.sum(solution.y[:, 0] * self.species.charge)
                charge_f         = np.sum(solution.y[:, -1] * self.species.charge)
                charge_change    = abs(charge_f - charge_i)
                charge_tolerance = 1e-6

                logging.info(f"   ► Charge conservation")
                logging.info(f"      • Initial total charge : {charge_i:.3e}")
                logging.info(f"      • Final total charge   : {charge_f:.3e}")
                
                if charge_change < charge_tolerance:
                    logging.info(f"      [OK] Charge conservation satisfied")
                else:
                    logging.warning(f"      [WARNING] Charge conservation violated")
                    logging.warning(f"      [WARNING] Difference exceeds charge_tolerance by {charge_change/charge_tolerance:.1e}x")


                logging.info("\n")
                logging.info("┏" + ("━" * 68) + "┓")
                logging.info("┃" + "".center(68) + "┃")
                logging.info("┃" + "Simulation Completed".center(68) + "┃")
                logging.info("┃" + "".center(68) + "┃")
                logging.info("┗" + ("━" * 68) + "┛")
                logging.info("\n")

                return {
                    'time': solution.t,
                    'abundances': solution.y,
                    'rates': self.rate_history,
                    'success': True,
                    'message': solution.message,
                    'species': self.species.name,
                    'reaction_labels': self.reactions.labels
                }
            else:
                return {
                    'success': False,
                    'message': solution.message
                }
                
        except TimeoutError as e:
            # Log the timeout error
            logging.warning(f"Integration timed out after {self.parameters.timeout_duration} seconds")
            
            # Log integration progress if any
            if len(self._all_times) > 0:
                last_t = self._all_times[-1]
                progress_percent = (last_t - t0) / (tf - t0) * 100
                logging.warning(f"Integration reached {progress_percent:.1f}% completion before timeout")
                logging.warning(f"Last successful time: {last_t:.2e} seconds")
            
            # Return failure status with timeout message
            return {
                'success': False,
                'message': str(e),
                'timeout': True,
                'last_time': self._all_times[-1] if len(self._all_times) > 0 else None
            }
            
        except Exception as e:
            # Log the error
            logging.error(f"Solver failed: {str(e)}")
            
            # Log integration progress if any
            if len(self._all_times) > 0:
                last_t = self._all_times[-1]
                progress_percent = (last_t - t0) / (tf - t0) * 100
                logging.error(f"Integration reached {progress_percent:.1f}% completion")
                logging.error(f"Last successful time: {last_t:.2e} seconds")
            
            # Return failure status
            return {
                'success': False,
                'message': str(e),
                'last_time': self._all_times[-1] if len(self._all_times) > 0 else None
            }
            
        finally:
            # Ensure progress bar is closed
            if self._progress_bar is not None:
                self._progress_bar.close()
