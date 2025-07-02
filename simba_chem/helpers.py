"""
SIMBA Helpers Module

This module provides utility functions for the SIMBA chemical network solver.
It includes functionalities for reading input files, creating readable 
reaction labels, formatting logger messages, and extracting/displaying DALI model parameters.

Key Features:
- Input file parsing and creation
- Chemical network file parsing and creation
- Reaction label and message formatting
- Logger formatting
- DALI model cell parameter extraction and display

Main Components:
- Chemical network reader with support for species and reaction parsing
- Reaction label generator for human-readable output
- DALI model parameter extractor and formatter

Dependencies:
    numpy: Required for numerical operations
"""

import numpy as np
import logging
import os
import shutil
try:
    # Python 3.9+
    from importlib.resources import files
except ImportError:
    # Python 3.7-3.8 fallback
    from importlib_resources import files


###################
# READ INPUT FILE #
###################

def read_input_file(file_path):
    """
    Reads a configuration file with key-value pairs and returns a dictionary.

    Args:
        file_path (str): Path to the input file.

    Returns:
        dict: A dictionary containing the parsed key-value pairs.
    """
    data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespace
            line = line.strip()
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
            
            # Split key and value at the equals sign
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                # Convert value to the appropriate type
                if value.lower() in ('true', 'false'):
                    value = value.lower() == 'true'  # Convert to boolean
                else:
                    try:
                        # Try to convert to a float or int
                        if '.' in value or 'e' in value.lower():
                            value = float(value)
                        else:
                            value = int(value)
                    except ValueError:
                        # Leave as string if conversion fails
                        value = value.strip('"').strip("'")

                data[key] = value

    return data


##############################
# READ CHEMICAL NETWORK FILE #
##############################

def read_chemnet(network_file):
    """
    Reads and parses a DALI-format chemical network file.

    This function reads a formatted chemical network file containing species
    information and reaction data. The file must follow a specific format with
    sections for elements, species, and reactions.

    Parameters:
        network_file (str): Path to the chemical network file

    Returns:
        tuple: A tuple containing network components:
            - n_elements (int): Number of elements
            - elements_name (list): Element names
            - n_species (int): Number of species
            - species_name (list): Species names
            - species_abu (numpy.ndarray): Initial abundances
            - species_mass (numpy.ndarray): Species masses
            - species_charge (numpy.ndarray): Species charges
            - n_reactions (int): Number of reactions
            - reactions_educts (list): Reactants for each reaction
            - reactions_products (list): Products for each reaction
            - reactions_reaction_id (numpy.ndarray): Reaction IDs
            - reactions_itype (numpy.ndarray): Reaction type identifiers
            - reactions_a (numpy.ndarray): Rate coefficient parameter a
            - reactions_b (numpy.ndarray): Rate coefficient parameter b
            - reactions_c (numpy.ndarray): Rate coefficient parameter c
            - reactions_temp_min (numpy.ndarray): Minimum valid temperatures
            - reactions_temp_max (numpy.ndarray): Maximum valid temperatures
            - reactions_pd_data (list): Photodissociation data (NOT USED BY THIS CODE)

    Notes:
        - File format must follow specific structure with header sections
        - File format is the same as that used by the DALI thermochemical code
        - Species data includes 4 columns: name, abundance, mass, charge
        - Reaction data includes educts, products, and rate parameters
    """

    with open(network_file, 'r') as file:
        for _ in range(4):
            file.readline()
        
        n_elements = int(file.readline().strip())

        file.readline() 
        
        elements_name = []
        for _ in range(n_elements):
            elem = file.readline().strip()
            elements_name.append(elem)
        
        for _ in range(4):
            file.readline()
                    
        n_species = int(file.readline().strip())

        file.readline()
        file.readline()
        
        # Initialize species arrays
        species_name   = []
        species_abu    = np.zeros(n_species)
        species_mass   = np.zeros(n_species)
        species_charge = np.zeros(n_species)
        
        # Extract species
        for i in range(n_species):
            line = file.readline().split()
            species_name.append(line[0])
            species_abu[i] = float(line[1])
            species_mass[i] = float(line[2])
            species_charge[i] = float(line[3])
            
        for _ in range(4):
            file.readline()
        
        n_reactions = int(file.readline().strip())
        
        file.readline()
        file.readline()
        
        # Initialize reactions arrays
        reactions_educts      = []
        reactions_products    = []
        reactions_reaction_id = np.zeros(n_reactions)
        reactions_itype       = np.zeros(n_reactions)
        reactions_a           = np.zeros(n_reactions)
        reactions_b           = np.zeros(n_reactions)
        reactions_c           = np.zeros(n_reactions)
        reactions_temp_min    = np.zeros(n_reactions)
        reactions_temp_max    = np.zeros(n_reactions)
        reactions_pd_data     = []
        
        # Extract reactions
        for i in range(0, n_reactions):
            line     = file.readline()
            educts   = [line[i:i+14].strip() for i in range(0, 3*14, 14)]
            products = [line[i:i+14].strip() for i in range(3*14, 8*14, 14)]
            
            reactions_educts.append(educts)
            reactions_products.append(products)

            remainder = line[8*14:].split()
            
            reactions_reaction_id[i] = remainder[0]
            reactions_itype[i] = remainder[1]
            reactions_a[i] = remainder[2]
            reactions_b[i] = remainder[3]
            reactions_c[i] = remainder[4]
            reactions_temp_min[i] = remainder[5]
            reactions_temp_max[i] = remainder[6]
            reactions_pd_data.append(remainder[7])
    
    return (n_elements, elements_name, n_species, species_name, species_abu, species_mass, species_charge, n_reactions, 
            reactions_educts, reactions_products, reactions_reaction_id, reactions_itype, reactions_a, reactions_b, reactions_c, 
            reactions_temp_min, reactions_temp_max, reactions_pd_data)


##########################
# CREATE REACTION LABELS #
##########################

def create_reaction_labels(n_reactions, educts, products):
    """
    Creates human-readable labels for chemical reactions.

    This function generates formatted strings representing chemical reactions,
    showing reactants and products connected by an arrow.

    Parameters:
        n_reactions (int): Number of reactions to process
        educts (list): List of lists containing reactant species names
        products (list): List of lists containing product species names

    Returns:
        list: List of formatted strings, each representing a reaction in the
              format "A + B -> C + D"

    Notes:
        - Empty species strings are filtered out
    """
    labels = []
    for i in range(0, n_reactions):
        educt_str   = ' + '.join(j for j in educts[i] if j != '')
        product_str = ' + '.join(j for j in products[i] if j != '')
        labels.append(f"{educt_str} -> {product_str}")
    return labels


##########################
# FORMAT LOGGER MESSAGES #
##########################

def log_section(title):
    logging.info("▓" * 24 + title.upper().center(22) + "▓" * 24 +"\n")


def log_param(name, value, unit=""):
    name_pad  = 25  
    value_pad = 15 
    if unit:
        logging.info(f"{name:<{name_pad}}{value:>{value_pad}} {unit}")
    else:
        logging.info(f"{name:<{name_pad}}{value:>{value_pad}}")

def format_scientific(number):
    """Format a number in scientific notation with proper symbols."""
    # For small numbers close to zero, just return the formatted number
    if abs(number) < 0.1 and abs(number) > 0.0001:
        return f"{number:.4f}"
    
    # For numbers that don't need scientific notation
    if abs(number) >= 0.1 and abs(number) < 1000:
        # Format with appropriate decimal places
        if abs(number) >= 100:
            return f"{number:.1f}"
        elif abs(number) >= 10:
            return f"{number:.2f}"
        else:
            return f"{number:.3f}"
    
    # Convert to scientific notation
    sci_notation = f"{number:.2e}"
    parts = sci_notation.split('e')
    mantissa = float(parts[0])
    exponent = int(parts[1])
    
    # Format exponent with superscripts
    exponent_map = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', 
                    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', 
                    '-': '⁻', '+': ''}
    exponent_str = ''.join(exponent_map[char] for char in str(exponent))
    
    # Return formatted string
    return f"{mantissa:.2f} × 10{exponent_str}"
    
def log_table_header(title):
    """Log a table section header."""
    logging.info(f" {title}")
    logging.info(" ┌────────────────────┬───────────────┬──────────┐")
    logging.info(" │ Parameter          │ Value         │ Unit     │")
    logging.info(" ├────────────────────┼───────────────┼──────────┤")

def log_table_footer():
    """Log the table footer."""
    logging.info(" └────────────────────┴───────────────┴──────────┘\n")

def log_table_row(parameter, value, unit=""):
    """Format and log a parameter as a table row."""
    formatted_value = format_scientific(value) if isinstance(value, float) else value
    # Handle boolean values
    if isinstance(value, bool):
        formatted_value = str(value)
    
    if unit:
        # Convert standard units to proper Unicode
        unit = (unit.replace("^-1", "⁻¹")
                   .replace("^-2", "⁻²")
                   .replace("^-3", "⁻³")
                   .replace("^2", "²")
                   .replace("^3", "³"))
        logging.info(f" │ {parameter:<18} │ {formatted_value:<13} │ {unit:<8} │")
    else:
        logging.info(f" │ {parameter:<18} │ {formatted_value:<13} │          │")

def log_table_header_analysis(title, col1, col2):
    """Log a table section header."""
    logging.info(f" {title}")
    logging.info(" ┌────────────────────────┬──────────────────────┐")
    logging.info(f" │ {col1:<12}           │ {col2:<15}      │")
    logging.info(" ├────────────────────────┼──────────────────────┤")

def log_table_row_species(species, abundance):
    formatted_value = format_scientific(abundance) if isinstance(abundance, float) else abundance
    logging.info(f" | {species:<8}               |  {formatted_value:<18}  |")

def log_table_row_reactions(reaction, rate):
    formatted_value = format_scientific(rate) if isinstance(rate, float) else rate
    logging.info(f" | {reaction:<23}|  {formatted_value:<13}       |")

def log_table_footer_analysis():
    """Log the table footer."""
    logging.info(" └────────────────────────┴──────────────────────┘\n")


###################
# PRINT DALI CELL #
###################

def print_dali_cell(dali_model_outdat_path, r, z):
    """
    Prints physical parameters from a specified cell in a DALI model.

    This function extracts and displays key physical parameters from a DALI
    model output file ('out.dat') for a specific grid cell location.

    Outputs are printed such that they can be copied directly into SIMBA input file.

    Parameters:
        dali_model_outdat_path (str): Path to the DALI out.dat file
        r (int): Radial index of the cell
        z (int): Vertical index of the cell

    Prints:
        Formatted output of cell parameters including:
            - Gas density and temperature
            - Dust density and temperature
            - Gas-to-dust ratio
            - Visual extinction (Av)
            - UV field strength (in Draine units)
            - X-ray ionization rate
            - H2 column density
    """
    
    def read_outdat(fname):
        '''
        Read DALI out.dat file
        Argument:
        - fname (full path to out.dat file)
        Notes:
        - this function written by Simon Bruderer
        '''

        data={}

        # read file
        lines=open(fname,"r").readlines()
        
        # read header
        data['n_r'], data['n_z']=map(int,lines[1].split())	
        head=lines[3].split()	
        data['n_col']=len(head)

        # set up array   
        for h in head:
            data[h]=[[0.0 for i_z in range(data['n_z'])] for i_r in range(data['n_r'])]
        
        # read data
        lin=5
        for i_r in range(data['n_r']):
            for i_z in range(data['n_z']):
                splt = list(map(float, lines[lin].split()))		
                for i_c in range(data['n_col']):
                    data[head[i_c]][i_r][i_z]=splt[i_c]
                lin+=1
                
        # find wavelength
        wave=[]
        for d in data.keys():
            if d[0:2]=="J=":			
                wave.append(float(d[2:]))
        wave.sort()
        data['wave_grid']=wave
        
        
        # define midpoint of the cell (used for interpolation)
        data['r_mid']=[[0.0 for i_z in range(data['n_z'])] for i_r in range(data['n_r'])]
        data['z_mid']=[[0.0 for i_z in range(data['n_z'])] for i_r in range(data['n_r'])]
        
        rmax=0.0
        zmax=0.0	
        for i_r in range(data['n_r']):
            for i_z in range(data['n_z']):
                data['r_mid'][i_r][i_z]=0.5*(data['ra'][i_r][i_z]+data['rb'][i_r][i_z])
                data['z_mid'][i_r][i_z]=0.5*(data['za'][i_r][i_z]+data['zb'][i_r][i_z])	
                zmax=max(zmax,data['zb'][i_r][i_z])
                rmax=max(rmax,data['rb'][i_r][i_z])				
        data['r_max']=rmax
        data['z_max']=zmax					
        return data
    
    dali = read_outdat(dali_model_outdat_path)

    rr = np.array(dali['ra'])[r, z]
    zz = np.array(dali['za'])[r, z]
    ngas = np.array(dali['n_gas'])[r, z]
    ndust = np.array(dali['n_dust'])[r, z]
    tgas = np.array(dali['t_gas'])[r, z]
    tdust = np.array(dali['t_dust'])[r, z]
    g0 = np.array(dali['G_0'])[r, z]
    av = np.array(dali['A_V'])[r, z]
    zetax = np.array(dali['Zeta_X'])[r, z]
    dg100 = np.array(dali['dg100'])[r, z]
    gtd = ngas/ndust
    pah_ism = np.array(dali['pah_abu'])[r, z]
    h2_col = np.array(dali['H2_COL'])[r, z]
    
    
    print('DALI MODEL')
    print('----------')
    print(f'n_gas          = {ngas:.10e}')
    print(f'n_dust         = {ndust:.10e}')
    print(f't_gas          = {tgas:.1f}')
    print(f't_dust         = {tdust:.1f}')
    print(f'gtd            = {gtd:.10e}')
    print(f'Av             = {av:.10e}')
    print(f'G_0            = {g0:.10e}')
    print(f'Zeta_X         = {zetax:.10e}')
    print(f'h2_col         = {h2_col:.10e}')
    
    
#####################
# CREATE INPUT FILE #
#####################

def create_input(path_to_output):
    """
    Creates a template SIMBA input file at the specified path with the content of the template file.
    
    Args:
        path_to_output (str): The path where the new Python file should be created.
    """
    file_content = '''
#############################################################################
# SIMBA Input File                                                          #
#                                                                           #
# Physical Parameters:                                                      #
# - Gas Properties:                                                         #
#    n_gas: Gas number density (cm^-3)                                      #
#    t_gas: Gas temperature (K)                                             #
#    gtd: Gas-to-dust ratio                                                 #
#    h2_col: H column density (ie. H+H2) (cm^-2)                            #
#                                                                           #
# - Dust Properties:                                                        #
#    n_dust: Dust number density (cm^-3)                                    #
#    t_dust: Dust temperature (K)                                           #
#                                                                           #
# - Radiation Field:                                                        #
#    Av: Visual extinction (mag)                                            #
#    G_0: Local FUV field strength                                          #
#     (normalised to Draine field ~2.7e-3 erg/s/cm^2 between 911-2067A)     #
#    Zeta_X: X-ray ionization rate (s^-1)                                   #
#    Zeta_CR: Cosmic ray ionization rate (s^-1)                             #
#                                                                           #
# - Chemistry Settings:                                                     #
#    pah_ism: PAH abundance relative to ISM                                 #
#    t_chem: Chemical evolution time (years)                                #
#    self_shielding: Enable self-shielding (H2, N2, CO, C)                  #
#    column: User-specified H2 column for self-shielding                    #
#        (otherwise approximate based on Av)                                #
#    network: Full path to chemical network file (str)                      #
#############################################################################


n_gas          = 4.6780000000e+09
n_dust         = 9.8240000000e+06
t_gas          = 56.4
t_dust         = 56.4
gtd            = 4.7618078176e+02
Av             = 1.1930000000e+01
G_0            = 2.6390000000e-11
Zeta_X         = 6.6870000000e-16
h2_col         = 6.3010000000e+22
Zeta_CR        = 5e-17
self_shielding = True
column         = True
pah_ism        = 0.1
t_chem         = 1e6
network        = '/path/to/network/network.dat'
'''

    with open(path_to_output, 'w') as file:
        file.write(file_content)


################################
# CREATE CHEMICAL NETWORK FILE #
################################

def create_network(output_dir):
    """
    Copy the chemnet_example_keyte2023.dat file to a user-specified directory.
    
    Args:
    output_dir : str
        Path to the directory where the data file should be copied
    """
    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the path to the data file using modern importlib.resources
    try:
        # Try the modern way first (Python 3.9+)
        data_files = files('simba_chem') / 'data'
        data_file = data_files / 'chemnet_example_keyte2023.dat'
        
        # For Python 3.9+, we can use as_posix() to get the path
        if hasattr(data_file, 'as_posix'):
            data_file_path = str(data_file)
        else:
            # Fallback for older versions - use context manager
            with files('simba_chem').joinpath('data/chemnet_example_keyte2023.dat').open('rb') as f:
                data_file_path = f.name
                
    except (ImportError, AttributeError):
        # Fallback to pkg_resources if importlib.resources fails
        import pkg_resources
        data_file_path = pkg_resources.resource_filename('simba_chem', 'data/chemnet_example_keyte2023.dat')
    
    # Define output path
    output_path = os.path.join(output_dir, 'chemnet_example_keyte2023.dat')
    
    # Copy the file
    if hasattr(data_file, 'read_bytes'):
        # Modern approach - read data and write to output
        data_content = data_file.read_bytes()
        with open(output_path, 'wb') as f:
            f.write(data_content)
    else:
        # Fallback approach
        shutil.copy2(data_file_path, output_path)
    
    print(f"Network data copied to: {output_path}")
    return output_path


#####################
# SAVE EXPONENTIALS #
#####################

def safe_exp(value):
    """Safely compute exponential, clipping to prevent underflow/overflow."""
    return np.exp(np.clip(value, -700.0, 700.0))