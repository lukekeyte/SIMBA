"""
Data files for SIMBA including self-shielding data for various molecules.
"""

import os

# Get the directory where the data files are stored
DATA_DIR = os.path.dirname(os.path.abspath(__file__))

# Define paths to specific data files
CO_SELFSHIELDING_FILE = os.path.join(DATA_DIR, 'data_selfshielding_co.dat')
N2_SELFSHIELDING_FILE = os.path.join(DATA_DIR, 'data_selfshielding_n2.dat')

# Export the file paths
__all__ = ['CO_SELFSHIELDING_FILE', 'N2_SELFSHIELDING_FILE']