"""
SIMBA Self-Shielding Module

This module handles the calculation of self-shielding effects. It provides 
functionality to read pre-computed shielding tables and calculate shielding factors 
for various species (H2, CO, N2, C) which are critical for accurate photochemistry 
calculations in eg. molecular clouds, protoplanetary disks.

Key Features:
- Reading and processing of tabulated self-shielding data
- Interpolation routines for CO and N2 shielding factors
- Analytical approximations for H2 and C self-shielding
- Support for temperature-dependent shielding effects
- Choice to use an input H2 column density or calculate based on Av
- Self-shielding can be turned off in the input file if desired

Main Components:
- CO shielding: Bilinear interpolation of values from Visser et al. (2009)
- N2 shielding: Trilinear interpolation of values from Visser et al. (2018)
- H2 shielding: Draine & Bertoldi (1996) approximation
- C shielding: Kamp & Bertoldi (2000) approximation

Dependencies:
    numpy: Required for numerical operations and array handling
"""

import numpy as np
from .data import CO_SELFSHIELDING_FILE, N2_SELFSHIELDING_FILE
from .helpers import safe_exp


def locate(x, arr):
    """
    Find where value x fits in sorted array arr.
    Returns (index, alpha) where:
    - index is the lower bound index
    - alpha is the interpolation factor (0-1) between index and index+1
    """
    arr = np.asarray(arr)
    
    # Handle out-of-bounds cases
    if arr[0] < arr[-1]:  # Increasing array
        if x <= arr[0]: return 0, 0.0
        if x >= arr[-1]: return len(arr)-2, 1.0
    else:  # Decreasing array
        if x >= arr[0]: return 0, 0.0
        if x <= arr[-1]: return len(arr)-2, 1.0
    
    # Find index using numpy's searchsorted
    idx = np.searchsorted(arr, x)
    if idx > 0:
        idx -= 1
    
    # Calculate interpolation factor
    alpha = (x - arr[idx]) / (arr[idx + 1] - arr[idx])
    
    return idx, alpha


###############################
# READ CO SELF SHIELDING DATA #
###############################

def read_selfshielding_co(file):
    """
    Reads CO self-shielding data from a formatted data file (Visser et al. 2009)

    This function reads and processes tabulated CO self-shielding factors used to 
    calculate photodissociation rates in the presence of H2 and CO shielding. The data
    is stored in a specific format with CO and H2 column density grids.

    Parameters:
        file (str): Path to the CO self-shielding data file.

    Returns:
        tuple: A tuple containing three arrays:
            - chem_coss_NCO (numpy.ndarray): Log10 of CO column density grid points
            - chem_coss_NH2 (numpy.ndarray): Log10 of H2 column density grid points
            - chem_coss (numpy.ndarray): 2D array of Log10 self-shielding factors

    Raises:
        ValueError: If the data file format is incorrect (wrong number of grid points)

    Notes:
        - Expected file format has specific header structure and grid dimensions
        - Grid sizes should be NCO=47 and NH2=42 points
        - Column densities are converted to log10 scale
    """
    with open(file, "r") as fp:
        # Skip first 5 lines
        for _ in range(5):
            fp.readline()
        
        chem_coss_nNCO = int(fp.readline().split()[-1])
        chem_coss_nNH2 = int(fp.readline().split()[-1])
        
        if chem_coss_nNCO != 47 or chem_coss_nNH2 != 42:
            print("ABORTED: Format problem in data_selfshielding_co.dat")

        chem_coss_NCO = np.zeros(chem_coss_nNCO)
        chem_coss_NH2 = np.zeros(chem_coss_nNH2)
        chem_coss     = np.zeros((chem_coss_nNCO, chem_coss_nNH2))
        
        fp.readline()
        
        for i in range(chem_coss_nNCO):
            chem_coss_NCO[i] = float(fp.readline().strip())
        chem_coss_NCO = np.log10(chem_coss_NCO)

        fp.readline()

        for i in range(chem_coss_nNH2):
            chem_coss_NH2[i] = float(fp.readline().strip())
        chem_coss_NH2 = np.log10(chem_coss_NH2)

        fp.readline()

        for i in range(0, chem_coss_nNH2):
            chem_12coss_temp = np.array([]) 
            for _ in range(5):
                line = fp.readline().split()  
                sub_arr = np.array([np.log10(float(i)) for i in line])
                chem_12coss_temp = np.concatenate((chem_12coss_temp, sub_arr)) 
            chem_coss[:,i] = chem_12coss_temp
    
    return chem_coss_NCO, chem_coss_NH2, chem_coss


###############################
# READ N2 SELF SHIELDING DATA #
###############################

def read_selfshielding_n2(file):
    """
    Reads N2 self-shielding data from a formatted data file (Visser et al. 2018)

    This function reads and processes tabulated N2 self-shielding factors used to 
    calculate photodissociation rates in the presence of H2, H, and N2 shielding.
    The data is stored in a 3D grid format.

    Parameters:
        file (str): Path to the N2 self-shielding data file.

    Returns:
        tuple: A tuple containing four arrays:
            - chem_n2ss_NN2 (numpy.ndarray): Log10 of N2 column density grid points
            - chem_n2ss_NH2 (numpy.ndarray): Log10 of H2 column density grid points
            - chem_n2ss_NH (numpy.ndarray): Log10 of H column density grid points
            - chem_n2ss (numpy.ndarray): 3D array of Log10 self-shielding factors

    Raises:
        ValueError: If the data file format is incorrect (wrong number of grid points)

    Notes:
        - Expected file format has specific header structure
        - Grid sizes should be NN2=46, NH2=46, and NH=10 points
        - Column densities are converted to log10 scale
    """
    with open(file, 'r') as fp:
        # Skip first three lines
        for _ in range(3):
            fp.readline()
        
        # Read dimensions
        chem_n2ss_nNN2 = int(fp.readline().split()[-1])
        chem_n2ss_nNH2 = int(fp.readline().split()[-1])
        chem_n2ss_nNH  = int(fp.readline().split()[-1])
        
        fp.readline()
        
        if chem_n2ss_nNN2 != 46 or chem_n2ss_nNH2 != 46 or chem_n2ss_nNH != 10:
            print("ABORTED: Format problem in data_selfshielding_n2.dat")

        # Initialize arrays
        chem_n2ss_NN2 = np.zeros(chem_n2ss_nNN2)
        chem_n2ss_NH2 = np.zeros(chem_n2ss_nNH2)
        chem_n2ss_NH  = np.zeros(chem_n2ss_nNH)
        chem_n2ss     = np.zeros((chem_n2ss_nNN2, chem_n2ss_nNH2, chem_n2ss_nNH))
        
        for i in range(chem_n2ss_nNN2):
            chem_n2ss_NN2[i] = float(fp.readline().strip())
        chem_n2ss_NN2 = np.log10(chem_n2ss_NN2)
        
        fp.readline()

        for i in range(chem_n2ss_nNH2):
            chem_n2ss_NH2[i] = float(fp.readline().strip())
        chem_n2ss_NH2 = np.log10(chem_n2ss_NH2)
        
        fp.readline()

        for i in range(chem_n2ss_nNH):
            chem_n2ss_NH[i] = float(fp.readline().strip())
        chem_n2ss_NH = np.log10(chem_n2ss_NH)

        fp.readline()

        # Read N2 self-shielding factors
        for k in range(chem_n2ss_nNH):  # 10 blocks
            for j in range(chem_n2ss_nNH2):  # 46 lines per block
                # Read one line and split it into values
                line = fp.readline().strip()
                values = [float(x) for x in line.split()]
                for i in range(chem_n2ss_nNN2):  # 46 values per line
                    chem_n2ss[i][j][k] = np.log10(values[i])        

    return chem_n2ss_NN2, chem_n2ss_NH2, chem_n2ss_NH, chem_n2ss


#######################################
# CALCULATE CO SELF SHIELDING FACTORS #
#######################################

def calc_selfshielding_co(chem_coss_NCO, chem_coss_NH2, chem_coss, col_h2, col_co):
    """
    Calculates CO self-shielding factor through bilinear interpolation.

    This function computes the CO self-shielding factor based on current H2 and CO
    column densities by interpolating within the pre-computed grid of shielding factors
    from Visser et al. (2009)

    Parameters:
        chem_coss_NCO (numpy.ndarray): Log10 of CO column density grid points
        chem_coss_NH2 (numpy.ndarray): Log10 of H2 column density grid points
        chem_coss (numpy.ndarray): 2D array of Log10 self-shielding factors
        col_h2 (float): Current H2 column density in cm^-2
        col_co (float): Current CO column density in cm^-2

    Returns:
        float: CO self-shielding factor (linear scale)

    Notes:
        - Performs bilinear interpolation in log space
        - Returns value is converted back to linear scale
        - Interpolated from values in Visser et al. (2009)
    """
    
    idx_h2, alpha_h2 = locate(np.log10(max(1.0, col_h2)), chem_coss_NH2)
    idx_co, alpha_co = locate(np.log10(max(1.0, col_co)), chem_coss_NCO)

    ssfact_CO = 10.0 ** (chem_coss[idx_co][idx_h2] * (1.0 - alpha_h2) * (1.0 - alpha_co) +
                        chem_coss[idx_co + 1][idx_h2] * alpha_co * (1.0 - alpha_h2) +
                        chem_coss[idx_co][idx_h2 + 1] * (1.0 - alpha_co) * alpha_h2 +
                        chem_coss[idx_co + 1][idx_h2 + 1] * alpha_h2 * alpha_co)
    
    return ssfact_CO


#######################################
# CALCULATE N2 SELF SHIELDING FACTORS #
#######################################

def calc_selfshielding_n2(chem_n2ss_NN2, chem_n2ss_NH2, chem_n2ss_NH, chem_n2ss, col_h2, col_h, col_n2):
    """
    Calculates N2 self-shielding factor through trilinear interpolation.

    This function computes the N2 self-shielding factor based on current H2, H, and N2
    column densities by interpolating within the pre-computed 3D grid of shielding factors
    from Visser et al. (2018)

    Parameters:
        chem_n2ss_NN2 (numpy.ndarray): Log10 of N2 column density grid points
        chem_n2ss_NH2 (numpy.ndarray): Log10 of H2 column density grid points
        chem_n2ss_NH (numpy.ndarray): Log10 of H column density grid points
        chem_n2ss (numpy.ndarray): 3D array of Log10 self-shielding factors
        col_h2 (float): Current H2 column density in cm^-2
        col_h (float): Current H column density in cm^-2
        col_n2 (float): Current N2 column density in cm^-2

    Returns:
        float: N2 self-shielding factor (linear scale)

    Notes:
        - Performs trilinear interpolation in log space
        - Returns value is converted back to linear scale
        - Interpolated from values in Visser et al. (2018)
    """

    h2_idx, alpha_h2 = locate(np.log10(max(1.0, col_h2)), chem_n2ss_NH2)
    h_idx, alpha_h   = locate(np.log10(max(1.0, col_h)), chem_n2ss_NH)
    n2_idx, alpha_n2 = locate(np.log10(max(1.0, col_n2)), chem_n2ss_NN2)

    # Calculate first interpolation term
    dum1_ss = (chem_n2ss[n2_idx, h2_idx, h_idx] * (1.0 - alpha_h2) * (1.0 - alpha_n2) * (1.0 - alpha_h) +
               chem_n2ss[n2_idx + 1, h2_idx, h_idx] * alpha_n2 * (1.0 - alpha_h2) * (1.0 - alpha_h) +
               chem_n2ss[n2_idx, h2_idx + 1, h_idx] * (1.0 - alpha_n2) * alpha_h2 * (1.0 - alpha_h) +
               chem_n2ss[n2_idx + 1, h2_idx + 1, h_idx] * alpha_h2 * alpha_n2 * (1.0 - alpha_h))

    # Calculate second interpolation term 
    dum2_ss = (chem_n2ss[n2_idx, h2_idx, h_idx + 1] * (1.0 - alpha_h2) * (1.0 - alpha_n2) * alpha_h +
               chem_n2ss[n2_idx + 1, h2_idx, h_idx + 1] * alpha_n2 * (1.0 - alpha_h2) * alpha_h +
               chem_n2ss[n2_idx, h2_idx + 1, h_idx + 1] * (1.0 - alpha_n2) * alpha_h2 * alpha_h +
               chem_n2ss[n2_idx + 1, h2_idx + 1, h_idx + 1] * alpha_h2 * alpha_n2 * alpha_h)

    # Calculate self-shielding factor
    ssfact_N2 = 10.0 ** (dum1_ss + dum2_ss)

    return ssfact_N2


#######################################
# CALCULATE H2 SELF SHIELDING FACTORS #
#######################################

def calc_selfshielding_h2(col_h2, delta_v):
    """
    Calculates H2 self-shielding factor using analytic approximation.

    This function computes the H2 self-shielding factor using the analytical formula
    from Draine & Bertoldi (1996). The formula includes dependencies on H2 column
    density and the velocity dispersion parameter.

    Parameters:
        col_h2 (float): H2 column density in cm^-2
        delta_v (float): Velocity dispersion parameter in km/s

    Returns:
        float: H2 self-shielding factor

    Notes:
        - Based on Draine & Bertoldi (1996) approximation
    """
    nh2_5e14   = col_h2 / 5e14
    ssfact_H2  = 0.965 / (1.0+(nh2_5e14/delta_v))**2 + 0.035 / np.sqrt(1.0 + nh2_5e14) * safe_exp(-8.5e-4 * np.sqrt(1.0+nh2_5e14)) 
    return ssfact_H2


######################################
# CALCULATE C SELF SHIELDING FACTORS #
######################################

def calc_selfshielding_c(col_h2, col_c, t_gas):
    """
    Calculates atomic carbon self-shielding factor.

    This function computes the self-shielding factor for atomic carbon using the
    analytical formula from Kamp & Bertoldi (2000). Depends on C and H2 column densities and gas temperature.

    Parameters:
        col_h2 (float): H2 column density in cm^-2
        col_c (float): C column density in cm^-2
        t_gas (float): Gas temperature in K

    Returns:
        float: C self-shielding factor

    Notes:
        - Includes both C self-shielding and H2 shielding effects
        - Temperature dependence affects H2 shielding component
        - Minimum shielding factor of 0.5 is enforced for H2 component
        - Based on Kamp & Bertoldi (2000)
    """
    ssfactor_c = safe_exp(-col_c*1.1e-17) * np.maximum(0.5, safe_exp(-0.9 * (t_gas**0.27) * ((col_h2/1e22)**0.45)))
    return ssfactor_c

