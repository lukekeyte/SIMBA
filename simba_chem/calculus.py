"""
SIMBA Calculus Module

This module implements the core numerical routines for solving chemical reaction networks
in the SIMBA solver. It provides optimized functions for calculating reaction rates,
species derivatives, and the system Jacobian matrix using Numba acceleration.

Key Features:
- Computation of chemical reaction rates and species derivatives
- Generation of sparse Jacobian matrices for ODE integration
- Performance optimization using Numba JIT compilation
- Support for reactions with up to 3 reactants and 5 products
- Efficient handling of sparse chemical networks

Main Components:
- calculate_derivatives: Computes formation/destruction rates and net changes
- calculate_jacobian_dense: Generates full Jacobian matrix
- calculate_jacobian: Converts dense Jacobian to sparse format for solver

Performance Notes:
- Uses Numba JIT compilation for compute-intensive operations
- Implements sparse matrix operations for memory efficiency

Dependencies:
   numpy: Required for numerical operations
   numba: Used for JIT compilation
   scipy.sparse: Used for sparse matrix operations
"""


import numpy as np
import numba
from scipy.sparse import lil_matrix
        
        
#########################
# CALCULATE DERIVATIVES #
#########################

@numba.jit(nopython=True, cache=True)
def calculate_derivatives(y, k, idx, ydot, nr):
    """
    Calculates the derivatives and reaction rates for a chemical network.

    This function computes the formation and destruction rates of chemical 
    species based on a set of reaction rate coefficients, reactant-product 
    relationships, and species concentrations. It uses the Numba JIT compiler 
    for performance optimization.

    Args:
        y (numpy.ndarray): Array of concentrations of chemical species.
        k (numpy.ndarray): Array of reaction rate coefficients.
        idx (numpy.ndarray): Array of indices defining reactants (educts) and 
            products for each reaction. The format is:
                idx[i,0:3] - indices of up to 3 reactants (educts).
                idx[i,3:8] - indices of up to 5 products.
        ydot (numpy.ndarray): Placeholder array for derivatives of concentrations.
        nr (int): Total number of reactions.

    Returns:
        tuple:
            numpy.ndarray: Array of net changes (formation - destruction) for 
                each species.
            numpy.ndarray: Array of reaction rates for each reaction.

    Notes:
        - The input `idx` is expected to have negative values for missing reactants 
          or products.
        - The function assumes that `y` and `k` are appropriately sized to 
          correspond to the number of reactions and species.
    """
    formation = np.zeros_like(ydot)
    destruction = np.zeros_like(ydot)
    rates = np.zeros(nr)
    
    # Loop over all reactions
    for i in range(nr):
        # Get indices for educts and products
        ir1 = idx[i,0]
        ir2 = idx[i,1]
        ir3 = idx[i,2]
        ip1 = idx[i,3]
        ip2 = idx[i,4]
        ip3 = idx[i,5]
        ip4 = idx[i,6]
        ip5 = idx[i,7]
        
        # Calculate reaction term
        term = k[i]
        
        # Multiply by concentrations of educts
        if ir1 >= 0:
            term *= y[ir1]
        if ir2 >= 0:
            term *= y[ir2]
        if ir3 >= 0:
            term *= y[ir3]
        
        # Save for the rate history 
        rates[i] = term
            
        # Add to formation rates for products
        if ip1 >= 0:
            formation[ip1] += term
        if ip2 >= 0:
            formation[ip2] += term
        if ip3 >= 0:
            formation[ip3] += term
        if ip4 >= 0:
            formation[ip4] += term
        if ip5 >= 0:
            formation[ip5] += term
            
        # Add to destruction rates for educts
        if ir1 >= 0:
            destruction[ir1] += term
        if ir2 >= 0:
            destruction[ir2] += term
        if ir3 >= 0:
            destruction[ir3] += term

    return (formation - destruction), rates


######################
# CALCULATE JACOBIAN #
######################

@numba.jit(nopython=True, cache=True)
def calculate_jacobian_dense(y, k, idx, nr):
    """
    Computes the dense Jacobian matrix for a chemical reaction network.

    This function calculates the Jacobian matrix, representing the partial derivatives 
    of reaction rates with respect to species concentrations. It accounts for 
    single-reactant and two-reactant reactions, optimizing the computation using 
    Numba's JIT compilation.

    Args:
        y (numpy.ndarray): Array of concentrations of chemical species.
        k (numpy.ndarray): Array of reaction rate coefficients.
        idx (numpy.ndarray): Array of indices defining reactants (educts) and 
            products for each reaction. The format is:
                idx[i, 0:3] - indices of up to 3 reactants (educts).
                idx[i, 3:8] - indices of up to 5 products.
        nr (int): Total number of reactions.

    Returns:
        numpy.ndarray: Dense Jacobian matrix of shape `(ns, ns)`, where `ns` 
            is the number of chemical species. Each entry `jac[j, i]` represents 
            the partial derivative of the rate of change of species `j` with 
            respect to species `i`.

    Notes:
        - Single-reactant reactions only affect the row and column of that reactant.
        - Two-reactant reactions account for the interaction between two species 
          and their contributions to products.
        - Unused indices in `idx` are expected to have a value of `-1`.
    """
    ns = len(y)
    jac = np.zeros((ns, ns))  # Use a dense NumPy array for calculation
    
    for i in range(nr):
        ir1 = idx[i, 0]
        ir2 = idx[i, 1]
        ir3 = idx[i, 2]
        ip1 = idx[i, 3]
        ip2 = idx[i, 4]
        ip3 = idx[i, 5]
        ip4 = idx[i, 6]
        ip5 = idx[i, 7]

        # Single reactant reactions
        if ir2 == -1 and ir3 == -1:
            if ir1 >= 0:
                # Effect on reactant
                jac[ir1, ir1] -= k[i]
                # Effect on products
                if ip1 >= 0:
                    jac[ip1, ir1] += k[i]
                if ip2 >= 0:
                    jac[ip2, ir1] += k[i]
                if ip3 >= 0:
                    jac[ip3, ir1] += k[i]
                if ip4 >= 0:
                    jac[ip4, ir1] += k[i]
                if ip5 >= 0:
                    jac[ip5, ir1] += k[i]

        # Two reactant reactions
        elif ir3 == -1:
            if ir1 >= 0 and ir2 >= 0:
                # Effect on first reactant
                term1 = k[i] * y[ir2]
                jac[ir1, ir1] -= term1
                jac[ir2, ir1] -= term1
                if ip1 >= 0:
                    jac[ip1, ir1] += term1
                if ip2 >= 0:
                    jac[ip2, ir1] += term1
                if ip3 >= 0:
                    jac[ip3, ir1] += term1
                if ip4 >= 0:
                    jac[ip4, ir1] += term1
                if ip5 >= 0:
                    jac[ip5, ir1] += term1

                # Effect on second reactant
                term2 = k[i] * y[ir1]
                jac[ir1, ir2] -= term2
                jac[ir2, ir2] -= term2
                if ip1 >= 0:
                    jac[ip1, ir2] += term2
                if ip2 >= 0:
                    jac[ip2, ir2] += term2
                if ip3 >= 0:
                    jac[ip3, ir2] += term2
                if ip4 >= 0:
                    jac[ip4, ir2] += term2
                if ip5 >= 0:
                    jac[ip5, ir2] += term2

    return jac


def calculate_jacobian(y, k, idx, nr):
    """
    Convert the dense Jacobian to sparse format after calculation
    """
    jac_dense = calculate_jacobian_dense(y, k, idx, nr)
    jac_sparse = lil_matrix(jac_dense)  
    return jac_sparse