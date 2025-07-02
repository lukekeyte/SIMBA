"""
SIMBA: Astrophysical Chemical Network Solver

A comprehensive solver for chemical reaction networks in astrophysical environments.
Designed to model and simulate complex chemical processes in various cosmic settings
such as the ISM, molecular clouds, and protoplanetary disks.

Author: Luke Keyte
"""

from .core import Simba
from .model_classes import Elements, Species, Gas, Dust, Environment, Reactions, Parameters
from .analysis import Analysis
from .helpers import create_input, create_network
from .gui_server import launch_gui, launch_gui_dev


__version__ = "1.0.2"
__author__ = "Luke Keyte"

# Export main classes
__all__ = [
    'Simba',
    'Elements',
    'Species',
    'Gas',
    'Dust',
    'Environment',
    'Reactions',
    'Parameters',
    'Analysis',
    'launch_gui',
    'launch_gui_dev'
]