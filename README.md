# SIMBA: Solver for Inferring Molecular aBundances in Astrophysical environments

SIMBA is a comprehensive solver for chemical reaction networks in astrophysical environments. It is designed to model and simulate complex chemical processes in various cosmic settings such as the interstellar medium (ISM), molecular clouds, and protoplanetary disks.

## Features

- Initialization of chemical species, reactions, and environmental parameters
- Efficient solving of stiff ODEs representing chemical reactions
- Support for various reaction types including:
  - Gas-phase reactions
  - Grain-surface chemistry
  - Photochemistry
- Integration of self-shielding factors for specific molecules (H2, CO, N2, C)
- Optimization using Numba JIT compilation for performance-critical functions
- Comprehensive logging and progress tracking

## Installation

You can install SIMBA using pip:

```bash
pip install simba_chem
```

## Quick Start

1. First, create an input file for your simulation:
```python
import simba_chem as simba

# Create default input file
simba.create_input('my_input.dat')
```

&nbsp;&nbsp;&nbsp;&nbsp; You will also need to have a correctly formatted chemical network file available. A template can be generated using:
```python
simba.create_network("directory/to/save/network/") 
```

&nbsp;&nbsp;&nbsp;&nbsp; Don't forget to specify the path to the chemical network file in your input file!


2. Modify the input parameters in `my_input.py` according to your needs

3. Run your simulation:
```python
# Initialize the network
network = simba.Simba()
network.init_simba('my_input.dat')

# Solve the network
result = network.solve_network()
```





## Dependencies

- numpy
- scipy
- matplotlib
- numba
- tqdm
- pandas

## Model Components

SIMBA consists of several key components:
- Elements: Handling of chemical elements
- Species: Management of atomic and molecular species
- Gas: Gas phase parameters and properties
- Dust: Dust grain properties and interactions
- Environment: Environmental conditions (UV field, cosmic rays, etc.)
- Reactions: Chemical reaction network and rates
- Parameters: System parameters and constants

## Example Usage

```python
import simba_chem as simba
import matplotlib.pyplot as plt

# Initialize and run
network = simba.Simba()
network.init_simba('my_input.dat')
result = network.solve_network()

# Plot results
time = result['time']
abundances = result['abundances']
plt.loglog(time, abundances[0,:])  # Plot first species abundance
plt.xlabel('Time (s)')
plt.ylabel('Number Density (cm^-3)')
plt.show()
```


## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Luke Keyte

## Contributing

Contributions are welcome! Please feel free to get in touch: l.keyte@qmul.ac.uk
