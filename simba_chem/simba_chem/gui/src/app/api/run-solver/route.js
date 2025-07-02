// /src/app/api/run-solver/route.js
const { spawn } = require('child_process');

export async function POST(req) {
  const params = await req.json();
  const fs = require('fs');
  const path = require('path');
  
  const lukenetPath = path.join(process.cwd(), '..', 'simba_input.dat');
  const lukenetDir = path.join(process.cwd(), '..');
  
  // Helper function to format values to scientific notation with 10 decimal places
  const formatScientific = (value) => {
    if (value === undefined || value === '') {
      return value; // Return original value if empty or undefined
    }
    
    try {
      const num = parseFloat(value);
      return num.toExponential(10);
    } catch (error) {
      console.error(`Error parsing value: ${value}`);
      return value; // Return original value in case of error
    }
  };
  
  // Helper function to convert ON/OFF to Python boolean
  const convertToPythonBool = (value) => {
    return value === 'ON' ? 'True' : 'False';
  };
  
  // Format the specified parameters to scientific notation with 10 decimal places
  const n_gasFormatted = formatScientific(params.n_gas);
  const n_dustFormatted = formatScientific(params.n_dust);
  const Zeta_XFormatted = formatScientific(params.Zeta_X);
  const Zeta_CRFormatted = formatScientific(params.Zeta_CR);
  const h2_colFormatted = formatScientific(params.h2_col);
  
  // Convert self_shielding and column to Python booleans
  const selfShieldingFormatted = convertToPythonBool(params.self_shielding);
  const columnFormatted = convertToPythonBool(params.column);
  
  // Create input file content
  const pythonContent = `
n_gas          = ${n_gasFormatted}
n_dust         = ${n_dustFormatted}
t_gas          = ${params.t_gas}
t_dust         = ${params.t_dust}
gtd            = ${params.gtd}
Av             = ${params.Av}
G_0            = ${params.G_0}
Zeta_X         = ${Zeta_XFormatted}
h2_col         = ${h2_colFormatted}
Zeta_CR        = ${Zeta_CRFormatted}
self_shielding = ${selfShieldingFormatted}
column         = ${columnFormatted}
pah_ism        = ${params.pah_ism}
t_chem         = ${params.t_chem}
network        = ${params.network}
`;
  // Write the input file
  fs.writeFileSync(lukenetPath, pythonContent);
  
  // Initialize global progress
  global.solverProgress = 0;

  return new Promise((resolve) => {
    // Modified Python code to capture and return results
    const python = spawn('python', ['-c', `
import simba_chem as simba
import json
import numpy as np

# Create and run network
network = simba.Simba()
network.init_simba("${lukenetPath.replace(/\\/g, '\\\\')}")
result = network.solve_network()

# Prepare results for JSON serialization
output_data = {
    'time': (result['time']/network.parameters.yr_sec ).tolist(),
    'abundances': (result['abundances'] / network.gas.n_gas).tolist(),
    'rates': result['rates'].tolist(),
    'species_names': network.species.name,
    'reaction_labels': network.reactions.labels
}

# Write results to a temporary file
with open('lukenet_results.json', 'w') as f:
    json.dump(output_data, f)
    `], { cwd: lukenetDir });

    // Handle progress updates from stderr (tqdm outputs to stderr)
    python.stderr.on('data', (data) => {
      const output = data.toString();
      const match = output.match(/(\d+)%/);
      if (match) {
        global.solverProgress = parseInt(match[1]);
        // console.log('Progress:', global.solverProgress);
      }
    });

    // Handle process completion
    python.on('close', async (code) => {
      if (code === 0) {
        try {
          // Read results from temporary file
          const resultsPath = path.join(lukenetDir, 'lukenet_results.json');
          const results = JSON.parse(fs.readFileSync(resultsPath, 'utf8'));
          
          // Clean up temporary file
          fs.unlinkSync(resultsPath);
          
          resolve(Response.json({ 
            success: true,
            progress: 100,
            data: results
          }));
        } catch (error) {
          resolve(Response.json({ 
            success: false,
            error: 'Failed to read results',
            progress: global.solverProgress 
          }));
        }
      } else {
        resolve(Response.json({ 
          success: false,
          error: 'Solver failed',
          progress: global.solverProgress 
        }));
      }
    });
  });
}