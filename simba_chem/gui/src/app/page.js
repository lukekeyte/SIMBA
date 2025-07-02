'use client';

import { useState, useEffect, useMemo, useRef } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '../components/ui/card';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import Logo from '/public/lukenet_logo.png'
import { Settings, Play, Save, Upload, RefreshCw, LineChart as ChartIcon, ChevronDown, ChevronUp, Check } from 'lucide-react';
import PathwayDiagram from '../components/PathwayDiagram';
import html2canvas from 'html2canvas';

// SpeciesSelector Component
const SpeciesSelector = ({ 
  species = [], 
  selectedSpecies = [], 
  onSpeciesChange,
  disabled = false
}) => {
  const [isOpen, setIsOpen] = useState(false);
  const [searchTerm, setSearchTerm] = useState("");

  const filteredSpecies = species.filter(name =>
    name.toLowerCase().includes(searchTerm.toLowerCase())
  );

  const handleSpeciesToggle = (speciesName) => {
    const newSelection = selectedSpecies.includes(speciesName)
      ? selectedSpecies.filter(name => name !== speciesName)
      : [...selectedSpecies, speciesName];
    onSpeciesChange(newSelection);
  };

  return (
    <div className="relative">
      <button
        onClick={() => setIsOpen(!isOpen)}
        disabled={disabled}
        className="w-full flex items-center justify-between rounded-md border border-gray-300 px-2 py-1 text-sm bg-white disabled:bg-gray-100"
      >
        <span className="truncate">
          {selectedSpecies.length 
            ? `${selectedSpecies.length} species selected`
            : "Select species"}
        </span>
        {isOpen ? <ChevronUp size={16} /> : <ChevronDown size={16} />}
      </button>

      {isOpen && (
        <div className="absolute z-50 mt-1 w-full rounded-md border border-gray-300 bg-white shadow-lg">
          <div className="p-2">
            <input
              type="text"
              placeholder="Search species..."
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
            />
          </div>
          <div className="max-h-48 overflow-auto">
            {filteredSpecies.map((name, index) => (
              <label
                key={`${name}-${index}`}
                className="flex items-center gap-2 px-2 py-1 hover:bg-gray-100 cursor-pointer"
                onClick={() => handleSpeciesToggle(name)}
              >
                <div className="flex items-center justify-center w-4 h-4 border border-gray-300 rounded">
                  {selectedSpecies.includes(name) && (
                    <Check size={12} className="text-indigo-600" />
                  )}
                </div>
                <span className="text-sm">{name}</span>
              </label>
            ))}
          </div>
        </div>
      )}
    </div>
  );
};

// ColorPickerModal Component
const ColorPickerModal = ({ 
  isOpen, 
  onClose, 
  items, 
  customColors,
  onColorChange 
}) => {
  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-white rounded-lg p-6 max-w-md w-full">
        <div className="flex justify-between items-center mb-4">
          <h3 className="text-lg font-semibold">Custom Colors</h3>
          <button 
            onClick={onClose}
            className="text-gray-500 hover:text-gray-700"
          >
            ×
          </button>
        </div>
        <div className="space-y-4">
          {items.map((item, index) => (
            <div key={item} className="flex items-center justify-between gap-4">
              <span className="text-sm font-medium">{item}</span>
              <input
                type="color"
                value={customColors[item] || '#000000'}
                onChange={(e) => onColorChange(item, e.target.value)}
                className="h-8 w-16 cursor-pointer rounded border border-gray-200"
              />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};


const LukeNetGUI = () => {
  const [plotType, setPlotType] = useState('Abundances');
  const [progress, setProgress] = useState(0);
  const [isRunning, setIsRunning] = useState(false);
  const [xAxisScale, setXAxisScale] = useState('log');
  const [yAxisScale, setYAxisScale] = useState('log');
  const [xAxisMin, setXAxisMin] = useState(1e-7);  
  const [xAxisMax, setXAxisMax] = useState(1e6);  
  const [yAxisMin, setYAxisMin] = useState(1e-8); 
  const [yAxisMax, setYAxisMax] = useState(3);    
  const [selfShielding, setSelfShielding] = useState('OFF');
  const [column, setColumn] = useState('OFF');
  const fileInputRef = useRef(null);
  const [numSpecies, setNumSpecies] = useState('8');  
  const [numRates, setNumRates] = useState('5');     
  const [numPathways, setNumPathways] = useState('');
  const [solverHasRun, setSolverHasRun] = useState(false);
  const [solverStarted, setSolverStarted] = useState(false);
  // Add these near the top with other state declarations
  const [solverData, setSolverData] = useState(null);
  const [chartData, setChartData] = useState([]);
  const [selectedSpecies, setSelectedSpecies] = useState([]);
  const [selectedReactions, setSelectedReactions] = useState([]);
  const [selectionMode, setSelectionMode] = useState('top-n'); // Controls whether we use top N species or manual selection
  const [manuallySelectedSpecies, setManuallySelectedSpecies] = useState([]); // Stores manually selected species

  const [rateSelectionMode, setRateSelectionMode] = useState('top-n');
  const [manuallySelectedRates, setManuallySelectedRates] = useState([]);
  const [lineOpacity, setLineOpacity] = useState(1);
  const [lineThickness, setLineThickness] = useState(1);
  const [showGrid, setShowGrid] = useState(true);
  const [lineStyle, setLineStyle] = useState('solid');
  const [showLines, setShowLines] = useState(true);
  const [showScatter, setShowScatter] = useState(true);
  const [colorScheme, setColorScheme] = useState('category10');  // Changed from 'default'

  const [customColors, setCustomColors] = useState({});
  const [isColorPickerOpen, setIsColorPickerOpen] = useState(false);
  const [lastColorScheme, setLastColorScheme] = useState('category10');

  const [maxConnections, setMaxConnections] = useState(5);
  const [centralNodeColor, setCentralNodeColor] = useState('#e6f2ff'); // Default blue color



  const COLOR_SCHEMES = {
    category10: index => ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'][index % 10],
    viridis: index => ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'][index % 6],
    spectral: index => ['#d53e4f', '#f46d43', '#fdae61', '#66c2a5', '#3288bd', '#5e4fa2'][index % 6],
    google: index => ['#4285F4', '#DB4437', '#F4B400', '#0F9D58', '#FF6D01', '#AA46BB'][index % 6],
    rainbow: index => ['#6e40aa', '#417de0', '#23abd8', '#29c172', '#b1c330', '#fb7e44', '#e6396c'][index % 7],
    accent: index => ['#7fc97f', '#beaed4', '#fdc086', '#fde725', '#386cb0', '#f0027f', '#bf5b17'][index % 7],
    terrain: index => ['#191980', '#0066cc', '#00b3e6', '#33cc33', '#b3d11f', '#cc9900', '#b35900', '#802b00'][index % 8],
    gates_of_tir_nanog: index => ['#408c6c', '#6eb9b4', '#65b1db', '#3d78c2', '#4f3794', '#indigo', '#00004d'][index % 7],
    luke: index => ['#2f7d8c', '#c45f5f', '#4a8435', '#8b4c7d', '#707eb5', '#c68930', '#456c8e', '#a05157'][index % 8],
    tol: index => ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"][index % 8],
    office: index => ["#4472c4", "#ed7d31", "#a5a5a5", "#ffc000", "#5b9bd5", "#70ad47", "#264478", "#9e480e", "#636363", "#997300"][index % 10],
    tableau20: index => ["#4e79a7", "#a0cbe8", "#f28e2b", "#ffbe7d", "#59a14f", "#8cd17d", "#b6992d", "#f1ce63", "#499894", "#86bcb6", 
                  "#e15759", "#ff9d9a", "#79706e", "#bab0ac", "#d37295", "#fabfd2", "#b07aa1", "#d4a6c8", "#9d7660", "#d7b5a6"][index % 20],
    custom: index => {
      const items = plotType === 'Abundances' ? selectedSpecies : selectedReactions;
      const item = items[index];
      return customColors[item] || COLOR_SCHEMES.category10(index);
    }
  };  


  // Handle custom color changes
  const handleCustomColorChange = (item, color) => {
    setCustomColors(prev => ({
      ...prev,
      [item]: color
    }));
  };

  // Get current plot items
  const currentPlotItems = useMemo(() => {
    return plotType === 'Abundances' ? selectedSpecies : selectedReactions;
  }, [plotType, selectedSpecies, selectedReactions]);



  // Add this function to get current colors before opening the modal
  const getCurrentColors = () => {
    const colors = {};
    const items = plotType === 'Abundances' ? selectedSpecies : selectedReactions;
    
    items.forEach((item, index) => {
      // If there's already a custom color, keep it
      if (customColors[item]) {
        colors[item] = customColors[item];
      } else {
        // Otherwise get the color from the current scheme
        colors[item] = COLOR_SCHEMES[colorScheme](index);
      }
    });
    
    return colors;
  };




  const [parameters, setParameters] = useState({
    n_gas: '',
    n_dust: '',
    t_gas: '',
    t_dust: '',
    gtd: '',
    Av: '',
    G_0: '',
    Zeta_X: '',
    h2_col: '',
    Zeta_CR: '',
    pah_ism: '',
    t_chem: '',
    network: ''
  });

  const parameterUnits = {
    n_gas: 'cm⁻³',
    n_dust: 'cm⁻³',
    t_gas: 'K',
    t_dust: 'K',
    gtd: 'ratio',  
    Av: 'mag',
    G_0: 'Draine',
    Zeta_X: 's⁻¹',
    Zeta_CR: 's⁻¹',
    pah_ism: '0.0 - 1.0', 
    h2_col: 'cm⁻²',
    t_chem: 'years',
    network: '.dat file'
  };

  const handleSaveData = () => {
    if (!chartData || !chartData.length) return;
    
    let content = '';
    // Get all data keys except 'time'
    const dataKeys = Object.keys(chartData[0]).filter(key => key !== 'time');
    
    // Create data rows for each data key
    dataKeys.forEach(key => {
      content += `\n# ${key}\n`; // Add header for each species/reaction
      chartData.forEach(point => {
        content += `${point.time.toExponential(3)} ${point[key].toExponential(3)}\n`;
      });
    });
    
    // Create and trigger download
    const blob = new Blob([content], { type: 'text/plain' });
    const downloadLink = document.createElement('a');
    downloadLink.href = URL.createObjectURL(blob);
    downloadLink.download = `lukenet-data-${Date.now()}.txt`;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
    URL.revokeObjectURL(downloadLink.href);
  };


  // Add ref for the chart container
  const chartRef = useRef(null);

  const handleSavePlot = () => {
    try {
      // Get the SVG element from the chart container
      const svgElement = chartRef.current?.querySelector('svg');
      if (!svgElement) {
        console.error('SVG element not found');
        return;
      }

      // Create a copy of the SVG element
      const svgClone = svgElement.cloneNode(true);
      
      // Set proper dimensions
      const bbox = svgElement.getBoundingClientRect();
      svgClone.setAttribute('width', bbox.width);
      svgClone.setAttribute('height', bbox.height);
      
      // Add font style to SVG
      const styleElement = document.createElement('style');
      styleElement.textContent = `
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&display=swap');
        
        * {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        }
        
        .recharts-legend-item-text {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif !important;
          font-size: 12px;
        }
        
        .recharts-cartesian-axis-tick-value {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif !important;
          font-size: 12px;
        }
        
        .recharts-label {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif !important;
          font-size: 12px;
        }
      `;
      svgClone.insertBefore(styleElement, svgClone.firstChild);

      // Convert SVG to string with XML declaration
      const finalSvgString = '<?xml version="1.0" encoding="UTF-8"?>\n' + 
                           new XMLSerializer().serializeToString(svgClone);
      
      // Create a Blob from the SVG string
      const svgBlob = new Blob([finalSvgString], { type: 'image/svg+xml;charset=utf-8' });
      
      // Create download link
      const downloadLink = document.createElement('a');
      downloadLink.href = URL.createObjectURL(svgBlob);
      downloadLink.download = `lukenet-plot-${Date.now()}.svg`;
      
      // Trigger download
      document.body.appendChild(downloadLink);
      downloadLink.click();
      document.body.removeChild(downloadLink);
      
      // Clean up the URL object
      URL.revokeObjectURL(downloadLink.href);
    } catch (error) {
      console.error('Error saving plot:', error);
    }
  };

  const handleSavePathway = () => {
    try {
      const diagramContainer = document.querySelector('.pathway-diagram-container');
      const svgElement = diagramContainer?.querySelector('svg');
      
      if (!svgElement) {
        alert('Could not find the pathway diagram to save. Please ensure a diagram is displayed.');
        return;
      }

      // Deep clone the SVG
      const svgClone = svgElement.cloneNode(true);
      
      // Set dimensions
      const bbox = svgElement.getBoundingClientRect();
      svgClone.setAttribute('width', bbox.width);
      svgClone.setAttribute('height', bbox.height);
      
      // Get computed styles from actual elements and apply them
      const centralNode = svgElement.querySelector('.central');
      if (centralNode) {
        const computedStyle = window.getComputedStyle(centralNode);
        // Apply the exact computed styles to the cloned central node
        const clonedCentralNode = svgClone.querySelector('.central');
        if (clonedCentralNode) {
          clonedCentralNode.style.fill = computedStyle.fill;
          clonedCentralNode.style.stroke = computedStyle.stroke;
          clonedCentralNode.style.strokeWidth = computedStyle.strokeWidth;
        }
      }

      // Ensure all Mermaid styles are included
      const styleElement = document.createElement('style');
      Array.from(svgElement.getElementsByTagName('style')).forEach(style => {
        styleElement.textContent += style.textContent;
      });
      
      // Add additional style overrides
      styleElement.textContent += `
        .node.central rect { 
          fill: #e6f2ff !important;
          stroke: #333 !important;
          stroke-width: 4px !important;
        }
        .node rect { 
          fill: white;
          stroke: #666;
          stroke-width: 2px;
        }
        .edgePath path {
          stroke-width: 1.5px;
        }
      `;
      
      svgClone.insertBefore(styleElement, svgClone.firstChild);

      // Convert to string with XML declaration
      const svgString = '<?xml version="1.0" encoding="UTF-8"?>\n' + 
                       new XMLSerializer().serializeToString(svgClone);
      
      // Create blob and download
      const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
      const downloadLink = document.createElement('a');
      downloadLink.href = URL.createObjectURL(blob);
      downloadLink.download = `lukenet-pathway-${Date.now()}.svg`;
      
      document.body.appendChild(downloadLink);
      downloadLink.click();
      document.body.removeChild(downloadLink);
      URL.revokeObjectURL(downloadLink.href);

    } catch (error) {
      alert('Error saving pathway diagram. Please try again.');
      console.error('Error saving pathway diagram:', error);
    }
  };


  // Add this handler function
  const handleFileLoad = async (event) => {
    const file = event.target.files[0];
    if (!file) return;

    const text = await file.text();
    
    // Parse the file content
    const lines = text.split('\n');
    const newParameters = { ...parameters };
    
    lines.forEach(line => {
      if (line.trim() === '') return;
      
      // Special handling for self_shielding and column
      if (line.includes('self_shielding')) {
        const selfShieldValue = line.includes('True') ? 'ON' : 'OFF';
        setSelfShielding(selfShieldValue);
        return;
      }

      if (line.includes('column')) {
        const columnValue = line.includes('True') ? 'ON' : 'OFF';
        setColumn(columnValue);
        return;
      }

      // Handle all other parameters
      if (line.includes('=')) {
        const [key, valueStr] = line.split('=').map(part => part.trim());
        const cleanKey = key.trim();
        
        if (cleanKey === 'network') {
          // Remove quotes from network path
          newParameters[cleanKey] = valueStr.replace(/['"]/g, '');
        } else {
          // Convert scientific notation to number
          try {
            newParameters[cleanKey] = parseFloat(valueStr);
          } catch (error) {
            console.error(`Error parsing value for ${cleanKey}`);
          }
        }
      }
    });

    setParameters(newParameters);
    
    // Reset the file input value so the same file can be selected again
    event.target.value = '';
  };


  const handleParameterChange = (key, value) => {
    setParameters(prev => ({
      ...prev,
      [key]: value
    }));
  };

  const sampleData = [
    { time: 1, H: 0.8, H2: 0.003, C: 0.4, N: 0.3, O: 0.5, CO: 0.2 },
    { time: 10, H: 0.7, H2: 0.1, C: 0.45, N: 0.35, O: 0.004, CO: 0.05 },
    { time: 100, H: 0.75, H2: 0.2, C: 0.5, N: 0.1, O: 0.02, CO: 0.006 },
    { time: 1000, H: 0.1, H2: 0.8, C: 0.55, N: 0.03, O: 0.2, CO: 0.002 },
    { time: 10000, H: 0.05, H2: 0.85, C: 0.6, N: 0.02, O: 0.35, CO: 0.001 },
    { time: 100000, H: 0.03, H2: 0.9, C: 0.1, N: 0.01, O: 0.3, CO: 0.0005 },
    { time: 1000000, H: 0.02, H2: 0.95, C: 0.07, N: 0.005, O: 0.15, CO: 0.0001 }
  ];







  // SOLVER

  const checkProgress = async () => {
    try {
      const response = await fetch('/api/check-progress');
      const data = await response.json();
      console.log('Progress check response:', data.progress); // Debug log
      return data.progress;
    } catch (error) {
      console.error('Error checking progress:', error);
      return null;
    }
  };

  const handleRunSolver = async () => {
    setSolverStarted(true);
    console.log('Starting solver...');
    setIsRunning(true);
    setProgress(0);
    
    const pollInterval = setInterval(async () => {
      console.log('Polling for progress...');
      const currentProgress = await checkProgress();
      if (currentProgress !== null) {
        console.log('Setting new progress:', currentProgress);
        setProgress(currentProgress);
      }
    }, 500);
    
    try {
      console.log('Making solver API call...');
      const response = await fetch('/api/run-solver', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          ...parameters,
          self_shielding: selfShielding,  // Add this line
          column: column                  // Add this line
        })
      });

      if (!response.ok) throw new Error('Failed to run solver');
      
      const data = await response.json();
      console.log('Solver API response:', data);

      if (data.success) {
        setProgress(100);
        setSolverHasRun(true);
        
        // Store the solver data
        setSolverData(data.data);
        
        // Set initial selected species/reactions
        if (data.data.species_names?.length > 0) {
          // Select first 6 species by default
          setSelectedSpecies(data.data.species_names.slice(0, 6));
        }
        if (data.data.reaction_labels?.length > 0) {
          // Select first 6 reactions by default
          setSelectedReactions(data.data.reaction_labels.slice(0, 6));
        }
        
        // Format initial chart data
        formatChartData(data.data, 'Abundances');
        
        console.log('Solver completed successfully');
      } else {
        console.error('Solver failed');
      }

    } catch (error) {
      console.error('Error:', error);
    } finally {
      console.log('Clearing poll interval');
      clearInterval(pollInterval);
      setIsRunning(false);
    }
  };


  const formatChartData = (data, type = plotType) => {
    if (!data?.time || !data?.abundances || !data?.rates) return;

    // Early return for Pathways type
    if (type === 'Pathways') {
      setChartData([]); // Clear chart data
      setSelectedSpecies([]); // Clear selected species
      setSelectedReactions([]); // Clear selected reactions
      return;
    }

    // Get indices to plot based on selection mode and plot type
    let indicesToPlot = [];  // Initialize with empty array

    if (type === 'Abundances') {
      if (selectionMode === 'top-n') {
        const numToPlot = Math.min(Math.max(parseInt(numSpecies) || 1, 1), data.species_names.length);
        const finalValues = data.abundances.map(speciesData => speciesData[speciesData.length - 1]);
        indicesToPlot = Array.from(Array(finalValues.length).keys())
          .sort((a, b) => finalValues[b] - finalValues[a])
          .slice(0, numToPlot);
      } else {
        indicesToPlot = manuallySelectedSpecies
          .map(name => data.species_names.indexOf(name))
          .filter(index => index !== -1);
        
        if (indicesToPlot.length === 0) {
          const finalValues = data.abundances.map(speciesData => speciesData[speciesData.length - 1]);
          indicesToPlot = Array.from(Array(finalValues.length).keys())
            .sort((a, b) => finalValues[b] - finalValues[a])
            .slice(0, 5);
        }
      }
    } else if (type === 'Rates') {
      if (rateSelectionMode === 'top-n') {
        const numToPlot = Math.min(parseInt(numRates) || 5, data.reaction_labels.length);
        const finalValues = data.rates[data.rates.length - 1];
        indicesToPlot = Array.from(Array(finalValues.length).keys())
          .sort((a, b) => finalValues[b] - finalValues[a])
          .slice(0, numToPlot);
      } else {
        indicesToPlot = manuallySelectedRates
          .map(label => data.reaction_labels.indexOf(label))
          .filter(index => index !== -1);
        
        if (indicesToPlot.length === 0) {
          const finalValues = data.rates[data.rates.length - 1];
          indicesToPlot = Array.from(Array(finalValues.length).keys())
            .sort((a, b) => finalValues[b] - finalValues[a])
            .slice(0, 5);
        }
      }
    }

    // Generate evenly spaced points in log space
    const numPoints = 27;
    const timeMin = Math.log10(Math.min(...data.time));
    const timeMax = Math.log10(Math.max(...data.time));
    const step = (timeMax - timeMin) / (numPoints - 1);
    
    const logSpacedTimes = Array.from({length: numPoints}, (_, i) => {
      return Math.pow(10, timeMin + i * step);
    });

    // Linear interpolation function
    const interpolateValue = (x, x1, x2, y1, y2) => {
      if (x1 === x2) return y1;
      const logX = Math.log10(x);
      const logX1 = Math.log10(x1);
      const logX2 = Math.log10(x2);
      const logY1 = Math.log10(Math.max(y1, 1e-40)); // Prevent log(0)
      const logY2 = Math.log10(Math.max(y2, 1e-40));
      
      // Linear interpolation in log space
      const logY = logY1 + (logX - logX1) * (logY2 - logY1) / (logX2 - logX1);
      return Math.pow(10, logY);
    };

    // Function to find nearest indices for interpolation
    const findNearestIndices = (time) => {
      let leftIndex = 0;
      while (leftIndex < data.time.length - 1 && data.time[leftIndex + 1] <= time) {
        leftIndex++;
      }
      return {
        left: leftIndex,
        right: Math.min(leftIndex + 1, data.time.length - 1)
      };
    };

    // Format the data with interpolation
    const formattedData = logSpacedTimes.map(t => {
      const point = { time: t };
      const { left, right } = findNearestIndices(t);

      if (type === 'Abundances') {
        indicesToPlot.forEach(speciesIndex => {
          const species = data.species_names[speciesIndex];
          const y1 = data.abundances[speciesIndex][left];
          const y2 = data.abundances[speciesIndex][right];
          point[species] = interpolateValue(t, data.time[left], data.time[right], y1, y2);
        });
      } else if (type === 'Rates') {
        indicesToPlot.forEach(reactionIndex => {
          const reaction = data.reaction_labels[reactionIndex];
          const y1 = data.rates[left][reactionIndex];
          const y2 = data.rates[right][reactionIndex];
          point[reaction] = interpolateValue(t, data.time[left], data.time[right], y1, y2);
        });
      }

      return point;
    });

    setChartData(formattedData);

    // Update selected species/reactions based on indices
    if (type === 'Abundances') {
      setSelectedSpecies(indicesToPlot.map(index => data.species_names[index]));
      setSelectedReactions([]); // Clear reactions when showing abundances
    } else if (type === 'Rates') {
      setSelectedSpecies([]); // Clear species when showing rates
      setSelectedReactions(indicesToPlot.map(index => data.reaction_labels[index]));
    }
  };



  useEffect(() => {
    if (solverData) {
      formatChartData(solverData, plotType);
    }
  }, [
    plotType, 
    numSpecies, 
    numRates, 
    manuallySelectedSpecies, 
    selectionMode,
    manuallySelectedRates,  
    rateSelectionMode,      
    solverData             
  ]);






  // Custom Tooltip Component
  const CustomTooltip = ({ active, payload, label }) => {
    if (!active || !payload || !payload.length) return null;

    const formatScientific = (value) => {
      if (value === 0) return '0.00e+0';
      const exp = Math.floor(Math.log10(Math.abs(value)));
      const coef = value / Math.pow(10, exp);
      return `${coef.toFixed(2)}e${exp >= 0 ? '+' : ''}${exp}`;
    };

    return (
      <div className="bg-white/90 backdrop-blur-sm border border-gray-200 rounded-lg p-3 shadow-lg">
        <p className="text-sm font-medium text-gray-600 mb-2">
          Time: {formatScientific(label)} years
        </p>
        {payload.map((entry, index) => (
          <p
            key={entry.name}
            className="text-sm"
            style={{ color: entry.color }}
          >
            {entry.name}: {formatScientific(entry.value)}
          </p>
        ))}
      </div>
    );
  };







  const getAxisDomain = (scale, min, max) => {
    if (scale === 'log') {
      return [
        min !== undefined ? Math.max(min, 1e-40) : 'auto',
        max !== undefined ? max : 'auto'
      ];
    } else {
      return [
        min !== undefined ? min : 'auto',
        max !== undefined ? max : 'auto'
      ];
    }
  };

  const generateTicks = (min, max, scale) => {
    if (min === undefined || max === undefined) return [];
    
    const ticks = [];
    let current = Math.floor(Math.log10(min));
    const end = Math.ceil(Math.log10(max));

    while (current <= end) {
      const tick = Math.pow(10, current);
      if (tick >= min && tick <= max) {
        ticks.push(tick);
      }
      current++;
    }

    return ticks;
  };

  const formatTickLabel = (value) => {
    // Convert to string in exponential format
    const expString = value.toExponential(0);
    const [base, exponent] = expString.split('e');
    
    // Remove '+' from positive exponents and convert to number to remove leading zeros
    const cleanExponent = Number(exponent).toString();
    
    // Convert the exponent to superscript numbers using Unicode
    const superscriptMap = {
      '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
      '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
      '-': '⁻'
    };
    
    const superscriptExp = cleanExponent
      .split('')
      .map(char => superscriptMap[char])
      .join('');
    
    // If the base is 1, just return 10 with superscript
    if (base === '1') {
      return `10${superscriptExp}`;
    }
    // For other cases, show the base multiplied by 10 with superscript
    return `${base}×10${superscriptExp}`;
  };


  const xTicks = useMemo(() => generateTicks(xAxisMin, xAxisMax, xAxisScale), [xAxisMin, xAxisMax, xAxisScale]);
  const yTicks = useMemo(() => generateTicks(yAxisMin, yAxisMax, yAxisScale), [yAxisMin, yAxisMax, yAxisScale]);

  return (
    <div className="w-[1440px] h-[920px] bg-gray-50 p-4 overflow-auto mx-auto">
      <div className="grid grid-cols-12 gap-4" style={{ width: '1400px' }}>
        {/* Left Panel - Parameters */}
        <Card className="col-span-3 bg-white h-full">
          <CardHeader className="pb-6">
            <CardTitle className="flex items-center gap-2 text-lg">
              <Settings size={18} />
              Input Parameters
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              {/* Gas/Dust Section */}
              <div className="space-y-1 pb-3">
                <h3 className="text-sm font-semibold text-gray-700 border-b pb-0.5">Gas/Dust</h3>
                {['n_gas', 'n_dust', 't_gas', 't_dust', 'gtd'].map((key) => (
                  <div key={key} className="flex items-center gap-1">
                    <label
                      className="text-sm font-medium text-gray-700 w-28"
                      dangerouslySetInnerHTML={{
                        __html: key
                          .replace(/n_gas/, 'n<sub>gas</sub>')
                          .replace(/n_dust/, 'n<sub>dust</sub>')
                          .replace(/t_gas/, 'T<sub>gas</sub>')
                          .replace(/t_dust/, 'T<sub>dust</sub>')
                          .replace(/gtd/, '&#916;<sub>gas/dust</sub>')
                      }}
                    />
                    <input
                      type="text"
                      className="flex-1 rounded-md border border-gray-300 px-2 py-0.5 text-sm"
                      value={parameters[key]}
                      onChange={(e) => handleParameterChange(key, e.target.value)}
                      placeholder={parameterUnits[key]}
                    />
                  </div>
                ))}
              </div>

              {/* Environment Section */}
              <div className="space-y-1 pb-3">
                <h3 className="text-sm font-semibold text-gray-700 border-b pb0.5">Environment</h3>
                
                {/* Regular input fields */}
                {['Av', 'G_0', 'Zeta_X', 'Zeta_CR', 'pah_ism'].map((key) => (
                  <div key={key} className="flex items-center gap-1">
                    <label
                      className="text-sm font-medium text-gray-700 w-28"
                      dangerouslySetInnerHTML={{
                        __html: key
                          .replace(/Av/, 'A<sub>V</sub>')
                          .replace(/G_0/, 'G<sub>0</sub>')
                          .replace(/Zeta_X/, '&zeta;<sub>X</sub>')
                          .replace(/Zeta_CR/, '&zeta;<sub>CR</sub>')
                          .replace(/pah_ism/, 'PAH<sub>ISM</sub>')
                      }}
                    />
                    <input
                      type="text"
                      className="flex-1 rounded-md border border-gray-300 px-2 py-0.5 text-sm"
                      value={parameters[key]}
                      onChange={(e) => handleParameterChange(key, e.target.value)}
                      placeholder={parameterUnits[key]}
                    />
                  </div>
                ))}

                {/* Self Shielding with nested structure */}
                <div className="space-y-1 pt-1">
                  {/* Self Shielding */}
                  <div className="flex items-center gap-1">
                    <label className="text-sm font-medium text-gray-700 w-28">
                      Self Shielding
                    </label>
                    <div className="flex gap-3">
                      <label className="flex items-center gap-1">
                        <input
                          type="radio"
                          name="selfShielding"
                          value="ON"
                          checked={selfShielding === 'ON'}
                          onChange={(e) => setSelfShielding(e.target.value)}
                          className="w-3 h-3"
                        />
                        <span className="text-sm">ON</span>
                      </label>
                      <label className="flex items-center gap-1">
                        <input
                          type="radio"
                          name="selfShielding"
                          value="OFF"
                          checked={selfShielding === 'OFF'}
                          onChange={(e) => {
                            setSelfShielding(e.target.value);
                            setColumn('OFF');
                          }}
                          className="w-3 h-3"
                        />
                        <span className="text-sm">OFF</span>
                      </label>
                    </div>
                  </div>

                  {/* Column (indented) */}
                  <div className={`flex items-center gap-1 ${selfShielding === 'OFF' ? 'opacity-50' : ''}`}>
                    <div className="w-28 flex items-center">
                      <span className="text-gray-400 ml-3">└─</span>
                      <label className="text-sm font-medium text-gray-700 ml-1">
                        Column
                      </label>
                    </div>
                    <div className="flex gap-3">
                      <label className="flex items-center gap-1">
                        <input
                          type="radio"
                          name="column"
                          value="ON"
                          checked={column === 'ON'}
                          onChange={(e) => setColumn(e.target.value)}
                          className="w-3 h-3"
                          disabled={selfShielding === 'OFF'}
                        />
                        <span className="text-sm">ON</span>
                      </label>
                      <label className="flex items-center gap-1">
                        <input
                          type="radio"
                          name="column"
                          value="OFF"
                          checked={column === 'OFF'}
                          onChange={(e) => setColumn(e.target.value)}
                          className="w-3 h-3"
                          disabled={selfShielding === 'OFF'}
                        />
                        <span className="text-sm">OFF</span>
                      </label>
                    </div>
                  </div>

                  {/* H2 Column (double indented) */}
                  <div className={`flex items-center gap-1 ${selfShielding === 'OFF' || column === 'OFF' ? 'opacity-50' : ''}`}>
                    <div className="w-28 flex items-center">
                      <span className="text-gray-400 ml-7">└─</span>
                      <label className="text-sm font-medium text-gray-700 ml-1" dangerouslySetInnerHTML={{__html: 'N(H<sub>2</sub>)'}}>
                      </label>
                    </div>
                    <input
                      type="text"
                      className="flex-1 rounded-md border border-gray-300 px-2 py-0.5 text-sm"
                      value={parameters['h2_col']}
                      onChange={(e) => handleParameterChange('h2_col', e.target.value)}
                      disabled={selfShielding === 'OFF' || column === 'OFF'}
                      placeholder={parameterUnits['h2_col']}
                    />
                  </div>
                </div>
              </div>


              {/* Network Section */}
              <div className="space-y-1 pb-3">
                <h3 className="text-sm font-semibold text-gray-700 border-b pb-0.5">Network</h3>
                {['t_chem', 'network'].map((key) => (
                  <div key={key} className="flex items-center gap-1">
                    <label
                      className="text-sm font-medium text-gray-700 w-28"
                      dangerouslySetInnerHTML={{
                        __html: key
                          .replace(/t_chem/, '<i>t</i><sub>chem</sub>')
                          .replace(/network/, 'Network path')
                      }}
                    />
                    <input
                      type="text"
                      className="flex-1 rounded-md border border-gray-300 px-2 py-0.5 text-sm"
                      value={parameters[key]}
                      onChange={(e) => handleParameterChange(key, e.target.value)}
                      placeholder={parameterUnits[key]}
                    />
                  </div>
                ))}
              </div>


              {/* Add this hidden file input */}
              <input 
                type="file"
                ref={fileInputRef}
                className="hidden"
                accept=".dat"
                onChange={handleFileLoad}
              />
              
              <div className="space-y-1 pt-4">
                <button 
                  className="flex w-full items-center justify-center gap-1 rounded-md bg-white px-3 py-1 border border-slate-400 text-sm text-slate-600 hover:bg-gray-100"
                  onClick={() => fileInputRef.current.click()}
                >
                  <Upload size={12} />
                  Load Config
                </button>

                <button 
                  className="flex w-full items-center justify-center gap-1 rounded-md bg-white px-3 py-1 border border-slate-400 text-sm font-bold text-slate-600 hover:bg-gray-100"
                  onClick={handleRunSolver}
                  disabled={isRunning}
                >
                  <Play size={12} strokeWidth={2.5} />
                  Run Solver
                </button>
              </div>

              {progress > -1 && (
                <div className="space-y-1 pt-4">
                  <div className="text-base text-gray-600">
                    {solverStarted ? `Progress: ${progress}%` : <span className="italic">Waiting</span>}
                  </div>
                  <div className="h-2.5 w-full rounded-full bg-gray-200">
                    <div
                      className="h-2.5 rounded-full bg-gray-400 transition-all duration-300"
                      style={{ width: `${progress}%` }}
                    />
                  </div>
                </div>

              )}



            </div>
          </CardContent>


        </Card>


        {/* Center Panel - Plot */}
        <Card className="col-span-6 bg-white h-full">
          <CardHeader className="pb-8">
            <div className="flex flex-col items-center">
              {solverHasRun && (
                <>
                  <CardTitle className="flex items-center gap-3 text-4xl text-gray-800">
                    <ChartIcon size={32} />
                    SIMBA
                  </CardTitle>
                  <p className="text-base text-gray-500 mt-1 mb-3">
                    
                  </p>
                  <div className="text-lg font-semibold text-gray-700 mt-2">
                    {plotType === 'Abundances' ? 'Abundances vs. Time' : 
                     plotType === 'Rates' ? 'Reaction Rates vs. Time' :
                     'Chemical Pathways'}
                  </div>
                </>
              )}
            </div>
          </CardHeader>
          <CardContent>
            {solverHasRun ? (
              <>
                {plotType === 'Pathways' ? (
                  <div className="h-[450px] w-full">
                    <PathwayDiagram
                      selectedSpecies={manuallySelectedSpecies[0]}
                      reactionData={solverData}
                      maxConnections={maxConnections}
                    />
                  </div>
                ) : (
                  <div ref={chartRef}>
                    <ResponsiveContainer width="100%" height={450}>
                      <LineChart 
                        data={chartData}
                        margin={{ top: 30, right: 40, left: -10, bottom: -55 }}
                      >
                        <CartesianGrid strokeDasharray="3 3" opacity={showGrid ? 1 : 0} />
                        <XAxis 
                          dataKey="time" 
                          label={{ value: 'Time (years)', position: 'bottom', offset: 0 }}
                          scale={xAxisScale}
                          type="number"
                          domain={getAxisDomain(xAxisScale, xAxisMin, xAxisMax)}
                          ticks={xTicks}
                          tickFormatter={formatTickLabel}
                          allowDataOverflow={true}
                        />
                        <YAxis 
                          label={{ 
                            value: plotType, 
                            angle: -90, 
                            position: 'insideLeft',
                            offset: -20,
                            dx: 55,
                            dy: 1,
                            textAnchor: 'middle',
                            style: {
                              textAnchor: 'middle',
                              dominantBaseline: 'middle',
                            }
                          }} 
                          scale={yAxisScale}
                          domain={getAxisDomain(yAxisScale, yAxisMin, yAxisMax)}
                          ticks={yTicks}
                          tickFormatter={formatTickLabel}
                          allowDataOverflow={true}
                          width={100}
                        />
                        <Tooltip content={<CustomTooltip />} />
                        <Legend 
                          layout="horizontal"
                          verticalAlign="bottom"
                          align="center"
                          wrapperStyle={{
                            paddingTop: '50px',
                            width: '100%',
                            display: 'flex',
                            justifyContent: 'space-around'
                          }}
                        />
                        {plotType === 'Abundances' ? 
                          selectedSpecies.map((species, index) => (
                            <Line 
                              key={species}
                              type="monotone"
                              dataKey={species}
                              stroke={COLOR_SCHEMES[colorScheme](index)}
                              strokeOpacity={showLines ? lineOpacity : 0}
                              strokeWidth={lineThickness}
                              strokeDasharray={lineStyle === 'dashed' ? '5 5' : lineStyle === 'dotted' ? '2 2' : '0'}
                              dot={showScatter ? {
                                fill: `hsl(${(index * 360) / selectedSpecies.length}, 70%, 50%)`,
                                r: 3,
                                fill: "white",
                                strokeOpacity: lineOpacity,  
                                fillOpacity: lineOpacity    
                              } : false}
                              name={species}
                            />
                          ))
                          :
                          selectedReactions.map((reaction, index) => (
                            <Line 
                              key={reaction}
                              type="monotone"
                              dataKey={reaction}
                              stroke={COLOR_SCHEMES[colorScheme](index)}
                              strokeOpacity={showLines ? lineOpacity : 0}
                              strokeWidth={lineThickness}
                              strokeDasharray={lineStyle === 'dashed' ? '5 5' : lineStyle === 'dotted' ? '2 2' : '0'}
                              dot={showScatter ? {
                                fill: `hsl(${(index * 360) / selectedReactions.length}, 70%, 50%)`,
                                r: 3,
                                fill: "white",
                                strokeOpacity: lineOpacity,  
                                fillOpacity: lineOpacity    
                              } : false}
                              name={reaction}
                            />
                          ))
                        }
                      </LineChart>
                    </ResponsiveContainer>
                  </div>
                )}
              </>
            ) : (
              <div className="flex items-center justify-center h-[700px] text-gray-500">
                <div className="flex items-center justify-center w-full h-full">
                  <img 
                    src={Logo.src}
                    alt="LukeNet Logo" 
                    className="w-30 h-auto max-w-full object-contain"
                    style={{ opacity: 1 }}  // Set your desired opacity value here (0-1)
                  />
                </div>
              </div>
            )}
          </CardContent>
        </Card>


        {/* Right Panel - Controls */}
        <Card className="col-span-3 bg-white h-full">
          <CardHeader className="pb-6">
            <CardTitle className="flex items-center gap-2 text-lg">
              <RefreshCw size={18} />
              Control Panel
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              {/* Data Selection Section */}
              <div className="space-y-1 pb-6">
                <h3 className="text-sm font-semibold text-gray-700 border-b pb-0.5">Data Selection</h3>
                
                {/* Plot Type */}
                <div className="space-y-1 pt-2">
                  <label className="text-sm font-medium text-gray-700">Plot Type</label>
                  <select 
                    className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                    value={plotType}
                    onChange={(e) => setPlotType(e.target.value)}
                  >
                    <option value="Abundances">Abundances</option>
                    <option value="Rates">Rates</option>
                    <option value="Pathways">Pathways</option>
                  </select>
                </div>

                {/* Species Selection Section */}
                {plotType === 'Abundances' && (
                  <div className="space-y-3 pt-2">
                    <div className="space-y-1">
                      <label className="text-sm font-medium text-gray-700">Selection Mode</label>
                      <div className="flex gap-4">
                        <label className="flex items-center gap-1">
                          <input
                            type="radio"
                            name="selectionMode"
                            value="top-n"
                            checked={selectionMode === 'top-n'}
                            onChange={(e) => setSelectionMode(e.target.value)}
                            className="w-3 h-3"
                          />
                          <span className="text-sm">Top N Species</span>
                        </label>
                        <label className="flex items-center gap-1">
                          <input
                            type="radio"
                            name="selectionMode"
                            value="manual"
                            checked={selectionMode === 'manual'}
                            onChange={(e) => setSelectionMode(e.target.value)}
                            className="w-3 h-3"
                          />
                          <span className="text-sm">Manual Selection</span>
                        </label>
                      </div>
                    </div>

                    {selectionMode === 'top-n' ? (
                      <div className="space-y-1">
                        <label className="text-sm font-medium text-gray-700">
                          Number of Species
                        </label>
                        <input
                          type="number"
                          className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                          value={numSpecies}
                          onChange={(e) => setNumSpecies(e.target.value)}
                          placeholder="Enter number"
                          min="1"
                        />
                      </div>
                    ) : (
                      <div className="space-y-1">
                        <label className="text-sm font-medium text-gray-700">
                          Select Species
                        </label>
                        <SpeciesSelector
                          species={solverData?.species_names || []}
                          selectedSpecies={manuallySelectedSpecies}
                          onSpeciesChange={setManuallySelectedSpecies}
                          disabled={!solverData}
                        />
                      </div>
                    )}
                  </div>
                )}

                {/* Rates Selection Section */}
                {plotType === 'Rates' && (
                  <div className="space-y-3 pt-2">
                    <div className="space-y-1">
                      <label className="text-sm font-medium text-gray-700">Selection Mode</label>
                      <div className="flex gap-4">
                        <label className="flex items-center gap-1">
                          <input
                            type="radio"
                            name="rateSelectionMode"
                            value="top-n"
                            checked={rateSelectionMode === 'top-n'}
                            onChange={(e) => setRateSelectionMode(e.target.value)}
                            className="w-3 h-3"
                          />
                          <span className="text-sm">Top N Rates</span>
                        </label>
                        <label className="flex items-center gap-1">
                          <input
                            type="radio"
                            name="rateSelectionMode"
                            value="manual"
                            checked={rateSelectionMode === 'manual'}
                            onChange={(e) => setRateSelectionMode(e.target.value)}
                            className="w-3 h-3"
                          />
                          <span className="text-sm">Manual Selection</span>
                        </label>
                      </div>
                    </div>

                    {rateSelectionMode === 'top-n' ? (
                      <div className="space-y-1">
                        <label className="text-sm font-medium text-gray-700">
                          Number of Rates
                        </label>
                        <input
                          type="number"
                          className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                          value={numRates}
                          onChange={(e) => setNumRates(e.target.value)}
                          placeholder="Enter number"
                          min="1"
                        />
                      </div>
                    ) : (
                      <div className="space-y-1">
                        <label className="text-sm font-medium text-gray-700">
                          Select Rates
                        </label>
                        <SpeciesSelector
                          species={solverData?.reaction_labels || []}
                          selectedSpecies={manuallySelectedRates}
                          onSpeciesChange={setManuallySelectedRates}
                          disabled={!solverData}
                        />
                      </div>
                    )}
                  </div>
                )}

                {/* Pathways Selection Section */}
                {plotType === 'Pathways' && (
                  <div className="space-y-3 pt-2">
                    <div className="space-y-1">
                      <label className="text-sm font-medium text-gray-700">
                        Select Species
                      </label>
                      <SpeciesSelector
                        species={solverData?.species_names || []}
                        selectedSpecies={manuallySelectedSpecies}
                        onSpeciesChange={(selected) => {
                          // Only take the last selected species
                          setManuallySelectedSpecies(selected.length > 0 ? [selected[selected.length - 1]] : []);
                        }}
                        disabled={!solverData}
                      />
                    </div>
                    <div className="space-y-1">
                      <label className="text-sm font-medium text-gray-700">
                        Number of Connections
                      </label>
                      <input
                        type="number"
                        className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                        value={maxConnections}
                        onChange={(e) => setMaxConnections(Math.max(1, parseInt(e.target.value) || 1))}
                        min="1"
                        max="10"
                        placeholder="Enter number (1-10)"
                      />
                    </div>
                  </div>
                )}
              </div>

              {/* Plot Settings Section - Only show for Abundances and Rates */}
              {plotType !== 'Pathways' && (
                <div className="space-y-3 pb-3">
                  <h3 className="text-sm font-semibold text-gray-700 border-b pb-0.5">Plot Settings</h3>
                  
                  {/* Line Appearance Group */}
                  <div className="space-y-2 pt-1">
                    {/* Plot Style Toggle Controls */}
                    <div className="grid grid-cols-3 gap-4">
                      <div className="flex items-center gap-2">
                        <label className="text-xs font-medium text-gray-700">Show lines</label>
                        <input
                          type="checkbox"
                          checked={showLines}
                          onChange={(e) => setShowLines(e.target.checked)}
                          className="rounded border-gray-300"
                        />
                      </div>
                      <div className="flex items-center gap-2">
                        <label className="text-xs font-medium text-gray-700">Show dots</label>
                        <input
                          type="checkbox"
                          checked={showScatter}
                          onChange={(e) => setShowScatter(e.target.checked)}
                          className="rounded border-gray-300"
                        />
                      </div>
                      <div className="flex items-center gap-2">
                        <label className="text-xs font-medium text-gray-700">Show grid</label>
                        <input
                          type="checkbox"
                          checked={showGrid}
                          onChange={(e) => setShowGrid(e.target.checked)}
                          className="rounded border-gray-300"
                        />
                      </div>
                    </div>

                    {/* Line Style and Color Scheme */}
                    <div className="grid grid-cols-2 gap-4">
                      <div className="space-y-1">
                        <label className={`text-xs font-medium ${showLines ? 'text-gray-700' : 'text-gray-400'}`}>Line Style</label>
                        <select
                          className={`w-full rounded-md border px-2 py-1 text-sm ${
                            showLines 
                              ? 'border-gray-300 bg-white text-gray-900' 
                              : 'border-gray-200 bg-gray-100 text-gray-400'
                          }`}
                          value={lineStyle}
                          onChange={(e) => setLineStyle(e.target.value)}
                          disabled={!showLines}
                        >
                          <option value="solid">Solid</option>
                          <option value="dashed">Dashed</option>
                          <option value="dotted">Dotted</option>
                        </select>
                      </div>

                      <div className="space-y-1">
                        <label className="text-xs font-medium text-gray-700">Color Scheme</label>
                        <select
                          className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                          value={colorScheme}
                          onChange={(e) => {
                            const value = e.target.value;
                            setColorScheme(value);
                            if (value !== 'custom') {
                              setLastColorScheme(value);
                            }
                          }}
                        >
                          <option value="category10">Default</option>
                          <option value="accent">Accent</option>
                          <option value="gates_of_tir_nanog">Gates of Tir-Nanog</option>
                          <option value="google">Google</option>
                          <option value="luke">Luke</option>
                          <option value="office">Office</option>
                          <option value="rainbow">Rainbow</option>
                          <option value="spectral">Spectral</option>
                          <option value="tableau20">Tableau 20</option>
                          <option value="terrain">Terrain</option>
                          <option value="tol">Tol</option>
                          <option value="viridis">Viridis</option>
                        </select>
                        <button
                          onClick={() => {
                            // Get colors from the last used scheme
                            const colors = {};
                            const items = plotType === 'Abundances' ? selectedSpecies : selectedReactions;
                            items.forEach((item, index) => {
                              colors[item] = COLOR_SCHEMES[lastColorScheme](index);
                            });
                            
                            setCustomColors(colors);
                            setColorScheme('custom');
                            setIsColorPickerOpen(true);
                          }}
                          className={`w-full text-left px-2 py-1 text-sm border rounded-md ${
                            colorScheme === 'custom' 
                              ? 'border-indigo-500 bg-indigo-50 text-indigo-600' 
                              : 'border-gray-300 hover:bg-gray-50'
                          }`}
                        >
                          Custom colors...
                        </button>
                      </div>
                    </div>

                    <div className="space-y-2 pt-1">
                      <div className="grid grid-cols-2 gap-4">
                        <div>
                          <label className="text-xs font-medium text-gray-700">Opacity</label>
                          <div className="flex items-center gap-2">
                            <input
                              type="range"
                              min="0"
                              max="1"
                              step="0.1"
                              value={lineOpacity} 
                              onChange={(e) => setLineOpacity(Number(e.target.value))}
                              className="w-full h-1 appearance-none bg-gray-300 rounded-full [&::-webkit-slider-thumb]:w-3 [&::-webkit-slider-thumb]:h-3 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:bg-blue-600 [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:cursor-pointer"
                            />
                            <span className="text-xs text-gray-500 min-w-[24px]">{lineOpacity.toFixed(1)}</span>
                          </div>
                        </div>

                        <div>
                          <label className="text-xs font-medium text-gray-700">Linewidth</label>
                          <div className="flex items-center gap-2">
                            <input
                              type="range" 
                              min="0.5"
                              max="5"
                              step="0.5"
                              value={lineThickness}
                              onChange={(e) => setLineThickness(Number(e.target.value))}
                              className="w-full h-1 appearance-none bg-gray-300 rounded-full [&::-webkit-slider-thumb]:w-3 [&::-webkit-slider-thumb]:h-3 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:bg-blue-600 [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:cursor-pointer"
                            />
                            <span className="text-xs text-gray-500 min-w-[24px]">{lineThickness.toFixed(1)}</span>
                          </div>
                        </div>
                      </div>
                    </div>
                  </div>

                  {/* Axis Settings Group */}
                  <div className="space-y-2 border-t pt-2">
                    <div className="grid grid-cols-2 gap-4">
                      <div className="space-y-1">
                        <label className="text-xs font-medium text-gray-700">X-Axis Scale</label>
                        <select
                          className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                          value={xAxisScale}
                          onChange={(e) => setXAxisScale(e.target.value)}
                        >
                          <option                   value="linear">Linear</option>
                          <option value="log">Logarithmic</option>
                        </select>
                      </div>

                      <div className="space-y-1">
                        <label className="text-xs font-medium text-gray-700">Y-Axis Scale</label>
                        <select
                          className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                          value={yAxisScale}
                          onChange={(e) => setYAxisScale(e.target.value)}
                        >
                          <option value="linear">Linear</option>
                          <option value="log">Logarithmic</option>
                        </select>
                      </div>
                    </div>

                    <div className="space-y-2">
                      <div className="space-y-1">
                        <label className="text-xs font-medium text-gray-700">X-Axis Limits</label>
                        <div className="flex gap-2">
                          <input
                            type="number"
                            placeholder="Min"
                            className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                            onChange={(e) => setXAxisMin(e.target.value ? Number(e.target.value) : undefined)}
                          />
                          <input
                            type="number"
                            placeholder="Max"
                            className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                            onChange={(e) => setXAxisMax(e.target.value ? Number(e.target.value) : undefined)}
                          />
                        </div>
                      </div>

                      <div className="space-y-1">
                        <label className="text-xs font-medium text-gray-700">Y-Axis Limits</label>
                        <div className="flex gap-2">
                          <input
                            type="number"
                            placeholder="Min"
                            className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                            onChange={(e) => setYAxisMin(e.target.value ? Number(e.target.value) : undefined)}
                          />
                          <input
                            type="number"
                            placeholder="Max"
                            className="w-full rounded-md border border-gray-300 px-2 py-1 text-sm"
                            onChange={(e) => setYAxisMax(e.target.value ? Number(e.target.value) : undefined)}
                          />
                        </div>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Pathway Diagram Settings - Only show for Pathways */}
              {plotType === 'Pathways' && (
                <div className="space-y-3 pb-3">
                  <h3 className="text-sm font-semibold text-gray-700 border-b pb-0.5">Diagram Settings</h3>
                  <div className="space-y-2">
                    <div className="space-y-1">
                      <label className="text-sm font-medium text-gray-700">
                        Central Node Color
                      </label>
                      <div className="flex items-center gap-2">
                        <input
                          type="color"
                          value={centralNodeColor}
                          onChange={(e) => setCentralNodeColor(e.target.value)}
                          className="h-8 w-16 cursor-pointer rounded border border-gray-200"
                        />
                        <input 
                          type="text"
                          value={centralNodeColor}
                          onChange={(e) => setCentralNodeColor(e.target.value)}
                          className="flex-1 rounded-md border border-gray-300 px-2 py-1 text-sm"
                          placeholder="#e6f2ff"
                        />
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Save/Export Buttons - Always visible but contextual */}
              <div className="flex gap-2 pt-1">
                <button 
                  className={`flex ${plotType === 'Pathways' ? 'w-full' : 'w-1/2'} items-center justify-center gap-2 rounded-md bg-white border border-slate-400 text-slate-600 px-3 py-1.5 text-sm hover:bg-gray-100`}
                  onClick={plotType === 'Pathways' ? handleSavePathway : handleSavePlot}
                >
                  <Save size={14} />
                  Save {plotType === 'Pathways' ? 'Diagram' : 'Plot'}
                </button>
                {plotType !== 'Pathways' && (
                  <button 
                    className="flex w-1/2 items-center justify-center gap-2 rounded-md bg-white border border-slate-400 text-slate-600 px-3 py-1.5 text-sm hover:bg-gray-100"
                    onClick={handleSaveData}
                  >
                    <ChartIcon size={14} />
                    Export Data
                  </button>
                )}
              </div>
            </div>
          </CardContent>
        </Card>

      </div>

      <ColorPickerModal
        isOpen={isColorPickerOpen}
        onClose={() => setIsColorPickerOpen(false)}
        items={currentPlotItems}
        customColors={customColors}
        onColorChange={handleCustomColorChange}
      />

    </div>
  );
};

export default LukeNetGUI;

