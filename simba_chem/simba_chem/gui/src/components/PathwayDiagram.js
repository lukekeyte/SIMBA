import React, { useEffect, useRef } from 'react';
import mermaid from 'mermaid';

// Initialize mermaid config
mermaid.initialize({
  startOnLoad: true,
  theme: 'default',
  securityLevel: 'loose',
  flowchart: {
    curve: 'basis',
    padding: 7,
    nodeSpacing: 50,
    rankSpacing: 50,
    htmlLabels: true,
    useMaxWidth: true
  }
});

// Helper function to sanitize node IDs for Mermaid
const sanitizeNodeId = (id) => {
  // Handle empty strings
  if (!id || id === '') return 'empty';
  
  // Replace special characters and spaces
  return id
    .replace(/[\s\[\](){}'"]/g, '_')  // Replace spaces and brackets with underscore
    .replace(/[+\-*/\\=,;.!@#$%^&*]/g, '_')  // Replace special characters
    .replace(/^[0-9]/, 'n$&');  // Prefix numbers with 'n'
};

// Helper function to sanitize node labels
const sanitizeLabel = (label) => {
  // Handle empty strings
  if (!label || label === '') return 'empty';
  
  // Escape special characters for display
  return label.replace(/["\\]/g, '\\$&');
};

const PathwayDiagram = ({ 
  selectedSpecies,
  reactionData,
  maxConnections = 5,
  centralNodeColor = '#e6f2ff' 
}) => {
  const containerRef = useRef(null);
  
  useEffect(() => {
    if (!containerRef.current) return;

    // Clear previous diagram
    containerRef.current.innerHTML = '';

    // If no species selected, show placeholder message
    if (!selectedSpecies || !reactionData) {
      const placeholder = document.createElement('div');
      placeholder.className = 'text-gray-500 text-center mt-20';
      placeholder.textContent = 'Select a species to view its reaction pathways';
      containerRef.current.appendChild(placeholder);
      return;
    }

    try {
      // Update Mermaid configuration for each render
      mermaid.initialize({
        startOnLoad: true,
        theme: 'default',
        securityLevel: 'loose',
        flowchart: {
          curve: 'basis',
          padding: 8,
          nodeSpacing: 50,
          rankSpacing: 50,
          defaultRenderer: 'dagre-d3'
        }
      });

      // Find reactions involving the selected species
      const relevantReactions = findRelevantReactions(selectedSpecies, reactionData);
      
      // Generate Mermaid diagram definition
      const mermaidDef = generateMermaidDefinition(
        selectedSpecies, 
        relevantReactions, 
        maxConnections,
        centralNodeColor
      );
      
      // Create new diagram container with unique ID
      const diagramId = `pathway-diagram-${Date.now()}`;
      
      // Render new diagram
      mermaid.render(diagramId, mermaidDef)
        .then(({ svg }) => {
          if (containerRef.current) {
            containerRef.current.innerHTML = svg;
            
            // Force style update after render
            const centralNode = containerRef.current.querySelector('.node.central rect');
            if (centralNode) {
              centralNode.style.fill = centralNodeColor;
            }
          }
        })
        .catch(error => {
          console.error('Mermaid rendering failed:', error);
          if (containerRef.current) {
            const errorMsg = document.createElement('div');
            errorMsg.className = 'text-red-500 text-center mt-20';
            errorMsg.textContent = 'Error rendering pathway diagram';
            containerRef.current.appendChild(errorMsg);
          }
        });
    } catch (error) {
      console.error('Error processing reaction data:', error);
      const errorMsg = document.createElement('div');
      errorMsg.className = 'text-red-500 text-center mt-20';
      errorMsg.textContent = 'Error processing reaction data';
      containerRef.current.appendChild(errorMsg);
    }
  }, [selectedSpecies, reactionData, maxConnections, centralNodeColor]);

  return (
    <div className="flex flex-col items-center">
      <div 
        ref={containerRef} 
        className="w-full h-full flex items-center justify-center overflow-auto pathway-diagram-container"
        style={{ minHeight: '450px' }}
      />
    </div>
  );
};

// Helper function to find relevant reactions
const findRelevantReactions = (species, reactionData) => {
  const { reaction_labels = [], rates = [] } = reactionData;
  
  // Get the final timestep rates
  const finalRates = rates[rates.length - 1] || [];
  
  // Helper function to check if a species is actually involved in a reaction
  const isSpeciesInvolved = (reaction, targetSpecies) => {
    const [reactants, products] = reaction.split('->').map(s => s.trim());
    const allSpecies = [...reactants.split(' + '), ...products.split(' + ')].map(s => s.trim());
    
    // Check if the exact species exists in the reaction
    return allSpecies.some(s => s === targetSpecies);
  };
  
  // Map reactions to their rates and filter for ones involving the species
  return reaction_labels
    .map((reaction, index) => ({
      reaction,
      rate: finalRates[index] || 0
    }))
    .filter(({ reaction }) => isSpeciesInvolved(reaction, species))
    .sort((a, b) => b.rate - a.rate); // Sort by rate descending
};

// Helper function to generate Mermaid diagram definition
// Only showing the modified generateMermaidDefinition function
const generateMermaidDefinition = (centralSpecies, reactions, maxConnections, centralNodeColor) => {
  let def = 'graph TD;\n';
  
  // Style definitions first (before any nodes)
  def += `
    classDef default fill:#fff,stroke:#666,stroke-width:2px,rx:10,ry:10;
    classDef central fill:${centralNodeColor},stroke:#333,stroke-width:4px;
  \n`;
  
  // Style for the central species
  const centralId = sanitizeNodeId(centralSpecies);
  def += `${centralId}(["${sanitizeLabel(centralSpecies)}"]):::central;\n`;
  
  // Add reactions and connections
  const addedConnections = new Set();
  let connectionCount = 0;
  
  reactions.forEach(({ reaction }) => {
    if (connectionCount >= maxConnections) return;
    
    const [reactants, products] = reaction.split('->').map(s => s.trim());
    const reactantList = reactants.split(' + ').map(s => s.trim());
    const productList = products.split(' + ').map(s => s.trim());
    
    const reactionText = reaction.replace('->', '→');
    
    // Check if central species is a reactant or product
    const isCentralReactant = reactantList.includes(centralSpecies);
    const isCentralProduct = productList.includes(centralSpecies);
    
    if (isCentralReactant) {
      // If central species is a reactant, only connect to products
      productList.forEach(product => {
        if (product !== centralSpecies && connectionCount < maxConnections) {
          const productId = sanitizeNodeId(product);
          const connectionKey = `${centralId}-${productId}`;
          if (!addedConnections.has(connectionKey)) {
            def += `${centralId} -->|${sanitizeLabel(reactionText)}| ${productId}["${sanitizeLabel(product)}"]:::default;\n`;
            addedConnections.add(connectionKey);
            connectionCount++;
          }
        }
      });
    }
    
    if (isCentralProduct) {
      // If central species is a product, only connect from reactants
      reactantList.forEach(reactant => {
        if (reactant !== centralSpecies && connectionCount < maxConnections) {
          const reactantId = sanitizeNodeId(reactant);
          const connectionKey = `${reactantId}-${centralId}`;
          if (!addedConnections.has(connectionKey)) {
            def += `${reactantId}["${sanitizeLabel(reactant)}"]:::default -->|${sanitizeLabel(reactionText)}| ${centralId};\n`;
            addedConnections.add(connectionKey);
            connectionCount++;
          }
        }
      });
    }
  });
  
  return def;
};

export default PathwayDiagram;