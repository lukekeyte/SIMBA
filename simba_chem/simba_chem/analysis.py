import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

class Analysis:
    def __init__(self, network):
        
        self.network = network

    ###################
    # PLOT ABUNDANCES #
    ###################
    
    def plot_abundance(self, species_list, 
                   figsize=(8, 6), 
                   xlabel="Time (yr)", ylabel="Abundance (X/H)", 
                   title="Species Abundance Over Time", 
                   xscale='log', yscale='log', 
                   xlim=None, ylim=None, 
                   fontsize_labels=12, fontsize_title=12, fontsize_legend=10, 
                   grid=True, grid_style='--', grid_width=0.5, alpha_lines=1, alpha_grid=1):
        """
        Plot the time evolution of species abundances.
        Parameters:
            species_list (list of str): List of species to plot.
            figsize (tuple): Figure size in inches (width, height).
            xlabel (str): Label for the x-axis.
            ylabel (str): Label for the y-axis.
            title (str): Title of the plot.
            xscale (str): Scale of the x-axis ('linear' or 'log').
            yscale (str): Scale of the y-axis ('linear' or 'log').
            xlim (tuple): Limits for the x-axis (xmin, xmax).
            ylim (tuple): Limits for the y-axis (ymin, ymax).
            fontsize_labels (int): Font size for x and y axis labels.
            fontsize_title (int): Font size for the plot title.
            fontsize_legend (int): Font size for the legend.
            grid (bool): Whether to display a grid.
            grid_style (str): Line style for the grid.
            grid_width (float): Line width for the grid.
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        times_yr = self.network.time_points / self.network.parameters.yr_sec
        species  = self.network.species.name
        for sp in species_list:
            if sp in species:
                idx = species.index(sp)
                abundance = self.network.abundance_history[:, idx] / self.network.gas.n_gas 
                ax.plot(times_yr, abundance, label=sp, alpha=alpha_lines)
            else:
                print(f"Warning: {sp} not found in species list.")

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        ax.set_xlabel(xlabel, fontsize=fontsize_labels)
        ax.set_ylabel(ylabel, fontsize=fontsize_labels)
        ax.set_title(title, fontsize=fontsize_title, pad=10)
        ax.legend(fontsize=fontsize_legend)
        
        if grid:
            ax.grid(True, which='both', linestyle=grid_style, linewidth=grid_width, alpha=alpha_grid)
        
        return fig



    #######################
    # PLOT REACTION RATES #
    #######################
    
    def plot_reaction_rate(self, reaction_id, 
                       figsize=(8, 6), 
                       xlabel="Time (yr)", ylabel="Reaction Rate (s$^{-1}$)", 
                       title="Reaction Rate Over Time", 
                       xscale='log', yscale='log', 
                       xlim=None, ylim=None, 
                       fontsize_labels=12, fontsize_title=14, fontsize_legend=10, 
                       grid=True, grid_style='--', grid_width=0.5, alpha_lines=1, alpha_grid=1):
        """
        Plot the time evolution of one or more reaction rates.

        Parameters:
            reaction_id (int or list of int): ID(s) of the reaction(s) to plot.
            figsize (tuple): Figure size in inches (width, height).
            xlabel (str): Label for the x-axis.
            ylabel (str): Label for the y-axis.
            title (str): Title of the plot.
            xscale (str): Scale of the x-axis ('linear' or 'log').
            yscale (str): Scale of the y-axis ('linear' or 'log').
            xlim (tuple): Limits for the x-axis (xmin, xmax).
            ylim (tuple): Limits for the y-axis (ymin, ymax).
            fontsize_labels (int): Font size for x and y axis labels.
            fontsize_title (int): Font size for the plot title.
            fontsize_legend (int): Font size for the legend.
            grid (bool): Whether to display a grid.
            grid_style (str): Line style for the grid.
            grid_width (float): Line width for the grid.
        """

        times_yr = self.network.time_points / self.network.parameters.yr_sec
        rates    = self.network.rate_history
        labels   = self.network.reactions.labels


        # Handle single reaction ID or list of reaction IDs
        if isinstance(reaction_id, int):
            reaction_id = [reaction_id]  # Convert to list for uniform handling

        # Validate reaction IDs
        for rid in reaction_id:
            if rid < 0 or rid >= len(labels):
                raise ValueError(f"Reaction ID {rid} is out of range. Available range: 0-{len(labels) - 1}")

        # Create the plot
        plt.figure(figsize=figsize)
        for rid in reaction_id:
            reaction_rate = rates[:, rid]
            reaction_label = labels[rid]
            plt.plot(times_yr, reaction_rate, label=reaction_label, alpha=alpha_lines)

        # Apply user customizations
        plt.xscale(xscale)
        plt.yscale(yscale)
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.xlabel(xlabel, fontsize=fontsize_labels)
        plt.ylabel(ylabel, fontsize=fontsize_labels)
        plt.title(title, fontsize=fontsize_title)
        plt.legend(fontsize=fontsize_legend)

        if grid:
            plt.grid(True, which='both', linestyle=grid_style, linewidth=grid_width, alpha=alpha_grid)

        plt.show()



    ##################
    # PLOT DASHBOARD #
    ##################
    
    def plot_dashboard(self, output_file=None):
        """
        Create a dashboard for SIMBA chemical network results.
        
        Parameters:
            ln (SIMBA): A SIMBA instance that has already run solve_network()
            output_file (str, optional): Path to save the dashboard figure
            
        Returns:
            matplotlib.figure.Figure: The dashboard figure object
        """
        
        # # Check if solve_network has been run
        # if not hasattr(self.network, 'abundance_history') or len(self.network.abundance_history) == 0:
        #     raise ValueError("No results found. Run solve_network() first.")
        
        # Create figure
        fig = plt.figure(figsize=(20, 13))
        gs = GridSpec(3, 3, figure=fig, height_ratios=[3, 3.25, 2], hspace=0.3, wspace=0.3)
        
        # 1. Parameters Panel 
        ax_params = fig.add_subplot(gs[0, :2])
        self.dash_parameters_panel(ax_params, self.network)
        
        # 2. Quasi Steady-State Analysis
        ax_qss = fig.add_subplot(gs[0, 2])
        self.dash_quasi_steady_state_analysis(ax_qss, self.network)
        
        # 3. Abundance Evolution
        ax_abund = fig.add_subplot(gs[1, 0])
        self.dash_abundance_evolution(ax_abund, self.network)
        
        # 4. Reaction Rates
        ax_rates = fig.add_subplot(gs[1, 1])
        self.dash_reaction_rates(ax_rates, self.network)
        
        # 5. Species Correlation Matrix
        ax_corr = fig.add_subplot(gs[1, 2])
        self.dash_species_correlation(ax_corr, self.network)
        
        # 6. Final abundances
        ax_abu_final = fig.add_subplot(gs[2, :])
        self.dash_final_abundances(ax_abu_final, self.network)
        
        plt.tight_layout()
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')

        return fig


    ###########################
    # Functions for dashboard #
    ###########################

    def dash_parameters_panel(self, ax, ln):
        """Plot the model input parameters."""
        ax.axis('on')
        
        # Column 1: Model Parameters

        params = [
            ("Gas density", f"{ln.gas.n_gas:.2e} cm^-3"),
            ("Gas temperature", f"{ln.gas.t_gas:.1f} K"),
            ("Dust density", f"{ln.dust.n_dust:.2e} cm^-3"),
            ("Dust temperature", f"{ln.dust.t_dust:.1f} K"),
            ("Gas/dust ratio", f"{ln.environment.gtd:.1e}"),
            ("Visual extinction", f"{ln.environment.Av:.2f} mag"),
            ("UV field", f"{ln.environment.G_0:.1e} G₀"),
            ("X-ray ionization rate", f"{ln.environment.Zeta_X:.2e} s^-1"),
            ("Cosmic ray ionization rate", f"{ln.environment.Zeta_CR:.2e} s^-1"),
            ("Chemical evolution timescale", f"{ln.parameters.time_final/ln.parameters.yr_sec:.1e} years"),
            ("Self-shielding", str(ln.parameters.self_shielding))
        ]
        
        ax.text(0.02, 0.925, "Model Parameters", ha='left', va='top', fontsize=10, fontweight='bold')
        
        y_pos = 0.8
        for name, value in params:
            ax.text(0.02, y_pos, name, ha='left', va='center', fontsize=9)
            ax.text(0.35, y_pos, value, ha='right', va='center', fontsize=9)
            y_pos -= 0.07
        
        
        # Column 2: Solver Performance
        
        ax.text(0.57, 0.925, "Solver Performance", ha='left', va='top', fontsize=10, fontweight='bold')
        
        # Check if we have solution data
        if hasattr(ln, 'time_points') and len(ln.time_points) > 0:
            # Extract performance data
            num_timesteps = len(ln.time_points)
            time_range = (ln.time_points[0], ln.time_points[-1])
            time_span_years = (time_range[1] - time_range[0]) / ln.parameters.yr_sec
            
            # Calculate step statistics
            if num_timesteps > 1:
                step_sizes = np.diff(ln.time_points)
                min_step = np.min(step_sizes) / ln.parameters.yr_sec  # in years
                max_step = np.max(step_sizes) / ln.parameters.yr_sec  # in years
                median_step = np.median(step_sizes) / ln.parameters.yr_sec  # in years
            else:
                min_step = max_step = median_step = 0
        else:
            metrics = [("No solver data available", "")]
            
        metrics = [
        ("Time Steps", f"{num_timesteps}"),
        ("Time Range", f"{time_range[0]/ln.parameters.yr_sec:.2e} to {time_range[1]/ln.parameters.yr_sec:.2e} years"),
        ("Total Evolution", f"{time_span_years:.2e} years"),
        ("Min Step Size", f"{min_step:.2e} years"),
        ("Max Step Size", f"{max_step:.2e} years"),
        ("Median Step Size", f"{median_step:.2e} years")]
        
        y_pos = 0.8
        for name, value in metrics:
            ax.text(0.57, y_pos, name, ha='left', va='center', fontsize=9)
            ax.text(0.87, y_pos, value, ha='right', va='center', fontsize=9)
            y_pos -= 0.07

        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xticks([])
        ax.set_yticks([])
        


    def dash_abundance_evolution(self, ax, ln):
        """
        Plot abundance evolution of top 10 species.
        """

        times_yr    = ln.time_points / ln.parameters.yr_sec
        top_indices = np.argsort(ln.abundance_history[-1])[-10:][::-1]
        
        for i, idx in enumerate(top_indices):
            species_name = ln.species.name[idx]
            abundance = ln.abundance_history[:, idx] / ln.gas.n_gas  
            ax.loglog(times_yr, abundance, label=species_name, alpha=0.8)
        
        ax.set_xlabel('Time (years)', fontsize=9)
        ax.set_ylabel('Abundance (X/H)', fontsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.set_title('Abundance Evolution of Top 10 Species', fontsize=10)
        ax.legend(loc='lower right', fontsize=8, ncol=2)
        ax.grid(True, which='both', linestyle='--', alpha=0.6)


    def dash_reaction_rates(self, ax, ln):
        """
        Plot reaction rates vs. time for ten most efficient reactions.
        """

        times_yr    = ln.time_points / ln.parameters.yr_sec
        top_indices = np.argsort(ln.rate_history[-1])[-30:][::-1]
        
        for idx in top_indices:
            reaction_label = ln.reactions.labels[idx]
            if len(reaction_label) > 30:
                reaction_label = reaction_label[:27] + '...'
            
            rate = ln.rate_history[:, idx]
            ax.loglog(times_yr, rate, label=reaction_label)
        
        ax.set_xlabel('Time (years)', fontsize=9)
        ax.set_ylabel('Reaction Rate (cm$^{-3}$ s$^{-1}$)', fontsize=9)
        ax.set_title('Top 10 Most Efficient Reactions', fontsize=10)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.legend(loc='lower right', fontsize=8, ncol=1)
        ax.grid(True, which='both', linestyle='--', alpha=0.6)


    def dash_species_correlation(self, ax, ln):
        """Plot species correlation matrix for top species."""
        
        # Get final abundances for top 15 species
        top_indices   = np.argsort(ln.abundance_history[-1])[-15:][::-1]
        species_names = [ln.species.name[i] for i in top_indices]
        abundances    = np.array([ln.abundance_history[:, i] for i in top_indices])
        
        # Compute correlation matrix
        corr_matrix    = np.corrcoef(abundances)
        display_matrix = corr_matrix.copy() # copy of to modify for display
        
        # Plotting
        im = ax.imshow(display_matrix, cmap='PRGn', vmin=-1, vmax=1, alpha=0.6)
        
        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.95, aspect=30)
        cbar.set_label('Correlation Coefficient', fontsize=10)
        cbar.ax.tick_params(labelsize=9)
        
        # Add correlation values as text
        for i in range(len(species_names)):
            for j in range(len(species_names)):
                if i != j:
                    text_color = 'white' if abs(corr_matrix[i, j]) > 0.7 else 'black'
                    ax.text(j, i, f"{corr_matrix[i, j]:.2f}", 
                            ha="center", va="center", color=text_color, fontsize=7)
                    
        # Add black squares for diagonal elements
        for i in range(len(species_names)):
            rect = plt.Rectangle((i-0.5, i-0.5), 1, 1, fill=True, color='black', alpha=1.0)
            ax.add_patch(rect)
        
        ax.set_xticks(np.arange(len(species_names)))
        ax.set_yticks(np.arange(len(species_names)))
        ax.set_xticklabels(species_names, rotation=90, fontsize=9)
        ax.set_yticklabels(species_names, fontsize=9)
        ax.set_title('Species Abundance Correlation Matrix', fontsize=10, pad=5)
        
        
    def dash_final_abundances(self, ax, ln):
        """Bar chart showing final abundance of all species (X/H)."""
        
        species    = ln.species.name
        n_gas      = ln.gas.n_gas
        abundances = ln.abundance_history[-1, :]
        
        ax.bar(species, abundances/n_gas, color='darkslateblue', alpha=0.7)
        
        ax.set_xticks(np.arange(len(species)))
        ax.set_xticklabels(species, rotation=90, fontsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.set_xlim(-1.2, len(species) - 0.2)
        ax.set_ylim(1e-20, 9)
        ax.set_yscale('log')
        ax.set_title('Time to Quasi-Steady State', fontsize=10)
        ax.set_xlabel('Time (years)', fontsize=9)
        ax.set_ylabel('Abundance (X/H)', fontsize=9)
        ax.set_title('Species Final Abundances', pad=7, fontsize=10)
        
        
    def dash_quasi_steady_state_analysis(self, ax, ln, threshold=0.1, min_points=3, skip_points=1, colormap='viridis_r', num_species=10):
        """
        Plot a visual timeline of when key species reach quasi-steady state.
        
        Parameters:
            ax (matplotlib.axes): The axis to plot on
            ln (SIMBA): A SIMBA instance with solved network
            threshold (float): Change threshold for QSS (fraction change per decade)
            min_points (int): Minimum consecutive points needed below threshold
            skip_points (int): Number of initial points to skip
            colormap (str or matplotlib.colors.Colormap): Colormap to use for timeline bars
                     Can be any matplotlib colormap name ('viridis', 'plasma', etc.) or a colormap object
            num_species (int): Number of most abundant species to show (default: 10)
        """

        ax.axis('off')
        
        if not hasattr(ln, 'abundance_history') or len(ln.abundance_history) == 0:
            ax.text(0.5, 0.5, "No solver data available", 
                    ha='center', va='center', fontsize=10)
            ax.set_title('Quasi-Steady State Analysis', fontsize=10)
            return
            
        # Get abundance history and time points
        abundances    = ln.abundance_history
        times         = ln.time_points / ln.parameters.yr_sec  # Convert to years
        species_names = ln.species.name
        max_time      = times[-1]
        
        # Top N species
        top_indices      = np.argsort(abundances[-1])[-(num_species):][::-1]  
        selected_species = [species_names[idx] for idx in top_indices]
        top_indices      = top_indices[::-1]      # reverse order
        selected_species = selected_species[::-1]
        
        # Calculate quasi steady-state times
        qss_times   = {}
        not_reached = []
        
        for i, idx in enumerate(top_indices):
            species   = species_names[idx]
            abundance = abundances[:, idx]
            
            # Skip trace species
            if np.max(abundance) < 1e-30:
                qss_times[species] = times[0] 
                continue
                
            rel_changes = np.zeros(len(abundance)-1)
            for i in range(skip_points, len(abundance)-1):
                if abundance[i] > 1e-20:  
                    # Calculate fractional change per decade in time
                    time_ratio = times[i+1]/times[i]
                    if time_ratio > 1:
                        decade_factor = np.log10(time_ratio)
                        rel_changes[i] = abs(abundance[i+1] - abundance[i]) / (abundance[i] * decade_factor)
            
            is_flat = rel_changes < threshold
            
            reached_qss = False
            for t in range(skip_points, len(rel_changes) - min_points + 1):
                if np.all(is_flat[t:t+min_points]):
                    qss_times[species] = times[t+1]
                    reached_qss = True
                    break
                    
            initial_flat = np.all(is_flat[skip_points:skip_points+min_points])
            if initial_flat and not reached_qss:
                qss_times[species] = times[skip_points]
                reached_qss = True
                
            if not reached_qss:
                qss_times[species] = max_time * 1.1  # Position beyond timeline
                not_reached.append(species)
        
        # Plotting
        
        bar_height = 0.7  
        
        if isinstance(colormap, str):
            cmap = plt.cm.get_cmap(colormap)
        else:
            cmap = colormap  
        
        log_min_time = np.log10(times[0])
        log_max_time = np.log10(max_time)
        log_time_range = log_max_time - log_min_time

        
        for i, species in enumerate(selected_species):
            qss_time = qss_times[species]
            
            ax.text(times[0]*0.6, i, species, ha='right', va='center', fontsize=10)

            if species not in not_reached:
                log_qss_time = np.log10(qss_time)
                norm_time = (log_qss_time - log_min_time) / log_time_range
                norm_time = np.clip(norm_time, 0, 1)
            else:
                norm_time = 1.0
                
            color = cmap(norm_time)
            
            ax.barh(i, qss_time, height=bar_height, left=times[0], 
                    color=color, alpha=0.7)
            
            if species not in not_reached:
                time_str = f"{qss_time:.1e} yr" if qss_time > times[0] else "initial"
                bar_width_log = np.log10(qss_time) - np.log10(times[0])
                
                if bar_width_log > 1.5:
                    ax.text(qss_time*0.6, i, time_str, ha='right', va='center', fontsize=9, color='white', zorder=10)
                else:
                    ax.text(qss_time*2, i, time_str, ha='left', va='center', fontsize=9, color='dimgrey', zorder=10)
            else:
                ax.text(max_time*1.2, i, "Not reached", ha='left', va='center', fontsize=8, style='italic')
        
        ax.set_xlabel('Time (years)', fontsize=9)
        ax.set_xscale('log')
        ax.set_xlim(times[0]*0.5, max_time*1.2)
        ax.set_ylim(-0.5, len(selected_species)-0.5)
        ax.set_title('Time to Quasi-Steady State', fontsize=10, pad=12)
        
        ax.text(0.5, 1.01, f"(<{threshold*100:.1f}% change per decade for {min_points} consecutive points)", 
                transform=ax.transAxes, fontsize=7, style='italic', ha='center')


    ###############
    # EXPORT DATA #
    ###############

    def export_abundance_data(self, file_name='abundance_data.csv'):
        """
        Export abundance data to a CSV file.
        
        Parameters:
            file_name (str): Name of the file to save data to.
        """
        
        time = self.results['time'] / 3.1556926e7  # convert time to years
        species = self.results['species']
        abundances = self.results['abundances']

        # Normalize abundances by n_gas
        n_gas = abundances[species.index('H')] + abundances[species.index('H2')]
        normalized_abundances = abundances / n_gas

        df = pd.DataFrame(data=normalized_abundances.T, columns=species)
        df.insert(0, 'Time (years)', time)

        # Format the 'Time (s)' column to scientific notation with 3 decimal places
        df['Time (years)'] = df['Time (years)'].apply(lambda x: f"{x:.3e}")
        
        # Save the dataframe to CSV
        df.to_csv(file_name, index=False)
        print(f"Abundance data saved to {file_name}")
