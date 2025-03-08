    def plot_combined_annotation_overlap(self, annotation1, annotation2):
        """
        Identifies spots that overlap with both specified annotations and plots them.
        Spots overlapping both annotations are colored red, all other spots are colored blue.
        
        Parameters:
        - annotation1: First annotation to check for overlap
        - annotation2: Second annotation to check for overlap
        
        Returns:
        - A matplotlib figure showing the spots with appropriate coloring
        """
        if self.overlay_computed is False:
            raise ValueError("Error: compute annotation overlay with Visium spots using compute_overlap_annotation() first")
        
        # Check if the annotations exist in the data
        available_annotations = self.data_exploded['Annotation'].unique()
        if annotation1 not in available_annotations or annotation2 not in available_annotations:
            raise ValueError(f"One or both annotations not found. Available annotations: {', '.join(available_annotations)}")
        
        # Create a copy of the spot data
        gdf = self.spot_geodata.copy()
        
        # Create a new column indicating whether a spot overlaps with both annotations
        overlap_col1 = str(annotation1) + '_overlap'
        overlap_col2 = str(annotation2) + '_overlap'
        
        # Check if the overlap columns exist
        if overlap_col1 not in gdf.columns or overlap_col2 not in gdf.columns:
            raise ValueError(f"Overlap data for one or both annotations not found. Please ensure compute_overlap_annotation() was called.")
        
        # Create a boolean column indicating spots that overlap with both annotations
        gdf['dual_overlap'] = (gdf[overlap_col1] > 0) & (gdf[overlap_col2] > 0)
        
        # Plot the results
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Plot spots that don't overlap with both annotations in blue
        non_overlapping = gdf[~gdf['dual_overlap']]
        if not non_overlapping.empty:
            non_overlapping.plot(ax=ax, color='blue', alpha=0.5, label='Non-overlapping spots')
        
        # Plot spots that overlap with both annotations in red
        overlapping = gdf[gdf['dual_overlap']]
        if not overlapping.empty:
            overlapping.plot(ax=ax, color='red', alpha=0.7, label=f'Spots overlapping {annotation1} & {annotation2}')
        
        # Add legend, title and labels
        ax.legend()
        ax.set_title(f'Spots Overlapping Both {annotation1} and {annotation2}', fontsize=15)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        
        # Add annotation statistics
        spot_count = len(gdf)
        overlap_count = len(overlapping)
        overlap_percentage = (overlap_count / spot_count) * 100 if spot_count > 0 else 0
        
        stats_text = (
            f"Total spots: {spot_count}\n"
            f"Spots overlapping both annotations: {overlap_count} ({overlap_percentage:.2f}%)"
        )
        
        # Place the statistics text box in the upper right corner
        ax.text(0.98, 0.98, stats_text, 
                transform=ax.transAxes, 
                horizontalalignment='right',
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.show()
        
        return fig, ax, overlapping