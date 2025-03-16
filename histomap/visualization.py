import geopandas as gpd
import pandas as pd
import ast
import gzip
import io
import zipfile
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
import spatialdata
import spatialdata_io
import geopandas as gpd
from shapely.geometry import Polygon, Point
from rtree import index
import numpy as np
from . import histomap_utils
import seaborn as sns 

def plot_cells(histomap, fill=True, contour='k', display_image=True):
    """
    Plots segmented cells based on a GeoDataFrame and optionally displays an image beneath the cells.

    Parameters:
    - histomap: HistoMap object containing the segmentation dataframe and spot geodata.
    - fill: False, True, or a list of colors.
        - False: No fill (only contours).
        - True: Uses a default colormap ('tab20') for fill.
        - List of colors: Specifies the fill color for each cell.
    - contour: Color or list of colors for the contours. If None, the default colormap is used.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    """
    # Check if 'n_cell' exists in spot_geodata
    if 'n_cell' not in histomap.spot_geodata.columns:
        raise ValueError("'n_cell' column does not exist in spot_geodata. Please ensure segmentation has been added correctly.")
    
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf = histomap.segmentation_dataframe
    
    # Default colormap
    cmap = plt.colormaps.get_cmap('tab20')
    unique_cells = len(gdf)

    # Handle fill colors
    if fill is True:
        fill = [cmap(idx % 20) for idx in range(unique_cells)]
    elif isinstance(fill, list):
        if len(fill) != unique_cells:
            raise ValueError("The length of the 'fill' list must match the number of cells.")
    elif fill is False:
        fill = [None] * unique_cells  # Ensure correct indexing

    # Handle contour colors
    if contour is None:
        contour = [cmap(idx % 20) for idx in range(unique_cells)]
    elif isinstance(contour, str):
        contour = [contour] * unique_cells  # Convert single color to list
    elif isinstance(contour, list) and len(contour) != unique_cells:
        raise ValueError("The length of the 'contour' list must match the number of cells.")

    # Display the image beneath the cells (if requested and available)
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    # Plot cells
    for idx, (geom, fill_color, contour_color) in enumerate(zip(gdf.geometry, fill, contour)):
        if geom.is_empty or not geom.is_valid:
            continue

        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            if fill_color:
                ax.fill(x, y, color=fill_color, edgecolor=contour_color)
            else:
                ax.plot(x, y, color=contour_color)
        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom.geoms:
                x, y = polygon.exterior.xy
                if fill_color:
                    ax.fill(x, y, color=fill_color, edgecolor=contour_color)
                else:
                    ax.plot(x, y, color=contour_color)

    ax.set_title('Segmented Cells')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', adjustable='box')

    # Adjust limits
    bounds = gdf.total_bounds
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()


def plot_combined_annotation_overlap(histomap, annotation1, annotation2):
    """
    Identifies spots that overlap with both specified annotations and plots them.
    Spots overlapping both annotations are colored red, all other spots are colored blue.
    
    Parameters:
    - annotation1: First annotation to check for overlap
    - annotation2: Second annotation to check for overlap
    
    Returns:
    - A matplotlib figure showing the spots with appropriate coloring
    """

    # Check if the annotations exist in the data
    available_annotations = histomap.data_exploded['Annotation'].unique()
    if annotation1 not in available_annotations or annotation2 not in available_annotations:
        raise ValueError(f"One or both annotations not found. Available annotations: {', '.join(available_annotations)}")
    
    # Create a copy of the spot data
    gdf = histomap.spot_geodata.copy()
    
    # Create a new column indicating whether a spot overlaps with both annotations
    overlap_col1 = str(annotation1) + '_overlap'
    overlap_col2 = str(annotation2) + '_overlap'
    
    # Check if the overlap columns exist
    if overlap_col1 not in gdf.columns or overlap_col2 not in gdf.columns:
        raise ValueError(f"Overlap data for one or both annotations not found. Please ensure compute_overlap_annotation() was called.")
    
    # Create a boolean column indicating spots that overlap with both annotations
    gdf['dual_overlap'] = (gdf[overlap_col1] > 0) & (gdf[overlap_col1] < 100) & (gdf[overlap_col2] > 0)& (gdf[overlap_col2] < 100)
    
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



def plot_annotation_overlay(histomap, annotation, disabled_annotations=None, display_image=True, max_cutoff=None):
    if isinstance(annotation, str):
        annotation = [annotation]
    
    # Handle disabled annotations
    if disabled_annotations is None:
        disabled_annotations = []
    
    # Filter out disabled annotations
    annotation = [ann for ann in annotation if ann not in disabled_annotations]
    
    # Check that all required annotation columns exist
    missing_columns = [ann + "_overlap" for ann in annotation if ann + "_overlap" not in histomap.spot_geodata.columns]
    
    if missing_columns:
        raise ValueError(f"Error: Missing columns {missing_columns}. Compute annotation overlay with Visium spot using compute_overlap_annotation() first.")
    
    # Create a subplot grid with one row and as many columns as there are annotations
    num_annotations = len(annotation)
    fig, axes = plt.subplots(1, num_annotations, figsize=(15, 10))

    if num_annotations == 1:
        axes = [axes]  # Ensure axes is always iterable if there's only one annotation

    for ax, ann in zip(axes, annotation):
        gdf = histomap.spot_geodata.copy()
        
        # Convert annotation overlap columns to float if necessary
        gdf[ann + "_overlap"] = gdf[ann + "_overlap"].astype(float)
        
        # Display the image beneath the annotations
        if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
            extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
            ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

        # If max_cutoff is provided, set vmax to it, otherwise, use the maximum value in the data
        vmin = gdf[ann + "_overlap"].min()
        vmax = max_cutoff if max_cutoff is not None else gdf[ann + "_overlap"].max()

        # Plot the polygons with a colormap based on the 'Positive_overlap'
        gdf.plot(ax=ax, column=ann + "_overlap", cmap='viridis', legend=True,
                 legend_kwds={'label': f"{ann} Positive Overlap (%)",
                              'orientation': "horizontal"},
                 vmin=vmin, vmax=vmax)  # Apply max_cutoff (vmax) to the color scale

        # Set title and axis labels
        ax.set_title(f'Spot Plot Colored by {ann} Positive Overlap', fontsize=15)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)

    plt.tight_layout()
    plt.show()




def plot_annotations(histomap, fill=False, contour=None, annotation=None, display_image=False, alpha=1):
    """Plots the annotations based on the DataFrame, respecting the plot order.
    
    Parameters:
    - histomap: HistoMap object containing annotation data.
    - fill: False, True, or a list of colors. 
            False means no fill (only contours), 
            True uses colors from histomap.annotation_colors, 
            a list of colors specifies the fill color for each annotation.
    - contour: a color or list of colors for the contours. 
               If None, uses colors from histomap.annotation_colors or default colormap.
    - annotation: a specific annotation (string) or a list of annotations to plot. If None, all annotations are plotted.
    - display_image: bool, whether to display `histomap.plotting_image` beneath the annotations.
    - alpha: float between 0 and 1, transparency of the fill color (0 is completely transparent, 1 is opaque).
             Only applies when fill is not False.
    """
    # Validate alpha parameter
    if not (0 <= alpha <= 1):
        raise ValueError("Alpha must be between 0 and 1")
    
    fig, ax = plt.subplots(figsize=(10, 10))

    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    cmap = plt.cm.get_cmap('tab20')
    unique_annotations = list(histomap.data_exploded['Annotation'].unique())

    # Validate annotation input
    if annotation is not None:
        if isinstance(annotation, str):
            annotation = [annotation]
        elif not isinstance(annotation, list) or not all(isinstance(a, str) for a in annotation):
            raise TypeError("Annotation must be a string, a list of strings, or None.")

        # Check for missing annotations
        missing_annotations = [ann for ann in annotation if ann not in unique_annotations]
        if missing_annotations:
            raise ValueError(
                f"Annotations {missing_annotations} not found. Available annotations are: {sorted(unique_annotations)}"
            )
        annotations_to_plot = annotation
    else:
        annotations_to_plot = unique_annotations

    annotations_to_plot = [ann for ann in annotations_to_plot if ann not in histomap.disabled_annotations]

    if histomap.disabled_annotations:
        print(f"Skipping disabled annotations: {', '.join(histomap.disabled_annotations)}")

    if len(histomap.activated_annotations) == 0:
        raise ValueError("No activated annotations to compute overlap.")

    # Process fill colors
    if fill is True and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        fill_dict = dict(zip(histomap.annotation_colors['annotation'], 
                            histomap.annotation_colors['color']))
        fill = [fill_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif fill is True and not hasattr(histomap, 'annotation_colors'):
        # Fall back to default colormap if annotation_colors doesn't exist
        fill = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(fill, list):
        if len(fill) != len(annotations_to_plot):
            raise ValueError("The length of the 'fill' list must match the number of unique annotations.")
    elif fill is False:
        fill = [None] * len(annotations_to_plot)
    
    # Process contour colors
    if contour is None and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        contour_dict = dict(zip(histomap.annotation_colors['annotation'], 
                               histomap.annotation_colors['color']))
        contour = [contour_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif contour is None:
        contour = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(contour, list):
        if len(contour) != len(annotations_to_plot):
            raise ValueError("The length of the 'contour' list must match the number of unique annotations.")
    elif isinstance(contour, str):
        contour = [contour] * len(annotations_to_plot)
        
    # Sort annotations based on plot_order (0 on top)
    data_sorted = histomap.data_exploded.sort_values(by="plot_order", ascending=False)
    
    # Plot annotations
    for idx, row in data_sorted.iterrows():
        geom = row['geometry']
        ann = row['Annotation']

        if ann not in annotations_to_plot or geom.is_empty or not geom.is_valid:
            continue
        
        ann_idx = annotations_to_plot.index(ann)
        fill_color = fill[ann_idx] if isinstance(fill, list) else None
        contour_color = contour[ann_idx]

        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            if fill_color:
                ax.fill(x, y, color=fill_color, edgecolor=contour_color, alpha=alpha, label=ann)
            else:
                ax.plot(x, y, color=contour_color, label=ann)
        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom:
                x, y = polygon.exterior.xy
                if fill_color:
                    ax.fill(x, y, color=fill_color, edgecolor=contour_color, alpha=alpha, label=ann)
                else:
                    ax.plot(x, y, color=contour_color, label=ann)

    ax.set_title('Annotations', fontsize=14)
    ax.set_ylim(ax.get_ylim()[::-1])  
    ax.set_aspect('auto')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 1))

    # Adjust plot limits
    bounds = histomap.data_exploded.total_bounds
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()





def plot_tissue_overlap(histomap, cmap="coolwarm", figsize=(10, 10), display_image=True):
    """
    Plot the tissue detection overlap using spot_geodata geometries.

    Parameters:
    - histomap (HistoMap): HistoMap object with the geospatial data containing the 'tissue_detection' column.
    - cmap (str): Colormap for visualization.
    - figsize (tuple): Figure size.

    Returns:
    - Matplotlib plot showing the tissue detection values.
    """

    if "tissue_detection" not in histomap.spot_geodata.columns:
        raise ValueError("Column 'tissue_detection' not found in spot_geodata. Run compute_tissue_overlap first.")



    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    histomap.spot_geodata.plot(column="tissue_detection", cmap=cmap, linewidth=0.1, edgecolor="black", legend=True, ax=ax)

    # Formatting
    ax.set_title("Tissue Detection Overlap", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    plt.show()

def violin_tissue_overlap(histomap, figsize=(8, 6)):
    """
    Plot a violin + boxplot of the tissue detection values, with individual spot points.

    Parameters:
    - histomap (HistoMap): histomap with the geospatial data containing the 'tissue_detection' column.
    - figsize (tuple): Figure size.

    Returns:
    - Violin + boxplot + strip plot visualization of tissue detection values.
    """

    if "tissue_detection" not in histomap.spot_geodata.columns:
        raise ValueError("Column 'tissue_detection' not found in spot_geodata. Run compute_tissue_overlap first.")

    # Set seaborn style
    sns.set_style("whitegrid")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Violin plot with refined colors and transparency
    sns.violinplot(
        y=histomap.spot_geodata["tissue_detection"], 
        inner=None, 
        color="lightblue", 
        linewidth=0.8, 
        alpha=0.7, 
        ax=ax
    )

    # Boxplot overlay with stronger visibility
    sns.boxplot(
        y=histomap.spot_geodata["tissue_detection"], 
        width=0.15, 
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 1.2}, 
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "black", "linewidth": 1.2},
        capprops={"color": "black", "linewidth": 1.2},
        flierprops={"marker": "o", "markerfacecolor": "red", "markeredgecolor": "black", "markersize": 4},
        ax=ax
    )

    # Strip plot to add individual data points
    sns.stripplot(
        y=histomap.spot_geodata["tissue_detection"], 
        color="black", 
        size=3, 
        alpha=0.5, 
        jitter=True, 
        ax=ax
    )

    # Formatting
    ax.set_title("Distribution of Tissue Detection Overlap", fontsize=14, fontweight="bold")
    ax.set_ylabel("Tissue Detection (%)", fontsize=12)
    ax.set_xlabel("")  # No x-label needed for a single variable
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.show()



def plot_positive_spots(histomap, annotation, display_image=True):
    """Plot spots colored by whether they are positive for the given annotation(s).
    
    Parameters:
    - histomap: HistoMap object containing the spot geodata with annotation_positive columns.
    - annotation: A string or list of strings representing the annotation(s) to plot.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    """
    
    # Ensure annotation is always a list
    if isinstance(annotation, str):
        annotation = [annotation]

    # Ensure the required columns exist
    missing_cols = []
    for ann in annotation:
        positive_col = ann + "_positive"
        if positive_col not in histomap.spot_geodata.columns:
            missing_cols.append(positive_col)
    
    if missing_cols:
        raise ValueError(f"Error: Missing columns {', '.join(missing_cols)}. Compute positive spots first using set_positive().")

    # Create figure and axis
    fig, axes = plt.subplots(1, len(annotation), figsize=(15, 10)) if len(annotation) > 1 else plt.subplots(1, 1, figsize=(10, 10))
    if len(annotation) == 1:
        axes = [axes]  # Ensure axes is iterable for a single plot

    # Get the custom color mapping if available
    color_dict = {}
    if hasattr(histomap, 'annotation_colors'):
        color_dict = dict(zip(histomap.annotation_colors['annotation'], 
                             histomap.annotation_colors['color']))

    for ax, ann in zip(axes, annotation):
        # Copy geodata for each annotation
        gdf = histomap.spot_geodata.copy()
        positive_col = ann + "_positive"
        
        # Display the image beneath the spots
        if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
            extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
            ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

        # Determine colors to use
        negative_color = 'lightgrey'
        positive_color = color_dict.get(ann, 'red') if ann in color_dict else 'red'
        
        # Create separate GeoDataFrames for positive and negative spots
        positive_spots = gdf[gdf[positive_col] == True].copy()
        negative_spots = gdf[gdf[positive_col] == False].copy()
        
        # Plot negative spots
        negative_spots.plot(ax=ax, color=negative_color, label='Negative', alpha=0.6)
        
        # Plot positive spots
        positive_spots.plot(ax=ax, color=positive_color, label=f'{ann} Positive', alpha=0.8)
        
        # Add legend with discrete values
        ax.legend(title=f'{ann}')
        
        # Set title and axis labels
        ax.set_title(f'Spot Plot for {ann} (Positive/Negative)', fontsize=15)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)

    plt.tight_layout()
    plt.show()

def plot_annotation_order(histomap, fill=False, contour=None, annotation=None, display_image=False):
    """Plots the annotations in 3D based on the DataFrame, respecting the plot order.
    
    Parameters:
    - histomap: HistoMap object containing annotation data.
    - fill: False, True, or a list of colors. 
            False means no fill (only contours), 
            True uses the default colormap for fill, 
            a list of colors specifies the fill color for each annotation.
            If None, uses colors from histomap.annotation_colors.
    - contour: a color or list of colors for the contours. 
               If None, uses colors from histomap.annotation_colors or default colormap.
    - annotation: a specific annotation (string) or a list of annotations to plot. If None, all annotations are plotted.
    - display_image: bool, whether to display `histomap.plotting_image` beneath the annotations.
    """
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    cmap = plt.cm.get_cmap('tab20')
    unique_annotations = list(histomap.data_exploded['Annotation'].unique())

    # Validate annotation input
    if annotation is not None:
        if isinstance(annotation, str):
            annotation = [annotation]
        elif not isinstance(annotation, list) or not all(isinstance(a, str) for a in annotation):
            raise TypeError("Annotation must be a string, a list of strings, or None.")

        # Check for missing annotations
        missing_annotations = [ann for ann in annotation if ann not in unique_annotations]
        if missing_annotations:
            raise ValueError(
                f"Annotations {missing_annotations} not found. Available annotations are: {sorted(unique_annotations)}"
            )
        annotations_to_plot = annotation
    else:
        annotations_to_plot = unique_annotations

    annotations_to_plot = [ann for ann in annotations_to_plot if ann not in histomap.disabled_annotations]

    if histomap.disabled_annotations:
        print(f"Skipping disabled annotations: {', '.join(histomap.disabled_annotations)}")

    if len(histomap.activated_annotations) == 0:
        raise ValueError("No activated annotations to compute overlap.")

    ####
    xmax, ymax = histomap.full_res_width, histomap.full_res_height
    # Add a rectable to the size of the image on each annotation 
    rectangle = Polygon([(0, 0), (xmax, 0), (xmax, ymax), (0, ymax)])
    tmp_data_exploded = histomap.data_exploded
    
    # subset only activated annotations
    tmp_data_exploded = tmp_data_exploded[tmp_data_exploded['Annotation'].isin(histomap.activated_annotations)]
    
    # Get unique annotations
    unique_annotations = tmp_data_exploded['Annotation'].unique()
    
    # Create a dictionary mapping annotations to their plot_order
    annotation_plot_order = tmp_data_exploded.groupby("Annotation")["plot_order"].first().to_dict()
    
    # Create a new GeoDataFrame with the rectangle for each unique annotation
    rectangle_gdf = gpd.GeoDataFrame({
        'geometry': [rectangle] * len(unique_annotations),  # List of rectangles
        'Annotation': unique_annotations,  # Corresponding annotations
        'plot_order': [annotation_plot_order[ann] for ann in unique_annotations]  # Matching plot_order
    })
    
    # Concatenate the new GeoDataFrame (rectangle_gdf) with the existing one (tmp_data_exploded)
    gdf = pd.concat([tmp_data_exploded, rectangle_gdf], ignore_index=True)
    
    # Process fill colors
    if fill is None and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        fill_dict = dict(zip(histomap.annotation_colors['annotation'], 
                            histomap.annotation_colors['color']))
        fill = [fill_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif fill is True:
        fill = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(fill, list):
        if len(fill) != len(annotations_to_plot):
            raise ValueError("The length of the 'fill' list must match the number of unique annotations.")
    elif fill is False:
        fill = [None] * len(annotations_to_plot)
    
    # Process contour colors
    if contour is None and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        contour_dict = dict(zip(histomap.annotation_colors['annotation'], 
                               histomap.annotation_colors['color']))
        contour = [contour_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif contour is None:
        contour = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(contour, list):
        if len(contour) != len(annotations_to_plot):
            raise ValueError("The length of the 'contour' list must match the number of unique annotations.")
    elif isinstance(contour, str):
        contour = [contour] * len(annotations_to_plot)

    # Plot annotations in 3D
    for idx, row in gdf.iterrows():
        geom = row['geometry']
        ann = row['Annotation']
        plot_order = row['plot_order']

        if ann not in annotations_to_plot or geom.is_empty or not geom.is_valid:
            continue
        
        ann_idx = annotations_to_plot.index(ann)
        fill_color = fill[ann_idx] if isinstance(fill, list) else None
        contour_color = contour[ann_idx]

        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            z = np.full_like(x, plot_order)  # Create a constant z based on plot_order
            # Swap x and z
            if fill_color:
                ax.fill(z, x, y, color=fill_color, edgecolor=contour_color, alpha=0.3, label=ann)
            else:
                ax.plot(z, x, y, color=contour_color, label=ann)
        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom:
                x, y = polygon.exterior.xy
                z = np.full_like(x, plot_order)  # Create a constant z based on plot_order
                # Swap x and z
                if fill_color:
                    ax.fill(z, x, y, color=fill_color, edgecolor=contour_color, alpha=0.3, label=ann)
                else:
                    ax.plot(z, x, y, color=contour_color, label=ann)

    ax.set_title('Annotations', fontsize=14)
    ax.set_xlabel('Plot Order')
    ax.set_ylabel('')
    ax.set_zlabel('')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 1))

    # Hide grid lines
    ax.grid(False)
    ax.set_zticks([])
    ax.set_yticks([])
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # Set the ticks on the z-axis to only show the actual plot_order values
    ax.set_xticks(sorted(set(gdf['plot_order'])))
    plt.tight_layout()
    plt.gca().invert_xaxis()
    plt.gca().invert_zaxis()
    plt.show()



def plot_annotation_map(histomap, display_image=True):
    """Plot spots colored by their assigned annotations in the Annotation_map.
    Uses colors from histomap.annotation_colors if available.
    
    Parameters:
    - histomap: HistoMap object containing the spot geodata with 'Annotation_map' column.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    """
    
    # Ensure the required column exists
    if 'Annotation_map' not in histomap.spot_geodata.columns:
        raise ValueError("Error: 'Annotation_map' column not found. Please generate the annotation map first.")
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Display the image beneath the spots
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    # Plot spots, coloring by the annotation in 'Annotation_map'
    gdf = histomap.spot_geodata.copy()
    
    # Check if there are any annotations and plot them
    unique_annotations = gdf['Annotation_map'].dropna().unique()
    if len(unique_annotations) == 0:
        raise ValueError("No annotations found in the 'Annotation_map' column.")
    
    # Get colors from annotation_colors if available, or use a default colormap
    color_dict = {}
    if hasattr(histomap, 'annotation_colors'):
        color_dict = dict(zip(histomap.annotation_colors['annotation'], 
                             histomap.annotation_colors['color']))
        
    # Default colormap as fallback
    cmap = plt.cm.get_cmap('tab20', len(unique_annotations))
    
    for i, annotation in enumerate(unique_annotations):
        # Filter spots for the current annotation
        annotation_spots = gdf[gdf['Annotation_map'] == annotation]
        
        # Extract x and y coordinates from the geometry column (assuming Point geometry)
        if annotation_spots.geometry.geom_type.iloc[0] == 'Point':  # Check if it's Point geometry
            x_coords = annotation_spots.geometry.x
            y_coords = annotation_spots.geometry.y
        else:
            # For non-Point geometries (e.g., Polygon), you might want to adjust this code
            # For example, you could take the centroid of the geometry
            x_coords = annotation_spots.geometry.centroid.x
            y_coords = annotation_spots.geometry.centroid.y

        # Use color from annotation_colors if available, otherwise use default colormap
        spot_color = color_dict.get(annotation, cmap(i))
        
        # Plot the spots for this annotation with a specific color
        ax.scatter(x_coords, y_coords, label=annotation, 
                   color=spot_color, s=20, edgecolors='k', alpha=0.6)
    
    # Set title and axis labels
    ax.set_title('Spots Colored by Annotations in Annotation_map', fontsize=15)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add a legend
    ax.legend(title='Annotations')

    plt.tight_layout()
    plt.show()

def plot_annotation_map_proportions(histomap):
    """
    Plots the proportion of each annotation present in histomap.spot_geodata['Annotation_map'] 
    as a stacked bar plot. Uses colors from histomap.annotation_colors if available.

    Parameters:
    - histomap: An instance of the histomap object.
    """
    # Ensure 'spot_geodata' and 'Annotation_map' exist in the histomap
    if 'Annotation_map' not in histomap.spot_geodata.columns:
        raise ValueError("The histomap object does not contain 'Annotation_map'. Run generate_annotation_map after computing overlap.")

    # Get the annotation data
    annotations = histomap.spot_geodata['Annotation_map']

    # Calculate the count of each annotation
    annotation_counts = annotations.value_counts()

    # Normalize to get proportions
    annotation_proportions = annotation_counts / annotation_counts.sum()
    
    # Create custom colors if annotation_colors is available
    colors = plt.cm.Paired.colors  # Default colors
    
    if hasattr(histomap, 'annotation_colors'):
        # Create a mapping from annotation to color
        color_dict = dict(zip(histomap.annotation_colors['annotation'], 
                              histomap.annotation_colors['color']))
        
        # Create a list of colors for each annotation in the proportion data
        custom_colors = []
        for annotation in annotation_proportions.index:
            # Use the custom color if available, otherwise use a default color
            custom_colors.append(color_dict.get(annotation, plt.cm.Paired(len(custom_colors) % 10)))
        
        # If we have custom colors, use them instead of the default
        if custom_colors:
            colors = custom_colors
    
    # Plotting
    plt.figure(figsize=(10, 6))
    ax = annotation_proportions.plot(kind='bar', color=colors)
    
    # Set plot labels and title
    plt.xlabel('Annotations')
    plt.ylabel('Proportion')
    plt.title('Proportions of Annotations in Annotation_map')
    plt.xticks(rotation=45)
    
    # Add percentage labels on top of each bar
    for i, v in enumerate(annotation_proportions):
        ax.text(i, v + 0.01, f"{v:.1%}", ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.show()

def plot_cell_density(histomap, display_image=True, max_cutoff=None):
    """
    Plots the cell density (n_cell) for each spot in the histomap, optionally displaying an image beneath the cells.

    Parameters:
    - histomap: HistoMap object containing the spot geodata with 'n_cell' column.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    - max_cutoff: Optional, the maximum value for the color scale. If None, the maximum value in 'n_cell' will be used.
    """
    # Check if 'n_cell' exists in spot_geodata
    if 'n_cell' not in histomap.spot_geodata.columns:
        raise ValueError("'n_cell' column does not exist in spot_geodata. Please ensure segmentation has been added correctly.")

    # Create a subplot grid for displaying the plot
    fig, ax = plt.subplots(figsize=(10, 10))

    # Copy the spot_geodata for plotting
    gdf = histomap.spot_geodata.copy()

    # Convert 'n_cell' to float if necessary
    gdf['n_cell'] = gdf['n_cell'].astype(float)

    # Display the image beneath the cells
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    # If max_cutoff is provided, set vmax to it, otherwise, use the maximum value in the data
    vmin = gdf['n_cell'].min()
    vmax = max_cutoff if max_cutoff is not None else gdf['n_cell'].max()

    # Plot the polygons with a colormap based on the 'n_cell' values
    gdf.plot(ax=ax, column='n_cell', cmap='viridis', legend=True,
             legend_kwds={'label': 'Number of Cells',
                          'orientation': "horizontal"},
             vmin=vmin, vmax=vmax)  # Apply max_cutoff (vmax) to the color scale

    # Set title and axis labels
    ax.set_title('Spot Plot Colored by Cell Density (n_cell)', fontsize=15)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)

    # Adjust layout
    plt.tight_layout()
    plt.show()


