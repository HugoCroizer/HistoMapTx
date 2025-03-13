# Utils 
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


# Define a function to calculate overlap percentage
def calculate_overlap(spot, annotation):
    # Calculate the intersection area
    intersection = spot.intersection(annotation)
    if intersection.is_empty:
        return 0
    else:
        # Calculate the area of overlap as a percentage of the spot's area
        return intersection.area / spot.area * 100

# Function to calculate overlap for each spot against annotations
def calculate_annotation_overlap(gdf, data_exploded, spot_idx):
    """ gdf = spot_geodata, data_exploded is a subset of histomap.data_exploded with only required annotation for computation"""
    n = 0
    # Loop over unique annotations in the 'Annotation' column of data_exploded
    for annotation_name in data_exploded['Annotation'].unique():
        n += 1
        print(f"{n}/{len(data_exploded['Annotation'].unique())} layer of annotation processed: {annotation_name}")
        
        # Filter the annotations that correspond to the current group
        annotation_group = data_exploded[data_exploded['Annotation'] == annotation_name]

        # Dynamically create a column name for this annotation
        column_name = f"{annotation_name}_overlap"
        
        # Initialize the new column 
        if column_name not in gdf.columns:
            gdf[column_name] = 0

        # Create a MultiPolygon for the current annotation group
        annotations = [row['geometry'] for _, row in annotation_group.iterrows()]
        from shapely.geometry import MultiPolygon
        multi_polygon = MultiPolygon(annotations)

        # Find all spots that intersect with this MultiPolygon using the spatial index
        possible_spots = list(spot_idx.intersection(multi_polygon.bounds))
        
        # Loop over possible spots and calculate overlap
        for spot_id in possible_spots:
            spot = gdf.loc[spot_id, 'geometry']
            overlap_percentage = calculate_overlap(spot, multi_polygon)

            # Update the overlap percentage for the spot in the corresponding annotation column
            gdf.at[spot_id, column_name] = overlap_percentage
    
    # Fill NaN values with 0 after processing all annotations
    gdf.fillna(0, inplace=True)

    return gdf



def filter_spatialdata(
    histomap, cell_index
    ):
    """
    Filter a SpatialData object based on cell indexes in place

    Returns:
        A SpatialData 
    """
    # Get and subset adata 
    adata = histomap.visium_spatialdata.tables['table']
    adata_subset = adata[cell_index]
    histomap.visium_spatialdata.tables['table'] = adata_subset
    # get and subset spot shapes
    spot_df_subset = histomap.visium_spatialdata.shapes[list(histomap.visium_spatialdata.shapes.keys())[0]].iloc[cell_index]
    histomap.visium_spatialdata.shapes[list(histomap.visium_spatialdata.shapes.keys())[0]] = spot_df_subset


def load_histomap(filename):
    """Load a HistoMap object from a pickle file."""
    try:
        import pickle
        with open(filename, 'rb') as f:
            loaded_obj = pickle.load(f)
        print(f"Object loaded from {filename}")
        return loaded_obj
    except Exception as e:
        print(f"Error loading object: {e}")
        return None