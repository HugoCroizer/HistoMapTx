    def plot_cells(gdf, fill=True, contour='k'):
        """
        Plots segmented cells based on a GeoDataFrame.

        Parameters:
        - gdf: GeoDataFrame containing segmented cell geometries.
        - fill: False, True, or a list of colors.
            - False: No fill (only contours).
            - True: Uses a default colormap ('tab20') for fill.
            - List of colors: Specifies the fill color for each cell.
        - contour: Color or list of colors for the contours. If None, the default colormap is used.
        """
        fig, ax = plt.subplots(figsize=(10, 10))

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
