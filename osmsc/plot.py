"""OSMsc plot module"""

import geopandas as gpd
import pydeck as pdk
from contextily import add_basemap
import matplotlib.pyplot as plt

from .utils import bbox_from_gdf


def viz_buildings(buildings_gdf, zoom=11, max_zoom=16, pitch=45, bearing=0, 
        html_name = "3D_buildings"):
    """
    Building DataFrame 3D Visualization.
    NOTE buildings_gdf should have "Building_height" column.

    Parameters
    ----------
    buildings_gdf : GeoDataFrame
    zoom : int
    max_zoom : int
    pitch : int
    bearing : int
    html_name : string
        filename to be saved in the current project location.

    Returns
    -------
    HTML
    """
    bbox = bbox_from_gdf(buildings_gdf,buffer_value=0.002)
    (minlat, minlon, maxlat, maxlon) = bbox
    land_cover = [[minlat,minlon], [minlat,maxlon], [maxlat,minlon], [maxlat,maxlon]]
    
    max_Height = max(buildings_gdf["Building_height"])

    INITIAL_VIEW_STATE = pdk.ViewState(latitude= (minlat + maxlat)/2, longitude=(minlon + maxlon)/2, 
                                        zoom = zoom, max_zoom=max_zoom, pitch=pitch, bearing=bearing)

    polygon = pdk.Layer(
        "PolygonLayer",
        land_cover,
        stroked=False,
        # processes the data as a flat longitude-latitude pair
        get_polygon="-",
        get_fill_color=[0, 0, 0, 20],
    )

    geojson = pdk.Layer(
        "GeoJsonLayer",
        buildings_gdf,
        opacity=0.8,
        stroked=False,
        filled=True,
        extruded=True,
        wireframe=True,
        get_elevation="Building_height",

        get_fill_color="[ Building_height/" + str(max_Height) + "*255, 250, 85]",
        get_line_color=[255, 255, 255]
    )

    r = pdk.Deck(layers=[polygon, geojson], initial_view_state=INITIAL_VIEW_STATE)

    return r.to_html(html_name + ".html")  


def gdf_with_basemap(gdf, figsize = (12, 6), column = None, color = None,markersize = None, 
                    cmap = None ,legend=False, categorical = False, edgecolor = None, linewidth = None, 
                    facecolor= None ,alpha = None, zoom = "auto", tile_source = None):
    """
    Plot gdf with basemap

    column : str, np.array, pd.Series (default None)
    The name of the dataframe column, np.array, or pd.Series to be plotted.
    If np.array or pd.Series are used then it must have same length as
    dataframe. Values are used to color the plot. Ignored if `color` is
    also set.

    cmap : str (default None)
    The name of a colormap recognized by matplotlib.

    color : str (default None)
    If specified, all objects will be colored uniformly.

    categorical : bool (default False)
        If False, cmap will reflect numerical values of the
        column being plotted.  For non-numerical columns, this
        will be set to True.

    legend : bool (default False)
        Plot a legend. Ignored if no `column` is given, or if `color` is given.

    markersize : str or float or sequence (default None)
        Only applies to point geometries within a frame.
        If a str, will use the values in the column of the frame specified
        by markersize to set the size of markers. Otherwise can be a value
        to apply to all points, or a sequence of the same length as the
        number of points.

    figsize : tuple of integers (default None)
        Size of the resulting matplotlib.figure.Figure. If the argument
        axes is given explicitly, figsize is ignored.


    edgecolor : str (default None)

    facecolor ： str (default None)

    linewidth ： str (default None)

    alpha ： float (default None)

    zoom : int or 'auto'
        [Optional. Default='auto'] Level of detail for the basemap. If 'auto',
        it is calculated automatically. Ignored if `source` is a local file.
        
    tile_source : contextily.providers object or str
        [Optional. Default: Stamen Terrain web tiles]
        The tile source: web tile provider or path to local file. The web tile
        provider can be in the form of a `contextily.providers` object or a
        URL. The placeholders for the XYZ in the URL need to be `{x}`, `{y}`,
        `{z}`, respectively. For local file paths, the file is read with
        `rasterio` and all bands are loaded into the basemap.
        IMPORTANT: tiles are assumed to be in the Spherical Mercator
        projection (EPSG:3857), unless the `crs` keyword is specified.

    Returns
    -------
    fig, ax
    """
    fig, ax = plt.subplots(figsize = figsize)
    gdf.plot( ax=ax, column = column, color = color,markersize = markersize, 
                    cmap = cmap ,legend =legend, categorical = categorical, edgecolor = edgecolor, 
                    linewidth = linewidth, facecolor= facecolor ,alpha = alpha)

    add_basemap(ax, crs = gdf.crs, zoom = "auto", source= tile_source)

    return fig, ax




