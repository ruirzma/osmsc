""" Data fusion for OSMsc objects"""


import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox

from .feature import add_gdf_elevation_features
from .grid import Grid



def add_building_height(OSMscBuilding_gdf, externalBuilding_gdf, extColName = None ):
    """
    Add building height features to OSM or building gdf.
    Please keep both gdf in the same crs

    Parameters
    ----------
    OSMscBuilding_gdf : GeoDataFrame
        OSM or local building dataset
    externalBuilding_gdf : GeoDataFrame
        external building dataset, including height info
    extColName : string 
        Height attr column

    Returns
    -------
    GeoDataFrame
    """
    print("Please wait a few minutes")

    # Spatial Join
    temp_gdf = gpd.sjoin(OSMscBuilding_gdf, externalBuilding_gdf, how="left", op="intersects")
    # For later quering
    uni_geom = temp_gdf['geometry'].unique()
    
    temp_gdf["Building_height"] = temp_gdf[extColName]

    building_height = []
    for index in range(len(uni_geom)):
        heights_list = list(temp_gdf[temp_gdf["geometry"] == uni_geom[index]].Building_height)
        # average height of all intersected buildings
        building_height.append(np.mean(heights_list))
    
    # Add attrs
    OSMscBuilding_gdf["Building_height"] = building_height
    
    return OSMscBuilding_gdf


def add_building_tags(OSMscBuilding_gdf, tagColName ,externalBuilding_gdf, extColName):

    """
    Add building tags to OSM or building gdf.
    Please keep both gdf in the same crs

    Parameters
    ----------
    OSMscBuilding_gdf : GeoDataFrame
        OSM or local building dataset
    tagColName : string
        Column name to be added in the building dataframe
    externalBuilding_gdf : GeoDataFrame
        external building dataset, including height info
    extColName : string 
        Original tag attr column

    Returns
    -------
    GeoDataFrame
    """

    intersection_area_gdf = gpd.overlay(OSMscBuilding_gdf, externalBuilding_gdf, how ='intersection')
    intersection_area_gdf["intersection_area"] = intersection_area_gdf["geometry"].area
    intersected_osmscID = list(intersection_area_gdf["osmscID"].unique())

    for ID in intersected_osmscID:
        # for one OSMscID object
        area_list = list(intersection_area_gdf[intersection_area_gdf["osmscID"] == ID].intersection_area)
        index_list = list(intersection_area_gdf[intersection_area_gdf["osmscID"] == ID].index)

        # Sometimes, one OSM building intersects with several external buildings.
        if len(index_list) != 1: 
            # Only keep one row with biggest intersection area
            # index_list for later delete indexing
            index_of_max_area = area_list.index(max(area_list))
            del index_list[index_of_max_area]

            for i in index_list:
                intersection_area_gdf = intersection_area_gdf.drop(index=i)    
    
    OSMscBuilding_gdf[tagColName] = ["Unknown" for i in range(len(OSMscBuilding_gdf))]
    
    for osmscID in intersected_osmscID:
        value = intersection_area_gdf[intersection_area_gdf["osmscID"] == osmscID][extColName].values
        index = OSMscBuilding_gdf[OSMscBuilding_gdf["osmscID"] == osmscID].index

        OSMscBuilding_gdf.loc[index, tagColName] = value 
    
    for osmscID in OSMscBuilding_gdf["osmscID"]:
        if osmscID not in intersected_osmscID:
            index = OSMscBuilding_gdf[OSMscBuilding_gdf["osmscID"] == osmscID].index
            OSMscBuilding_gdf.loc[index, tagColName] = "Unknown"
    
    # replace None with "Unknown" for CityJSON output function
    for i in range(len(OSMscBuilding_gdf[tagColName].values)):
        col_list = OSMscBuilding_gdf[tagColName].values
        if col_list[i] == None:
            col_list[i] = "Unknown"
    
    OSMscBuilding_gdf[tagColName] = col_list
    
    return OSMscBuilding_gdf


def add_elevation(cityobject_gdf, elevation_dataset = "aster30m", step = 5):
    """
    Mostly, add elevation to UrbanTile and Transportation GeoDataFrame

    Parameters
    ----------
    cityobject_gdf : GeoDataFrame
        UrbanTile or Transportation GeoDataFrame
    elevation_dataset : string
        Online elevation can be chose from https://www.opentopodata.org
    step : int
        Grid size

    Returns
    -------
    GeoDataFrame
    """

    # Create the grid GeoDataFrame (with elevation)
    grid_gdf = get_intersected_grid_gdf(cityobject_gdf, elevation_dataset = elevation_dataset, step = 5)
    grid_gdf["y_lat"] = grid_gdf["y_lat"].replace(np.nan, "None")
    grid_gdf["x_lng"] = grid_gdf["x_lng"].replace(np.nan, "None")
    grid_gdf["ground_elevation"] = grid_gdf["ground_elevation"].replace(np.nan, "None")  

    elevation_list = []
    cityobject_gdf["elevation"] = ["Unknown" for i in range(len(cityobject_gdf))]
    
    for i in range(len(cityobject_gdf)):
        _osmscid = list(cityobject_gdf["osmscID"])[i]

        # Extract all location in the grid.
        lat_list = list(grid_gdf[grid_gdf.osmscID == _osmscid].y_lat)
        lon_list = list(grid_gdf[grid_gdf.osmscID == _osmscid].x_lng)
        ele_list = list(grid_gdf[grid_gdf.osmscID == _osmscid].ground_elevation)

        # Make the mass point list
        massPoints_list = []
        for j in range(len(lat_list)):
            massPoints_list.append({'lat': lat_list[j], 'lon': lon_list[j],'ele': ele_list[j]})
            
        elevation_list.append(massPoints_list)
        
    cityobject_gdf["elevation"] = elevation_list

    return cityobject_gdf


def get_intersected_grid_gdf(cityobject_gdf, elevation_dataset = "aster30m", step = 5):
    """
    Mostly, add elevation to a grid

    Parameters
    ----------
    cityobject_gdf : GeoDataFrame
        UrbanTile or Transportation GeoDataFrame
    elevation_dataset : string
        Online elevation can be chose from https://www.opentopodata.org
    step : int
        Grid size

    Returns
    -------
    GeoDataFrame
    """

    # Uniform the crs to epsg:4326 and then project it.
    if cityobject_gdf.crs != "epsg:4326":
        cityobject_gdf = cityobject_gdf.to_crs("epsg:4326")
    cityobject_gdf_prj = ox.project_gdf(cityobject_gdf)

    # Make the gird dataframe
    _Grid = Grid(data_type= "Polygon", step = 5, bandwidth = 10)
    Grid_gdf_prj = _Grid.computed_grid(cityobject_gdf_prj, X_weights = None)
    Grid_gdf_prj.crs = cityobject_gdf_prj.crs

    # Keep Point geometry
    Grid_gdf_prj["Point"] = Grid_gdf_prj["geometry"]

    # In the projection crs, only keep grid points inside cityobject_gdf_prj
    intersected_grid_gdf_prj = gpd.sjoin(cityobject_gdf_prj, Grid_gdf_prj,how="left", op="intersects")
    intersected_grid_gdf_prj["geometry"] = intersected_grid_gdf_prj["Point"]

    # For elevaion query, change the crs into epsg:4326
    intersected_grid_gdf = intersected_grid_gdf_prj.to_crs("epsg:4326")
    intersected_grid_gdf_add_ele = add_gdf_elevation_features(gdf= intersected_grid_gdf, data_type="Point", 
                        elevation_dataset = elevation_dataset)
    intersected_grid_gdf_add_ele.crs = "epsg:4326"

    return intersected_grid_gdf_add_ele





