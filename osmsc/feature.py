"""Add Polygon-based features."""

import math
import requests
import networkx as nx
import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox

from .utils import create_regularity_gdf

def add_spatial_semantics_attr(left_gdf,right_gdf,semColName, how="left", op="intersects"):
    """
    Explore the spatial semantic relationship between two GeoDataFrame

    Parameters
    ----------
    left_gdf : GeoDataFrame
    right_gdf : GeoDataFrame
        # As for both gdf, more detail cna be found in 
        # https://automating-gis-processes.github.io/CSC18/lessons/L4/spatial-join.html

    semColName : string e.g. containsBuilding or within Tile
        New column name

    Returns
    -------
    GeoDataFrame
    """
    
    # Obtain the osmscID and geometry info
    temp_left_gdf = left_gdf[["osmscID","geometry"]]
    temp_right_gdf = right_gdf[["osmscID","geometry"]]
    # relationship stored in GeoDataFrame
    temp_relationship =  gpd.sjoin(temp_left_gdf, temp_right_gdf, how = how, op = op)
    
    relation_list = []
    for left_ID in temp_left_gdf["osmscID"]:

        # right_ID_list for each left objects
        try:
            right_ID_list = list(temp_relationship[temp_relationship["osmscID_left"] == left_ID].osmscID_right) 
        except:
            right_ID_list = None

        relation_list.append(right_ID_list)
    
    # Save the semantic info in the left_gdf
    left_gdf[semColName] = relation_list
    
    return left_gdf

#################### intra layer attrs ######################
################## For polygon cityobjects ##################

def add_minimum_rotated_rectangle_attr(polygon_gdf):
    """
    Obtain the minimum rotated rectangle and corresponding attributes for a GeoDataFrame.

    Parameters
    ----------
    polygon_gdf : GeoDataFrame
        Add minimum_rotated_rectangle_attr into this GeoDataFrame.
        Before the operation, the GeoDataFrame need to be projected.

    Returns
    -------
    GeoDataFrame
    """

    temp_mrr_gdf = gpd.GeoDataFrame()
    temp_mrr_gdf["geometry"] = polygon_gdf["geometry"]
    temp_mrr_gdf.crs = polygon_gdf.crs
    
    mrr_geom = [geom.minimum_rotated_rectangle for geom in list(temp_mrr_gdf["geometry"])]
    temp_mrr_gdf["geometry"] = mrr_geom
    temp_mrr_gdf["mrr_area"] = temp_mrr_gdf["geometry"].area
    
    polygon_gdf["mrr_geometry"] = temp_mrr_gdf["geometry"]
    polygon_gdf["mrr_area"] = temp_mrr_gdf["mrr_area"]
    
    return polygon_gdf


def add_minimum_circumscribed_circle_attr(polygon_gdf):
    """
    Obtain the minimum circumscribed circle and corresponding attributes for a GeoDataFrame.

    Parameters
    ----------
    polygon_gdf : GeoDataFrame
        Add minimum_circumscribed_circle_attr into this GeoDataFrame.
        Before the operation, the GeoDataFrame need to be projected.

    Returns
    -------
    GeoDataFrame
    """
    temp_mcc_gdf = create_regularity_gdf(polygon_gdf)
    # Add attr columns
    polygon_gdf["mcc_geometry"] = temp_mcc_gdf["geometry"]
    polygon_gdf["mcc_area"] = temp_mcc_gdf["geometry"].area
    
    return polygon_gdf


def add_shape_factor_attr(polygon_gdf):
    """
    Obtain the shape factor attribute for a GeoDataFrame.

    Parameters
    ----------
    polygon_gdf : GeoDataFrame
        Add shape factor attribute into this GeoDataFrame.
        Before the operation, the GeoDataFrame need to be projected.

    Returns
    -------
    GeoDataFrame
    """

    temp_gdf = polygon_gdf
    temp_gdf["Polygon_area"] = temp_gdf["geometry"].area

    # shape factor calculation is based on minimum circumscribed circle geometry
    if "mcc_area" not in temp_gdf.columns:
        temp_gdf = add_minimum_circumscribed_circle_attr(temp_gdf)

    # Add attr columns
    polygon_gdf["shape_factor"] = temp_gdf["Polygon_area"]/temp_gdf["mcc_area"]
    # delete unnecessary column
    polygon_gdf = polygon_gdf.drop("Polygon_area",axis =1)

    return polygon_gdf


def add_polygon_bearing_attr(polygon_gdf_prj, use_mrr = True):
    """
    Obtain the polygon bearing attribute for a GeoDataFrame.

    Parameters
    ----------
    polygon_gdf_prj : GeoDataFrame
        Add polygon bearing attribute into this GeoDataFrame.
    use_mrr : bool
        True, use the minimum rotated rectangle to represent the current polygon
        False, use the current polygon to calculate bearing attributes.

    Returns
    -------
    GeoDataFrame
    """
    # re-project the gdf_prj
    if polygon_gdf_prj.crs != "EPSG:4326":
        polygon_gdf = polygon_gdf_prj.to_crs("EPSG:4326")

    # Whether using minimum rotated rectangle
    if use_mrr:
        mbr_gdf = create_minimum_bounding_rectangle_gdf(polygon_gdf_prj)
    else:
        mbr_gdf = polygon_gdf

    polygon_brngs = []
    for i in range(len(mbr_gdf)):
        # Traverse each geometry
        polygon = mbr_gdf.geometry[i]
        # Obtain the direction of each side (azimuth angle +90)
        brngs_list = get_brngList_for_polygon(polygon)
        polygon_brngs.append(brngs_list)

    polygon_gdf_prj["bearings"] = polygon_brngs

    return polygon_gdf_prj


def create_minimum_bounding_rectangle_gdf(polygon_gdf):
    """
    Obtain the minimum bounding rectangle for a GeoDataFrame.
    More detail can be found in
    https://gis.stackexchange.com/questions/22895/finding-minimum-area-rectangle-for-given-points
    
    Parameters
    ----------
    polygon_gdf : GeoDataFrame
        Add minimum bounding rectangle into this GeoDataFrame.
        
    Returns
    -------
    GeoDataFrame
    """

    mbr_gdf = gpd.GeoDataFrame()
    mbr_gdf["geometry"] = [geo.minimum_rotated_rectangle for geo in polygon_gdf["geometry"]]
    mbr_gdf.crs = polygon_gdf.crs

    return mbr_gdf


def get_brng_btw_points(pointA, pointB):
    """
    Calculate the azimuth angle between two points
    More detail can be found in
    https://blog.csdn.net/qiannianlaoyao2010/article/details/102807883?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522162157829716780366528793%2522%252C%2522scm%2522%253A%252220140713.1
    
    Parameters
    ----------
    pointA : shapely.geometry.Point
    pointA : shapely.geometry.Point

    Returns
    -------
    Float
    """

    # Point (lon, lat) of GeoDataFrame
    lonA, latA = pointA
    lonB, latB = pointB
    
    radLatA = math.radians(latA)
    radLonA = math.radians(lonA)
    radLatB = math.radians(latB)
    radLonB = math.radians(lonB)
    dLon = radLonB - radLonA
    y = math.sin(dLon) * math.cos(radLatB)
    x = math.cos(radLatA) * math.sin(radLatB) - math.sin(radLatA) * math.cos(radLatB) * math.cos(dLon)
    brng = math.degrees(math.atan2(y, x))
    brng = round((brng + 360) % 360, 4)
    brng = int(brng)
    
    return brng


def get_brngList_for_polygon(polygon):
    """
    Calculate the azimuth angle for a certain polygon

    Parameters
    ----------
    polygon : shapely.geometry.Polygon

    Returns
    -------
    List
    """

    brngs_list = []

    coords = list(polygon.boundary.coords)
    for i in range(len(coords)):
        if i != len(coords)-1:
            pointA = coords[i]
            pointB = coords[i+1]
            
            # Because the polygon nodes are connected counterclockwise, 
            # the normal of the breath degree is rotated by 90
            brng = (get_brng_btw_points(pointA, pointB) + 90 ) % 360
            
            brngs_list.append(brng)
            
    return brngs_list



################## inter layer attrs ##################
def add_interlayer_building_attr(UrbanTile_gdf, Building_gdf):
    """
    Add interlayer building attrs to gdf, mostly to UrbanTile_gdf.
    The prerequisite is the spatial semantic relationship between two layer is known, namely,
    add_spatial_semantics_attr() has been run.

    Parameters
    ----------
    UrbanTile_gdf : GeoDataFrame
        Add interlayer building attrs into this GeoDataFrame.
    Building_gdf : GeoDataFrame

        
    Returns
    -------
    GeoDataFrame
    """

    buildingDensity_list = []
    avg_buildingHeight_list = []
    avg_buildingArea_list = []
    avg_buildingPerimeter_list = []

    for row_num in range(len(UrbanTile_gdf)):

        # For one UrbanTile object
        ########## Avergae height/area/perimeter of buildings ###########
        bldg_height_list = []
        bldg_area_list = []
        bldg_perimeter_list = []

        for bldg_ID in UrbanTile_gdf["containsBuilding"][row_num]:

            try:
                # Building_height 
                a_height = float(Building_gdf[Building_gdf["osmscID"] == bldg_ID].Building_height)
                bldg_height_list.append(a_height)
            except:
                # some buildings do not have height 
                bldg_height_list.append(0)

            try:    
                a_area = float(Building_gdf[Building_gdf["osmscID"] == bldg_ID].Building_area)
                bldg_area_list.append(a_area)

                a_perimeter = float(Building_gdf[Building_gdf["osmscID"] == bldg_ID].Building_perimeter)
                bldg_perimeter_list.append(a_perimeter)

            except:
                # Some buildings need to add more area and perimeter attrs.
                bldg_area_list.append(0)
                bldg_perimeter_list.append(0)

        avg_buildingHeight_list.append(np.mean(bldg_height_list)) 
        avg_buildingArea_list.append(np.mean(bldg_area_list)) 
        avg_buildingPerimeter_list.append(np.mean(bldg_perimeter_list))

        UrbanTile_area = float(UrbanTile_gdf["UrbanTile_area"].loc[row_num])
        buildingDensity_list.append(np.sum(bldg_area_list)/UrbanTile_area)

    # Add attrs columns
    UrbanTile_gdf["buildingDensity"] = buildingDensity_list
    UrbanTile_gdf["avg_buildingHeight"] = avg_buildingHeight_list
    UrbanTile_gdf["avg_buildingArea"] = avg_buildingArea_list
    UrbanTile_gdf["avg_buildingPerimeter"] = avg_buildingPerimeter_list
    
    return UrbanTile_gdf


def add_interlayer_vegetation_attr(UrbanTile_gdf, Vegetation_gdf):
    """
    Add interlayer vegetation attrs to gdf, mostly to UrbanTile_gdf.
    The prerequisite is the spatial semantic relationship between two layer is known, namely,
    add_spatial_semantics_attr() has been run.

    Parameters
    ----------
    UrbanTile_gdf : GeoDataFrame
        Add interlayer vegetation attrs into this GeoDataFrame.
    Building_gdf : GeoDataFrame

        
    Returns
    -------
    GeoDataFrame
    """
    vegetationDensity_list = []
    avg_vegetationArea_list = []
    avg_vegetationPerimeter_list = []

    for row_num in range(len(UrbanTile_gdf)):

        # For one UrbanTile object
        ########## Avergae area/perimeter of vegetation objects ###########

        veg_area_list = []
        veg_perimeter_list = []

        for veg_ID in UrbanTile_gdf["containsVegetation"][row_num]:

            try:

                a_area = float(Vegetation_gdf[Vegetation_gdf["osmscID"] == veg_ID].Vegetation_area)
                veg_area_list.append(a_area)

                a_perimeter = float(Vegetation_gdf[Vegetation_gdf["osmscID"] == veg_ID].Vegetation_perimeter)
                veg_perimeter_list.append(a_perimeter)

            except:
                veg_area_list.append(0)
                veg_perimeter_list.append(0)

        avg_vegetationArea_list.append(np.mean(veg_area_list)) 
        avg_vegetationPerimeter_list.append(np.mean(veg_perimeter_list))

        UrbanTile_area = float(UrbanTile_gdf["UrbanTile_area"].loc[row_num])
        vegetationDensity_list.append(np.sum(veg_area_list)/UrbanTile_area)

    # Add attrs columns
    UrbanTile_gdf["vegetationDensity"] = vegetationDensity_list
    UrbanTile_gdf["avg_vegetationArea"] = avg_vegetationArea_list
    UrbanTile_gdf["avg_vegetationPerimeter"] = avg_vegetationPerimeter_list
    
    return UrbanTile_gdf


def add_interlayer_waterbody_attr(UrbanTile_gdf, Waterbody_gdf):
    """
    Add interlayer waterbody attrs to gdf, mostly to UrbanTile_gdf.
    The prerequisite is the spatial semantic relationship between two layer is known, namely,
    add_spatial_semantics_attr() has been run.

    Parameters
    ----------
    UrbanTile_gdf : GeoDataFrame
        Add interlayer waterbody attrs into this GeoDataFrame.
    Building_gdf : GeoDataFrame

        
    Returns
    -------
    GeoDataFrame
    """
    waterbodyDensity_list = []
    avg_waterbodyArea_list = []
    avg_waterbodyPerimeter_list = []

    for row_num in range(len(UrbanTile_gdf)):

        # For one UrbanTile object
        ########## Avergae area/perimeter of waterbody objects ###########

        wat_area_list = []
        wat_perimeter_list = []

        for wat_ID in UrbanTile_gdf["containsWaterbody"][row_num]:

            try:

                a_area = float(Waterbody_gdf[Waterbody_gdf["osmscID"] == wat_ID].Waterbody_area)
                wat_area_list.append(a_area)

                a_perimeter = float(Waterbody_gdf[Waterbody_gdf["osmscID"] == wat_ID].Waterbody_perimeter)
                wat_perimeter_list.append(a_perimeter)

            except:
                wat_area_list.append(0)
                wat_perimeter_list.append(0)

        avg_waterbodyArea_list.append(np.mean(wat_area_list)) 
        avg_waterbodyPerimeter_list.append(np.mean(wat_perimeter_list))

        UrbanTile_area = float(UrbanTile_gdf["UrbanTile_area"].loc[row_num])
        waterbodyDensity_list.append(np.sum(wat_area_list)/UrbanTile_area)
    
    # Add attrs columns
    UrbanTile_gdf["waterbodyDensity"] = waterbodyDensity_list
    UrbanTile_gdf["avg_waterbodyArea"] = avg_waterbodyArea_list
    UrbanTile_gdf["avg_waterbodyPerimeter"] = avg_waterbodyPerimeter_list
    
    return UrbanTile_gdf

################## Elevation ##################

def add_graph_elevation_features(G, elevation_dataset= None, max_locations_per_batch=100):
    """
    Add graph elevation features to Graph object.

    Editted from osmnx https://github.com/gboeing/osmnx 
    Open Topo Data Dataset: https://www.opentopodata.org/#public-api

    Parameters
    ----------
    G : Graph
    elevation_dataset : string
        Online elevation can be chose from https://www.opentopodata.org
    max_locations_per_batch : int

    Returns
    -------
    Graph
    """

    # editted from osmnx
    # Open Topo Data Dataset: https://www.opentopodata.org/#public-api

    if elevation_dataset is None:
        url_template = 'https://api.opentopodata.org/v1/aster30m?locations={}'
    else:
        url_template = 'https://api.opentopodata.org/v1/{}?locations={}'.format(elevation_dataset,{})


    node_points = pd.Series({node:'{:.5f},{:.5f}'.format(data['y'], data['x']) for node, data in G.nodes(data=True)})

    # API format is locations=lat,lng|lat,lng|lat,lng|lat,lng...
    
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]
        locations = '|'.join(chunk)
        url = url_template.format(locations)

        try:
            # request the elevations from the API
            response = requests.get(url)
            response_json = response.json()
            
        except:
            print('Server responded with {}: {}'.format(response.status_code, response.reason))

        # append these elevation results to the list of all results
        results.extend(response_json['results'])

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G.nodes()) == len(node_points)):
        raise Exception('Graph has {} nodes but we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))
    else:
        print('Graph has {} nodes and we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))

    # add elevation as an attribute to the nodes
    df = pd.DataFrame(node_points, columns=['node_points'])
    df['elevation'] = [result['elevation'] for result in results]
    df['elevation'] = df['elevation'].round(3) # round to millimeter
    nx.set_node_attributes(G, name='elevation', values=df['elevation'].to_dict())
    
    # print('Added elevation data to all nodes.')

    return G


def add_gdf_elevation_features(gdf, elevation_dataset = None, data_type = None, max_locations_per_batch=100):
    """
    Add gdf elevation features to GeoDataFrame.
    Open Topo Data Dataset: https://www.opentopodata.org/#public-api

    Parameters
    ----------
    gdf : GeoDataFrame
    elevation_dataset : string
        Online elevation can be chose from https://www.opentopodata.org
    data_type : string
        "Polygon" or "Point" from shape.geometry
    max_locations_per_batch : int

    Returns
    -------
    GeoDataFrame
    """

    if elevation_dataset is None:
        url_template = 'https://api.opentopodata.org/v1/aster30m?locations={}'
    else:
        url_template = 'https://api.opentopodata.org/v1/{}?locations={}'.format(elevation_dataset,{})

    # For polygon object
    # Make a pandas series of all the nodes' coordinates as 'lat,lng'
    if data_type == "Polygon":

        gdf_prj = ox.project_gdf(gdf)
        centroid_prj = gdf_prj["geometry"].centroid
        gdf["centroid"] = centroid_prj.to_crs("EPSG:4326")
        gdf["centroid_y_lat"] = gdf["centroid"].y
        gdf["centroid_x_lng"] = gdf["centroid"].x

        points = pd.Series('{:.5f},{:.5f}'.format(gdf.iloc[i].centroid_y_lat, gdf.iloc[i].centroid_x_lng) for i in range(len(gdf)))

    elif data_type == "Point":
        gdf["y_lat"] = gdf.geometry.y
        gdf["x_lng"] = gdf.geometry.x
        points = pd.Series('{:.5f},{:.5f}'.format(gdf.iloc[i].y_lat, gdf.iloc[i].x_lng) for i in range(len(gdf)))
    
    elif data_type is None:
        print("Please input the data_type of geometry column")

    # Store the elevation results.
    results = []
    for i in range(0, len(points), max_locations_per_batch):
        chunk = points.iloc[i : i + max_locations_per_batch]
        locations = '|'.join(chunk)
        url = url_template.format(locations)
        response = requests.get(url)
        response_json = response.json()

        results.extend(response_json['results'])
    
    elevation = [results[i]['elevation'] for i in range(len(results))]
    
    gdf["ground_elevation"] = elevation
    
    return gdf








