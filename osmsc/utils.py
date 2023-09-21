"""General utility functions."""

import osmnx as ox
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import math
import requests
import random

from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.ops import polygonize


def download(overpass_query):
    """
    Download OSM via Overpass API

    Parameters
    ----------
    overpass_query : string
        Customized overpass API

    Returns
    -------
    JSON
    """
    overpass_url = "http://overpass-api.de/api/interpreter"
    response = requests.get(overpass_url, params={'data': overpass_query})
    osm_json = response.json()

    return osm_json 


def json_to_gdf(osm_json, data_type, tags = True, building_levels = False, height = False):
    """
    Convert the json data into GeoDataFrame.

    Parameters
    ----------
    osm_json : dict
        Query result from OSM
    data_type : string
        Adjust the GeoPandas geometry object with OSMsc object datatype
    tags : bool
        The final GeoDataFrame has "tags" column, that are downloaded from OSM
    building_levels : bool
        The final GeoDataFrame has "building_levels" column, that are downloaded from OSM
    height : bool
        The final GeoDataFrame has "height" column, that are downloaded from OSM

    Returns
    -------
    GeoDataFrame
    """

    # when osm_json is not empty
    if osm_json['elements']:       

        if data_type == "Polygon":
            # Collect the Polygons in osm_json
            coords_list = []
            osm_id = []
            tags_list = []
            building_levels_list = []
            building_height_list = []

            for element in osm_json['elements']:

                # Polygon id (from osm "way" id)
                # osm_id.append(element["id"])

                # Polygon nodes
                coords_one = []
                for node in element["geometry"]:
                    lon = node['lon']
                    lat = node['lat']

                    # in gpd geometry column, the order is lon, lat,
                    coords_one.append((lon, lat)) 
                # A Polygon has 3 points at least. 
                if len(coords_one) > 2:
                    coords_list.append(Polygon(coords_one))

                    # ensure that the two columns have the same number of columns
                    osm_id.append(element["id"])

                # Polygon tags
                if tags:

                    # element["tags"] --> dict
                    tags_list.append(element["tags"])
                
                # Polygon building:levels tag
                try:# unknown tags 
                    if building_levels:
                        if "building:levels" in element["tags"].keys():
                            building_levels_list.append(int(float(element["tags"]["building:levels"])))
                        else:
                            building_levels_list.append(1)
                except:
                        building_levels_list.append(1)
                try:# unknown tags
                    if height:
                        if "height" in element["tags"].keys():
                            # aviod "7 m"  "7.5 m"
                            height_string = element["tags"]["height"]
                            building_height_list.append([float(i) for i in height_string.split() if i.isdigit() or i.replace(".", '', 1).isdigit()][0])
                        else:
                            building_height_list.append(3)                            
                except:
                    building_height_list.append(3)
            
                # # if building level>1, height  = building level * 3
                if building_levels:
                    if building_levels_list[-1] !=1 and building_height_list[-1] == 3:
                        building_height_list[-1] = building_levels_list[-1] *3
                    if building_height_list[-1] == 0:
                        building_height_list[-1] = building_levels_list[-1] *3



        if data_type == "Point":
            coords_list = []
            osm_id = []
            tags_list = []

            for element in osm_json['elements']:

                # Polygon id (from osm "way" id)
                osm_id.append(element["id"])

                # Polygon nodes
                coords_one = []

                lon = element['lon']
                lat = element['lat']

                # in gpd, the order is lon, lat,

                coords_one.append((lon, lat)) 

                coords_list.append(Point(coords_one))

                # Polygon tags
                if tags:
                    # element["tags"] --> dict
                    tags_list.append(element["tags"])

        if data_type == "LineString":

            coords_list = []
            osm_id = []
            tags_list = []

            for element in osm_json['elements']:

                # Polygon id (from osm "way" id)
                osm_id.append(element["id"])

                # Polygon nodes
                coords_one = []
                for node in element["geometry"]:
                    lon = node['lon']
                    lat = node['lat']

                    # in gpd, the order is lon, lat,
                    coords_one.append([lon, lat]) 

                coords_list.append(LineString(coords_one))

                if tags:
                    # element["tags"] --> dict
                    tags_list.append(element["tags"])
                    
        # create a new GeoDataFrame to store the Polygons
        temp_gdf = gpd.GeoDataFrame()
        temp_gdf["osmid"] = osm_id
        temp_gdf['geometry'] = coords_list
        # set crs
        temp_gdf = temp_gdf.set_crs(crs='epsg:4326')
        #gpd.GeoDataFrame(temp_gdf,crs='epsg:4326')

        # Add "tags" column
        if tags:
            temp_gdf["tags"] = tags_list

        # Add "building_levels" column
        if building_levels:
            temp_gdf["building_levels"] = building_levels_list           
        if height:
            temp_gdf["Building_height"] = building_height_list

        return temp_gdf

    else:
        print("Your query returns nothing!")


def fill_nan_list_position(gdf,colName):
    """
    Fill the nan value of GeoDataFrame with "None"
    Designed for CityJSON object construction

    Parameters
    ----------
    gdf : GeoDataFrame
    colName : string
        Which column needs to be applied with this funciton.

    Returns
    -------
    GeoDataFrame
    """
    item_list = []
    for item in gdf[colName]:

        if item == [np.nan]:
            item_list.append("None")
        else:
            item_list.append(item)
    
    gdf[colName] = item_list
    
    return gdf


def tailor_gdf_with_bbox(target_gdf_prj, bbox):
    """
    Tailor the target GeoDataFrame, usually for vegetation or waterbody objects

    Parameters
    ----------
    target_gdf : GeoDataFrame
        vegetation or waterbody GeoDataFrames
    bbox : tuple
        bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)


    Returns
    -------
    GeoDataFrame
    """
    # create bbox GeoDataFrame
    bbox_gdf = get_bbox_gdf(bbox)
    # convert the crs
    bbox_gdf_prj = ox.project_gdf(bbox_gdf)

    # tailor the target_gdf
    tailored_gdf = gpd.overlay(target_gdf_prj, bbox_gdf_prj, how="intersection")
    
    return tailored_gdf


################# for bbox #############################
def bbox_from_gdf(gdf, buffer_value = 0.001):
    """
    Obtain bbox from the whole GeoDataFrame
    bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)

    Parameters
    ----------
    gdf : GeoDataFrame
    buffer_value : float
        buffer variable, see more in GeoDataFrame.Geometry.buffer()

    Returns
    -------
    Tuple
    """

    # the whole geometry
    gdf_boundary =gdf.cascaded_union.convex_hull.buffer(buffer_value)
    # get the bounds
    (West, South, East, North) =  gdf_boundary.bounds    
    # create the bbox
    bbox = (South, West, North, East)

    return bbox

def bbox_from_place_name(place_name):
    """
    Obtain bbox from place_name
    bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)
    Place name can be checked in https://nominatim.openstreetmap.org/ui/search.html
    
    Parameters
    ----------
    place_name : str

    Returns
    -------
    Tuple
    """
    # get place ploygon
    poly_gdf = ox.geocoder.geocode_to_gdf(place_name)
    # get the bounds of this ploygon
    bounds = poly_gdf.geometry.iloc[0].bounds
    min_lon, min_lat, max_lon, max_lat  = list(bounds)
    # create the bbox
    bbox = (min_lat, min_lon, max_lat, max_lon)

    return bbox

def get_bbox_gdf(bbox):
    """
    Create bbox GeoDataFrame

    Parameters
    ----------
    bbox : tuple
        bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)


    Returns
    -------
    GeoDataFrame
    """

    # create bbox geometry 
    min_lat, min_lon, max_lat, max_lon = bbox
    bbox_geometry = Polygon([(min_lon,min_lat),(min_lon,max_lat ),(max_lon, max_lat), (max_lon,min_lat)])

    # create bbox GeoDataFrame
    bbox_gdf = gpd.GeoDataFrame()
    bbox_gdf["geometry"] = [bbox_geometry]
    bbox_gdf.crs = "epsg:4326"

    return bbox_gdf

################# for osm_json #############################
def get_keys_from_osm_json(osm_json):
    """
    Obtain all keys from queried OSM results, and do not include repeated elements.

    Parameters
    ----------
    osm_json : JSON

    Returns
    -------
    List
    """

    osm_keys = []
    for element in osm_json['elements']:
        keys = list(element["tags"].keys())
    
        osm_keys = list(set(osm_keys + keys))
    
    return osm_keys

def get_values_from_osm_json(osm_json, osm_key):
    """
    Obtain all values for a specific key from queried OSM results

    Parameters
    ----------
    osm_json : JSON
    osm_key : string
        A specific key, for instance, "building"

    Returns
    -------
    List
    """

    osm_values = []
    for element in osm_json['elements']:
        
        # Not every element has this key
        try:
            values = element["tags"][str(osm_key)]
            if values not in osm_values:
                osm_values.append(values)
        except:
            pass
        
    return osm_values


################# for overpass_api #############################
def customize_query_content(osm_dict_k_v, bbox, union_set, intersection_set, Polygon, Graph, Point):
    """
    Customize a query statement

    Parameters
    ----------
    osm_dict_k_v : dict
        For example, osm_dict_k_v = {osm_key_1:osm_values_1, osm_key_2: osm_values_2}
        osm_key -- str
        osm_values -- list [str,str,str, ...]
    bbox : tuple
        bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)
    union_set : bool
        True, return the union set of all query results.
    intersection_set : bool
        True, return the intersection set of all query results
    Polygon : bool
        True, Polygon geometry saved in GeoPandas geometry colunmn
    Graph : bool
        True, Graph geometry saved in GeoPandas geometry colunmn
    Point : bool
        True, Point geometry saved in GeoPandas geometry colunmn

    Returns
    -------
    String
    """
    # overpass content is one row
    if intersection_set:
        conditions = ""
        for osm_key, osm_values in osm_dict_k_v.items():
        
            # Meet the conditions of different conditions at the same time. 
            # And the'|' of a condition is an "or" relationship
            conditions = conditions +  "[" + '"' + str(osm_key) + '"' + "~" +'"' +  '|'.join(osm_values) +'"' + "]"

        if Polygon or Graph:
            query_content = """way""" + conditions + str(bbox) +  """;"""
        elif Point:
            query_content = """node""" + conditions + str(bbox) +  """[~"."~"."];"""
    
    # overpass content is multiple rows
    if union_set:
        
        if Polygon or Graph:
            query_content = ""
            for osm_key, osm_values in osm_dict_k_v.items():
                query_content = query_content + "way[" + '"' + str(osm_key) + '"' + "~" +'"' +  '|'.join(osm_values) +'"' + "]"+ str(bbox) + ";"
        
        elif Point:
            query_content = ""
            for osm_key, osm_values in osm_dict_k_v.items():
                query_content = query_content + "node[" + '"' + str(osm_key) + '"' + "~" +'"' +  '|'.join(osm_values) +'"' + "]"+ str(bbox) + """[~"."~"."];"""

    return query_content

# TODO more overpass api settings
def create_overpass_query(osm_dict_k_v, bbox, union_set = True, intersection_set = False, 
                            Polygon = True, Graph = False, Point = False):

    """
    Create a overpass query according to k_v pairs

    Parameters
    ----------
    osm_dict_k_v : dict
        For example, osm_dict_k_v = {osm_key_1:osm_values_1, osm_key_2: osm_values_2}
        osm_key -- str
        osm_values -- list [str,str,str, ...]
    bbox : tuple
        bbox --> (South, West, North, East) or (# minlat, minlon, maxlat, maxlon)
    union_set : bool
        True, return the union set of all query results.
    intersection_set : bool
        True, return the intersection set of all query results
    Polygon : bool
        True, Polygon geometry saved in GeoPandas geometry colunmn
    Graph : bool
        True, Graph geometry saved in GeoPandas geometry colunmn
    Point : bool
        True, Point geometry saved in GeoPandas geometry colunmn

    Returns
    -------
    String
    """

    # Build main query content
    query_content = customize_query_content(osm_dict_k_v, bbox, union_set , intersection_set  ,Polygon , Graph , Point )

    if query_content:   
        # Complete overpass query         
        overpass_query = """
                    [out:json][timeout:50];
                    ( """ + query_content +  """ );
                    out geom;
                    """ 
    else:
        print("please enter your target data_type!") 
        return None

    return overpass_query 


###################### osmnx utils ################################
# TODO rewrite osmnx package parts for osmsc
def graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs=None):
    """
    Updated from osmnx https://github.com/gboeing/osmnx
    Convert node and edge GeoDataFrames to a MultiDiGraph.

    This function is the inverse of `graph_to_gdfs` and is designed to work in
    conjunction with it. However, you can convert arbitrary node and edge
    GeoDataFrames as long as gdf_nodes is uniquely indexed by `osmid` and
    gdf_edges is uniquely multi-indexed by `u`, `v`, `key` (following normal
    MultiDiGraph structure). This allows you to load any node/edge shapefiles
    or GeoPackage layers as GeoDataFrames then convert them to a MultiDiGraph
    for graph analysis.

    Parameters
    ----------
    gdf_nodes : geopandas.GeoDataFrame
        GeoDataFrame of graph nodes uniquely indexed by osmid
    gdf_edges : geopandas.GeoDataFrame
        GeoDataFrame of graph edges uniquely multi-indexed by u, v, key
    graph_attrs : dict
        the new G.graph attribute dict. if None, use crs from gdf_edges as the
        only graph-level attribute (gdf_edges must have crs attribute set)

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    if graph_attrs is None:
        graph_attrs = {"crs": gdf_edges.crs}
    G = nx.MultiDiGraph(**graph_attrs)

    # add edges and their attributes to graph, but filter out null attribute
    # values so that edges only get attributes with non-null values

    # attr_names = gdf_edges.columns.to_list()
    # for (u, v, k), attr_vals in zip(gdf_edges.index, gdf_edges.values):

    for (u, v, k) in gdf_edges.index:
        # data_all = zip(attr_names, attr_vals)
        # data = {name: val for name, val in data_all if isinstance(val, list) or pd.notnull(val)}
        # G.add_edge(u, v, key=k, **data)
        G.add_edge(u, v, key=k)

    # add nodes' attributes to graph
    for col in gdf_nodes.columns:
        nx.set_node_attributes(G, name=col, values=gdf_nodes[col].dropna())

    return G



##################### osmuf utils ################################
# TODO rewrite osmuf package parts for osmsc
def street_graph_from_gdf(gdf, network_type='drive'):
    """
    Updated from osmuf https://github.com/AtelierLibre/osmuf
    Download streets within a convex hull around a GeoDataFrame.

    Used to ensure that all streets around city blocks are downloaded, not just
    those inside an arbitrary bounding box.

    Parameters
    ----------
    gdf : geodataframe
        currently accepts a projected gdf

    network_type : string
        network_type as defined in osmnx

    Returns
    -------
    networkx multidigraph
    """
    # generate convex hull around the gdf
    boundary = gdf_convex_hull(gdf)

    # download the highway network within this boundary
    street_graph = ox.graph_from_polygon(boundary, network_type,
                                        simplify=True, retain_all=True,
                                        truncate_by_edge=True, clean_periphery=False)
    # remove duplicates which make polygonization fail
    street_graph = ox.get_undirected(street_graph)

    return street_graph

def streets_from_street_graph(street_graph):
    """
    Updated from osmuf https://github.com/AtelierLibre/osmuf
    Convert a networkx multidigraph to a GeoDataFrame.

    Primarily here to allow future filtering of streets data for osmuf purposes

    Parameters
    ----------
    street_graph : networkx multidigraph

    Returns
    -------
    GeoDataFrame
    """

    # convert to gdf
    streets = ox.graph_to_gdfs(street_graph, nodes=False)
    # write index into a column
    streets['street_id'] = streets.index

    # insert filtering/processing here for OSMuf purposes

    return streets

def gdf_convex_hull(gdf):
    
    """
    Updated from osmuf https://github.com/AtelierLibre/osmuf
    Creates a convex hull around the total extent of a GeoDataFrame.

    Used to define a polygon for retrieving geometries within. When calculating
    densities for urban blocks we need to retrieve the full extent of e.g.
    buildings within the blocks, not crop them to an arbitrary bounding box.

    Parameters
    ----------
    gdf : geodataframe
        currently accepts a projected gdf

    Returns
    -------
    shapely polygon
    """
    ### INSERT CHECK FOR CRS HERE?

    # project gdf back to geographic coordinates as footprints_from_polygon
    # requires it
    gdf_temp = ox.projection.project_gdf(gdf, to_latlong=True)
    # determine the boundary polygon to fetch buildings within
    # buffer originally 0.000225, buffer actually needs to go whole block away
    # to get complete highways therefor trying 0.001
    boundary=gdf_temp.cascaded_union.convex_hull.buffer(0.001)
    # NOTE - maybe more efficient to generate boundary first then reproject second?

    return boundary

def graph_to_polygons(G, node_geometry=True, fill_edge_geometry=True):
    """
    Updated from osmuf https://github.com/AtelierLibre/osmuf
    Convert the edges of a graph into a GeoDataFrame of polygons.

    Parameters
    ----------
    G : networkx multidigraph

    node_geometry : bool
        if True, create a geometry column from node x and y data

    fill_edge_geometry : bool
        if True, fill in missing edge geometry fields using origin and
        destination nodes
    Returns
    -------
    GeoDataFrame
        gdf_polygons
    """

    # create a list to hold edges
    edges = []
    # loop through the edges in the graph
    for u, v, key, data in G.edges(keys=True, data=True):

        # for each edge, add key and all attributes in data dict to the
        # edge_details
        edge_details = {'u':u, 'v':v, 'key':key}
        for attr_key in data:
            edge_details[attr_key] = data[attr_key]

        # if edge doesn't already have a geometry attribute, create one now
        # if fill_edge_geometry==True
        if 'geometry' not in data:
            if fill_edge_geometry:
                point_u = Point((G.nodes[u]['x'], G.nodes[u]['y']))
                point_v = Point((G.nodes[v]['x'], G.nodes[v]['y']))
                edge_details['geometry'] = LineString([point_u, point_v])
            else:
                edge_details['geometry'] = np.nan

        edges.append(edge_details)

    # extract the edge geometries from the list of edge dictionaries
    edge_geometry = []

    for edge in edges:
        edge_geometry.append(edge['geometry'])

    # create a list to hold polygons
    polygons = []

    polygons = list(polygonize(edge_geometry))

    # Create a GeoDataFrame from the list of polygons and set the CRS

    # an option here is to feed it a list of dictionaries with a 'geometry' key
    # this would be one step e.g.
    # gdf_polys = gpd.GeoDataFrame(polygons)

    # Greate an empty GeoDataFrame
    gdf_polygons = gpd.GeoDataFrame()

    # Create a new column called 'geometry' to the GeoDataFrame
    gdf_polygons['geometry'] = None

    # Assign the list of polygons to the geometry column
    gdf_polygons.geometry = polygons

    # Set the crs
    gdf_polygons.crs = G.graph['crs']
    gdf_polygons.gdf_name = '{}_polygons'.format('name') #G.graph['name'])

    return gdf_polygons

def gen_regularity(gdf):
    """
    Updated from osmuf https://github.com/AtelierLibre/osmuf
    Returns a geodataframe of  smallest enclosing 'circles' (shapely polygons)
    generated from the input geodataframe of polygons.

    It includes columns that contain the area of the original polygon and the
    circles and the 'form factor' ratio of the area of the polygon to the area
    of the enclosing circle.

    Parameters
    ----------
    poly_gdf : geodataframe
        a geodataframe containing polygons.

    Returns
    -------
    GeoDataFrame
    """
    #gdf_regularity = gdf[['morpho_Tile_id', 'geometry']].copy()
    gdf_regularity = gdf[['geometry']].copy()

    # write the area of each polygon into a column
    gdf_regularity['poly_area_m2'] = gdf.area.round(decimals=1)

    # replace the polygon geometry with the smallest enclosing circle
    gdf_regularity['geometry'] = gdf_regularity['geometry'].apply(circlizer)

    # calculate the area of the smallest enclosing circles
    gdf_regularity['SEC_area_m2'] = gdf_regularity.area.round(decimals=1)

    # calculate 'regularity' as "the ratio between the area of the polygon and
    # the area of the circumscribed circle C" Barthelemy M. and Louf R., (2014)
    gdf_regularity['regularity'] = gdf_regularity['poly_area_m2']/gdf_regularity['SEC_area_m2']

    return gdf_regularity

def circlizer(x):
    """Updated from osmuf https://github.com/AtelierLibre/osmuf"""
    # takes a shapely polygon or multipolygon and returns a shapely circle polygon
    SEC = sec()
    (cx, cy, buff) = SEC.make_circle(extract_poly_coords(x))
    donut = Point(cx, cy).buffer(buff)

    return donut

def extract_poly_coords(geom):
    """Updated from osmuf https://github.com/AtelierLibre/osmuf"""
    # extract the coordinates of shapely polygons and multipolygons
    # as a list of tuples
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
        interior_coords = []
        for interior in geom.interiors:
            interior_coords += interior.coords[:]
    elif geom.type == 'MultiPolygon':
        exterior_coords = []
        interior_coords = []
        for part in geom:
            epc = extract_poly_coords(part)  # Recursive call
            exterior_coords += epc['exterior_coords']
            interior_coords += epc['interior_coords']
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))

    return exterior_coords + interior_coords

def create_regularity_gdf(morpho_Tiles_gdf):
    """Updated from osmuf https://github.com/AtelierLibre/osmuf"""

    # poly_area_m2 --> block area 
    # SEC_area_m2 --> smallest enclosing circle
    # geometry --> smallest enclosing circle geometry
    # calculate "regularity" as "the ratio between the area of the polygon and
    # the area of the circumscribed circle C" Barthelemy M. and Louf R., (2014)
    # regularity  --> poly_area_m2 / SEC_area_m2
    regularity_gdf = gen_regularity(morpho_Tiles_gdf)

    # Add the radius of SEC_R_m circumscribed circle
    regularity_gdf['SEC_area_m2'] = regularity_gdf['SEC_area_m2'].astype('float')
    regularity_gdf['SEC_R_m'] = np.sqrt(regularity_gdf['SEC_area_m2']/math.pi)

    # Add the position of the center of the block circumscribed circle
    regularity_gdf['centroid'] = regularity_gdf['geometry'].centroid   

    return regularity_gdf 

##################### Smallest enclosing circle ###########################

class sec(object):
    """ 
    Smallest enclosing circle - Library (Python)
    https://www.nayuki.io/page/smallest-enclosing-circle
    
    Copyright (c) 2018 Project Nayuki
    https://www.nayuki.io/page/smallest-enclosing-circle
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with this program (see COPYING.txt and COPYING.LESSER.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Data conventions: A point is a pair of floats (x, y). A circle is a triple of floats (center x, center y, radius).

    Returns the smallest circle that encloses all the given points. Runs in expected O(n) time, randomized.
    Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
    Output: A triple of floats representing a circle.
    Note: If 0 points are given, None is returned. If 1 point is given, a circle of radius 0 is returned.
    """
    # Initially: No boundary points known
    def make_circle(self, points):
        # Convert to float and randomize order
        shuffled = [(float(x), float(y)) for (x, y) in points]
        random.shuffle(shuffled)
        
        # Progressively add points to circle or recompute circle
        c = None
        for (i, p) in enumerate(shuffled):
            if c is None or not self.is_in_circle(c, p):
                c = self._make_circle_one_point(shuffled[ : i + 1], p)
        return c


    # One boundary point known
    def _make_circle_one_point(self, points, p):
        c = (p[0], p[1], 0.0)
        for (i, q) in enumerate(points):
            if not self.is_in_circle(c, q):
                if c[2] == 0.0:
                    c = self.make_diameter(p, q)
                else:
                    c = self._make_circle_two_points(points[ : i + 1], p, q)
        return c


    # Two boundary points known
    def _make_circle_two_points(self, points, p, q):
        circ = self.make_diameter(p, q)
        left  = None
        right = None
        px, py = p
        qx, qy = q
        
        # For each point not in the two-point circle
        for r in points:
            if self.is_in_circle(circ, r):
                continue
            
            # Form a circumcircle and classify it on left or right side
            cross = self._cross_product(px, py, qx, qy, r[0], r[1])
            c = self.make_circumcircle(p, q, r)
            if c is None :
                continue
            elif cross > 0.0 and (left is None or self._cross_product(px, py, qx, qy, c[0], c[1]) > self._cross_product(px, py, qx, qy, left[0], left[1])):
                left = c
            elif cross < 0.0 and (right is None or self._cross_product(px, py, qx, qy, c[0], c[1]) < self._cross_product(px, py, qx, qy, right[0], right[1])):
                right = c
        
        # Select which circle to return
        if left is None and right is None:
            return circ
        elif left is None:
            return right
        elif right is None:
            return left
        else:
            return left if (left[2] <= right[2]) else right


    def make_diameter(self, a, b):
        cx = (a[0] + b[0]) / 2.0
        cy = (a[1] + b[1]) / 2.0
        r0 = math.hypot(cx - a[0], cy - a[1])
        r1 = math.hypot(cx - b[0], cy - b[1])
        return (cx, cy, max(r0, r1))


    def make_circumcircle(self, a, b, c):
        # Mathematical algorithm from Wikipedia: Circumscribed circle
        ox = (min(a[0], b[0], c[0]) + max(a[0], b[0], c[0])) / 2.0
        oy = (min(a[1], b[1], c[1]) + max(a[1], b[1], c[1])) / 2.0
        ax = a[0] - ox;  ay = a[1] - oy
        bx = b[0] - ox;  by = b[1] - oy
        cx = c[0] - ox;  cy = c[1] - oy
        d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
        if d == 0.0:
            return None
        x = ox + ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
        y = oy + ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d
        ra = math.hypot(x - a[0], y - a[1])
        rb = math.hypot(x - b[0], y - b[1])
        rc = math.hypot(x - c[0], y - c[1])
        return (x, y, max(ra, rb, rc))


    def is_in_circle(self, c, p):
        _MULTIPLICATIVE_EPSILON = 1 + 1e-14
        return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * _MULTIPLICATIVE_EPSILON


    # Returns twice the signed area of the triangle defined by (x0, y0), (x1, y1), (x2, y2).
    def _cross_product(self, x0, y0, x1, y1, x2, y2):
        return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)

