""" Basic geogroups in OSMsc. """

from .utils import json_to_gdf


class polygon_group(object):
    """
    polygon_group class for urban polygon objects
    """
    def __init__(self, bbox = None , file_path = None, overpass_query = None,
                trans_type = None, place_name = None):

        # bbox limits the study area
        # bbox --> (South, West, North, East)
        self.bbox = bbox
        self.overpass_query = overpass_query
        self.file_path = file_path
        self.data_type = "Polygon"
        self.trans_type = trans_type
        self.place_name = place_name
    
class point_group(object):
    """
    polygon_group class for urban point objects, like any tagged OSM points.
    """

    def __init__(self, bbox = None , overpass_query = None):
        self.bbox = bbox 
        self.overpass_query = overpass_query
        self.data_type = "Point"     

    def query(self):
        # At present, use subclass query() function
        return None

    def get_gdf(self):
        """
        Obtain OSM data and save as GeoDataFrame.

        Returns
        -------
        GeoDataFrame
        """
        return json_to_gdf(osm_json= self.query(), data_type= self.data_type)
    
class line_string_group(object):
    """
    line_string_group for urban line objects, like street, railways, etc.
    Possibly made for street graph or individual street
    Three way to retrive street graph: 
        1. input bbox into osmnx.graph_from_bbox
        2. self.scope to limit the study area and input it in osmuf
        3. user-defined overpass API
    """

    def __init__(self, bbox = None, BuildingGroup_gdf = None, overpass_query = None):
        self.bbox = bbox
        self.overpass_query = overpass_query
        # LineString or Graph from BuildingGroup_gdf scope
        self.scope = BuildingGroup_gdf
        self.data_type = "LineString"   

    def query(self):
        # At present, use subclass query() function
        return None

    def get_gdf(self, tags = False):
        """
        Obtain OSM data and save as GeoDataFrame.

        Parameters
        ----------
        tags : bool
            if False, the GeoDataFrame won't add OSM "tags" column.
            if True, need to extract tag info into current GeoDataFrame


        Returns
        -------
        GeoDataFrame
        """

        return json_to_gdf(osm_json= self.query(), data_type= self.data_type, tags= tags)




