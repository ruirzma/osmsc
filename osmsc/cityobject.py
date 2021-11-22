"""Construct OSMsc objects"""

import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
from shapely.ops import polygonize

from .utils import download, street_graph_from_gdf, graph_from_gdfs, json_to_gdf, streets_from_street_graph
from .geogroup import point_group, polygon_group, line_string_group

#TODO logging 


#################################################################################

class building_group(polygon_group): 
    """
    Construct OSMsc building objects
    """
    def query(self):
        """
        Query objects from OpenStreetMap via Overpass API

        Returns
        -------
        JSON
        """
        # Default query
        # query buildings with all kinds of tags
        # key is "building"
        if not self.overpass_query: 

            # print("The default setting is to query all buildings, if necessary, please enter your BuildingGroup overpass api query!!")
            # not including underground structure

            self.overpass_query = """
                [out:json][timeout:50];
                ( way["building"][!"building:levels:underground"]""" + str(self.bbox) +  """; 
                );
                out geom;
                """        
        
        # return data in JSON
        return download(self.overpass_query) 

    def read_file(self):
        """
        Read the local building datasets (if any)

        Returns
        -------
        GeoDataFrame
        """

        return gpd.GeoDataFrame.from_file(self.file_path)

    def get_gdf(self, tags = True, building_levels = False):
        """
        Obtain OSM data and save as GeoDataFrame.

        Parameters
        ----------
        tags : bool
            if True, need to extract tag info into current GeoDataFrame
            if False, the GeoDataFrame won't add OSM "tags" column.
            
        building_levels : bool
            if True, need to download building level into current GeoDataFrame
            if False, the GeoDataFrame won't add OSM "building_levels" column.

        Returns
        -------
        GeoDataFrame
        """

        temp_gdf = json_to_gdf(osm_json= self.query(), data_type= self.data_type, 
                            tags= tags, building_levels = building_levels) 
        # projection
        temp_gdf_prj = ox.project_gdf(temp_gdf)

        # Add GeoDataFrame columns
        temp_gdf["osmscID"] = ["Building_"+ str(i) for i in temp_gdf["osmid"]]
        temp_gdf["Building_area"] = temp_gdf_prj["geometry"].area 
        temp_gdf["Building_perimeter"] = temp_gdf_prj["geometry"].length 

        return temp_gdf

class vegetation_group(polygon_group):
    """
    Construct OSMsc vegetation objects
    """
    def query(self):
        """
        Query objects from OpenStreetMap via Overpass API

        Returns
        -------
        JSON
        """
        if not self.overpass_query:  
            #print("If necessary, please enter your VegetationGroup overpass api query!!")
            self.overpass_query = """
                [out:json][timeout:50];
                ( way["leisure"~"park"]""" + str(self.bbox) +  """; 
                way["landuse"~"grass"]""" + str(self.bbox) +  """; 
                way["leisure"~"pitch"]""" + str(self.bbox) +  """; 
                way["leisure"~"garden"]""" + str(self.bbox) +  """; 
                way["natural"~"scrub"]""" + str(self.bbox) +  """; 
                way["landuse"~"farmyard"]""" +str(self.bbox) +  """; 
                way["landuse"~"recreation_ground"]""" + str(self.bbox) +  """; 
                way["leisure"~"playground"]""" +str(self.bbox) +  """; 
                way["natural"~"wetland"]""" +str(self.bbox) +  """; 
                way["natural"~"wood"]""" +str(self.bbox) +  """; 
                way["landuse"~"meadow"]""" +str(self.bbox) +  """; 
                way["landuse"~"cemetery"]""" +str(self.bbox) +  """; 
                
                );
                out geom;
                """

        # return data in JSON
        return download(self.overpass_query) 

    def get_gdf(self, tags = True):
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

        temp_gdf = json_to_gdf(osm_json= self.query(), data_type= self.data_type, 
                            tags= tags) 

        temp_gdf_prj = ox.project_gdf(temp_gdf)

        temp_gdf["osmscID"] = ["Vegetation_"+ str(i) for i in temp_gdf["osmid"]]
        temp_gdf["Vegetation_area"] = temp_gdf_prj["geometry"].area
        temp_gdf["Vegetation_perimeter"] = temp_gdf_prj["geometry"].length

        return temp_gdf

class waterbody_group(polygon_group):
    """
    Construct OSMsc waterbody objects
    """
    def query(self):
        """
        Query objects from OpenStreetMap via Overpass API

        Returns
        -------
        JSON
        """
        if not self.overpass_query:  
            #print("If necessary, please enter your WaterbodyGroup overpass api query!!")
            self.overpass_query = """
                [out:json][timeout:50];
                ( way["leisure"~"ice_rink"]""" + str(self.bbox) +  """; 
                way["landuse"~"swimming_pool"]""" + str(self.bbox) +  """; 
                way["man_made"~"water_tower"]""" + str(self.bbox) +  """; 
                way["natural"~"water"]""" +str(self.bbox) +  """; 
                way["reservoir"~"water_storage"]""" +str(self.bbox) +  """; 

                );
                out geom;
                """
        # return data in JSON        
        return download(self.overpass_query)

    def get_gdf(self, tags = True):
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

        temp_gdf = json_to_gdf(osm_json= self.query(), data_type= self.data_type, 
                            tags= tags) 

        temp_gdf_prj = ox.project_gdf(temp_gdf)

        temp_gdf["osmscID"] = ["WaterBody_"+ str(i) for i in temp_gdf["osmid"]]
        temp_gdf["Waterbody_area"] = temp_gdf_prj["geometry"].area
        temp_gdf["Waterbody_perimeter"] = temp_gdf_prj["geometry"].length

        return temp_gdf

class transportation_group(polygon_group):
    """
    Construct OSMsc transportation objects
    """

    def get_gdf_prj(self, street_width = 3):
        """
        Obtain OSM data and save as GeoDataFrame.
        To define the street width, this function has to use the 
        projected coordinate system.

        Parameters
        ----------
        street_width : float or int
            Assumed street width for OSM street

        Returns
        -------
        GeoDataFrame
        """

        building_gdf = building_group(bbox= self.bbox).get_gdf()
        # 
        street_graph = street_graph_from_gdf(building_gdf)
        # Street geometry is LineString       
        street_temp_gdf = streets_from_street_graph(street_graph)

        # only keep user-defined trans_type
        if self.trans_type is not None: 
            street_type_1_gdf = street_temp_gdf[ street_temp_gdf["highway"]== list(self.trans_type)[0] ]
            
            concat_gdf = street_type_1_gdf
            for i in range(len(self.trans_type)-1):
                street_type = list(self.trans_type)[i+1]
                street_type_2_gdf = street_temp_gdf[ street_temp_gdf["highway"]== street_type ]

                concat_gdf = pd.concat([concat_gdf,street_type_2_gdf])
        # filtered street_temp_gdf
        street_temp_gdf = concat_gdf
        # set crs
        street_temp_gdf = gpd.GeoDataFrame(street_temp_gdf,crs='epsg:4326')


        # Empty street objects
        street_gdf_prj = gpd.GeoDataFrame()
        # re-project geometries to a projected CRS before buffer function
        street_temp_gdf_prj = ox.project_gdf(street_temp_gdf)
        # Assume street width is 3m and construct street outlines
        street_gdf_prj["geometry"] = street_temp_gdf_prj["geometry"].buffer(street_width, cap_style = 3).exterior 
        # Polygon street objects
        street_gdf_prj["geometry"] = list(polygonize(street_gdf_prj["geometry"]))
        # set crs
        street_gdf_prj.crs = street_temp_gdf_prj.crs
        # Add osmscID column
        street_gdf_prj["osmscID"] = ["Transportation_"+str(i) for i in range(len(street_gdf_prj))]

        return street_gdf_prj

class urban_patch_group(polygon_group):
    """
    Construct OSMsc urban patch objects
    """
    def get_gdf_prj(self):
        """
        Construct urban patch objects
        To define the street width, this function has to use the 
        projected coordinate system.

        Returns
        -------
        GeoDataFrame
        """
        # Download the street dataframe
        street_gdf_prj = transportation_group(bbox= self.bbox, trans_type= self.trans_type).get_gdf_prj()
        # Merged streets are the basis for generating blocks
        street_dis_geom = street_gdf_prj.dissolve().iloc[0].geometry       
        
        # Create urban_patch_gdf
        urban_patch_gdf_prj = gpd.GeoDataFrame()
        # All streets are merged into a polygon
        urban_patch_gdf_prj["geometry"] = list(polygonize(street_dis_geom.boundary))
        # delete the whole street polygon
        urban_patch_gdf_prj = urban_patch_gdf_prj.drop(index=0)
        # re-index
        urban_patch_gdf_prj.index = urban_patch_gdf_prj.index - 1
        # set crs
        urban_patch_gdf_prj.crs = street_gdf_prj.crs

        # Add osmscID column
        urban_patch_gdf_prj["osmscID"] = ["UrbanPatch_"+str(i) for i in range(len(urban_patch_gdf_prj))]
        # Add other attr columns       
        urban_patch_gdf_prj["UrbanPatch_area"] = urban_patch_gdf_prj["geometry"].area
        urban_patch_gdf_prj["UrbanPatch_perimeter"] = urban_patch_gdf_prj["geometry"].length

        return urban_patch_gdf_prj


####### TODO More city objects in OSMsc, namely, networks and points ########
class street_network(line_string_group):
    """
    Construct street networks, rather than street polygons
    """
    def query(self):
        """
        Query objects from OpenStreetMap via Overpass API

        Returns
        -------
        JSON
        """
        # Default query
        if not self.overpass_query: 

            print("The default setting is to query all street, if necessary, please enter your StreetNetwork overpass api query!!")

            self.overpass_query = """
                [out:json][timeout:50];
                ( way["highway"]["area"!~"yes"]""" + str(self.bbox) +  """; 
                );
                out geom;
                """

        # return data in JSON format
        return download(self.overpass_query)

    def query_with_network_type(self, network_type):
        """
        Query streets of a specific street type from OpenStreetMap via OSMnx
        More info about OSMnx: https://github.com/gboeing/osmnx 

        Parameters
        ----------
        network_type : string
            which type of street network to get, chose from
            {"all_private", "all", "bike", "drive", "drive_service", "walk"}

        Returns
        -------
        Graph
        """

        if self.bbox is None:
            print("Need to input bbox for StreetNetwork")
            return None
        else:
            (south, west, north, east) = self.bbox
            self.street_graph = ox.graph_from_bbox(south= south, west=west, north=north, 
                                                    east=east, network_type = network_type)
            return self.street_graph

    def query_all_streets_from_buildings(self):
        """
        Query all_streets from OpenStreetMap via OSMnx
        More info about OSMnx: https://github.com/gboeing/osmnx 

        Returns
        -------
        Graph
        """

        if self.scope is None:
            print("Need to input BuildingGroup_gdf for StreetNetwork")
            return None
        else:
            self.street_graph = street_graph_from_gdf(self.scope)
            return self.street_graph

    def query_streets_with_tags_from_buildings(self, primary = True, secondary = True, tertiary = False, 
                                            residential = False, service = False):
        # TODO Needed to be modified later
        """
        Query streets of a specific tag from OpenStreetMap via OSMnx
        More info about OSMnx: https://github.com/gboeing/osmnx 

        Parameters
        ----------
        primary : bool
        secondary : bool
        tertiary : bool
        residential : bool
        service : bool
            Specific tags needed to be kept

        Returns
        -------
        Graph
        """    
        if self.scope is None:
            print("Need to input BuildingGroup_gdf for StreetNetwork")
            return None
        else:
            # Only based on the query result from OSMnx and "highway" tags
            self.street_graph = street_graph_from_gdf(self.scope)
            street_edges = ox.graph_to_gdfs(self.street_graph)[1]
            street_nodes = ox.graph_to_gdfs(self.street_graph)[0]

            if primary:
                temp_gdf_1 = street_edges[street_edges["highway"] == 'primary']
                temp_gdf_edges = temp_gdf_1
            if secondary:
                temp_gdf_2 = street_edges[street_edges["highway"] == 'secondary']
                temp_gdf_edges = pd.concat([temp_gdf_edges, temp_gdf_2])
            if tertiary:
                temp_gdf_3 = street_edges[street_edges["highway"] == 'tertiary']
                temp_gdf_edges = pd.concat([temp_gdf_edges, temp_gdf_3])        
            if residential:
                temp_gdf_4 = street_edges[street_edges["highway"] == 'residential']
                temp_gdf_edges = pd.concat([temp_gdf_edges, temp_gdf_4])    
            if service:
                temp_gdf_5 = street_edges[street_edges["highway"] == 'service']
                temp_gdf_edges = pd.concat([temp_gdf_edges, temp_gdf_5])  

            # u_list = list(temp_gdf_edges["u"])
            # v_list = list(temp_gdf_edges["v"])
            u_list = list(temp_gdf_edges["from"])
            v_list = list(temp_gdf_edges["to"])
            nodes_within_edges = u_list + v_list

            # Only keep the nodes located in the specific edges
            temp_gdf_nodes = gpd.GeoDataFrame()
            for i in range(len(list(street_nodes.index))):
                node_index = list(street_nodes.index)[i]
                if node_index in nodes_within_edges:
                    temp_gdf_nodes = pd.concat([temp_gdf_nodes,street_nodes.loc[[node_index]]])

            temp_gdf_nodes.gdf_name = "'unnamed_nodes'"

            # Generate graph from nodes and edges
            self.street_graph = graph_from_gdfs(temp_gdf_nodes,temp_gdf_edges,)

            return self.street_graph


class tagged_point_group(point_group):
    """
    Construct OSMsc tagged points.
    """
    def query(self):
        """
        Query objects from OpenStreetMap via Overpass API.
        Please input your query condition, if any.

        Returns
        -------
        JSON
        """
        if not self.overpass_query:  

            print("The default setting is to query tagged points, if necessary, please enter your TaggedPointGroup overpass api query!!")

            self.overpass_query = """
                [out:json][timeout:50];
                ( 
                node""" + str(self.bbox) +  """[~"."~"."]; 
                );
                out geom;
                """

        return download(self.overpass_query)  # json 



