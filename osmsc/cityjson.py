""" Construct and output CityJSON object """

import json
from cjio import cityjson
from cjio.models import CityObject, Geometry
from .utils import fill_nan_list_position


class city_json(object):
    """
    Create CityJSON-schema objects, including building, vegetation, waterbody
    transportation and urban patch objects.

    More info about CityJSON can be found in 
    https://www.cityjson.org/specs/1.0.1/
    """

    def __init__(self, building_gdf = None , vegetation_gdf = None, 
    waterbody_gdf = None, transportation_gdf = None, urban_patch_gdf = None):

        self.building_gdf = building_gdf
        self.vegetation_gdf = vegetation_gdf
        self.waterbody_gdf = waterbody_gdf
        self.transportation_gdf = transportation_gdf
        self.urban_patch_gdf = urban_patch_gdf

        # create an empty CityModel
        self.cm = cityjson.CityJSON()


    def create_building_object(self):
        """
        Create CityJSON objects for all buildings
        Main class variable: self.building_gdf, self.cm

        Returns
        -------
        cm : cityjson.CityJSON
        """
        
        # Create a CityJSON object for each building
        for b_index in range(len(self.building_gdf)):
            
            # Create empty building CityJSON object
            buildingObject = CityObject(id = str(self.building_gdf.iloc[b_index].osmscID))
            
            temp_gdf = self.building_gdf.set_index(self.building_gdf["osmscID"])
            attr_dict= json.loads(temp_gdf.to_json())

            b_attrs = attr_dict['features'][b_index]['properties']
            buildingObject.attributes = b_attrs

            #######################################Geometry############################
            b_geom = Geometry(type='Solid', lod=1)

            building_height = float(self.building_gdf.iloc[b_index].Building_height)
            building_poly = self.building_gdf.iloc[b_index].geometry
            
            # Checking if vertices of polygon are in counter-clockwise 是否是逆时针
            building_poly_ccw = building_poly.exterior.is_ccw

            # Extract the point values that define the perimeter of the polygon
            x, y = building_poly.exterior.coords.xy
            x_list = list(x)
            y_list = list(y)

            # bottom surface
            bottom_sur = []
            # top surface
            top_sur = []
            # side surfaces
            side_sur = []
                    
            for i in range(len(x_list)-1):

                # The coordinates of the polygon are counterclockwise
                if  building_poly_ccw: 
                    top_sur.append((x_list[i],y_list[i],building_height))
                    bottom_sur.append((x_list[-i-2],y_list[-i-2],0))

                    side_coor = [(x_list[i],y_list[i],0), (x_list[i+1],y_list[i+1],0),
                                (x_list[i+1],y_list[i+1],building_height), (x_list[i],y_list[i],building_height)]

                    side_sur.append([side_coor])

                # The coordinates of polygon are clockwise            
                else:   
                    top_sur.append((x_list[-i-2],y_list[-i-2],building_height))
                    bottom_sur.append((x_list[i],y_list[i],0))

                    side_coor = [(x_list[i],y_list[i],building_height),(x_list[i+1],y_list[i+1],building_height),
                                (x_list[i+1],y_list[i+1],0),(x_list[i],y_list[i],0)]

                    side_sur.append([side_coor])
                    
            # Building boundary and geometry
            b_bdry =  [[top_sur],[bottom_sur]] + side_sur 
            b_geom.boundaries.append(b_bdry)

            # Add propertries to CityJSON objects
            buildingObject.geometry.append(b_geom)
            buildingObject.type = "Building"

            self.cm.cityobjects[buildingObject.id] = buildingObject

        return self.cm

    def create_vegetation_object(self):
        """
        Create CityJSON objects for all vegetation objects
        Main class variable: self.vegetation_gdf, self.cm

        Returns
        -------
        cm : cityjson.CityJSON
        """
        # # Create a CityJSON object for each vegetation object
        for v_index in range(len(self.vegetation_gdf)):
            
            vegetationObject = CityObject(id = str(self.vegetation_gdf.iloc[v_index].osmscID))
            
            temp_gdf = self.vegetation_gdf.set_index(self.vegetation_gdf["osmscID"])
            attr_dict= json.loads(temp_gdf.to_json())

            v_attrs = attr_dict['features'][v_index]['properties']
            vegetationObject.attributes = v_attrs

            #######################################Geometry############################
            v_geom = Geometry(type='MultiSurface', lod=0)

            vegetation_poly = self.vegetation_gdf.iloc[v_index].geometry
            
            try:
                # Don't consider the content of MultiPolygon at present, 
                # and do it separately later
            
                # Checking if vertices of polygon are in counter-clockwise 是否是逆时针
                vegetation_poly_ccw = vegetation_poly.exterior.is_ccw

                # Extract the point values that define the perimeter of the polygon
                x, y = vegetation_poly.exterior.coords.xy
                x_list = list(x)
                y_list = list(y)
                top_sur = []

                for i in range(len(x_list)-1):
                    # The coordinates of the polygon are counterclockwise
                    if  vegetation_poly_ccw: 
                        top_sur.append((x_list[i],y_list[i],0))
                    # The coordinates of polygon are clockwise
                    else:   
                        top_sur.append((x_list[-i-2],y_list[-i-2],0))

                # Building boundary and geometry
                v_bdry =  [top_sur]
                v_geom.boundaries.append(v_bdry)

                # Add propertries to CityJSON objects
                vegetationObject.geometry.append(v_geom)
                vegetationObject.type = "PlantCover"

                self.cm.cityobjects[vegetationObject.id] = vegetationObject
                
            except:
                pass

        return self.cm

    def create_waterbody_object(self):
        """
        Create CityJSON objects for all waterbody objects
        Main class variable: self.waterbody_gdf, self.cm

        Returns
        -------
        cm : cityjson.CityJSON
        """

        # Create a CityJSON object for each waterbody
        for w_index in range(len(self.waterbody_gdf)):
            
            # Create empty waterbody CityJSON object
            waterbodyObject = CityObject(id= str(self.waterbody_gdf.iloc[w_index].osmscID))
            
            temp_gdf = self.waterbody_gdf.set_index(self.waterbody_gdf["osmscID"])
            attr_dict= json.loads(temp_gdf.to_json())        
            

            w_attrs = attr_dict['features'][w_index]['properties']
            waterbodyObject.attributes = w_attrs

            #######################################Geometry############################
            w_geom = Geometry(type='CompositeSurface', lod=0)

            waterbody_poly = self.waterbody_gdf.iloc[w_index].geometry
            
            # Checking if vertices of polygon are in counter-clockwise 是否是逆时针
            waterbody_poly_ccw = waterbody_poly.exterior.is_ccw

            # Extract the point values that define the perimeter of the polygon
            x, y = waterbody_poly.exterior.coords.xy
            x_list = list(x)
            y_list = list(y)
            
            top_sur = []
            
            for i in range(len(x_list)-1):
                # The coordinates of the polygon are counterclockwise
                if  waterbody_poly_ccw:
                    top_sur.append((x_list[i],y_list[i],0))
                # The coordinates of polygon are clockwise
                else:
                    top_sur.append((x_list[-i-2],y_list[-i-2],0))

            # Building boundary and geometry
            w_bdry =  [top_sur]
            w_geom.boundaries.append(w_bdry)

            # Add propertries to CityJSON objects
            waterbodyObject.geometry.append(w_geom)
            waterbodyObject.type = "WaterBody"

            self.cm.cityobjects[waterbodyObject.id] = waterbodyObject
        return self.cm

    def create_transportation_object(self):
        """
        Create CityJSON objects for all transportation objects
        Main class variable: self.transportation_gdf, self.cm

        Returns
        -------
        cm : cityjson.CityJSON
        """
        # Create a CityJSON object for each transportation object
        # transportation_object refers to the street currently
        for r_index in range(len(self.transportation_gdf)):

            # Create empty transportation CityJSON object           
            roadObject = CityObject(id = str(self.transportation_gdf.iloc[r_index].osmscID))
            
            temp_gdf = self.transportation_gdf.set_index(self.transportation_gdf["osmscID"])
            attr_dict= json.loads(temp_gdf.to_json())        
            
            r_attrs = attr_dict['features'][r_index]['properties']
            roadObject.attributes = r_attrs

            #######################################Geometry############################
            r_geom = Geometry(type='CompositeSurface', lod=0)

            road_poly = self.transportation_gdf.iloc[r_index].geometry
            
            try:
                # Don't consider the content of MultiPolygon at present, 
                # and do it separately later
            
                # Checking if vertices of polygon are in counter-clockwise 是否是逆时针
                road_poly_ccw = road_poly.exterior.is_ccw

                # Extract the point values that define the perimeter of the polygon
                x, y = road_poly.exterior.coords.xy
                x_list = list(x)
                y_list = list(y)

                top_sur = []

                for i in range(len(x_list)-1):
                    # The coordinates of the polygon are counterclockwise
                    if  road_poly_ccw: 
                        top_sur.append((x_list[i],y_list[i],0))
                    # The coordinates of polygon are clockwise
                    else: 
                        top_sur.append((x_list[-i-2],y_list[-i-2],0))

                # Building boundary and geometry
                r_bdry =  [top_sur]
                r_geom.boundaries.append(r_bdry)

                # Add propertries to CityJSON objects
                roadObject.geometry.append(r_geom)
                roadObject.type = "Road"

                self.cm.cityobjects[roadObject.id] = roadObject
                
            except:
                pass
        return self.cm

    def create_urban_patch_object(self):
        """
        Create CityJSON objects for all urban patch objects
        Main class variable: self.urban_patch_gdf, self.cm

        Returns
        -------
        cm : cityjson.CityJSON
        """
        # Create a CityJSON object for each urban patch object
        for u_index in range(len(self.urban_patch_gdf)):
            # Create empty urban patch CityJSON object 
            urbanPatchObject = CityObject(id=str(self.urban_patch_gdf.iloc[u_index].osmscID))

            temp_gdf = self.urban_patch_gdf.set_index(self.urban_patch_gdf["osmscID"])
            attr_dict= json.loads(temp_gdf.to_json())

            u_attrs = attr_dict['features'][u_index]['properties']
            urbanPatchObject.attributes = u_attrs

            #######################################Geometry############################
            u_geom = Geometry(type='CompositeSurface', lod=0)

            urbanPatch_poly = self.urban_patch_gdf.iloc[u_index].geometry
            
            try:
                # Don't consider the content of MultiPolygon at present, 
                # and do it separately later
            
                # Checking if vertices of polygon are in counter-clockwise 是否是逆时针
                urbanPatch_poly_ccw = urbanPatch_poly.exterior.is_ccw

                # Extract the point values that define the perimeter of the polygon
                x, y = urbanPatch_poly.exterior.coords.xy

                x_list = list(x)
                y_list = list(y)

                top_sur = []

                for i in range(len(x_list)-1):
                    # The coordinates of the polygon are counterclockwise
                    if  urbanPatch_poly_ccw:
                        top_sur.append((x_list[i],y_list[i],0))
                    # The coordinates of polygon are clockwise                    
                    else:
                        top_sur.append((x_list[-i-2],y_list[-i-2],0))

                # Building boundary and geometry
                u_bdry =  [top_sur]
                u_geom.boundaries.append(u_bdry)

                # Add propertries to CityJSON objects
                urbanPatchObject.geometry.append(u_geom)
                urbanPatchObject.type = "GenericCityObject"

                self.cm.cityobjects[urbanPatchObject.id] = urbanPatchObject
                
            except:
                pass

        return self.cm

    def drop_extra_geom(self, gdf):
        """
        Drop extra geometry columns in GeoDataframe
        As the GeoDataframe.to_json() function only works for the GeoDataframe that have 
        one geometry column

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            mostly building_gdf and urban_patch_gdf

        Returns
        -------
        gdf : geopandas.GeoDataFrame
        """

        if 'mrr_geometry' in gdf.columns:
            gdf = gdf.drop("mrr_geometry",axis =1)
        if 'mcc_geometry' in gdf.columns:
            gdf = gdf.drop("mcc_geometry",axis =1)

        return gdf

    def output_json(self,filename = "cityName"):
        """
        Output CityJSON object into a JSON file.
        Main class variable: self.building_gdf, self.vegetation_gdf, self.waterbody_gdf
                            self.transportation_gdf, self.urban_patch_gdf, self.cm

        Parameters
        ----------
        filename : str
            the filename to be saved in the current project location
        """


        print("Please wait a few minutes")

        # pre cleaning work
        # fill spatial semantic attrs
        # make each gdf has only one shapely column
        # Build city objects

        ################# Building #################
        if self.building_gdf is not None:
            self.building_gdf = fill_nan_list_position(self.building_gdf,"within_UrbanPatch")
            self.building_gdf = self.drop_extra_geom(self.building_gdf)
            self.cm = self.create_building_object()

        ################# Vegetation #################
        if self.vegetation_gdf is not None:
            self.vegetation_gdf = fill_nan_list_position(self.vegetation_gdf,"within_UrbanPatch")
            self.cm = self.create_vegetation_object()       

        ################# Waterbody #################
        if self.vegetation_gdf is not None:
            self.waterbody_gdf = fill_nan_list_position(self.waterbody_gdf,"within_UrbanPatch")
            self.cm = self.create_waterbody_object()


        ################# Transportation #################
        if self.transportation_gdf is not None:
            self.cm = self.create_transportation_object()

        ################# UrbanPatch #################
        if self.urban_patch_gdf is not None:
            self.urban_patch_gdf = fill_nan_list_position(self.urban_patch_gdf, "contains_Vegetation")
            self.urban_patch_gdf = fill_nan_list_position(self.urban_patch_gdf, "contains_Building")
            self.urban_patch_gdf = fill_nan_list_position(self.urban_patch_gdf, "contains_Waterbody")

            self.urban_patch_gdf = self.drop_extra_geom(self.urban_patch_gdf)

            self.cm = self.create_urban_patch_object()


        # reference_geometry
        cityobjects, vertex_lookup = self.cm.reference_geometry()    
        # multipoint multilinestring multisurface compositesurface 
        # solid multisolid compositesolid
        self.cm.add_to_j(cityobjects,vertex_lookup)
        # cm.j["CityObjects"]

        self.cm.update_bbox()

        # https://www.cityjson.org/tutorials/validation/
        self.cm.validate() 

        # save the CityJSON file
        cityjson.save(self.cm, filename + '.json')




