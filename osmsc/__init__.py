"""OSMsc init https://github.com/Revisedrzma/osmsc"""

from .__version__ import __version__
from .cityobject import building_group, vegetation_group, waterbody_group, transportation_group, urban_Tile_group
from .feature import add_spatial_semantics_attr
from .cityjson import city_json


# osmsc.utils, osmsc.feature, osmsc.fusion
# import osmsc.plot,  osmsc.cityjson

# TODO from .__api__ import *   

print ("OSMsc has been successfully imported!!!")


