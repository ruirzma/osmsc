"""OSMsc init https://github.com/ruirzma/osmsc"""

from .__version__ import __version__
from .cityobject import building_group, vegetation_group, waterbody_group, transportation_group, urban_patch_group
from .grid import Grid
# from .patch import Patch

import osmsc.utils, osmsc.feature, osmsc.fusion
import osmsc.plot,  osmsc.cityjson

# TODO from .__api__ import *   

print ("OSMsc has been successfully imported!!!")


