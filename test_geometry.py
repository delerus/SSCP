from fenics import *
from mshr import *

from Cardiac_model import left_ventrical_geometry
from matplotlib import pyplot as plt

from time import time
def test_set_mesh():
    """
    This test tests if the mesh is sat propertly
    """
    
    geo = left_ventrical_geometry()
    print(geo.mesh)



