from fenics import *
from mshr import *
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class Left_ventrical_geometry:
    """
    Sets up a left ventrical geometry, used to run simulations in SSCP
    """

    def __init__(self,param=None):

        
        if param is None:
            self.set_param(param)
        self.set_mesh()
        self.set_markers()
        self.set_applied_pressure()

    def set_param(self,param=None):
        """
        Setting the parameters of the model, if none parameters is introduced, it sets the default.
        """
        
        if param:
            for k in param:
                setattr(self,k,Constant(param[k]))
    
        else:
            param = {'k':0.02, 'D':10**-5, 'k1':0.5, 'k2':5, 'k3':10, 'beta12':0.02, 'beta23':0.05, 'R_12':1,'R_23':1}
            for k in param:
                setattr(self,k,Constant(param[k]))


    def set_mesh(self,mesh_path=None):
        """
        Create the mesh for the geometry, the mesh design to represent the left ventrical
        """

        self.mesh = Mesh()
        
        f = XDMFFile(self.mesh.mpi_comm(),"Files/pressure_mesh.xdmf")
        f.read(self.mesh)
        f.close()


    def set_markers(self,markers=None):
        """
        Creates the markers for the left ventrical, this is where the blood will be introduced.
        """

        if markers:
            self.markers = MeshFunction("size_t",self.mesh,markers)
        else: 
            self.markers = MeshFunction("size_t",self.mesh,"Files/pressure_markers.xml")

        def set_applied_pressure(self,path=None,unit=None):
        """
        Sets the pressure to be aplied as a boundry condition to the geometry
        """

        if path is None:
            path = "Files/coronary_pressure.csv"
            unit = 'mmHg'

        df = pd.read_csv(path,names=['ti','pre'])
        self.timesteps = np.array(df['ti'])

        if unit == 'mmHg':
            self.pressure = np.array(df['pre'])
            self.pressure = self.pressure * 0.133322368 #Converting to kPa

        elif unit =='kPa':
            self.pressure = np.array(df['pre'])

