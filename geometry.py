from fenics import *
from mshr import *
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class Left_ventrical_geometry:

    def __init__(self,param=None):

        
        if param is None:
            self.set_param(param)
        self.set_mesh()
        self.set_markers()
        self.set_applied_pressure()

    def set_param(self,param=None):
        """
        Setting the parameters of the model
        """
        
        if param:
            for k in param:
                setattr(self,k,Constant(param[k]))
    
        else:
            param = {'k':0.02, 'D':10**-5, 'k1':0.5, 'k2':5, 'k3':10, 'beta12':0.02, 'beta23':0.05, 'R_12':1,'R_23':1}
            for k in param:
                setattr(self,k,Constant(param[k]))


    def set_mesh(self,mesh_path=None):
        

        self.mesh = Mesh()

        if mesh_path:
            try:
                f = XDMFFile(self.mesh.mpi_comm(),mesh_path)
                f.read(self.mesh)
                f.close()
            except exception as e:
                print(e)
        else:
            f = XDMFFile(self.mesh.mpi_comm(),"Files/pressure_mesh.xdmf")
            f.read(self.mesh)
            f.close()


    def set_markers(self):

        self.markers = MeshFunction("size_t",self.mesh,"Files/pressure_markers.xml")

    def set_applied_pressure(self,path=None,unit = 'mmHg'):
        """
        Sets the pressure to be aplied as a boundry condition to or geometry
        """

        if path is None:
            path = "Files/coronary_pressure.csv"

        df = pd.read_csv(path,names=['ti','pre'])
        self.timesteps = np.array(df['ti'])

        if unit == 'mmHg':
            self.pressure = np.array(df['pre'])
            self.pressure = self.pressure * 0.133322368 #Converting to kPa

        elif unit =='kPa':
            self.pressure = np.array(df['pre'])

