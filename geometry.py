from fenics import *
from mshr import *
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d



class left_ventrical_geometry():

    def __init__(self,param=None,):

        if param is None:

            dt = 0.02
            diff_02 = 10**-5

            self.k1 = Constant(0.5) # Permability Constant in comp 1
            self.k2 = Constant(5) # Permability Constant in comp 2
            self.k3 = Constant(10) # Permability Constant in comp 3

            self.beta12 = Constant(0.02) #Intercompartment permability
            self.beta23 = Constant(0.05) #Intercompartment permability

            self.D = Constant(diff_02)
            self.k = Constant(dt)

            self.R_12 = Constant(1)
            self.R_23 = Constant(1)

        self.set_mesh()
        self.set_markers()
        self.set_applied_pressure()

    def set_param(defult=True):
        """
        Setting the parameters of the model
        """
        if default:
        
            dt = 0.02
            diff_02 = 10**-5

            self.k1 = Constant(0.5) # Permability Constant in comp 1
            self.k2 = Constant(5) # Permability Constant in comp 2
            self.k3 = Constant(10) # Permability Constant in comp 3

            self.beta12 = Constant(0.02) #Intercompartment permability
            self.beta23 = Constant(0.05) #Intercompartment permability

            self.D = Constant(diff_02)
            self.k = Constant(dt)

            R_12 = Constant(1)
            R_23 = Constant(1)


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

