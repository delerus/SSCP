from fenics import *
from mshr import *
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d



class left_ventrical_geometry():

    def __init__(self,param=None):

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

            R_12 = Constant(1)
            R_23 = Constant(1)

        self.mesh = self.set_mesh()

    def set_param(defult=True):

        if param is True:
        
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
