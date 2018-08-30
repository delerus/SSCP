"This module will take care of the simulations"
"""
Things i want this class to handle:
    1. number of simulations:
    2. Filesaving/naming
        - Both creating file, saving file, making a description file of parameters, closing files
    3. Solver

"""
from fenics import *
from mshr import *
from model import perfusion_model
from scipy.interpolate import interp1d
import numpy as np
import time

parameters['allow_extrapolation'] = True
WARNING = 30
set_log_level(WARNING)


class Simulate():

    def __init__(self,model='perfusion',n_step=50,save=False,n_cycles=1):

        self.model = model
        self.n_step = n_step 
        if model == 'perfusion':
            self.mod = perfusion_model()
            self.set_boundry_perfusion()
        self.save=save
        self.n_cycles=n_cycles

        self.set_timesteps()
        self.set_pressure()
        self.save = save

    def set_timesteps(self):
        """
        Sets timesteps that will be run trough by the simulation
        """
        org_time = self.mod.geo.timesteps
        self.timesteps = np.linspace(org_time[0],org_time.max(),self.n_step)

    def set_pressure(self):

        pressure_function = interp1d(self.mod.geo.timesteps,self.mod.geo.pressure)
        self.pressure = pressure_function(self.timesteps)

    def set_boundry_perfusion(self):
        self.pD = Expression("p",p=0.0,degree=2)
        self.bc = DirichletBC(self.mod.FS.sub(0),self.pD,self.mod.geo.markers,1)

    def open_save_files(self):
        timestamp = time.time()
        p1_name ='Results/Perfusion/'+str(timestamp)+'/p1.xdmf'
        p2_name ='Results/Perfusion/'+str(timestamp)+'/p2.xdmf'
        p3_name ='Results/Perfusion/'+str(timestamp)+'/p3.xdmf'

        self.xdmffile_p1 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p1_name)
        self.xdmffile_p2 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p2_name)
        self.xdmffile_p3 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p3_name)
    
    def save_files(self,p,t):

        p1,p2,p3 = p.split()
        self.xdmffile_p1.write(p1,t)
        self.xdmffile_p2.write(p2,t)
        self.xdmffile_p3.write(p3,t)

    def close_save_files(self):
        
        self.xdmffile_p1.close()
        self.xdmffile_p2.close()
        self.xdmffile_p3.close()

    def simulate(self):
        if self.save:
            self.open_save_files()
        if self.model == 'perfusion':
            for t,init_p in zip(self.timesteps,self.pressure):
                print(t)
                self.pD.p = init_p
                solve(self.mod.F==0, self.mod.p, self.bc)

                if self.save:
                    self.save_files(self.mod.p,t)
        if self.save:
            self.close_save_files()
