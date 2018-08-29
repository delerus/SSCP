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
from darcymodel import darcy_model
from scipy.interpolate import interp1d
import numpy as np

class Simulate():

    def __init__(self,model='perfusion',n_step=50,save=True,n_cycles=1):

        self.model = model
        self.n_step = n_step 
        if model == 'perfusion':
            self.mod = darcy_model()
            self.set_boundry_perfusion()
        self.save=save
        self.n_cycles=n_cycles

        self.set_timesteps()
        self.set_pressure()
        
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

    def simulate(self):
        if self.model == 'perfusion':
            for t,init_p in zip(self.timesteps,self.pressure):
                
                self.pD.p = init_p
                solve(self.mod.F==0, self.mod.p, self.bc)

