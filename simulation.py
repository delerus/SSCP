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
from model import Perfusion
from model import Advection_diffusion
from scipy.interpolate import interp1d
import numpy as np
import time

parameters['allow_extrapolation'] = True
WARNING = 30
set_log_level(WARNING)


class Simulate:

    def __init__(self,model='perfusion',n_step=50,save=False,n_cycles=1):

        self.model = model
        self.n_step = n_step 
        self.save=save
        self.m_cycles = n_cycles

        self.xdmffile_p1 = None
        self.xdmffile_v1 = None
        self.xdmffile_c1 = None

        if model == 'perfusion':
            print('perfusion er valgt')
            self.mod = Perfusion()
        elif model == 'advection_diffusion':
            print('advection_diffusion valgt')
            self.mod = Advection_diffusion()

        self.set_boundry_perfusion()
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

    def open_save_files(self,p=False,v=False,c=False):
        timestamp = time.time()
        if p:
            p1_name ='Results/Perfusion/'+str(timestamp)+'/p1.xdmf'
            p2_name ='Results/Perfusion/'+str(timestamp)+'/p2.xdmf'
            p3_name ='Results/Perfusion/'+str(timestamp)+'/p3.xdmf'

            self.xdmffile_p1 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p1_name)
            self.xdmffile_p2 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p2_name)
            self.xdmffile_p3 = XDMFFile(self.mod.geo.mesh.mpi_comm(), p3_name)
        if v:
            v1_name ='Results/Perfusion/'+str(timestamp)+'/v1.xdmf'
            v2_name ='Results/Perfusion/'+str(timestamp)+'/v2.xdmf'
            v3_name ='Results/Perfusion/'+str(timestamp)+'/v3.xdmf'

            self.xdmffile_v1 = XDMFFile(self.mod.geo.mesh.mpi_comm(), v1_name)
            self.xdmffile_v2 = XDMFFile(self.mod.geo.mesh.mpi_comm(), v2_name)
            self.xdmffile_v3 = XDMFFile(self.mod.geo.mesh.mpi_comm(), v3_name)
        if c:

            c1_name ='Results/Perfusion/'+str(timestamp)+'/c1.xdmf'
            c2_name ='Results/Perfusion/'+str(timestamp)+'/c2.xdmf'
            c3_name ='Results/Perfusion/'+str(timestamp)+'/c3.xdmf'

            self.xdmffile_c1 = XDMFFile(self.mod.geo.mesh.mpi_comm(), c1_name)
            self.xdmffile_c2 = XDMFFile(self.mod.geo.mesh.mpi_comm(), c2_name)
            self.xdmffile_c3 = XDMFFile(self.mod.geo.mesh.mpi_comm(), c3_name)

    def save_files(self,t,p=None,v=None,c=None):

        if p:
            self.mod.p1,self.mod.p2,self.mod.p3 = self.mod.p.split()
            self.xdmffile_p1.write(self.mod.p1,t)
            self.xdmffile_p2.write(self.mod.p2,t)
            self.xdmffile_p3.write(self.mod.p3,t)
        
        if v:
            self.xdmffile_v1.write(self.mod.vd1,t)
            self.xdmffile_v2.write(self.mod.vd2,t)
            self.xdmffile_v3.write(self.mod.vd3,t)

        if c:
            self.mod.c1,self.mod.c2,self.mod.c3 = self.mod.c.split()
            self.xdmffile_c1.write(self.mod.c1,t)
            self.xdmffile_c2.write(self.mod.c2,t)
            self.xdmffile_c3.write(self.mod.c3,t)


    def close_save_files(self):
        
        if self.xdmffile_p1:
            self.xdmffile_p1.close()
            self.xdmffile_p2.close()
            self.xdmffile_p3.close()
        
        if self.xdmffile_v1:
            self.xdmffile_v1.close()
            self.xdmffile_v2.close()
            self.xdmffile_v3.close()
        
        if self.xdmffile_c1:
            self.xdmffile_c1.close()
            self.xdmffile_c2.close()
            self.xdmffile_c3.close()

    def simulate_perfusion(self):
        if self.save:
            self.open_save_files(p=True)
        if self.model == 'perfusion':
            for t,init_p in zip(self.timesteps,self.pressure):
                print(t)
                self.pD.p = init_p
                solve(self.mod.F==0, self.mod.p, self.bc)

                if self.save:
                    self.save_files(self.mod.p,t)
        if self.save:
            self.close_save_files()

    def simulate_advection_diffusion(self):
        if self.save:
            self.open_save_files(p=True,v=True,c=True)
            c_0 = 1.0
            initc = DirichletBC(self.mod.FSC.sub(0),c_0,self.mod.geo.markers,1)
            initc.apply(self.mod.c_n.vector())
        if self.model == 'advection_diffusion':
            for t,init_p in zip(self.timesteps,self.pressure):
                print(t)
                self.pD.p = init_p
                solve(self.mod.F==0, self.mod.p, self.bc)
                self.mod.calculate_velocity()
#                solve(self.mod.F3==0, self.mod.c)
#                self.mod.c_n.assign(self.mod.c)
                solve(self.mod.F2+self.mod.F3==0, self.mod.c)
                self.mod.c_n.assign(self.mod.c)
                if self.save:
                    self.save_files(t,p=True,v=True,c=True)
        if self.save:
            self.close_save_files()
