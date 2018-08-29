"This module will take care of the simulations"
"""
Things i want this class to handle:
    1. number of simulations:
    2. Filesaving/naming
        - Both creating file, saving file, making a description file of parameters, closing files
    3. Solver

"""

from darcymodel import darcy_model
from scipy.interpolate import interp1d

class Simulate():

    def __init__(self,model='perfusion',tstep = 50, ):
        if model = 'perfusion':
            self.mod = darcy_model()
    
    self.set_timesteps(t_steps)
    
    def set_timesteps():
        """
        Sets timesteps that will be run trough by the simulation
        """
        org_time = self.model.geo.timesteps
        self.timesteps = np.linspace(org_time[0],org_time.max(),t_steps())
