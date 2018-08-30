from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
from geometry import left_ventrical_geometry

class perfusion_model():
    """
    This class is gonna include the equations (Maths) and also the definitions of functionspaces etc.
    """
    def __init__(self,geo=None):
        
        if geo is None:
            self.geo = left_ventrical_geometry()
        self.function_space()
        self.test_functions()
        self.funcs()
        self.multi_comp_darcy_equation()

    def function_space(self):

        self.el = tetrahedron
        self.P = FiniteElement('P',self.el,2)

        self.element = MixedElement([self.P,self.P,self.P])
        self.FS = FunctionSpace(self.geo.mesh,self.element)

    def test_functions(self):

        self.q1,self.q2,self.q3 = TestFunctions(self.FS)
    
    def funcs(self):

        self.p = Function(self.FS)
        self.c = Function(self.FS)
        self.c_n = Function(self.FS)

        self.p1, self.p2, self.p3 = split(self.p)
        self.c1, self.c2, self.c3 = split(self.c)
        self.c1_n, self.c2_n, self.c3_n = split(self.c_n)

    def multi_comp_darcy_equation(self):
        """
        The simplified multi compartment dary equation

        A1-3: Intracompartment perfusion in comp 1-3.
        B1: Intercompartment perfusion between comp 1 and 2.
        B2: Intercompartment perfusion between comp 2 and 3.

        S: Sink term
        """
        
        A1 = -self.geo.k1 * dot(grad(self.p1),grad(self.q1))*dx
        A2 = -self.geo.k2 * dot(grad(self.p2),grad(self.q2))*dx
        A3 = -self.geo.k3 * dot(grad(self.p3),grad(self.q3))*dx

        B1 = dot(self.geo.beta12*(self.p1-self.p2),self.q1)*dx + dot(self.geo.beta12*(self.p2-self.p1),self.q2)*dx
        B2 = dot(self.geo.beta23*(self.p2-self.p3),self.q2)*dx + dot(self.geo.beta23*(self.p3-self.p2),self.q3)*dx

        S = - Constant(0.1)*(self.p3-Constant(3.0))

        self.F = A1+A2+A3+B1+B2+dot(S,self.q3)*dx


