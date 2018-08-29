from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
from geometry import left_ventrical_geometry

class darcy_model():
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
        self.P = FiniteElement('P',el,2)

        self.element = MixedElement([P,P,P])
        self.FS = FunctionSpace(mesh,element)

    def test_functions(self):

        self.q1,self.q2,self.q3 = TestFunctions(FS)
    
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
        
        A1 = -self.geo.K1 * dot(grad(self.p1),grad(self.q1))*dx
        A2 = -self.geo.K2 * dot(grad(self.p2),grad(self.q2))*dx
        A3 = -self.geo.K3 * dot(grad(self.p3),grad(self.q3))*dx

        B1 = dot(self.geo.beta12*(p1-p2),q1)*dx + dot(self.geo.beta12*(p2-p1),q2)*dx
        B2 = dot(self.geo.beta23*(p2-p3),q2)*dx + dot(self.geo.beta23*(p3-p2),q3)*dx

        S = - Constant(0.1)*(p3-Constant(3.0))

        self.F = A1+A2+A3+B1+B2+S


