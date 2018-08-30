from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
from geometry import left_ventrical_geometry

class Perfusion_model:
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

class Advection_diffusion:
    
    def __init__(self,geo=None):

        if geo is None:
            self.geo = left_ventrical_geometry()
        self.function_space()
        self.test_functions()
        self.funcs()
        self.multi_comp_darcy_equation()
        self.intercompartment_consentrationflow()
        self.advection_diffusion_equation()

    def function_space(self):

        self.el = tetrahedron
        self.P = FiniteElement('P',self.el,2)

        self.element = MixedElement([self.P,self.P,self.P])
        self.FS = FunctionSpace(self.geo.mesh,self.element) #Functionspace Pressure
        self.FSC = FunctionSpace(self.geo.mesh,self.element) #Functionspace Consentration
    

    def test_functions(self):

        self.q1,self.q2,self.q3 = TestFunctions(self.FS)
        self.v1,self.v2,self.v3 = TestFunctions(self.FSC)

    
    def funcs(self):

        self.p = Function(self.FS)

        self.vd1 = Function(self.FS)i #Darcy Velocity
        self.vd2 = Function(self.FS)
        self.vd3 = Function(self.FS)

        self.c = Function(self.FSC) #Current consentration
        self.c_n = Function(self.FSC)#Initial consentration

        self.p1, self.p2, self.p3 = split(self.p)
        self.c1, self.c2, self.c3 = split(self.c)
        self.c_n1,self.c_n2,self.c_n3 = split(c_n)


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

    def intercompartment_consentrationflow(self):
        """
        This equation transport a ratio of all concentration in comp1 to comp2 etc.
        
        R1: Tranfers consentration from comp1->2
        R2: Tranfers consentration from comp2->3
        """
        R1 = geo.R_12*self.c1*(self.v1-self.v2)*dx
        R2 = geo.R_23*self.c2*(self.v2-self.v3)*dx

        self.F3 = R1+R2


    def advection_diffusion_equation(self):
        """
        The multicompartment advection_diffusion equation

        C1-3: Timedependent consentration
        D1-3: Intracompartment Advection
        E1-3: Intracompartment Diffusion
        OBS, k not defined(dt)
        """

        C1 = ((self.c1 - self.c_n1)/geo.k)*self.v1*dx
        C2 = ((self.c2 - self.c_n2)/geo.k)*self.v2*dx
        C3 = ((self.c3 - self.c_n3)/geo.k)*self.v3*dx

        D1 = dot(self.vd1,grad(self.c1))*self.v1*dx
        D2 = dot(self.vd2,grad(self.c2))*self.v2*dx
        D3 = dot(self.vd3,grad(self.c3))*self.v3*dx

        E1 = -geo.D*dot(grad(self.c1),grad(self.v1))*dx
        E2 = -geo.D*dot(grad(self.c2),grad(self.v2))*dx
        E3 = -geo.D*dot(grad(self.c3),grad(self.v3))*dx

        self.F2 = C1+C2+C3+D1+D2+D3+E1+E2+E3

