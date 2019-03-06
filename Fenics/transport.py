"""
FEniCS tutorial demo program: transport equation with Dirichlet conditions.
"""

from __future__ import print_function
from dolfin import *

# Create mesh and define function space
N = 50 # 100 200 400 800
mesh = RectangleMesh(Point(0., 0.), Point(1., 0.2), 50, 10, 'crossed') # h = 1/N
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('exp(-2.*x[0])*sin(pi*x[1]/0.2)', degree=4)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
D = Constant(1E-2)
a = Constant((1., 0.))
r = Constant(2.)
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('(pow(pi, 2.)/4. - 4E-2)*exp(-2.*x[0])*sin(pi*x[1]/0.2)', degree=2)
B = dot(D*grad(u), grad(v))*dx - dot(u*a, grad(v))*dx + r*u*v*dx
L = f*v*dx

# Compute solution
u = Function(V, name='concentration')
solve(B == L, u, bc)

# Save solution to file in VTK format
vtkfile = File('transport/solution.pvd')
vtkfile << u

# Compute error in L2-norm and H1-norm
error_L2 = errornorm(u_D, u, 'L2')
error_H1 = errornorm(u_D, u, 'H1')

# Print errors
print('error_L2  =', error_L2)
print('error_H1  =', error_H1)
