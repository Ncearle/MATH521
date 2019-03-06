# coding=utf-8
"""
FEniCS program: Goal-oriented mesh adaptation for Poisson's equation with the
dual weighted residual method.

The boundary values and the source term are chosen such that
    u(x,y) = x(1-x^a)y(1-y)
is the exact solution of this problem (with a parameter a >=1).
"""

from dolfin import *
import numpy as np

###############################################################################
# DATA
###############################################################################

# Parameters
tol = 1E-4 # tolerance up to which the quantity of interest should be computed
#tol = 1E-5 # Now error pollution ('diffusion of errors into R') becomes important! Compare the final meshes for these different tolerances in ParaView!
N = 64 # the PDE will initially be solved on an NxN grid (N must be even)
a = 50.
f = Expression('a*(a+1)*pow(x[0],a-1)*x[1]*(1-x[1]) + 2*x[0]*(1-pow(x[0],a))', degree=3, a=a) # source term
u_D = Constant(0.) # boundary values

# Create mesh and compute extra mesh data
mesh = UnitSquareMesh(N, N, 'crossed')

###############################################################################
# SOLUTION OF PROBLEM (P)
###############################################################################

# Function space and boundary conditions
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
B = dot(grad(u), grad(v))*dx
F = f*v*dx

u = Function(V, name='primal solution')

# Define quantity of interest
j = Expression('x[0]>0.5 && x[1]<0.5 ? 4. : 0.', degree=0)
J = j*u*dx

# Solve problem, compute error indicators, refine the mesh and repeat until the
# estimated error in J is < tol
problem = LinearVariationalProblem(B, F, u, bc)
solver = AdaptiveLinearVariationalSolver(problem, J)
solver.solve(tol)

# Export solution on finest mesh
File('adaptivePoisson.pvd') << u.leaf_node()
# To see the mesh in ParaView: switch from 'Surface' to "Surface with Edges' or
# 'Wireframe'.