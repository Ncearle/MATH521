# coding=utf-8
"""
FEniCS program: A posteriori error estimates for Poisson's equation with
Dirichlet boundary conditions on the unit square.

The boundary values and the source term are chosen such that
    u(x,y) = x(1-x^a)y(1-y)
is the exact solution of this problem (with a parameter a >=1).
"""

from dolfin import *
from mshr import *
import numpy as np

###############################################################################
# DATA
###############################################################################

# Parameters
N = 64 # the PDE will be solved on an NxN grid (N must be even)
a = 50.
f = Expression('a*(a+1)*pow(x[0],a-1)*x[1]*(1-x[1]) + 2*x[0]*(1-pow(x[0],a))', degree=3, a=a) # source term
# NB: For solving the PDE with the optimal convergence rate, a piecewise
# constant approximation (midpoint rule) would already suffice. However, we'll
# also need f later on to compute the (ideally exact) residuals, and hence
# we're using a higher-order, piecewise cubic, interpolation here.
u_D = Constant(0.) # boundary values

# Create mesh and compute extra mesh data
coarsemesh = UnitSquareMesh(int(N/2), int(N/2), 'crossed') # for post-processing only (cheap Strategy 2)
mesh = refine(coarsemesh) # for solving the PDEs.
h = CellDiameter(mesh) # longest edge of each triangle
n = FacetNormal(mesh) # unit normal vectors on each triangle edge

###############################################################################
# SOLUTION OF THE PRIMAL PROBLEM
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

# Compute solution
u = Function(V, name='primal solution')
solve(B == F, u, bc)

# Compute quantity of interest
j = Expression('x[0]>0.5 && x[1]<0.5 ? 4. : 0.', degree=0)
# NB: j restrcted to one triangle in this mesh is constant. j is discontinuous
# over vertices and edges of some triangles, and very undesirable things will
# happen if you try to evaluate a linear or higher-order interpolating
# polynomial in a point where j is discontinuous:
# Due to inevitable rounding errors, it will be pure coincidence whether a
# function value on a triangle edge where j is discontinuous is computed as 0
# or 4. As a result, on these critical triangles  you'll interpolate data
# points that oscillate between 0 and 4. The result: 
# https://en.wikipedia.org/wiki/Runge%27s_phenomenon

Jh = assemble(j*u*dx) # quantity of interest computed from numerical solution

###############################################################################
# A POSTERIORI ERROR ESTIMATION: QUANTITY OF INTEREST
###############################################################################

# All error indicators and cell-wise norms are constant on each triangle
W = FunctionSpace(mesh, 'DG', 0) # space of piecewise constant functions
w = TestFunction(W) # basis functions which are = 1 on one triangle and = 0 elsewhere

# EXPENSIVE STRATEGY 1 ########################################################

# Function space and boundary conditions
Z = FunctionSpace(mesh, 'P', 2)
bc = DirichletBC(Z, u_D, boundary)

# Define variational problem
zh = TrialFunction(Z)
v = TestFunction(Z)
B = dot(grad(zh), grad(v))*dx
J = j*v*dx
# NB: We have to re-assemble the entire system from scratch since we're now
# using quadratic and not linear elements. Cannot re-use anything from the
# primal problem. It even has more degrees of freedom -> bigger linear system.

# Compute solution
zh = Function(Z, name='dual solution')
solve(B == J, zh, bc) # numerical solution of the dual problem

# Piecewise quadratic approximation of the exact solution
z = zh
# Its piecewise linear interpolant
Iz = interpolate(z,V)

# NB: In FEniCS we can only add / subtract functions that belong to the same
# function space. We need z-Iz, however:
#    z belongs to Z = {piecewise quadratic functions on mesh}
#    Iz belongs to V = {piecewise linear functions on mesh}
# Therefore, we re-write the piecewise linear function as a piecewise quadratic
# function (with quadratic terms = 0):
Iz = interpolate(Iz,Z)

# Local error indicators
etaT_J1 = (-div(grad(u))-f)*(z-Iz)*w*dx + jump(grad(u), n)*(z-Iz)*avg(w)*dS

# Assemble the vector which contains the local error indicators
etaT_J1_vec = np.abs(assemble(etaT_J1))

# Expensive a posteriori error estimate of the J-error
eta_J1 = np.sum(etaT_J1_vec)

print('=== A POSTERIORI ERROR ESTIMATION: QUANTITY OF INTEREST (EXPENSIVE) ==')
print('J(uh) - J(u) ≈', eta_J1)

# CHEAP STRATEGY 2 ############################################################

# Function space and boundary conditions
Z = FunctionSpace(mesh, 'P', 1)
bc = DirichletBC(Z, u_D, boundary)

# Define variational problem
zh = TrialFunction(Z)
v = TestFunction(Z)
B = dot(grad(zh), grad(v))*dx
J = j*v*dx
# NB: If we hadn't overwritten the data from the primal problem for the
# expensive error estimator, we could have re-used it here. FEniCS even
# provides the command 'adjoint' which computes the left hand side of the dual
# problem automatically from the bilinear form B of the primal problem.

# Compute solution
zh = Function(Z, name='dual solution')
solve(B == J, zh, bc) # numerical solution of the dual problem

# Piecewise quadratic interpolant on the coarse mesh
Z2 = FunctionSpace(coarsemesh, 'P', 2)
z = interpolate(zh, Z2)
# Its piecewise linear interpolant = zh
Iz = zh

# Write both as piecewise quadratics on the fine mesh, so that we can compute
# their difference.
Z = FunctionSpace(mesh, 'P', 2)
z = interpolate(z, Z)
Iz = interpolate(Iz, Z)

# Local error indicators
etaT_J2 = (-div(grad(u))-f)*(z-Iz)*w*dx + jump(grad(u), n)*(z-Iz)*avg(w)*dS

# Assemble the vector which contains the local error indicators
etaT_J2_vec = np.abs(assemble(etaT_J2))

# Cheap a posteriori error estimate of the J-error
eta_J2 = np.sum(etaT_J2_vec)

print('=== A POSTERIORI ERROR ESTIMATION: QUANTITY OF INTEREST (CHEAP) ======')
print('J(uh) - J(u) ≈', eta_J2)

###############################################################################
# EXPORT DATA FOR PLOTTING
###############################################################################

# PRIMAL SOLUTION #############################################################
File('u.pvd') << u

# DUAL SOLUTION ###############################################################

File('z.pvd') << zh

# CELL RESIDUALS ##############################################################

# We have to compute a piecewise constant function with the cell residuals as
# its function values.

# Squares of cell residuals ( ||rT||_L²(T) )²
rT = (-div(grad(u))-f)**2*w*dx
# Assemble the corresponding vector 
rT_vec = np.sqrt(np.abs(assemble(rT)))
# Define the function with these function values
rT_fun = Function(W, name='cell residuals')
rT_fun.vector().set_local(rT_vec)

rT_fun = project(rT_fun, V)
residuals = Function(V, name = 'cell residuals')
residuals.assign(rT_fun)

File('rT.pvd') << residuals

# DUAL WEIGHTS ################################################################

# Squares of the dual weights ( ||wT||_L²(T) )²
wT = (z-Iz)**2*w*dx
# Assemble the corresponding vector 
wT_vec = np.sqrt(np.abs(assemble(wT)))
# Define the function with these function values
wT_fun = Function(W, name='dual weights')
wT_fun.vector().set_local(wT_vec)

wT_fun = project(wT_fun, V)
weights = Function(V, name = 'dual weights')
weights.assign(wT_fun)

File('wT.pvd') << weights

# ERROR INDICATORS ############################################################

etaT_J2_fun = Function(W, name='cheap error indicators')
etaT_J2_fun.vector().set_local(etaT_J2_vec)

etaT_J2_fun = project(etaT_J2_fun, V)
indicators = Function(V, name = 'error indicators')
indicators.assign(etaT_J2_fun)

File('etaT.pvd') << indicators