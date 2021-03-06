# coding=utf-8
"""
FEniCS program: Solution of the wave equation with homogeneous Dirichlet (reflective)
boundary conditions.
"""

from fenics import *
from mshr import *

# Create a mesh on the unit disk
disk = Circle(Point(0., 0.), 1.)
mesh = generate_mesh(disk, 100) # h ~ 1/50

# Function space and boundary condition
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant(0.), boundary)

# Problem data
t0 = 0. # initial time
T = 5. # final time
t = t0 # current time
c = 1. # propagation speed
u0 = interpolate(Expression('pow(x[0], 2)+pow(x[1], 2) < 1./16. ? pow(1.-16.*(pow(x[0], 2)+pow(x[1], 2)), 2) : 0.', degree=1), V) # initial displacement
v0 = interpolate(Constant(0.), V) # initial velocity

# Parameters of the time-stepping scheme
tsteps = 500 # number of time steps
dt = T/tsteps # time step size
theta = 0.5 # degree of implicitness

# Define the variational problem
w = TrialFunction(V) # w = u in the 1st equation and w = v in the 2nd equation
z = TestFunction(V)
B1 = w*z*dx + (Constant(theta)*dt*c)**2*dot(grad(w),grad(z))*dx # LHS of the 1st equation
B2 = w*z*dx # LHS of the 2nd equation

A1 = assemble(B1)
A2 = assemble(B2)

bc.apply(A1)
bc.apply(A2)

# Solvers

# directSolver1 = LUSolver(A1)
# directSolver2 = LUSolver(A2)

iterativeSolver1 = KrylovSolver(A1, 'gmres', 'jacobi')
iterativeSolver2 = KrylovSolver(A2, 'gmres', 'jacobi')

# Parameters for the solvers

# directSolver2.parameters['symmetric'] = True # True for Cholesky, False for LU
# directSolver1.parameters['symmetric'] = True # True for Cholesky, False for LU

abstol = 10 ** -9
reltol = 10 ** -6
maxiter = 1000
iterativeSolver1.parameters['nonzero_initial_guess'] = True
iterativeSolver2.parameters['nonzero_initial_guess'] = True

iterativeSolver1.parameters['absolute_tolerance'] = abstol
iterativeSolver2.parameters['absolute_tolerance'] = abstol
iterativeSolver1.parameters['relative_tolerance'] = reltol
iterativeSolver2.parameters['relative_tolerance'] = reltol
iterativeSolver1.parameters['maximum_iterations'] = maxiter
iterativeSolver2.parameters['maximum_iterations'] = maxiter

iterativeSolver1.parameters['monitor_convergence'] = False
iterativeSolver2.parameters['monitor_convergence'] = False
# -------------------------------------------------------
# I guess I just really want it to stop at some point :)
# -------------------------------------------------------


# Set initial data
u = Function(V, name='Displacement')
v = Function(V, name='Velocity')
u.assign(u0)
v.assign(v0)

# Write initial data to file
displacement = File('wave/thetaCN_IterativeGMRESjacobi.pvd')
displacement << (u, t)

# Time stepping
for k in range(tsteps):

    # Current time
    t = t0 + (k+1)*dt
    print('Step = ', k+1, '/', tsteps , 'Time =', t)

    # System for the displacement
    L1 = ((u0 + dt*v0)*z - Constant(theta*(1.-theta)*(dt*c)**2)*dot(grad(u0),grad(z)))*dx
    b1 = assemble(L1)
    bc.apply(b1)
    # directSolver1.solve(u.vector(), b1)
    iterativeSolver1.solve(u.vector(), b1)
    # directSolver1.set_operator(A1)
    # solve(B1 == L1, u, bc)

    # System for the velocity
    L2 = (v0*z - Constant(dt*c**2)*dot(grad(Constant(theta)*u + Constant(1.-theta)*u0), grad(z)))*dx
    b2 = assemble(L2)
    bc.apply(b2)
    # directSolver2.solve(v.vector(), b2)
    iterativeSolver2.solve(v.vector(), b2)
    # directSolver2.set_operator(A2)
    # solve(B2 == L2, v, bc)

    # Write data to file
    displacement << (u, t)

    # Update
    u0.assign(u)
    v0.assign(v)

# list_krylov_solver_preconditioners()
