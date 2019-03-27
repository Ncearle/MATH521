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
theta = 1. # degree of implicitness

# Define the variational problem
w = TrialFunction(V) # w = u in the 1st equation and w = v in the 2nd equation
z = TestFunction(V)
B1 = w*z*dx + (Constant(theta)*dt*c)**2*dot(grad(w),grad(z))*dx # LHS of the 1st equation
B2 = w*z*dx # LHS of the 2nd equation

# Set initial data
u = Function(V, name='Displacement')
v = Function(V, name='Velocity')
u.assign(u0)
v.assign(v0)

# Write initial data to file
displacement = File('wave/thetaBE.pvd')
displacement << (u, t)

# Time stepping
for k in range(tsteps):

    # Current time
    t = t0 + (k+1)*dt
    print('Step = ', k+1, '/', tsteps , 'Time =', t)

    # System for the displacement
    L1 = (u0 + dt*v0)*z*dx - (Constant(theta)*(1-Constant(theta))*(dt*c)**2)*dot(grad(u0),grad(z))*dx
    solve(B1 == L1, u, bc)

    # System for the velocity
    L2 = v0*z*dx - dt*c*c*dot(Constant(theta)*grad(u) + (1-Constant(theta))*grad(u0),grad(z))*dx
    solve(B2 == L2, v, bc)

    # Write data to file
    displacement << (u, t)

    # Update
    u0.assign(u)
    v0.assign(v)
