#!/usr/bin/python2

# Derive the jacobian matrix of the 1-D Euler equations.
# Of course, using the symbolic manipulator.

from sympy import *

# Define the functions.

rho  = Symbol('rho')
rhou = Symbol('rhou')
rhov = Symbol('rhou')
e    = Symbol('e')
g    = Symbol('g')

# Define the primitives as a function of the Q vector.

rho = rho
u   = rhou/rho
v   = rhov/rho
e   = e

# Pressure is a function of the primitives dumb ass !

p = (g - 1.0)*(e - 0.5*rho*(u**2.0 + v**2.0))

# Define the vectors.

Q = Matrix([rho, rho*u, rho*v, e])

E = Matrix([rho*u, rho*u**2 + p, rho*u*v, u*(e + p)])

F = Matrix([rho*v, rho*u*v, rho*v**2 + p, v*(e + p)])

# First line of the Jacobian.

elem_00 = diff(E[0],Q[0])
elem_01 = diff(E[0],Q[1])
elem_02 = diff(E[0],Q[2])
elem_03 = diff(E[0],Q[3])

# Second line of the Jacobian.

elem_10 = diff(E[1],Q[0])
elem_11 = diff(E[1],Q[1])
elem_12 = diff(E[1],Q[2])
elem_13 = diff(E[1],Q[3])

# Third line of the Jacobian.

elem_20 = diff(E[2],Q[0])
elem_21 = diff(E[2],Q[1])
elem_22 = diff(E[2],Q[2])
elem_23 = diff(E[2],Q[3])

# Fourth line of the Jacobian.

elem_30 = diff(E[3],Q[0])
elem_31 = diff(E[3],Q[1])
elem_32 = diff(E[3],Q[2])
elem_33 = diff(E[3],Q[3])

# Print in Latex for us to simplify.

print latex(elem_00)

