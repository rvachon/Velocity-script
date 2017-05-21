
# coding: utf-8
"""
Computation and visualisation of velocity field,
strain rate, deviatoric strain rate, 
and second strain rate invariant at each point of a grid.


Import numpy and matplotlib.plot library and change prefixe to np and
plt respectively. Using numpy library should improve computation time 
vs Matlab. The Matplotlib library is used for graphical scientific plot.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sc

mat = sc.loadmat('strain_rate_var.mat')['epsII']


"""
Declaration of constants:

define the length X and depth Y of the model [m]

"""
X=1.e+6
Y=1.5e+6

#define the x and y number of nodes that will define the grid.
Xi=31.
Yi=31.

#define the size of the step between two nodal points [m]
Xstep=X/(Xi-1.)
Ystep=Y/(Yi-1.)

#define initial velocities x and y components [m/s]
Vx0=1.1e-9*X/2/Y
Vy0=1.0e-9

#defining x and y coordinate of each nodal points and
#creating a mesh using the meshgrid numpy function.
x,y = np.meshgrid(np.linspace(0,X,Xi),np.linspace(0,Y,Yi))
print(x)
print(y)

#calculation of the horizontal and vertical component of the velocity at each
#nodes of coordinate x and y through the entire box.
Vx=-Vx0*np.sin(np.pi*x/X*2)*np.cos(np.pi*y/Y)
Vy=-Vy0*np.sin(np.pi*y/Y)*np.cos(np.pi*x/X*2);

plt.imshow(Vx, extent = (0,X,0,Y))
plt.colorbar()
plt.figure()
plt.imshow(Vy)
plt.colorbar()

# calculation of the Velocity magnitude.
Vt = np.sqrt(Vx**2+Vy**2)

plt.figure()
plt.imshow(Vt)
plt.colorbar()

"""
Computing partial derivatives 
dVx/dx
"""
dvxdx=-Vx0*np.pi/X*2*np.cos(np.pi*x/X*2)*np.cos(np.pi*y/Y);

"""
# dVx/dy
"""
dvxdy=Vx0*np.pi/Y*np.sin(np.pi*x/X*2)*np.sin(np.pi*y/Y);

"""
# dVy/dy
"""
dvydy=Vy0*np.pi/Y*np.cos(np.pi*y/Y)*np.cos(np.pi*x/X*2);

"""
# dVy/dx
"""
dvydx=-Vy0*np.pi/X*2*np.sin(np.pi*y/Y)*np.sin(np.pi*x/X*2);

"""
Computing EPSkk=dVx/dx+dVy/dy=div(v)
"""
epskk=dvxdx+dvydy;

"""        
Computing EPS'xx
"""
eps1xx=dvxdx-1/3*epskk;
"""
Computing EPS'yy
"""
eps1yy=dvydy-1/3*epskk;

"""
Computing EPSxy
"""
eps1xy=1/2*(dvxdy+dvydx);

"""
Computing EPSII
"""
def strainII(A,B,C):
   epsii = (1/2*(A**2+B**2)+C**2)**0.5;
   return epsii;
	
epsii= strainII(eps1xx,eps1yy,eps1xy);

def test_strainII():
	assert np.all(strainII(eps1xx,eps1yy,eps1xy) == mat)