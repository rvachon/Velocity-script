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

"""
Declaration of constants:

define the length X and depth Y of the model [m]

"""
X=1.e+6
Y=1.e+6
"""
#define the x and y number of nodes that will define the grid.
"""
Xi=11
Yi=11

"""
define the size of the step between two nodal points [m]
"""
Xstep=X/Xi
Ystep=Y/Yi

"""
define initial velocities x and y components [m/s]
"""
Vx0=1.1e-9*X/2/Y
Vy0=1.0e-9

"""
defining x and y coordinate of each nodal points and
creating a mesh using the meshgrid numpy function.
"""
x,y=np.meshgrid(np.linspace(0,X,Xi),np.linspace(0,Y,Yi))

"""
calculation of the horizontal and vertical component of the velocity at each
nodes of coordinate x and y through the entire box.
Setup velocity field corresponding to circulation with central upwelling.
"""
Vx=-Vx0*np.sin(np.pi*x/X*2)*np.cos(np.pi*y/Y)
Vy=-Vy0*np.sin(np.pi*y/Y)*np.cos(np.pi*x/X*2);


"""
calculation of the Velocity magnitude.
"""

Vt = np.sqrt(Vx**2+Vy**2)

"""
plot of the Vt output
"""

plt.figure(1)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Velocity magnitude')
plt.text(100000,900000, 'Xi= 3001 and Yi=3001')
plt.axis([0,X,0,Y])
plt.pcolormesh(x,y,Vt)
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

plt.figure(2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Strain rate 2nd Invariant')
plt.text(100000,900000, 'Xi= 3001 and Yi=3001')
plt.axis([0,X,0,Y])
plt.pcolormesh(x,y,epsii)
plt.colorbar()
plt.show()