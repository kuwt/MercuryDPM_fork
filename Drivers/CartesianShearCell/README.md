# The Cartesian Shear Cell Simulations

This system can be used to study shear bands. 

To execute the simulation, you need to run three codes in this order: 
 - CSCWalls
 - CSCInit
 - CSCRun

## CSCWalls
This code is used to define the particle and material properties and the geometry:

We use particles that are non-dimensionalised such that their diameter d=1 and mass m=1. Gravitational acceleration is also non-dimensionalised to g=(0,0,1); thus, the gravitational time-scale is t_g=sqrt(d/g)=1. We use a linear spring-damper contact law with collision time t_c=t_g/20, restitution r=0.88, and a sliding friction of mu=0.5. 

This code defines a box of size W=30 in x-, L=20 in y-, and H=30 in z-direction. The walls consist of flat walls with particles attached to the wall to make it more frictional.

## CSCInit

This code restarts from the box created in CSCWalls, and fills the box with particles

## CSCRun

This code restarts from the box created in CSCInit, and adds a shear velocity by moving the walls at a steady velocity in y-direction: The right wall and the right half of the bottom wall is moved at a velocity (0,v,0), the left wall and the left half of the bottom wall is moved at a velocity (0,-v,0).

## Author

This code was written by Thomas Weinhart