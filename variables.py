### the variables used in the simulation stored in a single file

## imports
import variables as var
import numpy as np

## constants
# dimensionality of the system
dimension = 3.

# number of particles in the system
numParticlesAxis = 10
numParticles = numParticlesAxis**3

# time step
deltaT = 1.

# length of the box side of the box
boxSize = 5.

## classes
class Particles(object):


    # initialize the particles
    def __init__(self):

        # initiate positions and velocities vectors
        positions = np.empty(numParticles, dimension)
        velocities = np.empty()

        #TODO implement raster distribution
        #positions = np.empty(var.numParticles, var.dimension)

        #TODO implement velocity distribution
        #velocities = np.empty(var.numParticles, var.dimension)
    def set_init_pos(self):
        #Using cubic lattice
        #volumeBox   = boxSize**dimension

        numAxis = float(numParticlesAxis)
        side        = boxSize/(numAxis-1)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,numParticles,increment):
            positions[j:j+increment, 0] = posAxis
            if np.mod(j,increment**2)==0:
                positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if np.mod(j,increment)==0:
                positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1
        #Using fcc lattice, we know density of one fcc cube is 14/a**3 units/m3
        #volumeBox   = boxSize**dimension
        #partDenisty = numParticles/volumeBox
        #sideFcc     = (14./partDenisty)**(1./3)



    # update the particles
    def update(self):
        pass
        #TODO implement update functions
        # updateParticles(self)
        # updateVelocities

#   def updateParticles(self)
#   def updateVelocities(self)
