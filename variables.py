### the variables used in the simulation stored in a single file

## imports
import variables as var
import numpy as np

## constants
# dimensionality of the system
dimension = 3

# number of particles in the system
numParticles = 1000

# time step
deltaT = 1

# length of the box side of the box
boxSize = 5

## classes
class Particles(object):


    # initialize the particles
    def __init__(self):

        # initiate positions and velocities vectors
        self.positions = np.empty()
        self.velocities = np.empty()

        #TODO implement raster distribution
        positions = np.empty(var.numParticles, var.dimension)

        #TODO implement velocity distribution
        velocities = np.empty(var.numParticles, var.dimension)

    # update the particles
    def update(self):
        pass
        #TODO implement update functions
        # updateParticles(self)
        # updateVelocities

#   def updateParticles(self)
#   def updateVelocities(self)
