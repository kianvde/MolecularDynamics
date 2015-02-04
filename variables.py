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

        #TODO implement raster distribution
        self.positions = np.zeros((var.numParticles, var.dimension))

        #TODO implement velocity distribution
        self.velocities = np.ones((var.numParticles, var.dimension))

    # update the particles
    def update(self, dT):

        #TODO implement update functions
        self.updateParticles(dT)
        # updateVelocities

    # update the particle positions
    def updateParticles(self, dT):
        self.positions += self.velocities * dT

#   def updateVelocities(self)
