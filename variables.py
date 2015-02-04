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
        self.updateVelocities(dT)

    # update the particle positions
    def updateParticles(self, dT):

        self.positions += self.velocities * dT

        # translate the particles outside of the box
        # +boxSize if positionComponent < 0, -boxSize if positionComponent > 5
        posTranslation = self.positions < 0
        negTranslation = self.positions > boxSize
        self.positions[posTranslation] += boxSize
        self.positions[negTranslation] -= boxSize


    def updateVelocities(self, dT):

        # TODO calculate forces on particles with the positions (Leo):
        # function:
        # in -> positionVectors (self.position)
        # out -> forceVectors (FORCE)
        #
        # both numParticles by dimension matrices
        FORCE = 0.1

        self.velocities += FORCE * (dT**2)
