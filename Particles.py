import numpy as np
import variables as var

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