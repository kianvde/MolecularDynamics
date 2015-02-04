### the variables used in the simulation stored in a single file

## imports
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

# Temperature (in Kelvin)
T = 300

# Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)
a = 1

# Lennard-Jones depth of potential well
eps = 1.0

# Lennard-Jones distance at which potential is minimal
rMin = 2.0**(1.0/6.0)

## classes
class Particles(object):


    # initialize the particles
    def __init__(self):

        #TODO implement raster distribution
        self.positions = np.zeros((numParticles, dimension))

        #TODO implement velocity distribution
        self.initVelocities()

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
        FORCE = 0.

        self.velocities += FORCE * (dT**2)


    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed
        # (i.e. Gaussian distribution with mean=0 and sigma(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))
