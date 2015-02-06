### the variables used in the simulation stored in a single file

## imports
import numpy as np
import Potentials as Pot

## constants
## Constants are given in nm, ps units (i.e distance 1 = 1 nm, time 1 = 1 ps)
# dimensionality of the system
dimension = 3

# number of particles in the system
numParticlesAxis = 10
numParticles = numParticlesAxis**3

# time step
deltaT = 1.*10**6 # 1 microsecond, time step in ps, rescale the rest to fit

# length of the box side of the box
boxSize = 50. # Box size in nm, rescale everything else to fit

# Temperature (in Kelvin)
T = 300

# Mass
m = 6.64648*10**(-27) # 6.64648*10**(-27) kg

# Boltzmann constant
kB = 1.3806488*10**(-29) # kB = 1.38*10^-5 [nm^2 kg ps^-2 K^-1]

# Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)
a = (3.0 * kB * T) / m

# Lennard-Jones depth of potential well
eps = 10.22 * 10**(-6) * kB # Helium Cyrogenics - Steven van Sciver, eps/kB = 10.22
# eps = [J] = [kg m^2 s^-2] = 10^-6 [kg nm^2 ps^-2]

# Lennard-Jones distance at which potential is minimal
rMin = 0.2869 # Helium Cyrogenics - Steven van Sciver, rMin = 0.2869 nm

## classes
class Particles(object):


    # initialize the particles
    def __init__(self):

        self.initPositions()
        self.initVelocities()

    def initPositions(self):
        #Using cubic lattice
        #volumeBox   = boxSize**dimension

        self.positions = np.zeros((numParticles,dimension))
        numAxis = float(numParticlesAxis)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,numParticles,increment):
            self.positions[j:j+increment, 0] = posAxis      #For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                            #Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:                         #Here add the 'z' coordinates after n**2
                self.positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if (j%increment)==0:                  #Add the 'y' coordinates, after n repetitions of x0->xn
                self.positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1
        #Using fcc lattice, we know density of one fcc cube is 14/a**3 units/m3
        #volumeBox   = boxSize**dimension
        #partDenisty = numParticles/volumeBox
        #sideFcc     = (14./partDenisty)**(1./3)

    # update the particles
    def update(self, dT):

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
        FORCE, Potential = Pot.Len_Jones(self.positions)


        self.velocities += FORCE / m * dT


    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed
        # (i.e. Gaussian distribution with mean=0 and sigma(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))
