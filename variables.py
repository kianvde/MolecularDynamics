### the variables used in the simulation stored in a single file

## imports
import numpy as np
from potentials import lennardJones

## constants
## Constants are given in nm, ps units (i.e distance 1 = 1 nm, time 1 = 1 ps)
# dimensionality of the system
dimension = 3

# number of particles in the system
numParticlesAxis = 3
numParticles = numParticlesAxis**3

# time step
deltaT = 0.1 # 1 microsecond, time step in ps, rescale the rest to fit

# length of the box side of the box
boxSize = 10.0 # Box size in nm, rescale everything else to fit

# interaction range for the particles
rCutoff = 0.2*boxSize

# Temperature (in Kelvin)
T = 1.0

# Mass
m = 6.64648*10**(-27) # 6.64648*10**(-27) kg

# Boltzmann constant
kB = 1.3806488*10**(-29) # kB = 1.38*10^-29 [nm^2 kg ps^-2 K^-1]

# Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)
a = ((kB * T) / m)**0.5

# Lennard-Jones depth of potential well
eps =  10.22 * kB # Helium Cyrogenics - Steven van Sciver, eps/kB = 10.22 ADDED FACTOR FOR INCREASED INTERACTION
# eps = [J] = [kg m^2 s^-2] = 10^-6 [kg nm^2 ps^-2]

# Lennard-Jones distance at which potential is minimal
rMin = 0.2869 # Helium Cyrogenics - Steven van Sciver, rMin = 0.2869 nm

## Particles
class Particles(object):


    ## initialization ##
    def __init__(self):

        self.initPositions()
        self.initVelocities()
        self.force = np.zeros(np.shape(self.positions))
        self.energy = 0
        self.potential = 0
        self.temperature = T
        self.initposs = self.initPositions()
        self.initvelocc = self.initVelocities()
    # initialize the particle positions
    def initPositions(self):
        # Using cubic lattice
        # volumeBox   = boxSize**dimension

        self.positions = np.zeros((numParticles,dimension))
        numAxis = float(numParticlesAxis)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,numParticles,increment):
            self.positions[j:j+increment, 0] = posAxis      # For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                            # Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:                         # Here add the 'z' coordinates after n**2
                self.positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if (j%increment)==0:                  #Add the 'y' coordinates, after n repetitions of x0->xn
                self.positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1
        # Using fcc lattice, we know density of one fcc cube is 14/a**3 units/m3
        # volumeBox   = boxSize**dimension
        # partDenisty = numParticles/volumeBox
        # sideFcc     = (14./partDenisty)**(1./3)
        return self.positions

    # initialize the particle velocities
    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed, in nm / ps
        # (i.e. Gaussian distribution with mean=0 and std(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))

    ## update functions ##
    def update(self, dT):

        # 4th-order symplectic RK4 integrator (Forest & Ruth, 1989)
        X = (1.0/6.0) * (2.0**(1.0/3.0) + 2.0**(-1.0/3.0) - 1)
        self.updateParticles(dT, X + 0.5)
        self.updateForce()
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, -X)
        self.updateForce()
        self.updateVelocities(dT, -4.0*X - 1.0)
        self.updateParticles(dT, -X)
        self.updateForce()
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, X + 0.5)
        self.updateForce()

        v2 = np.sum(self.velocities**2,axis=None)
        self.energy = (0.5 * m * v2 + self.potential) * 10**6 # Energy in Joule
        self.temperature = (2.0/3.0) * (0.5 * m * v2) / (numParticles * kB)  # Paper Verlet 1967

    # update the particle positions
    def updateParticles(self, dT, D):

        self.positions += D * self.velocities * dT # + 0.5 * (self.force / m) * (dT**2)

        # translate the particles outside of the box
        # +boxSize if positionComponent < 0, -boxSize if positionComponent > 5
        posTranslation = self.positions < 0
        negTranslation = self.positions > boxSize
        self.positions[posTranslation] += boxSize
        self.positions[negTranslation] -= boxSize

    # update the particle velocities and forces
    def updateVelocities(self, dT, C):

        self.velocities += C * (self.force / m) * dT

    # update the forces on the particles and calculate the potential energy
    def updateForce(self):

        self.force, self.potential = lennardJones(self.positions)



    ## UNUSED ##

    # get the particles within a distance imageSize from the box boundaries and translate them to outside
    # the box for force computation to simulate an infinite volume
    #
    # return -> a by d matrix with all the imaged particles
    def getTranslatedImages(self):

        p = self.positions                          # original particle position matrix
        t = np.zeros((np.shape(p)))                 # matrix containing translation values per coordinate
        t[p < rCutoff] = boxSize
        t[p > boxSize - rCutoff] = -boxSize
        tBool = t != 0
        numWalls = np.sum(tBool, axis=1)             # number of walls the particle is close to

        # select cases with one or more translations and add translations and translate
        # all coordinates close to the wall
        translatedImages = p[numWalls >= 1, :] + t[numWalls >= 1, :]

        # select cases with two or more translations and translate for single
        # coordinates close to the wall
        p2 = p[numWalls >= 2, :]
        t2 = t[numWalls >=2, :]

        for column in range(3):
            pCol = p2[:,column]
            tCol = t2[:,column]
            pConcat = p2[tCol != 0, :]
            pConcat[:, column] += tCol[tCol != 0]
            translatedImages = np.concatenate((translatedImages, pConcat), axis=0)

        # select cases with 3 translations and translate for two
        # coordinates close to the wall
        p3 = p[numWalls == 3, :]
        t3 = t [numWalls == 3, :]
        for column in range(3):
            t2Col = t3.copy()
            t2Col[:, column] = 0
            translatedImages = np.concatenate((translatedImages, p3 + t2Col), axis=0)

        return translatedImages
